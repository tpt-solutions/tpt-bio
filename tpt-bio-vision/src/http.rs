//! HTTP range reads over a hand-written HTTP/1.1 client.
//!
//! [`HttpSource`] implements [`ByteSource`] by issuing `Range: bytes=start-end`
//! requests against a single object URL. Its length is learned on construction.
//!
//! # Limitations
//!
//! - **Plain `http://` only.** This project has no TLS. An `https://` URL is rejected
//!   with [`io::ErrorKind::Unsupported`]; nothing ever falls back to `http://`.
//! - **Signing over plain HTTP is insecure.** When [`S3Signing`] is configured, the
//!   `SigV4` signature and the session token travel in cleartext. Anyone on the path
//!   can read them and replay the signed requests until the signature expires.
//!   Use this only on a trusted network or against a local gateway.
//! - Every request sends `Connection: close`. Connection reuse is not implemented.
//! - Chunked transfer-encoding is rejected. The server must send `Content-Length`.
//! - Redirects are not followed. Any non-2xx status other than 5xx is an error.
//! - If a server ignores `Range` and returns `200`, the reader discards the prefix
//!   before the requested range. This works but transfers the whole prefix.
//! - Userinfo (`user:pass@host`) in URLs is rejected. Use [`S3Signing`] instead.
//!
//! # Retries
//!
//! Connection failures, timeouts, 5xx responses and truncated bodies are retried
//! according to [`RetryPolicy`]. 4xx responses and protocol violations are not.
//! Range reads are idempotent, so a retry is always safe.
//!
//! # Presigned and public URLs
//!
//! A presigned URL is used as given: its query string is sent byte for byte and is
//! not re-signed. A public-bucket URL is sent with no `Authorization` header.
//! Passing a presigned URL together with [`S3Signing`] is an error.

use std::io::{self, BufRead, BufReader, Read, Write};
use std::net::{TcpStream, ToSocketAddrs};
use std::thread;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use crate::sigv4::{self, Credentials, Request, Scope, SignedHeaders};
use crate::source::{self, ByteSource};

/// Upper bound on the size of a response head (status line and headers).
const MAX_HEAD_BYTES: usize = 64 * 1024;

/// Service name used in `SigV4` scopes for S3 requests.
const S3_SERVICE: &str = "s3";

/// Retry behaviour for idempotent range reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RetryPolicy {
    /// Retries after the first attempt. Zero disables retries.
    pub max_retries: u32,
    /// Delay before the first retry. Doubles on each retry.
    pub initial_backoff: Duration,
    /// Upper bound on any single delay.
    pub max_backoff: Duration,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        Self {
            max_retries: 3,
            initial_backoff: Duration::from_millis(100),
            max_backoff: Duration::from_secs(2),
        }
    }
}

impl RetryPolicy {
    /// Delay before retry number `retry` (zero-based): `initial * 2^retry`, capped at
    /// `max_backoff`.
    #[must_use]
    pub fn backoff(&self, retry: u32) -> Duration {
        let factor = 1u32.checked_shl(retry.min(31)).unwrap_or(u32::MAX);
        self.initial_backoff.saturating_mul(factor).min(self.max_backoff)
    }
}

/// Static credentials for signing requests to an S3 object over plain HTTP.
///
/// Both path-style (`http://host/bucket/key`) and virtual-hosted-style
/// (`http://bucket.host/key`) URLs are signed the same way: the signer uses the
/// path and host of the URL exactly as given.
#[derive(Clone)]
pub struct S3Signing {
    /// Access key ID.
    pub access_key: String,
    /// Secret access key.
    pub secret_key: String,
    /// Session token for temporary credentials.
    pub session_token: Option<String>,
    /// Region the bucket lives in, for example `us-east-1`.
    pub region: String,
}

impl std::fmt::Debug for S3Signing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("S3Signing")
            .field("access_key", &self.access_key)
            .field("secret_key", &"<redacted>")
            .field("session_token", &self.session_token.as_ref().map(|_| "<redacted>"))
            .field("region", &self.region)
            .finish()
    }
}

/// Options for [`HttpSource::open_with`].
#[derive(Debug, Clone)]
pub struct HttpOptions {
    /// Retry behaviour for range reads and the initial length request.
    pub retry: RetryPolicy,
    /// Connect, read and write timeout for each request. Must be non-zero.
    pub timeout: Duration,
    /// When set, each request is signed with `SigV4` for the `s3` service.
    pub s3: Option<S3Signing>,
}

impl Default for HttpOptions {
    fn default() -> Self {
        Self { retry: RetryPolicy::default(), timeout: Duration::from_secs(30), s3: None }
    }
}

/// A parsed `http://` URL.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HttpUrl {
    /// Host name or IP address, without IPv6 brackets.
    pub host: String,
    /// TCP port (80 when the URL has none).
    pub port: u16,
    /// Authority as written in the URL, used for the `Host` header.
    pub authority: String,
    /// Request target: path plus query, starting with `/`. Sent verbatim.
    pub target: String,
}

/// Parses an `http://` URL.
///
/// # Errors
/// Returns [`io::ErrorKind::Unsupported`] with the message
/// `https is not supported: no audited pure-Rust TLS stack is available` for `https://`
/// URLs. Returns [`io::ErrorKind::InvalidInput`] for other schemes, userinfo, an empty
/// host or an invalid port.
pub fn parse_url(url: &str) -> io::Result<HttpUrl> {
    let (scheme, rest) =
        url.split_once("://").ok_or_else(|| invalid("URL must start with http://"))?;
    if scheme.eq_ignore_ascii_case("https") {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "https is not supported: no audited pure-Rust TLS stack is available",
        ));
    }
    if !scheme.eq_ignore_ascii_case("http") {
        return Err(invalid(format!("unsupported URL scheme `{scheme}`")));
    }
    // The fragment is never sent to the server.
    let rest = rest.split('#').next().unwrap_or_default();
    let split = rest.find(['/', '?']).unwrap_or(rest.len());
    let (authority, target) = rest.split_at(split);
    if authority.contains('@') {
        return Err(invalid("userinfo in URLs is not supported; use S3Signing for credentials"));
    }
    let (host, port) = split_host_port(authority)?;
    let target = if target.is_empty() {
        "/".to_string()
    } else if target.starts_with('?') {
        format!("/{target}")
    } else {
        target.to_string()
    };
    Ok(HttpUrl { host, port, authority: authority.to_string(), target })
}

fn split_host_port(authority: &str) -> io::Result<(String, u16)> {
    let (host, port) = if let Some(bracketed) = authority.strip_prefix('[') {
        let end = bracketed.find(']').ok_or_else(|| invalid("unclosed `[` in URL host"))?;
        let host = &bracketed[..end];
        match &bracketed[end + 1..] {
            "" => (host, None),
            after => match after.strip_prefix(':') {
                Some(port) => (host, Some(port)),
                None => return Err(invalid("malformed URL host")),
            },
        }
    } else {
        match authority.rsplit_once(':') {
            Some((host, port)) => (host, Some(port)),
            None => (authority, None),
        }
    };
    if host.is_empty() {
        return Err(invalid("URL has an empty host"));
    }
    let port = match port {
        Some(p) => p.parse::<u16>().map_err(|_| invalid(format!("invalid URL port `{p}`")))?,
        None => 80,
    };
    Ok((host.to_string(), port))
}

/// A remote object read with HTTP range requests.
#[derive(Debug)]
pub struct HttpSource {
    url: HttpUrl,
    options: HttpOptions,
    len: u64,
}

/// Outcome of a failed attempt. `Retry` failures may be retried; `Fatal` ones may not.
enum Failure {
    Retry(io::Error),
    Fatal(io::Error),
}

impl HttpSource {
    /// Opens an unsigned (public or presigned) object with default options.
    ///
    /// # Errors
    /// See [`HttpSource::open_with`].
    pub fn open(url: &str) -> io::Result<Self> {
        Self::open_with(url, HttpOptions::default())
    }

    /// Opens an object and learns its length with a `Range: bytes=0-0` request.
    ///
    /// # Errors
    /// Returns [`io::ErrorKind::Unsupported`] for `https://` URLs. Returns
    /// [`io::ErrorKind::InvalidInput`] for an invalid URL, a zero timeout, or a
    /// presigned URL combined with `options.s3`. Returns the last failure once the
    /// retry budget is exhausted, or a fatal error for a 4xx response or an
    /// unusable length response.
    pub fn open_with(url: &str, options: HttpOptions) -> io::Result<Self> {
        let parsed = parse_url(url)?;
        if options.timeout.is_zero() {
            return Err(invalid("timeout must be non-zero"));
        }
        if options.s3.is_some() && parsed.target.to_ascii_lowercase().contains("x-amz-signature=") {
            return Err(invalid("a presigned URL cannot also be signed with S3Signing"));
        }
        let mut source = Self { url: parsed, options, len: 0 };
        source.len = source.with_retry(|| source.fetch_length())?;
        Ok(source)
    }

    /// The parsed URL this source reads from.
    #[must_use]
    pub fn url(&self) -> &HttpUrl {
        &self.url
    }

    fn with_retry<T>(&self, mut op: impl FnMut() -> Result<T, Failure>) -> io::Result<T> {
        let mut retry = 0;
        loop {
            match op() {
                Ok(value) => return Ok(value),
                Err(Failure::Fatal(err)) => return Err(err),
                Err(Failure::Retry(err)) => {
                    if retry >= self.options.retry.max_retries {
                        return Err(err);
                    }
                    let wait = self.options.retry.backoff(retry);
                    if !wait.is_zero() {
                        thread::sleep(wait);
                    }
                    retry += 1;
                }
            }
        }
    }

    fn fetch_length(&self) -> Result<u64, Failure> {
        let (head, _body) = self.perform(0, 0)?;
        match head.status {
            206 => {
                check_framing(&head)?;
                let (start, _end, total) = header_content_range(&head)?;
                if start != 0 {
                    return Err(protocol("Content-Range does not start at 0"));
                }
                Ok(total)
            }
            // An empty object cannot satisfy any range; the server reports its size.
            416 => match head.get("content-range").and_then(parse_unsatisfied_total) {
                Some(0) => Ok(0),
                _ => Err(status_failure(&head)),
            },
            200 => {
                check_framing(&head)?;
                content_length(&head)?.ok_or_else(|| {
                    Failure::Fatal(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "server sent neither Content-Range nor Content-Length",
                    ))
                })
            }
            _ => Err(status_failure(&head)),
        }
    }

    fn fetch_into(&self, start: u64, end: u64, buf: &mut [u8]) -> Result<(), Failure> {
        let (head, mut body) = self.perform(start, end)?;
        match head.status {
            206 => {
                check_framing(&head)?;
                let (got_start, got_end, total) = header_content_range(&head)?;
                if got_start != start || got_end != end || total != self.len {
                    return Err(protocol("Content-Range does not match the request"));
                }
                if content_length(&head)?.is_some_and(|n| n != buf.len() as u64) {
                    return Err(protocol("Content-Length does not match the requested range"));
                }
                body.read_exact(buf).map_err(Failure::Retry)
            }
            200 => {
                // The server ignored Range. Skip the prefix and read the wanted bytes.
                check_framing(&head)?;
                if content_length(&head)?.is_some_and(|n| n != self.len) {
                    return Err(protocol("object length changed during the read"));
                }
                let skipped = io::copy(&mut (&mut body).take(start), &mut io::sink())
                    .map_err(Failure::Retry)?;
                if skipped != start {
                    return Err(Failure::Retry(truncated()));
                }
                body.read_exact(buf).map_err(Failure::Retry)
            }
            _ => Err(status_failure(&head)),
        }
    }

    /// Sends one request for the inclusive byte range `start..=end` and reads the head.
    fn perform(&self, start: u64, end: u64) -> Result<(Head, BufReader<TcpStream>), Failure> {
        let request = self.build_request(start, end).map_err(Failure::Fatal)?;
        let mut stream = connect(&self.url, self.options.timeout).map_err(Failure::Retry)?;
        stream.set_read_timeout(Some(self.options.timeout)).map_err(Failure::Retry)?;
        stream.set_write_timeout(Some(self.options.timeout)).map_err(Failure::Retry)?;
        stream.write_all(request.as_bytes()).map_err(Failure::Retry)?;
        let mut reader = BufReader::new(stream);
        let head = read_head(&mut reader)?;
        Ok((head, reader))
    }

    fn build_request(&self, start: u64, end: u64) -> io::Result<String> {
        let range = format!("bytes={start}-{end}");
        let mut lines = vec![
            format!("GET {} HTTP/1.1", self.url.target),
            format!("Host: {}", self.url.authority),
            format!("Range: {range}"),
            "User-Agent: tpt-bio-vision".to_string(),
            "Accept-Encoding: identity".to_string(),
            "Connection: close".to_string(),
        ];
        if let Some(s3) = &self.options.s3 {
            let payload_hash = sigv4::payload_hash_hex(b"");
            let signed_headers = [
                ("host", self.url.authority.as_str()),
                ("range", range.as_str()),
                ("x-amz-content-sha256", payload_hash.as_str()),
            ];
            let signed =
                sign_s3(s3, &amz_timestamp(SystemTime::now()), &self.url.target, &signed_headers)?;
            lines.push(format!("x-amz-date: {}", signed.amz_date));
            lines.push(format!("x-amz-content-sha256: {payload_hash}"));
            if let Some(token) = &signed.security_token {
                lines.push(format!("x-amz-security-token: {token}"));
            }
            lines.push(format!("Authorization: {}", signed.authorization));
        }
        Ok(format!("{}\r\n\r\n", lines.join("\r\n")))
    }
}

impl ByteSource for HttpSource {
    fn len(&self) -> u64 {
        self.len
    }

    fn read_exact_at(&self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
        source::check_range(offset, buf.len(), self.len)?;
        if buf.is_empty() {
            return Ok(());
        }
        let end = offset + buf.len() as u64 - 1;
        self.with_retry(|| self.fetch_into(offset, end, &mut *buf))
    }
}

/// Signs a request to S3 with `SigV4` for the `s3` service.
fn sign_s3(
    s3: &S3Signing,
    timestamp: &str,
    target: &str,
    headers: &[(&str, &str)],
) -> io::Result<SignedHeaders> {
    let creds = Credentials {
        access_key: &s3.access_key,
        secret_key: &s3.secret_key,
        session_token: s3.session_token.as_deref(),
    };
    sign_target(&creds, &s3.region, S3_SERVICE, timestamp, "GET", target, headers)
}

/// Signs a `GET` to `target` (path plus raw query) with the given scope.
///
/// The query is split into decoded name and value pairs, which the signer then
/// re-encodes in canonical form.
fn sign_target(
    creds: &Credentials<'_>,
    region: &str,
    service: &str,
    timestamp: &str,
    method: &str,
    target: &str,
    headers: &[(&str, &str)],
) -> io::Result<SignedHeaders> {
    let (path, raw_query) = target.split_once('?').unwrap_or((target, ""));
    let query = parse_query(raw_query)?;
    let query_refs: Vec<(&str, &str)> =
        query.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
    let payload_hash = sigv4::payload_hash_hex(b"");
    let req = Request { method, path, query: &query_refs, headers, payload_hash: &payload_hash };
    Ok(sigv4::sign(creds, &Scope { region, service }, timestamp, &req))
}

fn parse_query(query: &str) -> io::Result<Vec<(String, String)>> {
    query
        .split('&')
        .filter(|pair| !pair.is_empty())
        .map(|pair| {
            let (name, value) = pair.split_once('=').unwrap_or((pair, ""));
            Ok((percent_decode(name)?, percent_decode(value)?))
        })
        .collect()
}

fn percent_decode(s: &str) -> io::Result<String> {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' {
            let value = bytes
                .get(i + 1..i + 3)
                .and_then(|hex| std::str::from_utf8(hex).ok())
                .and_then(|hex| u8::from_str_radix(hex, 16).ok())
                .ok_or_else(|| invalid("malformed percent-escape in query"))?;
            out.push(value);
            i += 3;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8(out).map_err(|_| invalid("query is not valid UTF-8 after decoding"))
}

/// Formats a Unix time as an ISO 8601 basic timestamp, `YYYYMMDDTHHMMSSZ`.
fn amz_timestamp(now: SystemTime) -> String {
    let secs = now.duration_since(UNIX_EPOCH).map_or(0, |d| d.as_secs());
    format_amz_timestamp(secs)
}

fn format_amz_timestamp(secs: u64) -> String {
    // Civil date from days since the epoch (Howard Hinnant's algorithm).
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!("{year:04}{month:02}{day:02}T{:02}{:02}{:02}Z", rem / 3600, rem % 3600 / 60, rem % 60)
}

/// A parsed response status line and headers.
struct Head {
    status: u16,
    headers: Vec<(String, String)>,
}

impl Head {
    fn get(&self, name: &str) -> Option<&str> {
        self.headers.iter().find(|(n, _)| n == name).map(|(_, v)| v.as_str())
    }
}

fn connect(url: &HttpUrl, timeout: Duration) -> io::Result<TcpStream> {
    let mut last_err = None;
    for addr in (url.host.as_str(), url.port).to_socket_addrs()? {
        match TcpStream::connect_timeout(&addr, timeout) {
            Ok(stream) => return Ok(stream),
            Err(err) => last_err = Some(err),
        }
    }
    Err(last_err
        .unwrap_or_else(|| io::Error::new(io::ErrorKind::NotFound, "host resolved to no address")))
}

fn read_head(reader: &mut BufReader<TcpStream>) -> Result<Head, Failure> {
    let mut lines: Vec<String> = Vec::new();
    let mut total = 0usize;
    loop {
        let mut line = Vec::new();
        let n = reader.read_until(b'\n', &mut line).map_err(Failure::Retry)?;
        if n == 0 {
            return Err(Failure::Retry(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "connection closed before the response head",
            )));
        }
        total += n;
        if total > MAX_HEAD_BYTES {
            return Err(protocol("response head is too large"));
        }
        let text = String::from_utf8_lossy(&line);
        let text = text.trim_end_matches(['\r', '\n']);
        if text.is_empty() {
            break;
        }
        lines.push(text.to_string());
    }
    let status_line = lines.first().ok_or_else(|| protocol("missing status line"))?;
    let mut parts = status_line.splitn(3, ' ');
    let version = parts.next().unwrap_or_default();
    if !version.starts_with("HTTP/1.") {
        return Err(protocol("unsupported HTTP version"));
    }
    let status = parts
        .next()
        .and_then(|s| s.parse::<u16>().ok())
        .ok_or_else(|| protocol("malformed status code"))?;
    let headers = lines[1..]
        .iter()
        .filter_map(|line| line.split_once(':'))
        .map(|(name, value)| (name.trim().to_ascii_lowercase(), value.trim().to_string()))
        .collect();
    Ok(Head { status, headers })
}

/// Rejects any transfer coding. Only identity framing (`Content-Length` or close) is read.
fn check_framing(head: &Head) -> Result<(), Failure> {
    match head.get("transfer-encoding") {
        Some(coding) if !coding.eq_ignore_ascii_case("identity") => {
            Err(Failure::Fatal(io::Error::new(
                io::ErrorKind::Unsupported,
                format!(
                    "transfer-encoding `{coding}` is not supported; the server must send Content-Length"
                ),
            )))
        }
        _ => Ok(()),
    }
}

fn content_length(head: &Head) -> Result<Option<u64>, Failure> {
    head.get("content-length")
        .map(|v| v.parse::<u64>().map_err(|_| protocol("malformed Content-Length")))
        .transpose()
}

fn header_content_range(head: &Head) -> Result<(u64, u64, u64), Failure> {
    head.get("content-range")
        .and_then(parse_content_range)
        .ok_or_else(|| protocol("missing or malformed Content-Range"))
}

/// Parses `bytes START-END/TOTAL`.
fn parse_content_range(value: &str) -> Option<(u64, u64, u64)> {
    let spec = value.strip_prefix("bytes ")?;
    let (range, total) = spec.split_once('/')?;
    let (start, end) = range.split_once('-')?;
    Some((start.parse().ok()?, end.parse().ok()?, total.parse().ok()?))
}

/// Parses `bytes */TOTAL`, the form a `416` response uses.
fn parse_unsatisfied_total(value: &str) -> Option<u64> {
    value.strip_prefix("bytes */")?.parse().ok()
}

fn status_failure(head: &Head) -> Failure {
    let code = head.status;
    let message = format!("HTTP {code} from server");
    if (500..600).contains(&code) {
        return Failure::Retry(io::Error::other(message));
    }
    let kind = match code {
        403 => io::ErrorKind::PermissionDenied,
        404 => io::ErrorKind::NotFound,
        416 => io::ErrorKind::UnexpectedEof,
        _ => io::ErrorKind::Other,
    };
    Failure::Fatal(io::Error::new(kind, message))
}

fn protocol(message: &str) -> Failure {
    Failure::Fatal(io::Error::new(io::ErrorKind::InvalidData, message.to_string()))
}

fn truncated() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "response body was truncated")
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use std::io::{BufRead, BufReader, Read, Write};
    use std::net::{SocketAddr, TcpListener, TcpStream};
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{Arc, Mutex};
    use std::thread;
    use std::time::Duration;

    use super::*;

    /// A fault injected into the next request.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Fault {
        /// Respond with this status and a small plain body.
        Status(u16),
        /// Close the connection without writing a response.
        Drop,
        /// Send the normal response head, then only this many body bytes.
        Truncate(usize),
        /// Send a chunked body.
        Chunked,
    }

    #[derive(Debug, Clone)]
    struct Recorded {
        target: String,
        headers: Vec<(String, String)>,
    }

    impl Recorded {
        fn header(&self, name: &str) -> Option<&str> {
            self.headers.iter().find(|(n, _)| n == name).map(|(_, v)| v.as_str())
        }
    }

    #[derive(Default)]
    struct Script {
        faults: VecDeque<Fault>,
        ignore_range: bool,
        requests: Vec<Recorded>,
    }

    /// A single-object HTTP server on 127.0.0.1 with an ephemeral port.
    struct TestServer {
        addr: SocketAddr,
        script: Arc<Mutex<Script>>,
        stop: Arc<AtomicBool>,
    }

    impl TestServer {
        fn start(body: Vec<u8>) -> Self {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let addr = listener.local_addr().unwrap();
            let script = Arc::new(Mutex::new(Script::default()));
            let stop = Arc::new(AtomicBool::new(false));
            let (thread_script, thread_stop) = (Arc::clone(&script), Arc::clone(&stop));
            thread::spawn(move || {
                for stream in listener.incoming() {
                    if thread_stop.load(Ordering::SeqCst) {
                        break;
                    }
                    if let Ok(stream) = stream {
                        handle(stream, &body, &thread_script);
                    }
                }
            });
            Self { addr, script, stop }
        }

        fn url(&self, target: &str) -> String {
            format!("http://{}{target}", self.addr)
        }

        fn push_fault(&self, fault: Fault) {
            self.script.lock().unwrap().faults.push_back(fault);
        }

        fn set_ignore_range(&self, ignore: bool) {
            self.script.lock().unwrap().ignore_range = ignore;
        }

        fn requests(&self) -> Vec<Recorded> {
            self.script.lock().unwrap().requests.clone()
        }
    }

    impl Drop for TestServer {
        fn drop(&mut self) {
            self.stop.store(true, Ordering::SeqCst);
            // Wake the blocked accept so the thread can exit.
            let _ = TcpStream::connect(self.addr);
        }
    }

    /// Parses `bytes=START-END` against a body of `len` bytes. `None` means unsatisfiable.
    fn parse_range_request(value: &str, len: u64) -> Option<(u64, u64)> {
        let spec = value.strip_prefix("bytes=")?;
        let (a, b) = spec.split_once('-')?;
        let start: u64 = a.parse().ok()?;
        let end: u64 = if b.is_empty() { len.checked_sub(1)? } else { b.parse().ok()? };
        if start >= len || end < start {
            return None;
        }
        Some((start, end.min(len - 1)))
    }

    fn handle(mut stream: TcpStream, body: &[u8], script: &Mutex<Script>) {
        let mut reader = BufReader::new(stream.try_clone().unwrap());
        let mut lines = Vec::new();
        loop {
            let mut line = String::new();
            if reader.read_line(&mut line).unwrap_or(0) == 0 {
                return;
            }
            let trimmed = line.trim_end().to_string();
            if trimmed.is_empty() {
                break;
            }
            lines.push(trimmed);
        }
        let target = lines[0].split(' ').nth(1).unwrap_or("/").to_string();
        let headers: Vec<(String, String)> = lines[1..]
            .iter()
            .filter_map(|l| l.split_once(':'))
            .map(|(n, v)| (n.trim().to_ascii_lowercase(), v.trim().to_string()))
            .collect();
        let (fault, ignore_range) = {
            let mut s = script.lock().unwrap();
            s.requests.push(Recorded { target, headers: headers.clone() });
            (s.faults.pop_front(), s.ignore_range)
        };
        let range = headers.iter().find(|(n, _)| n == "range").map(|(_, v)| v.clone());

        if fault == Some(Fault::Drop) {
            return;
        }
        if let Some(Fault::Status(code)) = fault {
            let _ = write_response(&mut stream, &format!("{code} Error"), "", b"fail", None);
            return;
        }
        if fault == Some(Fault::Chunked) {
            let _ = stream.write_all(
                b"HTTP/1.1 200 OK\r\nTransfer-Encoding: chunked\r\nConnection: close\r\n\r\n5\r\nhello\r\n0\r\n\r\n",
            );
            return;
        }

        let limit = match fault {
            Some(Fault::Truncate(n)) => Some(n),
            _ => None,
        };
        let len = body.len() as u64;
        match range.filter(|_| !ignore_range).map(|value| parse_range_request(&value, len)) {
            Some(Some((start, end))) => {
                let part = &body[start as usize..=end as usize];
                let extra = format!("Content-Range: bytes {start}-{end}/{len}\r\n");
                let _ = write_response(&mut stream, "206 Partial Content", &extra, part, limit);
            }
            Some(None) => {
                let extra = format!("Content-Range: bytes */{len}\r\n");
                let _ = write_response(&mut stream, "416 Range Not Satisfiable", &extra, b"", None);
            }
            None => {
                let _ = write_response(&mut stream, "200 OK", "", body, limit);
            }
        }
    }

    fn write_response(
        stream: &mut TcpStream,
        status: &str,
        extra_headers: &str,
        body: &[u8],
        limit: Option<usize>,
    ) -> io::Result<()> {
        let head = format!(
            "HTTP/1.1 {status}\r\n{extra_headers}Content-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        );
        stream.write_all(head.as_bytes())?;
        let cut = limit.map_or(body.len(), |n| n.min(body.len()));
        stream.write_all(&body[..cut])
    }

    fn fast_retry(max_retries: u32) -> HttpOptions {
        HttpOptions {
            retry: RetryPolicy {
                max_retries,
                initial_backoff: Duration::ZERO,
                max_backoff: Duration::ZERO,
            },
            timeout: Duration::from_secs(5),
            s3: None,
        }
    }

    fn sample_body() -> Vec<u8> {
        (0u8..=99).collect()
    }

    #[test]
    fn https_is_rejected_without_network() {
        let err = HttpSource::open("https://example.invalid/object").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Unsupported);
        assert_eq!(
            err.to_string(),
            "https is not supported: no audited pure-Rust TLS stack is available"
        );
        let err = parse_url("HTTPS://example.invalid/").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Unsupported);
    }

    #[test]
    fn url_parsing_covers_ports_queries_and_ipv6() {
        let u = parse_url("http://example.com/a/b?x=1#frag").unwrap();
        assert_eq!(u.host, "example.com");
        assert_eq!(u.port, 80);
        assert_eq!(u.target, "/a/b?x=1");
        let u = parse_url("http://127.0.0.1:9000").unwrap();
        assert_eq!((u.port, u.target.as_str()), (9000, "/"));
        let u = parse_url("http://[::1]:8080/k").unwrap();
        assert_eq!((u.host.as_str(), u.port, u.authority.as_str()), ("::1", 8080, "[::1]:8080"));
        assert!(parse_url("ftp://x/").is_err());
        assert!(parse_url("http://user:pw@x/").is_err());
        assert!(parse_url("http://x:notaport/").is_err());
        assert!(parse_url("http:///path").is_err());
    }

    #[test]
    fn reads_byte_ranges_and_learns_length() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj.bin"), fast_retry(0)).unwrap();
        assert_eq!(src.len(), 100);
        assert_eq!(src.read_vec_at(10, 5).unwrap(), [10, 11, 12, 13, 14]);
        assert_eq!(src.read_vec_at(95, 5).unwrap(), [95, 96, 97, 98, 99]);
        assert!(src.read_vec_at(96, 5).is_err());
        let reqs = server.requests();
        assert_eq!(reqs[0].header("range"), Some("bytes=0-0"));
        assert_eq!(reqs[1].header("range"), Some("bytes=10-14"));
        assert_eq!(reqs[1].header("connection"), Some("close"));
        assert_eq!(reqs[1].header("authorization"), None);
    }

    #[test]
    fn empty_object_has_zero_length() {
        let server = TestServer::start(Vec::new());
        let src = HttpSource::open_with(&server.url("/empty"), fast_retry(0)).unwrap();
        assert_eq!(src.len(), 0);
        assert!(src.is_empty());
    }

    #[test]
    fn server_ignoring_range_is_handled() {
        let server = TestServer::start(sample_body());
        server.set_ignore_range(true);
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(0)).unwrap();
        assert_eq!(src.len(), 100);
        assert_eq!(src.read_vec_at(20, 3).unwrap(), [20, 21, 22]);
    }

    #[test]
    fn chunked_responses_are_rejected() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(3)).unwrap();
        server.push_fault(Fault::Chunked);
        let before = server.requests().len();
        let err = src.read_vec_at(0, 4).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Unsupported);
        assert!(err.to_string().contains("chunked"));
        assert_eq!(server.requests().len() - before, 1, "chunked must not be retried");
    }

    #[test]
    fn retries_on_server_error() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(2)).unwrap();
        server.push_fault(Fault::Status(500));
        let before = server.requests().len();
        assert_eq!(src.read_vec_at(1, 2).unwrap(), [1, 2]);
        assert_eq!(server.requests().len() - before, 2);
    }

    #[test]
    fn retries_on_dropped_connection() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(2)).unwrap();
        server.push_fault(Fault::Drop);
        let before = server.requests().len();
        assert_eq!(src.read_vec_at(7, 3).unwrap(), [7, 8, 9]);
        assert_eq!(server.requests().len() - before, 2);
    }

    #[test]
    fn retries_on_truncated_body() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(2)).unwrap();
        server.push_fault(Fault::Truncate(2));
        let before = server.requests().len();
        assert_eq!(src.read_vec_at(30, 6).unwrap(), [30, 31, 32, 33, 34, 35]);
        assert_eq!(server.requests().len() - before, 2);
    }

    #[test]
    fn length_request_is_retried() {
        let server = TestServer::start(sample_body());
        server.push_fault(Fault::Status(503));
        server.push_fault(Fault::Drop);
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(3)).unwrap();
        assert_eq!(src.len(), 100);
        assert_eq!(server.requests().len(), 3);
    }

    #[test]
    fn gives_up_after_retry_budget() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(3)).unwrap();
        for _ in 0..4 {
            server.push_fault(Fault::Status(500));
        }
        let before = server.requests().len();
        let err = src.read_vec_at(0, 1).unwrap_err();
        assert!(err.to_string().contains("500"));
        assert_eq!(server.requests().len() - before, 4, "1 attempt plus 3 retries");
    }

    #[test]
    fn client_errors_are_not_retried() {
        let server = TestServer::start(sample_body());
        let src = HttpSource::open_with(&server.url("/obj"), fast_retry(3)).unwrap();
        server.push_fault(Fault::Status(404));
        let before = server.requests().len();
        let err = src.read_vec_at(0, 1).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        assert_eq!(server.requests().len() - before, 1);
    }

    #[test]
    fn backoff_doubles_and_is_capped() {
        let policy = RetryPolicy {
            max_retries: 5,
            initial_backoff: Duration::from_millis(10),
            max_backoff: Duration::from_millis(35),
        };
        assert_eq!(policy.backoff(0), Duration::from_millis(10));
        assert_eq!(policy.backoff(1), Duration::from_millis(20));
        assert_eq!(policy.backoff(2), Duration::from_millis(35));
        assert_eq!(policy.backoff(200), Duration::from_millis(35));
    }

    #[test]
    fn presigned_query_is_sent_unchanged() {
        let server = TestServer::start(sample_body());
        let query = "?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKID%2F20240101%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240101T000000Z&X-Amz-Expires=300&X-Amz-SignedHeaders=host&X-Amz-Signature=ab%2Fcd";
        let url = server.url(&format!("/bucket/key.bin{query}"));
        let src = HttpSource::open_with(&url, fast_retry(0)).unwrap();
        assert_eq!(src.read_vec_at(0, 2).unwrap(), [0, 1]);
        for req in server.requests() {
            assert_eq!(req.target, format!("/bucket/key.bin{query}"));
            assert_eq!(req.header("authorization"), None);
        }
    }

    #[test]
    fn presigned_url_with_s3_signing_is_refused() {
        let url = "http://127.0.0.1:1/k?X-Amz-Signature=00";
        let options = HttpOptions {
            s3: Some(S3Signing {
                access_key: "AKID".into(),
                secret_key: "secret".into(),
                session_token: None,
                region: "us-east-1".into(),
            }),
            ..HttpOptions::default()
        };
        let err = HttpSource::open_with(url, options).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn s3_signing_matches_recomputed_signature() {
        let server = TestServer::start(sample_body());
        let options = HttpOptions {
            retry: RetryPolicy::default(),
            timeout: Duration::from_secs(5),
            s3: Some(S3Signing {
                access_key: "AKIDEXAMPLE".into(),
                secret_key: "wJalrXUtnFEMI/K7MDENG+bPxRfiCYEXAMPLEKEY".into(),
                session_token: Some("SESSION".into()),
                region: "us-east-1".into(),
            }),
        };
        // Path-style URL: the bucket is part of the path.
        let src = HttpSource::open_with(&server.url("/bucket/dir/obj.bin"), options).unwrap();
        assert_eq!(src.read_vec_at(4, 2).unwrap(), [4, 5]);

        let req = server.requests().last().unwrap().clone();
        let date = req.header("x-amz-date").unwrap().to_string();
        assert!(date.len() == 16 && date.ends_with('Z'));
        assert_eq!(req.header("x-amz-security-token"), Some("SESSION"));
        let host = format!("127.0.0.1:{}", server.addr.port());
        let payload = sigv4::payload_hash_hex(b"");
        let creds = Credentials {
            access_key: "AKIDEXAMPLE",
            secret_key: "wJalrXUtnFEMI/K7MDENG+bPxRfiCYEXAMPLEKEY",
            session_token: Some("SESSION"),
        };
        let expected = sigv4::sign(
            &creds,
            &Scope { region: "us-east-1", service: "s3" },
            &date,
            &Request {
                method: "GET",
                path: "/bucket/dir/obj.bin",
                query: &[],
                headers: &[
                    ("host", host.as_str()),
                    ("range", "bytes=4-5"),
                    ("x-amz-content-sha256", payload.as_str()),
                ],
                payload_hash: &payload,
            },
        );
        assert_eq!(req.header("authorization"), Some(expected.authorization.as_str()));
        assert!(expected.authorization.contains(
            "SignedHeaders=host;range;x-amz-content-sha256;x-amz-date;x-amz-security-token"
        ));
    }

    #[test]
    fn signing_glue_reproduces_aws_documented_example() {
        // The IAM ListUsers example from the AWS SigV4 documentation, as in sigv4.rs.
        let creds = Credentials {
            access_key: "AKIDEXAMPLE",
            secret_key: "wJalrXUtnFEMI/K7MDENG+bPxRfiCYEXAMPLEKEY",
            session_token: None,
        };
        let signed = sign_target(
            &creds,
            "us-east-1",
            "iam",
            "20150830T123600Z",
            "GET",
            "/?Action=ListUsers&Version=2010-05-08",
            &[
                ("Content-Type", "application/x-www-form-urlencoded; charset=utf-8"),
                ("Host", "iam.amazonaws.com"),
            ],
        )
        .unwrap();
        assert_eq!(
            signed.authorization,
            "AWS4-HMAC-SHA256 Credential=AKIDEXAMPLE/20150830/us-east-1/iam/aws4_request, \
             SignedHeaders=content-type;host;x-amz-date, \
             Signature=5d672d79c15b13162d9279b0855cfba6789a8edb4c82c400e06b5924a6f2b5d7"
        );
    }

    #[test]
    fn timestamps_format_as_amz_basic_iso8601() {
        assert_eq!(format_amz_timestamp(0), "19700101T000000Z");
        assert_eq!(format_amz_timestamp(1_440_938_160), "20150830T123600Z");
        assert_eq!(format_amz_timestamp(1_709_164_800), "20240229T000000Z");
    }

    #[test]
    fn percent_decoding_and_query_parsing() {
        assert_eq!(percent_decode("a%20b%2Fc").unwrap(), "a b/c");
        assert!(percent_decode("bad%2").is_err());
        assert!(percent_decode("bad%zz").is_err());
        let pairs = parse_query("Action=ListUsers&flag").unwrap();
        assert_eq!(pairs[0], ("Action".to_string(), "ListUsers".to_string()));
        assert_eq!(pairs[1], ("flag".to_string(), String::new()));
    }

    #[test]
    fn test_server_rejects_bad_ranges_with_416() {
        let server = TestServer::start(sample_body());
        let mut stream = TcpStream::connect(server.addr).unwrap();
        stream
            .write_all(b"GET /x HTTP/1.1\r\nRange: bytes=200-300\r\nConnection: close\r\n\r\n")
            .unwrap();
        let mut response = String::new();
        stream.read_to_string(&mut response).unwrap();
        assert!(response.starts_with("HTTP/1.1 416"));
        assert!(response.contains("Content-Range: bytes */100"));
    }
}
