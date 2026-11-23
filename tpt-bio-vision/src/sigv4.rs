//! AWS Signature Version 4 request signing.
//!
//! SHA-256 and HMAC come from `tpt-crypto-hash`, so no other crypto library is
//! linked. This module only builds the canonical request, the string to sign, and
//! the derived signing key. It does not send requests.

use std::fmt::Write as _;

use tpt_crypto_hash::mac::hmac_sha256;
use tpt_crypto_hash::sha2::sha256;

/// The algorithm identifier used in the `Authorization` header.
pub const ALGORITHM: &str = "AWS4-HMAC-SHA256";

/// Payload hash to send when the body is not signed (for example, streamed uploads).
pub const UNSIGNED_PAYLOAD: &str = "UNSIGNED-PAYLOAD";

/// Access key pair, with an optional session token for temporary credentials.
#[derive(Debug, Clone, Copy)]
pub struct Credentials<'a> {
    /// Access key ID.
    pub access_key: &'a str,
    /// Secret access key.
    pub secret_key: &'a str,
    /// Session token for temporary credentials, sent as `x-amz-security-token`.
    pub session_token: Option<&'a str>,
}

/// A request to sign. Names in `headers` must be lowercase and must include `host`.
#[derive(Debug, Clone, Copy)]
pub struct Request<'a> {
    /// HTTP method, for example `GET`.
    pub method: &'a str,
    /// Canonical URI path, already percent-encoded. Use `/` for an empty path.
    pub path: &'a str,
    /// Query parameters as raw (unencoded) name and value pairs.
    pub query: &'a [(&'a str, &'a str)],
    /// Headers to sign, as raw name and value pairs.
    pub headers: &'a [(&'a str, &'a str)],
    /// Lowercase hex SHA-256 of the body, or [`UNSIGNED_PAYLOAD`].
    pub payload_hash: &'a str,
}

/// The region and service the signature is scoped to.
#[derive(Debug, Clone, Copy)]
pub struct Scope<'a> {
    /// Region, for example `us-east-1`.
    pub region: &'a str,
    /// Service, for example `s3`.
    pub service: &'a str,
}

/// Headers to add to the request so that it is authorized.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SignedHeaders {
    /// Value of the `Authorization` header.
    pub authorization: String,
    /// Value of the `x-amz-date` header (the request timestamp).
    pub amz_date: String,
    /// Value of `x-amz-security-token`, when a session token is used.
    pub security_token: Option<String>,
}

/// Percent-encodes `s` as `SigV4` requires: unreserved characters (`A-Z a-z 0-9 - _ . ~`)
/// are kept, everything else becomes `%XX` in uppercase. `/` is kept only when
/// `encode_slash` is false, which is the rule for URI paths.
#[must_use]
pub fn uri_encode(s: &str, encode_slash: bool) -> String {
    let mut out = String::with_capacity(s.len());
    for &b in s.as_bytes() {
        let unreserved = b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~');
        if unreserved || (!encode_slash && b == b'/') {
            out.push(char::from(b));
        } else {
            // Writing to a String cannot fail.
            let _ = write!(out, "%{b:02X}");
        }
    }
    out
}

/// Lowercase hexadecimal encoding.
#[must_use]
pub fn hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        out.push(char::from(DIGITS[usize::from(b >> 4)]));
        out.push(char::from(DIGITS[usize::from(b & 0xf)]));
    }
    out
}

/// Lowercase hex SHA-256 of `data`. Use it for the payload hash of a request body.
#[must_use]
pub fn payload_hash_hex(data: &[u8]) -> String {
    hex(&sha256(data))
}

/// Signs a request.
///
/// `timestamp` is the request time in `YYYYMMDDTHHMMSSZ` form. The date is its first
/// eight characters. The signer adds `x-amz-date` (and the session token, if any)
/// to the signed headers itself.
///
/// # Panics
/// Panics if `timestamp` is not at least eight characters long.
#[must_use]
pub fn sign(
    creds: &Credentials<'_>,
    scope: &Scope<'_>,
    timestamp: &str,
    req: &Request<'_>,
) -> SignedHeaders {
    assert!(timestamp.len() >= 8, "timestamp must be YYYYMMDDTHHMMSSZ");
    let date = &timestamp[..8];

    // Headers to sign: the caller's headers plus x-amz-date and the session token.
    let mut headers: Vec<(String, String)> =
        req.headers.iter().map(|(n, v)| (n.to_ascii_lowercase(), normalize_value(v))).collect();
    headers.push(("x-amz-date".to_string(), timestamp.to_string()));
    if let Some(token) = creds.session_token {
        headers.push(("x-amz-security-token".to_string(), token.to_string()));
    }
    headers.sort_by(|a, b| a.0.cmp(&b.0));

    let canonical_headers: String = headers.iter().fold(String::new(), |mut acc, (n, v)| {
        // Writing to a String cannot fail.
        let _ = writeln!(acc, "{n}:{v}");
        acc
    });
    let signed_headers: Vec<&str> = headers.iter().map(|(n, _)| n.as_str()).collect();
    let signed_headers = signed_headers.join(";");

    let canonical_request = format!(
        "{}\n{}\n{}\n{}\n{}\n{}",
        req.method,
        req.path,
        canonical_query(req.query),
        canonical_headers,
        signed_headers,
        req.payload_hash,
    );

    let credential_scope = format!("{date}/{}/{}/aws4_request", scope.region, scope.service);
    let string_to_sign = format!(
        "{ALGORITHM}\n{timestamp}\n{credential_scope}\n{}",
        hex(&sha256(canonical_request.as_bytes())),
    );

    let signing_key = derive_signing_key(creds.secret_key, date, scope);
    let signature = hex(&hmac_sha256(&signing_key, string_to_sign.as_bytes()));

    SignedHeaders {
        authorization: format!(
            "{ALGORITHM} Credential={}/{credential_scope}, SignedHeaders={signed_headers}, Signature={signature}",
            creds.access_key,
        ),
        amz_date: timestamp.to_string(),
        security_token: creds.session_token.map(str::to_string),
    }
}

/// The `SigV4` signing key: an HMAC chain over date, region, service and the terminator.
#[must_use]
pub fn derive_signing_key(secret_key: &str, date: &str, scope: &Scope<'_>) -> [u8; 32] {
    let k_secret = format!("AWS4{secret_key}");
    let k_date = hmac_sha256(k_secret.as_bytes(), date.as_bytes());
    let k_region = hmac_sha256(&k_date, scope.region.as_bytes());
    let k_service = hmac_sha256(&k_region, scope.service.as_bytes());
    hmac_sha256(&k_service, b"aws4_request")
}

/// Trims a header value and collapses internal runs of spaces to one, as `SigV4` requires.
fn normalize_value(v: &str) -> String {
    v.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Sorted, encoded query string.
fn canonical_query(query: &[(&str, &str)]) -> String {
    let mut pairs: Vec<(String, String)> =
        query.iter().map(|(k, v)| (uri_encode(k, true), uri_encode(v, true))).collect();
    pairs.sort();
    pairs.iter().map(|(k, v)| format!("{k}={v}")).collect::<Vec<_>>().join("&")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// AWS documentation example: IAM `ListUsers`, signed at `2015-08-30T12:36:00Z`.
    /// The expected signature is the one published in the AWS `SigV4` examples.
    #[test]
    fn matches_aws_documented_example() {
        let creds = Credentials {
            access_key: "AKIDEXAMPLE",
            secret_key: "wJalrXUtnFEMI/K7MDENG+bPxRfiCYEXAMPLEKEY",
            session_token: None,
        };
        let scope = Scope { region: "us-east-1", service: "iam" };
        let req = Request {
            method: "GET",
            path: "/",
            query: &[("Action", "ListUsers"), ("Version", "2010-05-08")],
            headers: &[
                ("Content-Type", "application/x-www-form-urlencoded; charset=utf-8"),
                ("Host", "iam.amazonaws.com"),
            ],
            payload_hash: &payload_hash_hex(b""),
        };
        let signed = sign(&creds, &scope, "20150830T123600Z", &req);
        assert_eq!(
            signed.authorization,
            "AWS4-HMAC-SHA256 Credential=AKIDEXAMPLE/20150830/us-east-1/iam/aws4_request, \
             SignedHeaders=content-type;host;x-amz-date, \
             Signature=5d672d79c15b13162d9279b0855cfba6789a8edb4c82c400e06b5924a6f2b5d7"
        );
    }

    #[test]
    fn session_token_is_signed_and_returned() {
        let creds =
            Credentials { access_key: "AKID", secret_key: "secret", session_token: Some("TOKEN") };
        let scope = Scope { region: "us-east-1", service: "s3" };
        let req = Request {
            method: "GET",
            path: "/bucket/key",
            query: &[],
            headers: &[("host", "example.s3.amazonaws.com")],
            payload_hash: UNSIGNED_PAYLOAD,
        };
        let signed = sign(&creds, &scope, "20240101T000000Z", &req);
        assert!(signed.authorization.contains("x-amz-security-token"));
        assert_eq!(signed.security_token.as_deref(), Some("TOKEN"));
    }

    #[test]
    fn uri_encoding_rules() {
        assert_eq!(uri_encode("a b/c~d", true), "a%20b%2Fc~d");
        assert_eq!(uri_encode("a b/c", false), "a%20b/c");
        assert_eq!(hex(&[0x00, 0xab, 0xff]), "00abff");
    }
}
