//! End-to-end tests that run the built `tptbio` binary.
//!
//! Uses `std::process::Command` and Cargo's `CARGO_BIN_EXE_tptbio`, so no test
//! dependency is needed.

use std::fs::{self, File};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use tpt_bio_seq::bcf::build_bcf_index;
use tpt_bio_seq::{
    BcfHeader, BcfWriter, BgzfWriter, IndexKind, TabixConfig, TabixIndex, VcfRecord,
};

const EXIT_OK: i32 = 0;
const EXIT_RUNTIME: i32 = 1;
const EXIT_USAGE: i32 = 2;

/// A unique scratch directory under the system temp dir, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Self {
        let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        let dir =
            std::env::temp_dir().join(format!("tptbio-cli-{tag}-{}-{nanos}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn write(&self, name: &str, contents: &str) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, contents).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn tptbio(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_tptbio")).args(args).output().unwrap()
}

fn code(out: &Output) -> i32 {
    out.status.code().expect("process terminated by signal")
}

fn stdout(out: &Output) -> String {
    String::from_utf8(out.stdout.clone()).unwrap()
}

fn stderr(out: &Output) -> String {
    String::from_utf8(out.stderr.clone()).unwrap()
}

/// Parses CSV output into its data rows as `f64` (skipping the header line).
fn numeric_rows(text: &str) -> Vec<Vec<f64>> {
    text.lines().skip(1).map(|l| l.split(',').map(|f| f.parse().unwrap()).collect()).collect()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn help_exits_zero() {
    assert_eq!(code(&tptbio(&["--help"])), EXIT_OK);
    assert_eq!(code(&tptbio(&["stats", "--help"])), EXIT_OK);
}

#[test]
fn usage_errors_exit_two() {
    assert_eq!(code(&tptbio(&["no-such-command"])), EXIT_USAGE);
    assert_eq!(code(&tptbio(&["-q", "-v", "stats", "pca", "x.csv"])), EXIT_USAGE);
    assert_eq!(code(&tptbio(&["stats", "normalize", "x.csv", "--method", "bogus"])), EXIT_USAGE);
}

#[test]
fn image_commands_without_required_arguments_exit_two() {
    for args in
        [&["vision", "segment"][..], &["vision", "features"], &["pipeline", "run"], &["view"]]
    {
        let out = tptbio(args);
        assert_eq!(code(&out), EXIT_USAGE, "{args:?}");
        assert!(stderr(&out).contains("required"), "{args:?}: {}", stderr(&out));
    }
}

#[test]
fn normalize_zscore_on_known_column() {
    let s = Scratch::new("norm-z");
    let input = s.write(
        "in.csv",
        "x
1
3
5
7
9
",
    );
    let out = tptbio(&["stats", "normalize", input.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let text = stdout(&out);
    assert!(
        text.starts_with(
            "x
"
        ),
        "{text}"
    );
    // Population variance of 1..9 in steps of 2 is 8, so the first value is -4 / sqrt(8).
    let rows = numeric_rows(&text);
    assert_eq!(rows.len(), 5);
    assert!(close(rows[0][0], -4.0 / f64::sqrt(8.0)), "{rows:?}");
    assert!(close(rows[2][0], 0.0));
}

#[test]
fn normalize_per_group_removes_plate_offset() {
    let s = Scratch::new("norm-group");
    let input = s.write("in.csv", "plate,x\nA,1\nA,3\nB,101\nB,103\n");
    let out = tptbio(&["stats", "normalize", input.to_str().unwrap(), "--group", "plate"]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let text = stdout(&out);
    assert_eq!(text, "plate,x\nA,-1\nA,1\nB,-1\nB,1\n");
}

#[test]
fn normalize_robust_writes_expected_values() {
    let s = Scratch::new("norm-robust");
    let input = s.write("in.csv", "x\n1\n2\n3\n4\n100\n");
    let out = tptbio(&["stats", "normalize", input.to_str().unwrap(), "-m", "robust"]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let rows = numeric_rows(&stdout(&out));
    // Median 3, MAD 1, so the scale is 1.4826.
    assert!(close(rows[2][0], 0.0));
    assert!(close(rows[0][0], -2.0 / 1.4826), "{:?}", rows[0]);
}

#[test]
fn normalize_unknown_group_column_is_runtime_error() {
    let s = Scratch::new("norm-missing-col");
    let input = s.write("in.csv", "plate,x\nA,1\nA,2\n");
    let out = tptbio(&["stats", "normalize", input.to_str().unwrap(), "--group", "well"]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("no column named"), "{}", stderr(&out));
}

#[test]
fn normalize_non_numeric_field_is_runtime_error() {
    let s = Scratch::new("norm-text");
    let input = s.write("in.csv", "x\n1\nabc\n");
    let out = tptbio(&["stats", "normalize", input.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("is not a number"), "{}", stderr(&out));
}

#[test]
fn missing_input_file_is_runtime_error() {
    let s = Scratch::new("missing");
    let missing = s.0.join("does-not-exist.csv");
    let out = tptbio(&["stats", "pca", missing.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_RUNTIME);
}

#[test]
fn pca_recovers_diagonal_direction() {
    let s = Scratch::new("pca");
    // Points on the line y = x, so the first component is (1,1)/sqrt2 and carries all variance.
    let input = s.write("in.csv", "x,y\n1,1\n2,2\n-1,-1\n-2,-2\n");
    let out = tptbio(&["stats", "pca", input.to_str().unwrap(), "-k", "1"]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let text = stdout(&out);
    assert!(text.starts_with("PC1\n"), "{text}");
    let rows = numeric_rows(&text);
    assert_eq!(rows.len(), 4);
    // Score of (1,1) is sqrt2 along the unit axis.
    assert!(close(rows[0][0], f64::sqrt(2.0)), "{rows:?}");
    assert!(close(rows[1][0], 2.0 * f64::sqrt(2.0)), "{rows:?}");
}

#[test]
fn pca_rejects_too_many_components() {
    let s = Scratch::new("pca-k");
    let input = s.write("in.csv", "x,y\n1,2\n3,5\n");
    let out = tptbio(&["stats", "pca", input.to_str().unwrap(), "-k", "3"]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("n_components"), "{}", stderr(&out));
}

#[test]
fn quiet_suppresses_log_lines_but_not_errors() {
    let s = Scratch::new("quiet");
    let input = s.write("in.csv", "x\n1\n3\n");
    let out = tptbio(&["-q", "stats", "normalize", input.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK);
    assert!(stderr(&out).is_empty(), "{}", stderr(&out));
    let out = tptbio(&["stats", "normalize", input.to_str().unwrap()]);
    assert!(stderr(&out).contains("normalized"), "{}", stderr(&out));
}

#[test]
fn seq_view_header_only_prints_sam_header() {
    let s = Scratch::new("sam");
    let sam: &Path = &s.write("in.sam", "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n");
    let out = tptbio(&["seq", "view", "-H", sam.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).contains("@SQ\tSN:chr1\tLN:100"), "{}", stdout(&out));
}

#[test]
fn seq_sort_replaces_sort_order_instead_of_duplicating_it() {
    let s = Scratch::new("sort-so");
    let sam: &Path = &s.write(
        "in.sam",
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\nr2\t0\tchr1\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let bam = s.0.join("out.bam");
    let out = tptbio(&["seq", "sort", sam.to_str().unwrap(), "-o", bam.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let out = tptbio(&["seq", "view", "-H", bam.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let hd = stdout(&out).lines().find(|l| l.starts_with("@HD")).unwrap().to_owned();
    assert_eq!(hd.matches("SO:").count(), 1, "{hd}");
    assert!(hd.ends_with("SO:coordinate"), "{hd}");
}

const VCF_TEXT: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##FILTER=<ID=q10,Description=\"Quality below 10\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\trs1\tA\tG\t50\tPASS\tDP=10\nchr1\t200\t.\tC\tT\t.\tPASS\tDP=20\nchr1\t300\trs3\tG\tA,C\t60\tq10\tDP=30\n";

const SAM_TEXT: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

const CSV_HEADER: &str = "CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO\n";

/// Writes `VCF_TEXT` as a bgzipped VCF with a tabix index beside it.
fn write_vcf_gz(s: &Scratch, name: &str) -> PathBuf {
    let path = s.0.join(name);
    let mut w = BgzfWriter::new(File::create(&path).unwrap());
    w.write_all(VCF_TEXT.as_bytes()).unwrap();
    w.finish().unwrap();
    let index =
        TabixIndex::build_text(File::open(&path).unwrap(), TabixConfig::vcf(), IndexKind::Tbi)
            .unwrap();
    index.write_to(File::create(s.0.join(format!("{name}.tbi"))).unwrap()).unwrap();
    path
}

/// Writes the same records as a BCF with a tabix-style index beside it.
fn write_bcf(s: &Scratch, name: &str) -> PathBuf {
    let path = s.0.join(name);
    let header_text: String =
        VCF_TEXT.lines().filter(|l| l.starts_with('#')).collect::<Vec<_>>().join("\n");
    let header = BcfHeader::parse(&header_text).unwrap();
    let mut w = BcfWriter::new_bgzf(File::create(&path).unwrap(), &header).unwrap();
    for (i, line) in VCF_TEXT.lines().filter(|l| !l.starts_with('#')).enumerate() {
        let rec = VcfRecord::parse(line, 0, i as u64 + 1).unwrap();
        w.write_record(&rec).unwrap();
    }
    w.finish().unwrap();
    // BCF uses CSI: TBI needs a text configuration, which a binary file lacks.
    let kind = IndexKind::Csi { min_shift: 14, depth: 5 };
    let index = build_bcf_index(File::open(&path).unwrap(), kind).unwrap();
    index.write_to(File::create(s.0.join(format!("{name}.csi"))).unwrap()).unwrap();
    path
}

/// Runs the binary with `dir` as the current directory, so `tptbio.toml` is found there.
fn tptbio_in(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_tptbio")).current_dir(dir).args(args).output().unwrap()
}

#[test]
fn completions_print_script_for_each_shell() {
    for shell in ["bash", "zsh", "fish", "powershell", "elvish"] {
        let out = tptbio(&["completions", shell]);
        assert_eq!(code(&out), EXIT_OK, "{shell}: {}", stderr(&out));
        assert!(stdout(&out).contains("tptbio"), "{shell}: {}", stdout(&out));
    }
}

#[test]
fn completions_reject_unknown_shell() {
    let out = tptbio(&["completions", "tcsh"]);
    assert_eq!(code(&out), EXIT_USAGE);
}

#[test]
fn man_writes_one_page_per_command() {
    let s = Scratch::new("man");
    let dir = s.0.join("man1");
    let out = tptbio(&["man", dir.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    for page in ["tptbio.1", "tptbio-seq.1", "tptbio-seq-view.1", "tptbio-stats-normalize.1"] {
        let text = fs::read_to_string(dir.join(page)).unwrap_or_else(|e| panic!("{page}: {e}"));
        assert!(text.contains("tptbio"), "{page}");
    }
}

#[test]
fn seq_view_vcf_header_and_native_records() {
    let s = Scratch::new("vcf-view");
    let input = s.write("in.vcf", VCF_TEXT);
    let out = tptbio(&["seq", "view", "-H", input.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).starts_with("##fileformat=VCFv4.2"), "{}", stdout(&out));
    assert!(!stdout(&out).contains("chr1\t100"));

    let out = tptbio(&["seq", "view", input.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert_eq!(stdout(&out), VCF_TEXT);
}

#[test]
fn seq_view_vcf_csv_jsonl_and_table_output() {
    let s = Scratch::new("vcf-formats");
    let input = s.write("in.vcf", VCF_TEXT);
    let path = input.to_str().unwrap();

    let out = tptbio(&["--output-format", "csv", "seq", "view", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert_eq!(
        stdout(&out),
        format!(
            "{CSV_HEADER}chr1,100,rs1,A,G,50,PASS,DP=10\nchr1,200,.,C,T,.,PASS,DP=20\nchr1,300,rs3,G,\"A,C\",60,q10,DP=30\n"
        )
    );

    let out = tptbio(&["--output-format", "jsonl", "seq", "view", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    let first = stdout(&out).lines().next().unwrap().to_owned();
    assert!(first.contains("\"POS\":100") && first.contains("\"QUAL\":50"), "{first}");

    let out = tptbio(&["--output-format", "table", "seq", "view", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    // Columns are padded to their widest value, so match the row by its fields.
    assert!(
        stdout(&out).lines().any(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            f.first() == Some(&"chr1") && f.get(1) == Some(&"100") && f.get(2) == Some(&"rs1")
        }),
        "{}",
        stdout(&out)
    );
}

#[test]
fn seq_view_bgzipped_vcf_region_uses_tabix() {
    let s = Scratch::new("vcf-gz");
    let input = write_vcf_gz(&s, "in.vcf.gz");
    let path = input.to_str().unwrap();
    let out = tptbio(&["--output-format", "csv", "seq", "view", "-r", "chr1:150-250", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert_eq!(stdout(&out), format!("{CSV_HEADER}chr1,200,.,C,T,.,PASS,DP=20\n"));

    let out = tptbio(&["seq", "view", "-H", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).contains("#CHROM"), "{}", stdout(&out));
}

#[test]
fn seq_view_bcf_header_records_and_region() {
    let s = Scratch::new("bcf");
    let input = write_bcf(&s, "in.bcf");
    let path = input.to_str().unwrap();

    let out = tptbio(&["seq", "view", "-H", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).contains("##contig=<ID=chr1"), "{}", stdout(&out));

    let out = tptbio(&["--output-format", "csv", "seq", "view", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).contains("chr1,300,rs3,G,\"A,C\",60,q10,DP=30\n"), "{}", stdout(&out));

    let out = tptbio(&["--output-format", "csv", "seq", "view", "-r", "chr1:1-150", path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert_eq!(stdout(&out), format!("{CSV_HEADER}chr1,100,rs1,A,G,50,PASS,DP=10\n"));
}

#[test]
fn seq_view_rejects_extension_that_contradicts_contents() {
    let s = Scratch::new("mislabelled");
    let bcf = write_bcf(&s, "data.bcf");
    let renamed = s.0.join("renamed.vcf");
    fs::copy(&bcf, &renamed).unwrap();
    let out = tptbio(&["seq", "view", "-H", renamed.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(
        stderr(&out).contains("extension says VCF, but the contents are BCF"),
        "{}",
        stderr(&out)
    );
}

#[test]
fn seq_view_variant_misuse_is_runtime_error() {
    let s = Scratch::new("variant-misuse");
    let plain = s.write("plain.vcf", VCF_TEXT);
    let out = tptbio(&["seq", "view", "-r", "chr1:1-10", plain.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("bgzipped VCF"), "{}", stderr(&out));

    let out = tptbio(&["seq", "view", "-b", plain.to_str().unwrap()]);
    assert_eq!(code(&out), EXIT_RUNTIME);
}

#[test]
fn export_format_follows_flag_then_config_then_default() {
    let s = Scratch::new("precedence-format");
    let sam = s.write("in.sam", SAM_TEXT);
    let cfg = s.write("cfg.toml", "output_format = \"jsonl\"\n");
    let sam_path = sam.to_str().unwrap();
    let cfg_path = cfg.to_str().unwrap();

    // Default: CSV.
    let out = tptbio(&["seq", "export", sam_path]);
    assert!(stdout(&out).starts_with("qname,"), "{}", stdout(&out));

    // The config file beats the default.
    let out = tptbio(&["--config", cfg_path, "seq", "export", sam_path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).starts_with('{'), "{}", stdout(&out));

    // The flag beats the config file.
    let out = tptbio(&["--config", cfg_path, "--output-format", "csv", "seq", "export", sam_path]);
    assert!(stdout(&out).starts_with("qname,"), "{}", stdout(&out));

    // A table request on export falls back to CSV with a warning.
    let out = tptbio(&["--output-format", "table", "seq", "export", sam_path]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stdout(&out).starts_with("qname,"), "{}", stdout(&out));
    assert!(stderr(&out).contains("no table layout"), "{}", stderr(&out));
}

#[test]
fn log_level_flag_beats_config_and_config_beats_default() {
    let s = Scratch::new("precedence-log");
    let input = s.write("in.csv", "x\n1\n3\n");
    let quiet_cfg = s.write("quiet.toml", "log_level = \"error\"\n");
    let input_path = input.to_str().unwrap();
    let quiet_path = quiet_cfg.to_str().unwrap();

    // The default level shows the info message.
    let out = tptbio(&["stats", "normalize", input_path]);
    assert!(stderr(&out).contains("normalized"), "{}", stderr(&out));

    // The config hides it.
    let out = tptbio(&["--config", quiet_path, "stats", "normalize", input_path]);
    assert_eq!(code(&out), EXIT_OK);
    assert!(stderr(&out).is_empty(), "{}", stderr(&out));

    // A verbose flag beats the config.
    let out = tptbio(&["--config", quiet_path, "-v", "stats", "normalize", input_path]);
    assert!(stderr(&out).contains("normalized"), "{}", stderr(&out));
}

#[test]
fn default_config_file_is_read_from_current_directory() {
    let s = Scratch::new("default-config");
    s.write("in.csv", "x\n1\n3\n");
    s.write("tptbio.toml", "log_level = \"error\"\n");
    let out = tptbio_in(&s.0, &["stats", "normalize", "in.csv"]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stderr(&out).is_empty(), "{}", stderr(&out));

    // An explicit --config path beats the file in the current directory.
    let other = s.write("other.toml", "log_level = \"info\"\n");
    let out =
        tptbio_in(&s.0, &["--config", other.to_str().unwrap(), "stats", "normalize", "in.csv"]);
    assert!(stderr(&out).contains("normalized"), "{}", stderr(&out));
}

#[test]
fn config_threads_is_recorded_at_debug_level() {
    let s = Scratch::new("threads");
    let input = s.write("in.csv", "x\n1\n3\n");
    let cfg = s.write("cfg.toml", "threads = 4\n");
    let out = tptbio(&[
        "--config",
        cfg.to_str().unwrap(),
        "-v",
        "stats",
        "normalize",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code(&out), EXIT_OK, "{}", stderr(&out));
    assert!(stderr(&out).contains("threads = 4"), "{}", stderr(&out));
}

#[test]
fn config_errors_are_runtime_errors() {
    let s = Scratch::new("config-errors");
    let input = s.write("in.csv", "x\n1\n3\n");
    let input_path = input.to_str().unwrap();

    let unknown = s.write("unknown.toml", "colour = \"blue\"\n");
    let out = tptbio(&["--config", unknown.to_str().unwrap(), "stats", "normalize", input_path]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("invalid config"), "{}", stderr(&out));

    let bad_level = s.write("level.toml", "log_level = \"trace\"\n");
    let out = tptbio(&["--config", bad_level.to_str().unwrap(), "stats", "normalize", input_path]);
    assert_eq!(code(&out), EXIT_RUNTIME);

    let zero = s.write("zero.toml", "threads = 0\n");
    let out = tptbio(&["--config", zero.to_str().unwrap(), "stats", "normalize", input_path]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("threads must be at least 1"), "{}", stderr(&out));

    let missing = s.0.join("missing.toml");
    let out = tptbio(&["--config", missing.to_str().unwrap(), "stats", "normalize", input_path]);
    assert_eq!(code(&out), EXIT_RUNTIME);
    assert!(stderr(&out).contains("cannot read config"), "{}", stderr(&out));
}
