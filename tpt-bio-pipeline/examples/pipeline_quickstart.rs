//! Pipeline quickstart: parse a TOML pipeline, run it, and drive a plain `Dag` directly.
//! Run with:
//!
//! ```text
//! cargo run -p tpt-bio-pipeline --example pipeline_quickstart
//! ```
//!
//! The pipeline only uses an in-memory table, so no input files are needed.

use std::convert::Infallible;
use std::error::Error;
use std::sync::atomic::AtomicBool;

use tpt_bio_pipeline::{Dag, NodeId, Pipeline, Value};

type Failure = Box<dyn Error>;

const STATS: &str = r#"
[[nodes]]
name = "table"
type = "table_source"
params = { values = [[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]] }

[[nodes]]
name = "scaled"
type = "normalize"
params = { method = "zscore" }
inputs = ["table"]

[[nodes]]
name = "reduced"
type = "pca"
params = { components = 1 }
inputs = ["scaled"]

[[nodes]]
name = "groups"
type = "kmeans"
params = { k = 2, seed = 1 }
inputs = ["reduced"]
"#;

fn main() -> Result<(), Failure> {
    // 1. Parse and validate the definition. Unknown kinds, bad parameters and cycles
    //    are reported here, before anything runs.
    let pipeline = Pipeline::from_toml(STATS)?;

    // 2. Run every node. Outputs come back in definition order: table, scaled, reduced,
    //    groups.
    for (i, value) in pipeline.run()?.iter().enumerate() {
        match value {
            Value::Table(m) => {
                println!("output {i}: table, {} rows x {} columns", m.rows(), m.cols());
            }
            Value::Labels(labels) => println!("output {i}: labels {labels:?}"),
            other => println!("output {i}: {other:?}"),
        }
    }

    // 3. The same graph without TOML: nodes carry payloads, edges say what runs after
    //    what, and the executor passes each node its predecessors' results.
    let mut dag: Dag<u64> = Dag::new();
    let a = dag.add_node(2);
    let b = dag.add_node(3);
    let c = dag.add_node(4);
    dag.add_edge(a, c)?;
    dag.add_edge(b, c)?;

    let cancel = AtomicBool::new(false);
    let sums = dag.execute_parallel(&cancel, |_: NodeId, payload: &u64, inputs: &[&u64]| {
        Ok::<u64, Infallible>(payload + inputs.iter().copied().sum::<u64>())
    })?;
    println!("dag results in node order: {sums:?}");
    Ok(())
}
