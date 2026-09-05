use gbz::{GBWT, GBZ, Orientation, Pos, ENDMARKER};
use gbz::bwt::BWT;
use gbz::support;
use rusqlite::{params, Connection};
use simple_sds::serialize;

use std::env;
use std::process;

const USAGE: &str = "Usage: gbz-haplotype-index [--interval BP] [--forward-only] graph.gbz graph.gbz.db

Walks every path in the GBZ in both orientations and writes a sample every
--interval bp (default 4096) into table HaplotypeSamples of the database,
plus the path lengths into HaplotypeLengths. The path start and end are
always sampled. Existing tables are replaced.
";

struct Args {
    gbz: String,
    db: String,
    interval: usize,
    forward_only: bool,
}

fn parse_args() -> Args {
    let mut positional = Vec::new();
    let mut interval = 4096;
    let mut forward_only = false;
    let mut iter = env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--interval" => {
                let value = iter.next().unwrap_or_default();
                interval = value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid --interval: {}", value);
                    process::exit(1);
                });
            }
            "--forward-only" => forward_only = true,
            "-h" | "--help" => {
                eprint!("{}", USAGE);
                process::exit(0);
            }
            _ => positional.push(arg),
        }
    }
    if positional.len() != 2 || interval == 0 {
        eprint!("{}", USAGE);
        process::exit(1);
    }
    Args { gbz: positional.remove(0), db: positional.remove(0), interval, forward_only }
}

struct Sample {
    node_handle: usize,
    node_offset: usize,
    path_handle: usize,
    orientation: usize,
    path_offset: usize,
}

fn walk(graph: &GBZ, index: &GBWT, path_handle: usize, orientation: Orientation, interval: usize) -> (Vec<Sample>, usize) {
    let sequence_id = support::encode_path(path_handle, orientation);
    let mut samples = Vec::new();
    let mut pos = index.start(sequence_id);
    let mut offset = 0;
    let mut next_sample = 0;
    let mut last: Option<(Pos, usize)> = None;
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let node_len = graph.sequence_len(support::node_id(current.node)).unwrap();
        if offset >= next_sample {
            samples.push(Sample {
                node_handle: current.node,
                node_offset: current.offset,
                path_handle,
                orientation: orientation as usize,
                path_offset: offset,
            });
            next_sample = offset + interval;
            last = None;
        } else {
            last = Some((current, offset));
        }
        offset += node_len;
        let bwt: &BWT = index.as_ref();
        let record = bwt.record(index.node_to_record(current.node)).unwrap();
        pos = record.lf(current.offset);
    }
    if let Some((end, end_offset)) = last {
        samples.push(Sample {
            node_handle: end.node,
            node_offset: end.offset,
            path_handle,
            orientation: orientation as usize,
            path_offset: end_offset,
        });
    }
    let length = offset;
    if orientation == Orientation::Reverse {
        for sample in samples.iter_mut() {
            let node_len = graph.sequence_len(support::node_id(sample.node_handle)).unwrap();
            sample.path_offset = length - sample.path_offset - node_len;
        }
    }
    (samples, length)
}

fn main() {
    let args = parse_args();
    let graph: GBZ = serialize::load_from(&args.gbz).unwrap_or_else(|e| {
        eprintln!("Cannot load {}: {}", args.gbz, e);
        process::exit(1);
    });
    let index: &GBWT = graph.as_ref();
    let metadata = graph.metadata().unwrap_or_else(|| {
        eprintln!("The GBZ has no path metadata");
        process::exit(1);
    });
    let mut connection = Connection::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", args.db, e);
        process::exit(1);
    });

    connection.execute_batch(
        "DROP TABLE IF EXISTS HaplotypeSamples;
         DROP TABLE IF EXISTS HaplotypeLengths;
         CREATE TABLE HaplotypeSamples (
             node_handle INTEGER NOT NULL,
             node_offset INTEGER NOT NULL,
             path_handle INTEGER NOT NULL,
             orientation INTEGER NOT NULL,
             path_offset INTEGER NOT NULL,
             PRIMARY KEY (node_handle, node_offset)
         ) STRICT;
         CREATE TABLE HaplotypeLengths (
             path_handle INTEGER PRIMARY KEY,
             length INTEGER NOT NULL
         ) STRICT;",
    ).unwrap();

    let transaction = connection.transaction().unwrap();
    let mut inserted = 0;
    {
        let mut insert_sample = transaction.prepare(
            "INSERT INTO HaplotypeSamples(node_handle, node_offset, path_handle, orientation, path_offset) VALUES (?1, ?2, ?3, ?4, ?5)"
        ).unwrap();
        let mut insert_length = transaction.prepare(
            "INSERT INTO HaplotypeLengths(path_handle, length) VALUES (?1, ?2)"
        ).unwrap();
        let mut insert_tag = transaction.prepare(
            "INSERT OR REPLACE INTO Tags(key, value) VALUES (?1, ?2)"
        ).unwrap();
        let orientations: &[Orientation] = if args.forward_only {
            &[Orientation::Forward]
        } else {
            &[Orientation::Forward, Orientation::Reverse]
        };
        for path_handle in 0..metadata.paths() {
            let mut length = 0;
            for &orientation in orientations {
                let (samples, walked) = walk(&graph, index, path_handle, orientation, args.interval);
                length = walked;
                for sample in samples {
                    insert_sample.execute(params![
                        sample.node_handle as i64, sample.node_offset as i64, sample.path_handle as i64,
                        sample.orientation as i64, sample.path_offset as i64
                    ]).unwrap();
                    inserted += 1;
                }
            }
            insert_length.execute(params![path_handle as i64, length as i64]).unwrap();
        }
        insert_tag.execute(params!["haplotype_index_interval", args.interval.to_string()]).unwrap();
        insert_tag.execute(params!["haplotype_index_orientations", if args.forward_only { "forward" } else { "both" }]).unwrap();
    }
    transaction.commit().unwrap();
    eprintln!("Inserted {} samples for {} paths", inserted, metadata.paths());
}
