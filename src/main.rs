use gbz::bwt::BWT;
use gbz::support;
use gbz::{GBWT, GBZ, Orientation, Pos, ENDMARKER};
use gbz_base::{GBZBase, GraphInterface};
use rusqlite::{params, Connection};
use simple_sds::serialize;

use std::env;
use std::process;

const USAGE: &str = "Usage: gbz-haplotype-index [--interval BP] [--forward-only] graph.gbz graph.gbz.db
       gbz-haplotype-index [--interval BP] [--forward-only] --from-db graph.gbz.db

Walks every path in both orientations and writes a sample every --interval bp
(default 4096) into table HaplotypeSamples of the database, plus the path
lengths into HaplotypeLengths. The path start and end are always sampled.
Existing tables are replaced. With --from-db the walk reads node records from
the database itself, so the GBZ is not needed.
";

struct Args {
    gbz: Option<String>,
    db: String,
    interval: usize,
    forward_only: bool,
}

fn parse_args() -> Args {
    let mut positional = Vec::new();
    let mut interval = 4096;
    let mut forward_only = false;
    let mut from_db = false;
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
            "--from-db" => from_db = true,
            "-h" | "--help" => {
                eprint!("{}", USAGE);
                process::exit(0);
            }
            _ => positional.push(arg),
        }
    }
    let expected = if from_db { 1 } else { 2 };
    if positional.len() != expected || interval == 0 {
        eprint!("{}", USAGE);
        process::exit(1);
    }
    let db = positional.pop().unwrap();
    let gbz = positional.pop();
    Args { gbz, db, interval, forward_only }
}

struct Sample {
    node_handle: usize,
    node_offset: usize,
    path_handle: usize,
    orientation: usize,
    path_offset: usize,
}

trait PathSource {
    fn path_count(&mut self) -> usize;
    fn start(&mut self, path_handle: usize, orientation: Orientation) -> Option<Pos>;
    fn step(&mut self, pos: Pos) -> (usize, Option<Pos>);
    fn node_len(&mut self, handle: usize) -> usize;
}

struct GbzSource {
    graph: GBZ,
}

impl PathSource for GbzSource {
    fn path_count(&mut self) -> usize {
        self.graph.metadata().map(|m| m.paths()).unwrap_or(0)
    }

    fn start(&mut self, path_handle: usize, orientation: Orientation) -> Option<Pos> {
        let index: &GBWT = self.graph.as_ref();
        index.start(support::encode_path(path_handle, orientation))
    }

    fn step(&mut self, pos: Pos) -> (usize, Option<Pos>) {
        let node_len = self.graph.sequence_len(support::node_id(pos.node)).unwrap();
        let index: &GBWT = self.graph.as_ref();
        let bwt: &BWT = index.as_ref();
        let record = bwt.record(index.node_to_record(pos.node)).unwrap();
        (node_len, record.lf(pos.offset))
    }

    fn node_len(&mut self, handle: usize) -> usize {
        self.graph.sequence_len(support::node_id(handle)).unwrap()
    }
}

struct DbSource<'a> {
    interface: GraphInterface<'a>,
    paths: usize,
}

impl<'a> PathSource for DbSource<'a> {
    fn path_count(&mut self) -> usize {
        self.paths
    }

    fn start(&mut self, path_handle: usize, orientation: Orientation) -> Option<Pos> {
        let path = self.interface.get_path(path_handle).unwrap()?;
        Some(if orientation == Orientation::Forward { path.fw_start } else { path.rev_start })
    }

    fn step(&mut self, pos: Pos) -> (usize, Option<Pos>) {
        let record = self.interface.get_record(pos.node).unwrap().unwrap();
        (record.sequence_len(), record.to_gbwt_record().lf(pos.offset))
    }

    fn node_len(&mut self, handle: usize) -> usize {
        self.interface.get_record(handle).unwrap().unwrap().sequence_len()
    }
}

fn walk(source: &mut dyn PathSource, path_handle: usize, orientation: Orientation, interval: usize) -> (Vec<Sample>, usize) {
    let mut samples = Vec::new();
    let mut pos = source.start(path_handle, orientation);
    let mut offset = 0;
    let mut next_sample = 0;
    let mut last: Option<(Pos, usize)> = None;
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let (node_len, next) = source.step(current);
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
        pos = next;
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
            let node_len = source.node_len(sample.node_handle);
            sample.path_offset = length - sample.path_offset - node_len;
        }
    }
    (samples, length)
}

fn create_tables(connection: &Connection) {
    connection
        .execute_batch(
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
        )
        .unwrap();
}

fn run(source: &mut dyn PathSource, connection: &mut Connection, args: &Args) {
    let orientations: Vec<Orientation> = if args.forward_only {
        vec![Orientation::Forward]
    } else {
        vec![Orientation::Forward, Orientation::Reverse]
    };
    let paths = source.path_count();
    let mut inserted = 0;
    let batch = 64;
    let mut handle = 0;
    while handle < paths {
        let transaction = connection.transaction().unwrap();
        {
            let mut insert_sample = transaction
                .prepare("INSERT INTO HaplotypeSamples(node_handle, node_offset, path_handle, orientation, path_offset) VALUES (?1, ?2, ?3, ?4, ?5)")
                .unwrap();
            let mut insert_length = transaction.prepare("INSERT INTO HaplotypeLengths(path_handle, length) VALUES (?1, ?2)").unwrap();
            for path_handle in handle..(handle + batch).min(paths) {
                let mut length = 0;
                for &orientation in orientations.iter() {
                    let (samples, walked) = walk(source, path_handle, orientation, args.interval);
                    length = walked;
                    for sample in samples {
                        insert_sample
                            .execute(params![
                                sample.node_handle as i64,
                                sample.node_offset as i64,
                                sample.path_handle as i64,
                                sample.orientation as i64,
                                sample.path_offset as i64
                            ])
                            .unwrap();
                        inserted += 1;
                    }
                }
                insert_length.execute(params![path_handle as i64, length as i64]).unwrap();
            }
        }
        transaction.commit().unwrap();
        handle += batch;
        eprintln!("{} / {} paths, {} samples", handle.min(paths), paths, inserted);
    }
    connection
        .execute(
            "INSERT OR REPLACE INTO Tags(key, value) VALUES ('haplotype_index_interval', ?1)",
            params![args.interval.to_string()],
        )
        .unwrap();
    connection
        .execute(
            "INSERT OR REPLACE INTO Tags(key, value) VALUES ('haplotype_index_orientations', ?1)",
            params![if args.forward_only { "forward" } else { "both" }],
        )
        .unwrap();
    eprintln!("Inserted {} samples for {} paths", inserted, paths);
}

fn path_count_from_db(db: &str) -> usize {
    let connection = Connection::open(db).unwrap();
    let value: String = connection
        .query_row("SELECT value FROM Tags WHERE key = 'paths'", [], |row| row.get(0))
        .unwrap_or_else(|_| "0".to_string());
    value.parse().unwrap_or(0)
}

fn main() {
    let args = parse_args();
    let mut connection = Connection::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", args.db, e);
        process::exit(1);
    });
    create_tables(&connection);
    match &args.gbz {
        Some(gbz) => {
            let graph: GBZ = serialize::load_from(gbz).unwrap_or_else(|e| {
                eprintln!("Cannot load {}: {}", gbz, e);
                process::exit(1);
            });
            if graph.metadata().is_none() {
                eprintln!("The GBZ has no path metadata");
                process::exit(1);
            }
            let mut source = GbzSource { graph };
            run(&mut source, &mut connection, &args);
        }
        None => {
            let paths = path_count_from_db(&args.db);
            let database = GBZBase::open(&args.db).unwrap_or_else(|e| {
                eprintln!("Cannot open {} as a GBZ-base: {}", args.db, e);
                process::exit(1);
            });
            let interface = GraphInterface::new(&database).unwrap();
            let mut source = DbSource { interface, paths };
            run(&mut source, &mut connection, &args);
        }
    }
}
