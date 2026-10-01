# gbz-haplotype-index

Writes a haplotype index for a [gbz-base](https://github.com/jltsiren/gbz-base)
database. gbz-base reports every walk but the query path as `unknown`. With the
index beside the database, [@gmod/gbz-base](https://github.com/GMOD/gbz-base-js)
names every walk in a query with its haplotype and coordinates.

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 graph.gbz graph.gbz.db graph.haplotype-index.db
gbz-haplotype-index --interval 16384 --from-db graph.gbz.db graph.haplotype-index.db
```

Give `graph.gbz.db` with the GBZ: the tool reads the top-level snarls from it,
which a query that keeps a few haplotypes and fills snarls needs.

The [documentation](docs/haplotype-index.md) describes the tables, how a query
uses them, and the proof that the keep route finds every piece.

## Developing

```bash
cargo test --release
```

The stray rows are complete only for the walks that @gmod/gbz-base's keep route
takes, so `src/strays.rs` here and
[`src/chosenPaths.ts`](https://github.com/GMOD/gbz-base-js/blob/main/src/chosenPaths.ts)
there state the same rule. A change to that rule on either side must bump the
`haplotype_index_stray_format` tag that `src/main.rs` writes and the reader
tests for; a reader that meets a format it does not know identifies every walk
instead of trusting the rows. The
[proof](docs/haplotype-index.md#proof-that-the-keep-route-finds-every-piece)
names the comparisons on both sides that the rule depends on.

Nothing here checks that the two sides agree. The fuzzer in gbz-base-js does:
build this crate, then run it from a gbz-base-js checkout with the binary as the
indexer.

```bash
cargo build --release
node test/fuzz/run.ts --seeds 0..150 --scale medium --jobs 4 \
  --indexer ../gbz-haplotype-index/target/release/gbz-haplotype-index
```

The `fuzz` job in gbz-base-js's `push.yml` lists the index options CI uses.
After a change to the walk rule, `test/data/build-indexes.sh` there rebuilds the
fixture indexes and reports the ones that change.

The test fixture `test/data/split-contig.gbz.db` is a copy of the one in
gbz-base-js.
