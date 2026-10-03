# Developing

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

The tables' layout is a separate contract, the `haplotype_index_format` tag
(3). `src/encode.rs` holds the blob encodings; the reader decodes them in its
`db.ts`. A change to a table or a blob bumps the tag, and the reader refuses a
format it does not know.

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
fixture indexes and reports the ones that change. The fixtures use
`--page-size 4096` to stay small; micb-kir3dl1 keeps the default 64 KiB pages,
and its format 2 index stays beside it for the old format's tests.

The test fixture `test/data/split-contig.gbz.db` is a copy of the one in
gbz-base-js.

`docs/CLAUDE.md` holds the writing rules for the pages in `docs/`.
