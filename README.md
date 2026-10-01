# gbz-haplotype-index

Writes a haplotype index for a [gbz-base](https://github.com/jltsiren/gbz-base)
database. gbz-base reports every walk but the query path as `unknown`. With the
index beside the database, @gmod/gbz-base names every walk in a query with its
haplotype and coordinates.

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 graph.gbz graph.gbz.db graph.haplotype-index.db
gbz-haplotype-index --interval 16384 --from-db graph.gbz.db graph.haplotype-index.db
```

Give `graph.gbz.db` with the GBZ: the tool reads the top-level snarls from it,
which a query that keeps a few haplotypes and fills snarls needs.

The stray rows are complete only for the walks that @gmod/gbz-base's keep route
takes, so `strays.rs` and `src/chosenPaths.ts` state the same rule. A change to
that rule on either side must bump the `haplotype_index_stray_format` tag that
`main.rs` writes and the reader tests for; a reader that meets a format it does
not know identifies every walk instead of trusting the rows. The
[proof](https://github.com/GMOD/gbz-base-js/blob/main/docs/haplotype-index.md#proof-that-the-keep-route-finds-every-piece)
that the keep route finds every piece names the comparisons on both sides that
the rule depends on.

[@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads the index beside the
database, locally or over HTTP. The
[documentation](https://github.com/GMOD/gbz-base-js/blob/main/docs/haplotype-index.md)
describes the tables and how a query uses them.
