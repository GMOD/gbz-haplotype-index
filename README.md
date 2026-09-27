# gbz-haplotype-index

Writes a haplotype index for a [gbz-base](https://github.com/jltsiren/gbz-base)
database. The index is a companion SQLite file that names every walk in a query,
with its haplotype and coordinates, where gbz-base reports all but the query
path as `unknown`.

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 --anchor-spacing 131072 graph.gbz graph.haplotype-index.db
gbz-haplotype-index --from-db graph.gbz.db graph.haplotype-index.db
```

[@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads the index beside the
database, locally or over HTTP. The
[documentation](https://github.com/GMOD/gbz-base-js/blob/main/docs/haplotype-index.md)
describes the tables and how a query uses them.
