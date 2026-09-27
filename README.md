# gbz-haplotype-index

Writes a haplotype index for a [gbz-base](https://github.com/jltsiren/gbz-base)
database. gbz-base reports every walk but the query path as `unknown`. With the
index beside the database, @gmod/gbz-base names every walk in a query with its
haplotype and coordinates.

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 --anchor-spacing 131072 graph.gbz graph.haplotype-index.db
gbz-haplotype-index --from-db graph.gbz.db graph.haplotype-index.db
```

[@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads the index beside the
database, locally or over HTTP. The
[documentation](https://github.com/GMOD/gbz-base-js/blob/main/docs/haplotype-index.md)
describes the tables and how a query uses them.
