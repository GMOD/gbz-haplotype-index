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
which a query that keeps a few haplotypes and fills snarls needs. The index also
carries an [overview](docs/haplotype-index.md#the-overview) of every haplotype
in bins along each reference path, for views of megabases or a whole
chromosome. The index is written in format 3, which @gmod/gbz-base reads from
version 6.1. Earlier versions read only the format 2 index that
gbz-haplotype-index 0.2 wrote, and fail on a format 3 file with "SQLite table
b-tree contains an index page", so keep a format 2 file hosted while clients
on 6.0 remain.

## Docs

- [The haplotype index](docs/haplotype-index.md)
  - [Samples and anchors](docs/haplotype-index.md#samples-and-anchors)
  - [Building the haplotype index](docs/haplotype-index.md#building-the-haplotype-index)
  - [Using the haplotype index](docs/haplotype-index.md#using-the-haplotype-index)
  - [Querying a subset of the haplotypes](docs/haplotype-index.md#querying-a-subset-of-the-haplotypes)
  - [How a query identifies walks](docs/haplotype-index.md#how-a-query-identifies-walks)
  - [Proof that the keep route finds every piece](docs/haplotype-index.md#proof-that-the-keep-route-finds-every-piece)
  - [The overview](docs/haplotype-index.md#the-overview)
  - [Measured on HPRC v2.1](docs/haplotype-index.md#measured-on-hprc-v21)
- [Formal proof that the keep route finds every piece](docs/proof.md):
  the same argument as definitions, lemmas and a theorem
- [DEVELOPING.md](DEVELOPING.md): tests, and keeping the indexer and the reader
  in agreement
