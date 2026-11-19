# gbz-haplotype-index

Names every walk in a [gbz-base](https://github.com/jltsiren/gbz-base) pangenome
query by its haplotype.

[@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads the gbz-base
database (`graph.gbz.db`) and the haplotype index (`graph.haplotype-index.db`)
and reports the sample, haplotype, contig and coordinate of every walk in a
window. With `graph.gbz.db` alone, a query labels each walk other than the query
path `unknown#N`. `gbz-haplotype-index` leaves `graph.gbz.db` unchanged, so an
index also works with a database someone else hosts, such as HPRC's.

To show haplotypes in a browser, we recommend
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix)'s walk files over this
index ([why](#gfa-to-tabix-or-gbz-haplotype-index)).

![graph.gbz feeds gbz-base construct, which writes graph.gbz.db, and gbz-haplotype-index, which writes graph.haplotype-index.db. @gmod/gbz-base reads both. With the database alone a query names the reference walk and prints the rest as unknown#N#chr6; with the haplotype index it names each one, such as HG002#1#chr6, with its coordinate](docs/img/overview.svg)

## Build an index

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 graph.gbz graph.gbz.db graph.haplotype-index.db

# without the GBZ, reading the paths from the database
gbz-haplotype-index --interval 16384 --from-db graph.gbz.db graph.haplotype-index.db
```

- The inputs are a GBZ graph (`graph.gbz`) and the `graph.gbz.db` that
  `gbz-base construct` built from it. Give both, because the tool reads the
  top-level snarls from `graph.gbz.db`. Without them, a query that uses the
  `keep` option and fills a snarl identifies every walk
- `--interval` sets the bp between samples along each path (default 4096). A
  larger interval gives a smaller index, and a query follows each walk further
  to reach a sample of its path
- The HPRC v2.1 GRCh38 graph (10 GB database, 464 haplotypes) gives a 5.1 GB
  index in 35 minutes on 22 threads; its chr22 alone gives 66 MB. The
  [recipe](docs/haplotype-index.md#building-hprc-release-21s-index) lists the
  downloads and options

## Query with the index

```ts
const db = await GBZBase.open(new RemoteFile(graphUrl), {
  haplotypeIndex: new RemoteFile(indexUrl),
})
const records = await db.getAlignmentsForRange('GRCh38#0#chr6', start, end, {
  keep: name => name.sample === 'HG002',
})
```

```bash
gbz-base-query https://host/graph.gbz.db \
  --haplotype-index https://host/graph.haplotype-index.db --keep HG002 ...
```

- Without `keep`, a query returns every haplotype in the window
- With `keep`, it returns the reference, the chosen haplotypes and the nodes
  they visit. It walks the chosen haplotypes from anchor to anchor, so its cost
  grows with the number chosen
- `haplotypeOverview` and `gbz-base-query --overview` read a per-bin summary of
  every haplotype, for views of megabases or a whole chromosome in about a dozen
  requests

## gfa-to-tabix or gbz-haplotype-index

[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) `--walks` writes every
haplotype's path through the graph, with its name, into one file per reference
genome, which a browser reads from any static web host. We recommend these walk
files over gbz-base with this index for showing haplotypes in a browser:

- A browser reads a 64 kb window from a walk file in a few byte ranges, where
  gbz-base with this index takes tens of round trips
- On HPRC v2.1, the GRCh38 walk file without bases is a fraction of the size
  of the gbz-base database, and smaller than this index alone
- gfa-to-tabix reads the GFA that Minigraph-Cactus and pggb write, with P or W
  lines, where gbz-base needs a GBZ

Use gbz-base with this index instead when you already serve a gbz-base
database, such as HPRC's, or need what walk files lack:

- coordinates on any path. gfa-to-tabix gives coordinates on the reference
  genomes chosen when it writes the files
- an overview of every haplotype across a whole chromosome
- a graph with many thousands of haplotypes, since a gbz-base database grows
  more slowly than the number of haplotypes
- each walk kept whole where it jumps to a paralog. gfa-to-tabix files the
  steps after such a jump under the paralog's position, which affects under
  0.1% of chr22's steps

Sizes, request counts and timings for both tools, and for gfaidx, are in
gfa-to-tabix's
[comparison](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/comparison.md)
and its
[details](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/comparison-details.md).

## Docs

- [The haplotype index](docs/haplotype-index.md): the tables, the two routes a
  query takes, the overview, and measurements on HPRC v2.1
- [Proof that the keep route finds every piece](docs/proof.md): a query that
  uses the `keep` option returns every piece of each chosen haplotype's walk
- [DEVELOPING.md](DEVELOPING.md): tests, and keeping the indexer and the reader
  in agreement

## See also

- [@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) - range-request reader
  for `.gbz.db` databases
- [jbrowse-plugin-graphgenomeviewer](https://github.com/GMOD/jbrowse-plugin-graphgenomeviewer) -
  JBrowse 2 plugin that browses these graphs by locus
- [BandageJS](https://github.com/cmdcolin/BandageJS) - browser viewer for GFA
  files, gfa-to-tabix walk files and gbz-base databases
- [gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) - indexes a GFA by genome
  coordinate

Tutorials on [jbrowse.org](https://jbrowse.org/jb2/docs/tutorials/)

- [Hosting your own graph](https://jbrowse.org/jb2/docs/tutorials/pangenome_prepare_graph/)
- [HPRC part 1: graph alleles and haplotypes](https://jbrowse.org/jb2/docs/tutorials/pangenome_hprc/)
