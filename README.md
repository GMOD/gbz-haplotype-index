# gbz-haplotype-index

Names every walk in a [gbz-base](https://github.com/jltsiren/gbz-base) pangenome
query by its haplotype.

[@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads the gbz-base
database (`graph.gbz.db`) and the haplotype index (`graph.haplotype-index.db`)
and reports the sample, haplotype, contig and coordinate of every walk in the
queried region. With `graph.gbz.db` alone, @gmod/gbz-base names the queried path
and labels every other walk `unknown#N`. `gbz-haplotype-index` leaves
`graph.gbz.db` unchanged, so you can build an index for a database someone else
hosts, such as HPRC's.

To show haplotypes in a browser, gfa-to-tabix's walk files take fewer requests
and less storage than gbz-base with this index, and we recommend them. Use this
index when you already use a gbz-base database or need its whole-chromosome
overview ([gfa-to-tabix or gbz-haplotype-index](#gfa-to-tabix-or-gbz-haplotype-index)).

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

[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) writes every haplotype's
path through the graph, with its name, into one walk file per reference genome
(`gfa-to-tabix graph.gfa --walks --refs GRCh38`). JBrowse's graph plugin reads
walk files from any static web host. We recommend them over gbz-base with this
index for showing haplotypes in a browser:

- A browser reads a 64 kb window of a walk file in one to four range requests,
  where gbz-base with this index takes tens of round trips
- On HPRC v2.1, the GRCh38 and CHM13 walk files, without bases, are together a
  third of the size of the gbz-base database, which holds the bases, and
  smaller than this index
- gfa-to-tabix reads the GFA that Minigraph-Cactus and pggb write. gbz-base
  needs a GBZ, which `vg gbwt -G graph.gfa --gbz-format -g graph.gbz` makes,
  and then `gbz-base construct`

Use gbz-base with this index instead when you want to use a gbz-base database
that you or someone else already hosts, such as HPRC's, or need what walk files
lack:

- an overview of every haplotype across a whole chromosome
- each walk kept whole where it jumps to a paralog. When a haplotype's walk
  steps onto the nodes of a distant paralog and back, gfa-to-tabix stores the
  steps after the jump at the paralog's position, so at the haplotype's own
  locus a browser shows the walk in several pieces, with a record of each jump.
  This affects under 0.1% of steps on HPRC v2.1 Minigraph-Cactus chr22

The gbz-base database grows far more slowly than the number of haplotypes, but
this index stores a sample every `--interval` bp along every path, so it grows
with them.

Sizes, request counts and timings for gfa-to-tabix, gbz-base and gfaidx are in
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
