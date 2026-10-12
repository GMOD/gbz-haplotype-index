# gbz-haplotype-index

Names every walk in a [gbz-base](https://github.com/jltsiren/gbz-base) pangenome
query by its haplotype.

- **In:** a GBZ graph and the `graph.gbz.db` that `gbz-base construct` built
  from it
- **Out:** `graph.haplotype-index.db`, a SQLite file to host beside the database

Then [@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads both files and
reports the sample, haplotype, contig and coordinate of every walk in a window.
With the database alone, gbz-base labels each walk other than the query path
`unknown#N`. The tool leaves `graph.gbz.db` as it is, so an index also works
with a database someone else hosts, such as HPRC's

![graph.gbz feeds gbz-base construct, which writes graph.gbz.db, and gbz-haplotype-index, which writes graph.haplotype-index.db. @gmod/gbz-base reads both. With the database alone a query names the reference walk and prints the rest as unknown#N#chr6; with the haplotype index it names each one, such as HG002#1#chr6, with its coordinate](docs/img/overview.svg)

## Build an index

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 graph.gbz graph.gbz.db graph.haplotype-index.db

# without the GBZ, reading the paths from the database
gbz-haplotype-index --interval 16384 --from-db graph.gbz.db graph.haplotype-index.db
```

- Give `graph.gbz.db` with the GBZ: the tool reads the top-level snarls from it,
  which a query that keeps a few haplotypes and fills snarls needs
- `--interval` trades size for speed: a query follows each walk up to that many
  bp to find the row that names it
- The HPRC v2.1 GRCh38 graph (10 GB database, 464 haplotypes) gives a 5.1 GB
  index in 35 minutes on 22 threads; its chr22 alone gives 66 MB. The
  [recipe](docs/haplotype-index.md#building-hprc-release-21s-index) lists the
  downloads and options

## Read it

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

## What the index holds

- **Samples:** a row every `--interval` bp along each path that maps a GBWT
  position to the path's name and coordinate, like the sampled suffix array of
  an FM-index
- **Anchors:** a reference node every `--anchor-spacing` bp with a row for each
  haplotype visit, so a query that uses the `keep` option walks the chosen
  haplotypes from one anchor to the next
- **Bins and stray rows:** the nodes near each 16 kb of reference, and the
  stretches of path there that a walk between anchors would miss
- **Overview:** each haplotype's class (absent, variant, partial or
  reference-like) per bin along each reference path, at several zoom levels

## Alternative: gfa-to-tabix

[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) `--walks` writes named,
Tabix-indexed BED files of each haplotype's walk. On 60 chr22 windows of HPRC
v2.1 it returned the same walks as gbz-base with this index in 56. It needs
fewer requests and takes any GFA. gbz-base grows sublinearly with haplotypes
and keeps a walk whole where it jumps to a paralog, which on chr22 affects
under 0.1% of steps.

See the [full comparison](docs/comparison.md).

## Docs

- [The haplotype index](docs/haplotype-index.md): overview of format and measurements on HPRC v2.1
- [Proof that the keep route finds every piece](docs/proof.md): mathematical proof that this actually is a valid approach
- [Comparison with gfa-to-tabix](docs/comparison.md): similarities and differences to our 'tabix based' gfa format
- [DEVELOPING.md](DEVELOPING.md): tests, and keeping the indexer and the reader
  in agreement

## See also

- [@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) - range-request reader
  for `.gbz.db` databases
- [jbrowse-plugin-graphgenomeviewer](https://github.com/GMOD/jbrowse-plugin-graphgenomeviewer) -
  JBrowse 2 plugin that browses these graphs by locus
- [BandageJS](https://github.com/cmdcolin/BandageJS) - standalone page for GFA
  and gbz-base graphs
- [gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) - indexes a GFA by genome
  coordinate

Tutorials on [jbrowse.org](https://jbrowse.org/jb2/docs/tutorials/)

- [Hosting your own graph](https://jbrowse.org/jb2/docs/tutorials/pangenome_prepare_graph/)
- [HPRC part 1: graph alleles and haplotypes](https://jbrowse.org/jb2/docs/tutorials/pangenome_hprc/)
