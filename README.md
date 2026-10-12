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

[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) `--walks` writes each
haplotype's path through the graph, named, as bgzip-compressed, Tabix-indexed
BED files that the JBrowse graph track reads. A window on those files and the
same window on gbz-base with this index return nearly the same walks:

- On 60 chr22 windows of HPRC v2.1 with eight haplotypes, both returned the same
  haplotypes and alignments in 56. The other four lie in the 22q11 repeats and
  the immunoglobulin lambda locus
- At the ABCA7 VNTR, all 462 walk rows match gbz-base's to the bp

The walk files have advantages:

- **Fewer requests:** at the HPRC tutorial's loci, JBrowse drew a window from
  the walk files in 5 to 7 requests and 2.4 to 7.9 s, and from gbz-base in 17 to
  57 requests and 6.4 to 13.3 s
- **Any GFA:** the tool reads the GFA that pggb, odgi, vg or Minigraph-Cactus
  write, so a graph needs no GBZ
- **Standard files:** htslib, `tabix` and any static web host handle BGZF and
  Tabix, and a reader needs no SQLite
- **Names in every row:** each row names its path as `sample#haplotype#contig`
  and the header lists every haplotype, so naming a walk takes no second index
- **More views:** `gfa-to-tabix build` also writes a coarse tier for whole
  chromosomes, an allele inventory and the graph's bubbles

gbz-base and this index have advantages too:

- **Repeats:** the walk files give each node one coordinate, so a pass over a
  collapsed repeat copy lands under the reference's first copy. gbz-base follows
  each walk through the GBWT, and the two cuts differ in places such as the
  22q11 repeats
- **Growth with haplotypes:** a per-haplotype walk store grows linearly with the
  haplotypes, an estimated 6 to 8 GB genome-wide at 464 and 50 to 70 GB at
  4,000. The GBWT in a `.gbz.db` grows sublinearly, because haplotypes share
  runs
- **One file for every reference:** one `.gbz.db` serves queries on GRCh38 and
  CHM13, where gfa-to-tabix writes a set of walk files per `--refs` sample
- **Proof:** the keep route has a [proof](docs/proof.md) and a fuzzer behind it
- **Input:** `--walks` needs a base-level graph with integer node ids and W or P
  lines, so a minigraph rGFA gets nodes and links but no walks

## Docs

- [The haplotype index](docs/haplotype-index.md): the tables, every build
  option, how a query names walks, the overview, and measurements on HPRC v2.1
- [Proof that the keep route finds every piece](docs/proof.md): first against
  the code on both sides, then as definitions, lemmas and a theorem
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
