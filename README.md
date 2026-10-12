# gbz-haplotype-index

Names every walk in a [gbz-base](https://github.com/jltsiren/gbz-base) pangenome
query by its haplotype.

- **In:** a GBZ graph and the `graph.gbz.db` that `gbz-base construct` built
  from it
- **Out:** `graph.haplotype-index.db`, a SQLite file to host beside the database
- **Then:** [@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reads both
  files and reports the sample, haplotype, contig and coordinate of every walk
  in a window. With the database alone, gbz-base labels each walk other than
  the query path `unknown#N`
- The tool leaves `graph.gbz.db` as it is, so an index also works with a
  database someone else hosts, such as HPRC's

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
- The tool writes format 3, the one format @gmod/gbz-base reads. Rebuild an
  index that gbz-haplotype-index 0.2 or earlier wrote

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
  every haplotype, for views of megabases or a whole chromosome in about a
  dozen requests

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
BED files that the JBrowse graph track reads. Compared with gbz-base and this
index, the walk files have advantages:

- **Any GFA:** it reads the GFA that pggb, odgi, vg or Minigraph-Cactus write,
  so a graph needs no GBZ
- **Standard files:** htslib, `tabix` and any static web host handle BGZF and
  Tabix, and a reader needs no SQLite
- **One set of files:** each row names its path as `sample#haplotype#contig`
  and each file's header lists every haplotype, so naming a walk takes no
  second index
- **More views:** `gfa-to-tabix build` also writes a coarse tier for whole
  chromosomes, an allele inventory and the graph's bubbles

They also have weaknesses:

- **Approximate placement:** a node takes one coordinate, from its first visit,
  so a haplotype that passes a collapsed repeat files that pass under the
  reference's first copy, and the graph track draws it in more fragments. On 45 chr22 windows with
  8 haplotypes, 41 matched gbz-base's cut
- **Paths off the reference:** a path that visits no node of a reference has
  no rows under it. gbz-base returns every walk through the window's nodes,
  and its keep route has a [proof](docs/proof.md) and a fuzzer behind it
- **One reference per set:** the tool writes three files per `--refs` sample,
  where one `.gbz.db` serves queries on GRCh38 and CHM13 alike
- **Size:** on HPRC v2.1 chr22 the GRCh38 and CHM13 sets total 408 MB without
  sequences, against 196 MB for the chr22 `.gbz.db`, which holds the
  sequences, and 66 MB for its haplotype index
- **Less from one query:** gbz-base expands `context` bp around the reference,
  fills snarls, and returns each walk's GBWT positions and coordinates; a
  Tabix query returns the rows filed under the window's chunks
- **Base-level graphs only:** `--walks` needs integer node ids and W or P
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
