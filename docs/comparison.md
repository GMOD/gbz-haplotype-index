# Comparison with gfa-to-tabix

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

- **Segmental duplications:** where a haplotype steps onto the reference nodes
  of a distant paralog and back, the walk files put the steps after the jump
  under the paralog's chunk, so a window at the haplotype's own locus shows it
  in more pieces and leaves out some of its inserted nodes. gbz-base reaches
  those nodes by graph distance from the window. On HPRC v2.1 chr22 the walk
  files place under 0.1% of steps this way, in at most 74 of 1,130 paths,
  nearly all in the 22q11 low-copy repeats
- **Growth with haplotypes:** a per-haplotype walk store grows linearly with the
  haplotypes, an estimated 6 to 8 GB genome-wide at 464 and 50 to 70 GB at
  4,000. The GBWT in a `.gbz.db` grows sublinearly, because haplotypes share
  runs
- **One file for every reference:** one `.gbz.db` serves queries on GRCh38 and
  CHM13, where gfa-to-tabix writes a set of walk files per `--refs` sample
- **Proof:** the keep route has a [proof](proof.md) and a fuzzer behind it
- **Input:** `--walks` needs a base-level graph with integer node ids and W or P
  lines, so a minigraph rGFA gets nodes and links but no walks
