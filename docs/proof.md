# Formal proof that the keep route finds every piece

A query that uses the `keep` option returns, for each kept haplotype, the
pieces of its walk: the maximal stretches whose nodes all lie in the window's
subgraph. It finds them by walking the haplotype out from the anchor rows and
along its stray rows. This page proves that whenever the query returns without
falling back to identifying every walk, it has found every piece. It states the
argument of [the proof on the haplotype index page](haplotype-index.md#proof-that-the-keep-route-finds-every-piece)
as definitions, lemmas and a theorem, without that page's code identifiers.
Both follow the code on both sides, [`src/strays.rs`](../src/strays.rs) here
and [`src/chosenPaths.ts`](https://github.com/GMOD/gbz-base-js/blob/main/src/chosenPaths.ts)
in gbz-base-js, and name the comparisons they rest on, so a change to either
program can be checked against them. The fuzzer in gbz-base-js tests the same
claim.

## Objects

Every node $x$ of the graph has a length $\ell(x) \ge 1$. A path $P$ is the
sequence of visits $v_0, \dots, v_{m-1}$ of its walk, taken in the forward
orientation of the path. Visit $v_i$ is to node $x_i$ at offset $o_i$ along
$P$, with $o_0 = 0$ and $o_{i+1} = o_i + \ell(x_i)$, and it may pass its node
in either orientation. Each visit has a GBWT position, and a position belongs
to one visit of one path in one orientation. Two consecutive visits of a path
are joined by an edge of the graph. For a set of nodes $S$, a **piece** of $P$
in $S$ is a maximal run of consecutive visits $v_i, \dots, v_j$ whose nodes all
lie in $S$.

The index is built with an anchor spacing $s$ (`--anchor-spacing`), a bin
length $\beta$ (`--stray-bin`), a reach $d$ for the walks (`--stray-bound`) and
a context $c$ (`--stray-context`). It covers a set $\mathcal{R}$ of reference
paths, those of one reference sample. The query names one of them, $R$, as its
reference, and $|R|$ is the length of $R$.

**Definition 1 (Anchors).** For each $R' \in \mathcal{R}$ and each integer $k$
with $0 \le k \le \lfloor |R'|/s \rfloor$, the indexer chooses one visit of
$R'$, the **anchor** of multiple $k$, and records its node $n_{R'}(k)$ and its
offset $a_{R'}(k)$ along $R'$. For $k = 0$ it is the first visit of $R'$. For
$k \ge 1$ it is a visit whose node overlaps the interval
$[ks - \lfloor s/2 \rfloor, ks)$ of offsets along $R'$: among those, the one
whose node has the most GBWT positions, and the last on a tie.

**Lemma 1 (Anchor order).** For every $R' \in \mathcal{R}$, writing $a$ for
$a_{R'}$: every multiple from $0$ to $\lfloor |R'|/s \rfloor$ has an anchor,
and $a(k) \le a(k+1)$, with equality only when the two anchors are the same
visit.

*Proof.* When $ks \le |R'|$, the interval $[ks - \lfloor s/2 \rfloor, ks)$ lies
within the offsets of $R'$, so some visit overlaps it. The visits of $R'$
occupy disjoint intervals of offsets, in order along $R'$. The anchor of $k$
starts before $ks$, since its interval overlaps one that ends at $ks$. The
anchor of $k+1$ ends after $(k+1)s - \lfloor s/2 \rfloor$, which exceeds $ks$.
If the anchor of $k+1$ were a different visit starting before the anchor of
$k$, it would end at or before $a(k) < ks$, a contradiction. So
$a(k+1) \ge a(k)$, and equal offsets mean the same visit. ∎

**Definition 2 (Cuts and sections).** A node is an **anchor node** of the
sample when it is $n_{R'}(k)$ for some $R' \in \mathcal{R}$ and $k$. It is
**shared** when it is so for more than one pair $(R', k)$, and otherwise
**simple**, with the label $(R', k)$. The **cuts** of a path $P$ are its visits
to anchor nodes, in either orientation. A **section** of $P$ on $R$ is a pair
of consecutive cuts, with no cut between them, that are simple with the labels
$(R, k)$ and $(R, k+1)$ for some $k$; the visits of the section are the two
cuts and every visit between them. For a section write $u$ for the offset
along $P$ of the cut labelled $k$ and $u'$ for that of the cut labelled $k+1$,
and $a = a_R(k)$ and $a' = a_R(k+1)$. By Lemma 1, $a < a'$: equality would
make the two anchors one visit and their node shared.

**Definition 3 (Bins and listed nodes).** Bin $b$ of $R$ is the interval
$[b_{\mathrm{lo}}, b_{\mathrm{hi}}) = [b\beta, (b+1)\beta)$ of offsets along
$R$. The index lists for it a set of nodes $N_R(b)$: the nodes within $c$ bp of
a node of $R$ that overlaps the bin, found by the search over node sides that
the query's context expansion runs. The proof uses $N_R(b)$ only through the
query's check in Definition 6, so it holds whatever that search finds.

**Definition 4 (Reached visits and stray rows).** Let $b$ be a bin of $R$ and
$v_t$ a visit of a section of $P$ on $R$. The section **reaches** $v_t$ for $b$
when one of these holds:

- (i) $a < b_{\mathrm{hi}}$ and $a' \ge b_{\mathrm{lo}}$;
- (ii) $a \ge b_{\mathrm{hi}}$, $a - b_{\mathrm{hi}} \le d$ and $|o_t - u| \le d$;
- (iii) $a' < b_{\mathrm{lo}}$, $b_{\mathrm{lo}} - a' \le d$ and $|o_t - u'| \le d$.

A visit $v_t$ of $P$ with $x_t \in N_R(b)$ is a **stray** of $(R, b, P)$ when
no section of $P$ on $R$ that contains $v_t$ reaches it for $b$. A visit that
is a cut lies in up to two sections, any other visit in up to one. The index
holds for $(R, b, P)$ a set of **stray rows**, each an interval $[p, q]$ of
offsets along $P$ that starts at a stray and holds that stray's GBWT position,
such that the offset of every stray of $(R, b, P)$ lies in a row. The indexer
closes a row and opens the next where two strays lie more than `--stray-gap`
apart; the proof does not depend on where.

**Definition 5 (Snarls).** A top-level snarl is given by a chain link of the
graph database between two boundary nodes $y$ and $z$. Its **region** is the
set of nodes a search reaches from the successors of $y$ and the predecessors
of $z$, following the edges on both sides of every node it reaches and never
entering $y$ or $z$. The indexer **models** the graph's snarls when no two
regions share a node, and leaves out a region of more than $2^{20}$ nodes. A
path $P$ whose nodes all lie in one region that the indexer kept has a **snarl
row** for every $R' \in \mathcal{R}$ and bin $b$ of $R'$ with
$\min(y, z) \in N_{R'}(b)$. The row names the snarl by $\min(y, z)$ and
$\max(y, z)$ and holds the GBWT position of $v_0$ and the interval
$[0, o_{m-1}]$.

## What the proof takes from the two programs

- **A1.** The index holds a sample at every visit of every path to an anchor
  node, in both orientations of the path, so the rows at an anchor node list
  every visit to it. Every sample, anchor and stray row records a visit's GBWT
  position together with the offset along its path at which the node starts,
  the quantity a walk tracks.
- **A2.** The indexer and the query fill a snarl with the search of
  Definition 5 from the same chain link. In a GBZ the graph's edges are the
  GBWT's, so both searches follow the same edges, and one of those edges joins
  any two consecutive visits of a path. The query uses snarl rows only when the
  index's tags say that it modelled snarls and record the same number of chain
  links as the database; the count does not compare the links themselves.
- **A3.** The index was built from the graph the query reads. The query
  compares the path and node counts of the two when it opens them, which would
  miss a different graph with the same counts.

## What the query does

**Definition 6 (The query).** The query has a window $W = [w, w')$ of offsets
along $R$, which touches the bins $F = \lfloor w/\beta \rfloor$ to
$L = \lfloor (w'-1)/\beta \rfloor$ of $R$. Write $\mathit{lo} = F\beta$ and
$\mathit{hi} = (L+1)\beta$, so that
$\mathit{lo} \le b_{\mathrm{lo}} < b_{\mathrm{hi}} \le \mathit{hi}$ for every
bin $b$ from $F$ to $L$. The query builds the subgraph $S$, a set of nodes: the
nodes of $R$'s walk through $W$, the nodes within its context of them, and the
regions of the snarls it fills. $I \subseteq S$ is the set of nodes the fills
added. The query then:

1. checks that its context is at most $c$; that every node of $S \setminus I$
   lies in $N_R(b)$ for some $b$ from $F$ to $L$; and that for every snarl it
   filled, the index models snarls, the region has at most $2^{20}$ nodes,
   $y \ne z$, and both $y$ and $z$ were in $S$ before any fill;
2. reads the anchors of $R$ for every multiple from $k_{\min}$ to $k_{\max}$,
   extending the range until either $k_{\min} = 0$ or
   $a_R(k_{\min}) < \mathit{lo} - d$, and until either multiple $k_{\max}$ has
   no anchor or $a_R(k_{\max}) > \mathit{hi} + d$; and reads the samples at
   both orientations of each anchor node it read, keeping those of kept paths
   in the forward orientation;
3. for each kept path $P$, takes the visits of step 2 in order of offset as
   its cuts, labelling each by the one multiple its node anchors among those
   read, or as shared when the node anchors two; and for each pair of
   consecutive cuts labelled $k$ and $k+1$, with $a = a_R(k)$,
   $a' = a_R(k+1)$ and $u$, $u'$ the offsets of the two cuts, plans
   - (i) a **section walk** over $[\min(u, u') - d,\ \max(u, u') + d]$ when
     $a < \mathit{hi}$ and $a' \ge \mathit{lo}$;
   - (ii) otherwise a walk over $[u - d, u + d]$ when $a \ge \mathit{hi}$ and
     $a - \mathit{hi} \le d$;
   - (iii) otherwise a walk over $[u' - d, u' + d]$ when $a' < \mathit{lo}$
     and $\mathit{lo} - a' \le d$;
4. reads the stray rows and snarl rows of $(R, b, P)$ for every $b$ from $F$
   to $L$ and every kept $P$, keeps the snarl rows of the snarls it filled,
   and plans a walk over each row's interval;
5. runs the planned walks (Definition 7), skipping one whose interval lies
   within an interval that an earlier walk covered;
6. falls back to identifying every walk when a check of step 1 fails, when
   more than 32 kept paths pass the anchors it read, when a walk runs farther
   outside $S$ than its cap, when a sample of a kept path on a node of $S$
   lies in no recorded piece, or when it cannot find the other orientation of
   a piece that needs one.

**Definition 7 (Walks).** A walk over an interval $[p, q]$ of offsets along a
kept path $P$ starts from a GBWT position of $P$ at a visit of offset
$p_0 \in [p, q]$, which an anchor row or a stray row supplies together with
$p_0$. Forward from there, while the current visit's offset is at most $q$,
the query records the piece of $P$ through the visit if its node lies in $S$,
continues from the last visit of that piece, and steps to the next visit.
Backward from the start, while the current visit's offset exceeds $p$, it
steps to the previous visit and, if that visit's node lies in $S$, records the
piece through it and continues from the first visit of the piece. Both
directions stop at an end of $P$. The query then records as covered the
interval $[p, q]$ widened to the offsets the walk reached.

## Proof

**Lemma 2 (Coverage).** After a walk over $[p, q]$, every visit of $P$ to a
node of $S$ whose offset lies in $[p, q]$ lies in a recorded piece. The same
holds for a planned walk the query skipped.

*Proof.* Forward, the walk processes every visit from its start through the
last visit whose offset is at most $q$, or through the end of $P$. Backward,
it processes every visit from its start down to the first whose offset is at
most $p$, or down to the start of $P$. Offsets increase along $P$, so every
visit with an offset in $[p, q]$ lies in that stretch, except when $P$ ends
first, and then no visit lies beyond the end. Each visit of the stretch the
walk either steps on, recording the piece through it when its node is in $S$,
or passes inside a piece it recorded. The covered interval the query records
is $[p, q]$ widened to the offsets the walk reached, and a visit in the
widening lies in a recorded piece by the same argument. A skipped walk's
interval lies within such an interval. ∎

**Lemma 3 (The query reads both anchors of a reaching section).** Let a
section of $P$ on $R$ labelled $k$ and $k+1$ reach some visit for a bin $b$
from $F$ to $L$. Then $k_{\min} \le k$ and $k + 1 \le k_{\max}$.

*Proof.* Each case of Definition 4 gives $a \le \mathit{hi} + d$ and
$a' \ge \mathit{lo} - d$. In case (i), $a < b_{\mathrm{hi}} \le \mathit{hi}$
and $a' \ge b_{\mathrm{lo}} \ge \mathit{lo}$. In case (ii),
$a \le b_{\mathrm{hi}} + d \le \mathit{hi} + d$ and
$a' > a \ge b_{\mathrm{hi}} > \mathit{lo}$. In case (iii),
$a' \ge b_{\mathrm{lo}} - d \ge \mathit{lo} - d$ and
$a < a' < b_{\mathrm{lo}} \le \mathit{hi}$.

Suppose $k < k_{\min}$. Then $k_{\min} \ge 1$, so
$a_R(k_{\min}) < \mathit{lo} - d$, and $k + 1 \le k_{\min}$ gives
$a' \le a_R(k_{\min})$ by Lemma 1, a contradiction. Suppose
$k + 1 > k_{\max}$, so $k \ge k_{\max}$. If multiple $k_{\max}$ has an anchor,
then $a \ge a_R(k_{\max}) > \mathit{hi} + d$, a contradiction. Otherwise,
since the multiples with an anchor run from $0$ without a gap (Lemma 1), no
multiple above $k_{\max}$ has one, so $k + 1$ has none, and the section does
not exist. ∎

**Lemma 4 (The query sees the section).** Let a section of $P$ on $R$
labelled $k$ and $k+1$ reach some visit for a bin $b$ from $F$ to $L$. Then
in step 3 the query finds the two cuts of the section consecutive among its
cuts of $P$, labelled $k$ and $k+1$, and plans for them by the three cases.

*Proof.* By Lemma 3 the query reads the anchors of $k$ and $k+1$, so their
nodes are among the anchor nodes it read, and by A1 the forward samples at
those nodes are exactly the visits of $P$ to them. So both cuts of the section
are cuts for the query. Every cut of the query is a visit to an anchor node of
$R$ and hence a cut for the indexer, and the section has no indexer cut
between its two, so the query has none between them either, and the two are
consecutive for the query. The node of the cut labelled $k$ is simple for the
indexer, so it anchors exactly one pair over all of $\mathcal{R}$, namely
$(R, k)$, and so exactly one multiple among those read, so the query labels it
$k$. The same argument labels the other cut $k+1$. ∎

**Lemma 5 (Every reached visit is covered).** Let a section of $P$ on $R$
reach the visit $v_t$ for a bin $b$ from $F$ to $L$, and let $x_t \in S$. Then
$v_t$ lies in a recorded piece.

*Proof.* By Lemma 4 the query plans for the section with the $a$, $a'$, $u$
and $u'$ of Definition 4, and $o_t$ lies between $u$ and $u'$, since $v_t$ is a
visit of the section.

In case (i), $a < b_{\mathrm{hi}} \le \mathit{hi}$ and
$a' \ge b_{\mathrm{lo}} \ge \mathit{lo}$, so the query plans the section walk,
whose interval $[\min(u, u') - d, \max(u, u') + d]$ contains $o_t$.

In case (ii), if $a < \mathit{hi}$ then
$a' > a \ge b_{\mathrm{hi}} > \mathit{lo}$ and the query plans the section
walk, whose interval contains every offset within $d$ of $u$. Otherwise
$a \ge \mathit{hi}$ and $a - \mathit{hi} \le a - b_{\mathrm{hi}} \le d$, so it
plans the walk over $[u - d, u + d]$. Either interval contains $o_t$, since
$|o_t - u| \le d$.

In case (iii), if $a' \ge \mathit{lo}$ then
$a < a' < b_{\mathrm{lo}} < \mathit{hi}$ and the query plans the section
walk. Otherwise $a' < \mathit{lo}$ and
$\mathit{lo} - a' \le b_{\mathrm{lo}} - a' \le d$, so it plans the walk over
$[u' - d, u' + d]$. Either interval contains $o_t$, since $|o_t - u'| \le d$.

In each case Lemma 2 puts $v_t$ in a recorded piece. ∎

**Theorem.** Suppose the query returns without falling back. Then for every
kept path $P$, every piece of $P$ in $S$ is a piece the query recorded.

*Proof.* Let $v_t$ be a visit of $P$ to a node $x_t \in S$. We show that $v_t$
lies in a recorded piece. The query records the piece through a visit as the
maximal run of $P$ in $S$ that contains the visit, so the recorded piece
containing $v_t$ is the piece of $P$ in $S$ that contains $v_t$, and every
piece of $P$ in $S$, which contains some visit, is then recorded.

Suppose $x_t \notin I$. The check of step 1 found $x_t \in N_R(b)$ for some
$b$ from $F$ to $L$. If a section of $P$ on $R$ containing $v_t$ reaches it
for $b$, Lemma 5 applies. Otherwise $v_t$ is a stray of $(R, b, P)$, so its
offset lies in a stray row of $(R, b, P)$, which the query read in step 4 and
walked or skipped in step 5, and Lemma 2 applies.

Suppose $x_t \in I$. Then $x_t$ lies in the region of a snarl the query
filled, with boundary nodes $y$ and $z$. By the check of step 1 the index
models snarls and the region has at most $2^{20}$ nodes, so by A2 the indexer
computed the same region and kept it, and $y, z \in S \setminus I$. Let
$v_i, \dots, v_j$ be the maximal run of consecutive visits of $P$ that
contains $v_t$ and whose nodes lie in the region. If $j < m - 1$, the edge from
$x_j$ to $x_{j+1}$ leaves a node of the region, so the search of Definition 5
followed it, and $x_{j+1}$ lies in the region or is $y$ or $z$; by maximality
it is $y$ or $z$. Likewise, if $i > 0$ then $x_{i-1}$ is $y$ or $z$, since the
search follows the edges on both sides of $x_i$. Such a visit $v_{j+1}$ or
$v_{i-1}$ is to a node of $S \setminus I$, so by the first case it lies in a
recorded piece, and since the nodes of $v_i, \dots, v_j$ lie in the region and
so in $S$, $v_t$ lies in the same piece. If $i = 0$ and $j = m - 1$, every node
of $P$ lies in the region, and the indexer wrote a snarl row for $P$ in every
bin $b$ of $R$ with $\min(y, z) \in N_R(b)$. The check of step 1 found
$\min(y, z) \in N_R(b)$ for some $b$ from $F$ to $L$, and the query filled this
snarl, so step 4 keeps the row and its walk over $[0, o_{m-1}]$ covers $v_t$
by Lemma 2. ∎

## What the proof leaves to the checks

The walks find each piece in the forward orientation of its path. The query
returns a piece in the orientation whose end nodes are canonical, so for a
piece kept in the other orientation it reads the reverse-orientation samples
along the path until one gives a position in that piece, and falls back when
none does. The proof uses neither the context check nor any agreement between
the indexer's search for $N_R(b)$ and the query's context expansion. The check
that every node of $S \setminus I$ is listed establishes what the proof needs,
and the agreement only keeps that check from failing. The sample check of
step 6 never fails while the theorem holds; it guards against an index built
for another graph (A3) and against a change to the rule on one side alone.

The proof rests on inclusive comparisons on both sides. Step 2 must keep
reading anchors while $a_R(k_{\min}) \ge \mathit{lo} - d$ and while
$a_R(k_{\max}) \le \mathit{hi} + d$; with a strict comparison in the second,
an anchor at exactly $\mathit{hi} + d$ stops the reading before the multiple
past it, and a section that reaches visits by case (ii) goes unplanned. That
is the omission the audit found, which
[`test/data/anchor-at-bound.gfa`](https://github.com/GMOD/gbz-base-js/blob/main/test/data/anchor-at-bound.gfa)
in gbz-base-js now exercises. The three cases of Definition 4 and of step 3
use the same comparisons, and the walks of Definition 7 process the visits at
both ends of their interval.
