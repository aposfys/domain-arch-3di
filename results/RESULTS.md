# Results

58 proteins across two clades (transporter 29, globular 29), 20 distinct domain architectures.

## Setup

| | |
| --- | --- |
| Transporter | IPR020846 — Major facilitator superfamily domain |
| Globular | IPR001254 — Serine protease, trypsin domain |
| Characters | amino acid vs 3Di, alignment-free 3-mer cosine, identical for both alphabets |
| Tree | neighbour joining (Biopython) |
| Masking | residues below pLDDT 70 masked; proteins below 50% confident excluded |
| Excluded | 2 proteins |

## The topology changes

Robinson-Foulds distance **80**, normalised **0.727**. Roughly three-quarters of the bipartitions in one tree are absent from the other, so switching from amino acids to 3Di is not a small perturbation of the phylogeny.

## And so do the events

Events are read off sister pairs (cherries), which is the one place two extant architectures can be compared without reconstructing an ancestor.

| | Sequence tree | Structure tree |
| --- | ---: | ---: |
| Cherries | 21 | 18 |
| complex | 7 | 3 |
| no change | 10 | 13 |
| internal deletion | 1 | 1 |
| terminal addition | 3 | 1 |

**The sequence tree implies 11 rearrangement events; the structure tree implies 5.** Only 9 sister pairs are shared between the two trees, a Jaccard of 0.30.

So the answer to the question this repository asks is that **both change**. It is not the case that 3Di reshuffles the tree while leaving the evolutionary story intact — the set of inferred rearrangements moves too, and in this dataset the structural characters imply substantially fewer of them.

A consistency check: of the 9 cherries both trees found, **0 were classified differently**. That number has to be zero — the same pair of architectures must classify the same way whichever tree produced it — and it is reported rather than assumed, because a non-zero value would mean the classifier depends on something other than the architectures.

## Limitations, stated plainly

- **Alignment-free distances are coarse.** Both alphabets go through an identical 3-mer cosine distance, which removes the confound of comparing two substitution matrices that were never calibrated against each other. The cost is resolution: these trees are weaker than a model-based inference would give, and the RF distance is correspondingly noisier.
- **Cherries are a small sample.** Around twenty sister pairs per tree is enough to show the event sets differ and not enough to estimate rates.
- **Two human clades, not a phylogeny.** These are paralogues within one species. The result is about how characters change an inference, not about the evolution of these families.
- **0 proteins carry no Pfam domain at all**, and contribute only identity events.
