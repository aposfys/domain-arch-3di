# Results

A pilot run. 58 proteins across two clades (transporter 29, globular 29), 20 distinct domain architectures.

## Setup

| | |
| --- | --- |
| Transporter | IPR020846, Major facilitator superfamily domain |
| Globular | IPR001254, Serine protease, trypsin domain |
| Characters | amino acid vs 3Di, alignment-free 3-mer cosine, identical for both alphabets |
| Tree | neighbour joining (Biopython) |
| Masking | residues below pLDDT 70 masked, proteins below 50% confident excluded |
| Excluded | 2 proteins |

## Topology

Robinson-Foulds distance **80**, normalised **0.727**. Most of the non-trivial splits in one tree are absent from the other. There is no bootstrap or resampling baseline, so how much two alignment-free neighbour-joining trees differ from noise alone is not known.

## Events

Events are read off sister pairs (cherries), the one place two extant architectures can be compared without reconstructing an ancestor. Sister paralogues have no direction without an outgroup, so an indel is not split into an addition and a deletion.

| | Sequence tree | Structure tree |
| --- | ---: | ---: |
| Cherries | 21 | 18 |
| no change | 10 | 13 |
| terminal indel | 3 | 1 |
| internal indel | 1 | 1 |
| complex | 7 | 3 |

The sequence tree's cherries show 11 architecture changes out of 21, the structure tree's 5 out of 18. The difference is not significant (two-sided Fisher exact p = 0.19, treating cherries as independent, which overstates the evidence because the two trees share 9 of them). Only 9 sister pairs are shared between the trees, a Jaccard of 0.30.

So the two alphabets give different trees and largely different sister pairs. Whether they imply different numbers of rearrangements is not settled by this run.

## Limitations

- **Complex includes Pfam family swaps.** Pfam splits the MFS clan (CL0015) into several families, for example PF07690 (MFS_1) and PF00083 (Sugar_tr). Two single-domain MFS proteins with different family calls count as complex although no domain was gained, lost or moved, so the complex counts may overstate real rearrangements.
- **Alignment-free distances are coarse.** Both alphabets go through an identical 3-mer cosine distance, which removes the confound of comparing two substitution matrices that were never calibrated against each other. The cost is resolution. These trees are weaker than a model-based inference would give, and the RF distance is correspondingly noisier.
- **Cherries are a small sample.** Around twenty sister pairs per tree is too few to estimate rates or to separate a modest difference from chance.
- **Two human clades, not a phylogeny.** These are paralogues within one species. The result is about how characters change an inference, not about the evolution of these families.
- **Inputs not recorded.** This run did not record its accession list or the UniProt, InterPro and Foldseek releases. The protein set is the first 30 UniProt hits per clade, which depends on the release, so a rerun can give different numbers.
