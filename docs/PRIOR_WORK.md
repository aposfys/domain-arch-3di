# Prior work, and the gap this repository addresses

Structural phylogenetics with Foldseek's 3Di alphabet (van Kempen et al., *Nature
Biotechnology* 42, 243-246, 2024, published online 2023,
[doi:10.1038/s41587-023-01773-0](https://doi.org/10.1038/s41587-023-01773-0)) is an active
field, and two preprints define its current practice.

- **Puente-Lelievre et al.** (bioRxiv preprint,
  [doi:10.1101/2023.12.12.571181](https://doi.org/10.1101/2023.12.12.571181)) treat 3Di
  letters as standard phylogenetic characters, with IQ-TREE maximum likelihood, a
  3Di-specific rate matrix, partitioning and ultrafast bootstrap. Combining amino acids
  with 3Di best matches a reference structural-distance tree and avoids long-branch
  attraction.
- **Fullmer et al.** (bioRxiv preprint 2025,
  [doi:10.1101/2025.06.30.662300](https://doi.org/10.1101/2025.06.30.662300)) find that 3Di
  combined with sequence resolves better than either alone, and that the gain is weaker in
  alpha-helical proteins, because high helical content reduces the information 3Di
  alignments carry.

Both ask whether 3Di improves phylogenetic resolution. Neither asks whether it changes the
downstream evolutionary events read off the tree, which is the question here.

Two caveats follow. Fullmer's alpha-helix result bears directly on the clades chosen here.
MFS transporters are almost entirely helical, so this is close to the regime where 3Di
carries least information. And the method used here, a shared alignment-free 3-mer cosine
distance with neighbour joining, is deliberately weaker than the field's. It was chosen so
both alphabets pass through an identical step and no uncalibrated substitution matrix is
compared against another. That control is worth having, but it means the event difference
in this repository cannot yet be separated from method noise. Redoing it with model-based
inference and bootstrap support is the next step, and until then the result is a
motivation for that work rather than a finding.
