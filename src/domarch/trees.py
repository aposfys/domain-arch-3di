"""Distance matrices and trees, from amino acids and from 3Di.

The comparison is only fair if the two alphabets go through an identical method, so both
use the same alignment-free *k*-mer distance and the same neighbour-joining step. Both
alphabets have twenty letters, so the feature spaces are the same size and neither is
advantaged by the representation.

Alignment-free is a deliberate limitation rather than a shortcut. A substitution-matrix
alignment would need a 3Di matrix and an amino acid matrix that were calibrated against
each other, and they are not -- any difference in the resulting trees would then be partly
a difference between two matrices. Sharing one method removes that confound at the cost of
resolution, and the cost is stated in the results.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

#: k-mer length. 3 gives 8,000 possible features over a 20-letter alphabet, which is
#: informative for sequences of a few hundred residues without being mostly zeros.
KMER_SIZE = 3


@dataclass
class DistanceMatrix:
    """A labelled symmetric distance matrix."""

    names: list[str]
    matrix: list[list[float]]

    def __len__(self) -> int:
        return len(self.names)


def kmer_profile(sequence: str, k: int = KMER_SIZE, *, skip: str = "X") -> Counter[str]:
    """Normalised *k*-mer counts, skipping windows containing a masked residue.

    A masked residue is not a letter -- it is the absence of trustworthy information -- so
    any window touching one is dropped rather than being counted as a distinct k-mer. The
    alternative silently turns low-confidence regions into their own signal, which is the
    artefact pLDDT masking exists to remove.
    """
    counts: Counter[str] = Counter()
    for i in range(len(sequence) - k + 1):
        window = sequence[i : i + k]
        if any(character in skip for character in window):
            continue
        counts[window] += 1
    total = sum(counts.values())
    if total == 0:
        return counts
    return Counter({kmer: count / total for kmer, count in counts.items()})


def cosine_distance(a: Counter[str], b: Counter[str]) -> float:
    """1 - cosine similarity between two k-mer profiles."""
    if not a or not b:
        return 1.0
    shared = set(a) & set(b)
    dot = sum(a[kmer] * b[kmer] for kmer in shared)
    norm_a = sum(value * value for value in a.values()) ** 0.5
    norm_b = sum(value * value for value in b.values()) ** 0.5
    if norm_a == 0 or norm_b == 0:
        return 1.0
    return max(0.0, 1.0 - dot / (norm_a * norm_b))


def distance_matrix(names: Sequence[str], sequences: Sequence[str]) -> DistanceMatrix:
    """Pairwise k-mer cosine distances."""
    if len(names) != len(sequences):
        raise ValueError("names and sequences must be the same length")
    profiles = [kmer_profile(sequence) for sequence in sequences]
    size = len(names)
    matrix = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1, size):
            distance = cosine_distance(profiles[i], profiles[j])
            matrix[i][j] = matrix[j][i] = distance
    return DistanceMatrix(names=list(names), matrix=matrix)


def neighbour_joining(distances: DistanceMatrix):
    """A neighbour-joining tree, via Biopython."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix as BioMatrix
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

    # Biopython wants the lower triangle including the diagonal.
    lower = [
        [distances.matrix[i][j] for j in range(i + 1)] for i in range(len(distances.names))
    ]
    return DistanceTreeConstructor().nj(BioMatrix(names=list(distances.names), matrix=lower))


def splits(tree) -> set[frozenset[str]]:
    """The set of non-trivial bipartitions a tree induces.

    Two unrooted trees have the same topology exactly when their split sets match, which is
    what makes this the basis of the Robinson-Foulds comparison.
    """
    leaves = {leaf.name for leaf in tree.get_terminals()}
    found: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        subset = frozenset(leaf.name for leaf in clade.get_terminals())
        # Trivial splits (everything, or a single leaf) carry no topological information.
        if 1 < len(subset) < len(leaves) - 1:
            found.add(subset)
    return found


def robinson_foulds(tree_a, tree_b) -> tuple[int, float]:
    """Robinson-Foulds distance, and the same normalised to [0, 1].

    Normalisation is by the total number of non-trivial splits in both trees, so a value
    of 1.0 means the two trees share no grouping at all.
    """
    splits_a, splits_b = splits(tree_a), splits(tree_b)
    symmetric_difference = len(splits_a ^ splits_b)
    total = len(splits_a) + len(splits_b)
    return symmetric_difference, (symmetric_difference / total if total else 0.0)


def cherries(tree) -> list[tuple[str, str]]:
    """Sister leaf pairs.

    A cherry is the one place in a tree where two extant architectures can be compared
    without reconstructing an ancestor. Restricting the event comparison to cherries keeps
    it free of an ancestral-state model whose assumptions would otherwise be doing part of
    the work.
    """
    pairs: list[tuple[str, str]] = []
    for clade in tree.get_nonterminals():
        children = clade.clades
        terminals = [child for child in children if child.is_terminal()]
        if len(terminals) == 2 and len(children) == 2:
            names = sorted(child.name for child in terminals)
            pairs.append((names[0], names[1]))
    return pairs
