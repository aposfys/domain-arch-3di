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
    """The set of non-trivial bipartitions a tree induces, read as an unrooted tree.

    Two unrooted trees have the same topology exactly when their split sets match, which is
    what makes this the basis of the Robinson-Foulds comparison. Neighbour joining places
    the root arbitrarily, so each split is stored as the side that excludes a fixed
    reference leaf. Without that, the same bipartition read from differently rooted copies
    of one tree could be stored as two different subsets.
    """
    leaves = frozenset(leaf.name for leaf in tree.get_terminals())
    if not leaves:
        return set()
    reference = min(leaves)
    found: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        subset = frozenset(leaf.name for leaf in clade.get_terminals())
        side = leaves - subset if reference in subset else subset
        # Trivial splits (a single leaf against the rest) carry no topological information.
        if 1 < len(side) < len(leaves) - 1:
            found.add(side)
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
    """Sister leaf pairs of the tree, read as an unrooted tree.

    A cherry is the one place in a tree where two extant architectures can be compared
    without reconstructing an ancestor. Restricting the event comparison to cherries keeps
    it free of an ancestral-state model whose assumptions would otherwise be doing part of
    the work.

    A cherry is an internal node with exactly two leaf neighbours. Neighbour joining
    returns a trifurcating root, so the root is treated like any other node, and a
    bifurcating root is suppressed by joining its two children.
    """
    root = tree.root
    neighbours: dict[int, list] = {}

    def link(a, b) -> None:
        neighbours.setdefault(id(a), []).append(b)
        neighbours.setdefault(id(b), []).append(a)

    root_children = list(root.clades)
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            if clade is root and len(root_children) == 2:
                continue
            link(clade, child)
    if len(root_children) == 2:
        link(root_children[0], root_children[1])

    pairs: list[tuple[str, str]] = []
    for clade in tree.get_nonterminals():
        if clade is root and len(root_children) == 2:
            continue
        adjacent = neighbours.get(id(clade), [])
        leaves = sorted(node.name for node in adjacent if node.is_terminal())
        # Three leaf neighbours only happens in a three-leaf tree, where no pair is special.
        if len(leaves) == 2:
            pairs.append((leaves[0], leaves[1]))
    return pairs
