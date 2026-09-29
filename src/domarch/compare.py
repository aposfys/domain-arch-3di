"""Event-level comparison between the sequence tree and the structure tree.

The question this repository asks is not whether the two trees differ -- they will -- but
whether the *rearrangement events* differ. A tree can be reshuffled substantially while
every inferred gain, loss and duplication stays the same, and that outcome would mean 3Di
changes the phylogeny without changing the evolutionary story. The opposite outcome is the
one that would matter.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from math import comb

from domarch.architecture import PairEvent, classify_pair, parse


@dataclass
class EventComparison:
    """Events inferred from each tree, and how much they overlap."""

    n_cherries_sequence: int
    n_cherries_structure: int
    shared_cherries: int
    #: Event counts from the cherries each tree found.
    events_sequence: dict[str, int] = field(default_factory=dict)
    events_structure: dict[str, int] = field(default_factory=dict)
    #: Cherries present in both trees where the two disagree on the event. Zero by
    #: construction, because the pair classifier is a function of the two architectures
    #: only. Kept as a self-test of the code, not as a finding.
    conflicting_shared_cherries: int = 0
    #: Events inferred from one tree's cherries but not the other's.
    events_only_in_sequence: dict[str, int] = field(default_factory=dict)
    events_only_in_structure: dict[str, int] = field(default_factory=dict)

    @property
    def cherry_jaccard(self) -> float:
        union = self.n_cherries_sequence + self.n_cherries_structure - self.shared_cherries
        return self.shared_cherries / union if union else 1.0


def events_from_cherries(
    pairs: list[tuple[str, str]], architectures: dict[str, str]
) -> dict[tuple[str, str], PairEvent]:
    """Classify the architecture change across each sister pair, without a direction."""
    classified: dict[tuple[str, str], PairEvent] = {}
    for left, right in pairs:
        if left not in architectures or right not in architectures:
            continue
        key = (left, right) if left <= right else (right, left)
        classified[key] = classify_pair(
            parse(architectures[left]), parse(architectures[right])
        )
    return classified


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for the 2x2 table ``[[a, b], [c, d]]``.

    Tables are summed when they are no more probable than the observed one, the same
    convention as ``scipy.stats.fisher_exact``.
    """
    if min(a, b, c, d) < 0:
        raise ValueError("counts must be non-negative")
    row, col, n = a + b, a + c, a + b + c + d
    if n == 0:
        return 1.0

    def probability(x: int) -> float:
        return comb(row, x) * comb(n - row, col - x) / comb(n, col)

    observed = probability(a)
    low, high = max(0, col - (n - row)), min(row, col)
    total = sum(
        p for p in (probability(x) for x in range(low, high + 1)) if p <= observed * (1 + 1e-7)
    )
    return min(1.0, total)


def compare_events(
    sequence_cherries: list[tuple[str, str]],
    structure_cherries: list[tuple[str, str]],
    architectures: dict[str, str],
) -> EventComparison:
    """Compare the events each tree's cherries imply."""
    from_sequence = events_from_cherries(sequence_cherries, architectures)
    from_structure = events_from_cherries(structure_cherries, architectures)

    shared = set(from_sequence) & set(from_structure)
    conflicting = sum(1 for pair in shared if from_sequence[pair] != from_structure[pair])

    counts_sequence = Counter(event.name for event in from_sequence.values())
    counts_structure = Counter(event.name for event in from_structure.values())

    only_sequence = Counter(
        from_sequence[pair].name for pair in set(from_sequence) - set(from_structure)
    )
    only_structure = Counter(
        from_structure[pair].name for pair in set(from_structure) - set(from_sequence)
    )

    return EventComparison(
        n_cherries_sequence=len(from_sequence),
        n_cherries_structure=len(from_structure),
        shared_cherries=len(shared),
        events_sequence=dict(sorted(counts_sequence.items())),
        events_structure=dict(sorted(counts_structure.items())),
        conflicting_shared_cherries=conflicting,
        events_only_in_sequence=dict(sorted(only_sequence.items())),
        events_only_in_structure=dict(sorted(only_structure.items())),
    )
