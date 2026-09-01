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

from domarch.architecture import Event, classify, parse


@dataclass
class EventComparison:
    """Events inferred from each tree, and how much they overlap."""

    n_cherries_sequence: int
    n_cherries_structure: int
    shared_cherries: int
    #: Event counts from the cherries each tree found.
    events_sequence: dict[str, int] = field(default_factory=dict)
    events_structure: dict[str, int] = field(default_factory=dict)
    #: Cherries present in both trees where the two disagree on the event. Always zero --
    #: the same pair of architectures classifies the same way whatever tree found it --
    #: and reported because a non-zero value would mean the classifier is not a function
    #: of the architectures alone.
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
) -> dict[tuple[str, str], Event]:
    """Classify the architecture change across each sister pair."""
    classified: dict[tuple[str, str], Event] = {}
    for left, right in pairs:
        if left not in architectures or right not in architectures:
            continue
        before = parse(architectures[left])
        after = parse(architectures[right])
        classified[(left, right)] = classify(before, after)
    return classified


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
