"""Domain architectures, and the rearrangement events that separate two of them.

An architecture is an ordered tuple of domain accessions, N-terminus first. Order matters:
the same domains in a different order is a rearrangement, not the same protein.

The event vocabulary follows the standard modular-evolution literature. Most architecture
change happens by terminal addition or deletion; repeat families expand mainly by internal
duplication. Classifying an observed pair into one of these is the unit of analysis for the
whole repo, so it lives here, in pure Python, and is tested exhaustively.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum

Architecture = tuple[str, ...]


class Event(Enum):
    """The minimal event explaining a change from one architecture to another."""

    IDENTITY = "identity"
    TERMINAL_ADDITION = "terminal_addition"
    TERMINAL_DELETION = "terminal_deletion"
    INTERNAL_DUPLICATION = "internal_duplication"
    INTERNAL_INSERTION = "internal_insertion"
    INTERNAL_DELETION = "internal_deletion"
    COMPLEX = "complex"


def classify(before: Architecture, after: Architecture) -> Event:
    """Classify the single most parsimonious event turning ``before`` into ``after``.

    Anything needing more than one event is ``COMPLEX``: it is reported rather than
    decomposed, because silently decomposing multi-step change into a chain of plausible
    single steps is how parsimony reconstructions manufacture the pattern they set out to
    find.
    """
    if before == after:
        return Event.IDENTITY

    if len(after) > len(before):
        if _is_prefix(before, after) or _is_suffix(before, after):
            return Event.TERMINAL_ADDITION
        added = _single_insertion(before, after)
        if added is not None:
            return Event.INTERNAL_DUPLICATION if added in before else Event.INTERNAL_INSERTION
        return Event.COMPLEX

    if len(after) < len(before):
        if _is_prefix(after, before) or _is_suffix(after, before):
            return Event.TERMINAL_DELETION
        if _single_insertion(after, before) is not None:
            return Event.INTERNAL_DELETION
        return Event.COMPLEX

    return Event.COMPLEX


def repeat_expansion(architecture: Architecture) -> dict[str, int]:
    """Count consecutive repeats per domain, which is how repeat families actually grow."""
    counts: dict[str, int] = {}
    run_domain: str | None = None
    run_length = 0
    for domain in architecture:
        if domain == run_domain:
            run_length += 1
        else:
            if run_domain is not None and run_length > 1:
                counts[run_domain] = max(counts.get(run_domain, 0), run_length)
            run_domain, run_length = domain, 1
    if run_domain is not None and run_length > 1:
        counts[run_domain] = max(counts.get(run_domain, 0), run_length)
    return counts


def parse(architecture: str, *, separator: str = "-") -> Architecture:
    """Parse ``PF00083-PF07690`` into a tuple, rejecting empty components."""
    parts = tuple(part.strip() for part in architecture.split(separator))
    if any(not part for part in parts):
        raise ValueError(f"empty domain in architecture {architecture!r}")
    return parts


def _is_prefix(short: Sequence[str], long: Sequence[str]) -> bool:
    return tuple(long[: len(short)]) == tuple(short)


def _is_suffix(short: Sequence[str], long: Sequence[str]) -> bool:
    return tuple(long[len(long) - len(short) :]) == tuple(short)


def _single_insertion(short: Architecture, long: Architecture) -> str | None:
    """Return the inserted domain if ``long`` is ``short`` with exactly one domain added."""
    if len(long) != len(short) + 1:
        return None
    for position in range(len(long)):
        if long[:position] + long[position + 1 :] == short:
            return long[position]
    return None
