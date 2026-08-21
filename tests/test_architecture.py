"""Event classification and masking invariants.

Rearrangement counts are the unit of analysis for this repo, so a misclassification here
would not crash anything -- it would quietly change the conclusion.
"""

from __future__ import annotations

import pytest

from domarch import structure
from domarch.architecture import Event, classify, parse, repeat_expansion


def test_identical_architectures_are_not_an_event() -> None:
    assert classify(("A", "B"), ("A", "B")) is Event.IDENTITY


@pytest.mark.parametrize(
    ("before", "after"),
    [(("A", "B"), ("A", "B", "C")), (("A", "B"), ("C", "A", "B"))],
)
def test_terminal_additions_at_either_end(
    before: tuple[str, ...], after: tuple[str, ...]
) -> None:
    """Most architecture change is terminal, at N- or C-terminus alike."""
    assert classify(before, after) is Event.TERMINAL_ADDITION


def test_terminal_deletion_is_the_mirror_of_addition() -> None:
    assert classify(("A", "B", "C"), ("A", "B")) is Event.TERMINAL_DELETION


def test_internal_duplication_is_distinguished_from_insertion() -> None:
    """Repeat families expand internally; that must not be scored as a novel domain."""
    assert classify(("A", "B", "C"), ("A", "B", "B", "C")) is Event.INTERNAL_DUPLICATION
    assert classify(("A", "B", "C"), ("A", "B", "Z", "C")) is Event.INTERNAL_INSERTION


def test_reordering_is_complex_not_silently_decomposed() -> None:
    """Multi-step change is reported, never invented as a chain of plausible single steps."""
    assert classify(("A", "B", "C"), ("C", "B", "A")) is Event.COMPLEX


def test_repeat_expansion_counts_consecutive_runs() -> None:
    assert repeat_expansion(("A", "B", "B", "B", "C", "B")) == {"B": 3}


def test_parse_rejects_empty_domains() -> None:
    assert parse("PF00083-PF07690") == ("PF00083", "PF07690")
    with pytest.raises(ValueError):
        parse("PF00083--PF07690")


def test_low_confidence_residues_are_masked_not_trusted() -> None:
    """The central guard: a guessed linker must never be scored as a structural match."""
    masked = structure.mask_by_plddt("ABCDE", [95.0, 30.0, 88.0, 12.0, 71.0])
    assert masked == "AXCXE"


def test_mask_character_is_outside_the_3di_alphabet() -> None:
    """If the mask collided with a real state, masked positions would score as matches."""
    assert structure.MASK_CHARACTER not in structure.THREE_DI_ALPHABET
    assert len(structure.THREE_DI_ALPHABET) == 20


def test_length_mismatch_between_structure_and_confidence_is_fatal() -> None:
    with pytest.raises(ValueError):
        structure.mask_by_plddt("ABC", [90.0, 90.0])


def test_confident_fraction_is_reported_per_structure() -> None:
    assert structure.confident_fraction([90.0, 30.0, 80.0, 20.0]) == 0.5
