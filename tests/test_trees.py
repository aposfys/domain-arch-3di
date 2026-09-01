"""Distance, topology and cherry extraction, none of which need a structure or a network."""

from __future__ import annotations

import pytest

from domarch import trees
from domarch.compare import compare_events


def test_kmer_profile_skips_windows_touching_a_mask():
    """A masked residue is absent information, not a 21st letter."""
    profile = trees.kmer_profile("AAXAA", k=3)
    # Only AAX, AXA, XAA exist as windows and all three touch the mask.
    assert profile == {}
    assert trees.kmer_profile("AAAA", k=3) == {"AAA": 1.0}


def test_identical_sequences_are_distance_zero():
    a = trees.kmer_profile("ACDEFGHIKLMNPQ")
    assert trees.cosine_distance(a, a) == pytest.approx(0.0)


def test_disjoint_sequences_are_distance_one():
    a = trees.kmer_profile("AAAAAAAA")
    b = trees.kmer_profile("CCCCCCCC")
    assert trees.cosine_distance(a, b) == pytest.approx(1.0)


def test_distance_matrix_is_symmetric_with_zero_diagonal():
    names = ["a", "b", "c"]
    matrix = trees.distance_matrix(names, ["AAAACCCC", "AAAAGGGG", "CCCCGGGG"])
    for i in range(3):
        assert matrix.matrix[i][i] == 0.0
        for j in range(3):
            assert matrix.matrix[i][j] == pytest.approx(matrix.matrix[j][i])


def test_distance_matrix_refuses_mismatched_inputs():
    with pytest.raises(ValueError, match="same length"):
        trees.distance_matrix(["a", "b"], ["AAAA"])


def _tree(newick: str):
    import io

    from Bio import Phylo

    return Phylo.read(io.StringIO(newick), "newick")


def test_robinson_foulds_is_zero_for_a_tree_against_itself():
    tree = _tree("(((a,b),(c,d)),(e,f));")
    assert trees.robinson_foulds(tree, tree) == (0, 0.0)


def test_robinson_foulds_detects_a_regrouping():
    a = _tree("(((a,b),(c,d)),(e,f));")
    b = _tree("(((a,c),(b,d)),(e,f));")
    absolute, normalised = trees.robinson_foulds(a, b)
    assert absolute > 0
    assert 0.0 < normalised <= 1.0


def test_cherries_finds_sister_leaf_pairs():
    tree = _tree("(((a,b),(c,d)),(e,f));")
    assert sorted(trees.cherries(tree)) == [("a", "b"), ("c", "d"), ("e", "f")]


def test_cherries_ignores_a_node_with_a_non_leaf_child():
    tree = _tree("((a,(b,c)),d);")
    assert trees.cherries(tree) == [("b", "c")]


def test_the_same_pair_classifies_the_same_way_in_either_tree():
    """If this ever fails, the event classifier depends on more than the architectures."""
    architectures = {"a": "PF1-PF2", "b": "PF1-PF2-PF3", "c": "PF1", "d": "PF1"}
    comparison = compare_events(
        [("a", "b"), ("c", "d")], [("a", "b"), ("c", "d")], architectures
    )
    assert comparison.conflicting_shared_cherries == 0
    assert comparison.shared_cherries == 2
    assert comparison.cherry_jaccard == pytest.approx(1.0)


def test_disjoint_cherry_sets_share_nothing():
    architectures = {"a": "PF1", "b": "PF1", "c": "PF2", "d": "PF2"}
    comparison = compare_events([("a", "b")], [("c", "d")], architectures)
    assert comparison.shared_cherries == 0
    assert comparison.cherry_jaccard == pytest.approx(0.0)


def test_cherries_over_absent_architectures_are_skipped_not_guessed():
    comparison = compare_events([("a", "zz")], [], {"a": "PF1"})
    assert comparison.n_cherries_sequence == 0
