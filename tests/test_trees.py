"""Distance, topology and cherry extraction, none of which need a structure or a network."""

from __future__ import annotations

import pytest

from domarch import trees
from domarch.compare import compare_events, fisher_exact_two_sided


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


def test_robinson_foulds_ignores_where_the_root_sits():
    """Neighbour joining roots arbitrarily, so a re-rooted copy is the same tree."""
    rooted_one_way = _tree("(((a,b),(c,d)),(e,f),(g,h));")
    rooted_another = _tree("((a,b),((c,d),((e,f),(g,h))));")
    rooted_at_a_leaf = _tree("(a,b,((c,d),((e,f),(g,h))));")
    assert trees.robinson_foulds(rooted_one_way, rooted_another) == (0, 0.0)
    assert trees.robinson_foulds(rooted_one_way, rooted_at_a_leaf) == (0, 0.0)
    assert len(trees.splits(rooted_one_way)) == 8 - 3


def test_cherries_ignores_a_node_with_a_non_leaf_child():
    tree = _tree("((a,(b,c)),(d,e));")
    assert sorted(trees.cherries(tree)) == [("b", "c"), ("d", "e")]


def test_cherries_include_a_pair_on_a_trifurcating_root():
    """Biopython's neighbour joining leaves three children on the root."""
    tree = _tree("(a,b,((c,d),e));")
    assert sorted(trees.cherries(tree)) == [("a", "b"), ("c", "d")]


def test_cherries_do_not_depend_on_where_the_root_sits():
    one = _tree("(((a,b),(c,d)),(e,f),(g,h));")
    other = _tree("(a,b,((c,d),((e,f),(g,h))));")
    assert sorted(trees.cherries(one)) == sorted(trees.cherries(other))


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


def test_proteins_without_domains_are_compared_not_rejected():
    architectures = {"a": "", "b": "", "c": "PF1", "d": ""}
    comparison = compare_events([("a", "b"), ("c", "d")], [], architectures)
    assert comparison.events_sequence == {"IDENTITY": 1, "TERMINAL_INDEL": 1}


def test_a_cherry_is_the_same_pair_in_either_order():
    architectures = {"a": "PF1", "b": "PF1-PF2"}
    comparison = compare_events([("a", "b")], [("b", "a")], architectures)
    assert comparison.shared_cherries == 1
    assert comparison.events_sequence == comparison.events_structure


def test_fisher_exact_matches_known_values():
    # Fisher's tea-tasting table.
    assert fisher_exact_two_sided(3, 1, 1, 3) == pytest.approx(0.4857, abs=1e-4)
    # The committed run, events against no change per tree.
    assert fisher_exact_two_sided(11, 10, 5, 13) == pytest.approx(0.1923, abs=1e-4)
    assert fisher_exact_two_sided(0, 0, 0, 0) == 1.0
    with pytest.raises(ValueError):
        fisher_exact_two_sided(-1, 0, 0, 0)


def test_cherries_over_absent_architectures_are_skipped_not_guessed():
    comparison = compare_events([("a", "zz")], [], {"a": "PF1"})
    assert comparison.n_cherries_sequence == 0
