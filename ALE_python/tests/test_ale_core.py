"""Tests for core ALE data structures (ApproxPosterior, newick, etc.)."""

import os
import tempfile

import pytest

from ALE_python.ale import ApproxPosterior
from ALE_python.ale_util import (
    observe_ALE_from_strings,
    load_ALE_from_file,
    save_ALE_to_file,
)
from ALE_python.newick import parse_newick, get_leaves


class TestConstructSimpleTree:

    def test_construct_simple_tree(self):
        """Construct from '(A:1,B:1,C:1);', verify Gamma_size=3 and leaf_ids."""
        ale = ApproxPosterior("(A:1,B:1,C:1);")

        assert ale.Gamma_size == 3, f"Expected Gamma_size=3, got {ale.Gamma_size}"
        assert len(ale.leaf_ids) == 3, (
            f"Expected 3 leaf_ids, got {len(ale.leaf_ids)}"
        )
        assert set(ale.leaf_ids.keys()) == {"A", "B", "C"}

        # IDs should be 1-based and assigned alphabetically
        assert ale.leaf_ids["A"] == 1
        assert ale.leaf_ids["B"] == 2
        assert ale.leaf_ids["C"] == 3

        # Reverse mapping
        assert ale.id_leaves[1] == "A"
        assert ale.id_leaves[2] == "B"
        assert ale.id_leaves[3] == "C"


class TestDecomposeAndObserve:

    def test_decompose_and_observe(self):
        """Observe multiple trees, verify Bip_counts are populated."""
        trees = [
            "((A:1,B:1):1,C:1);",
            "((A:1,C:1):1,B:1);",
            "((A:1,B:1):1,C:1);",
            "((B:1,C:1):1,A:1);",
        ]

        ale = observe_ALE_from_strings(trees)

        assert ale.observations == pytest.approx(4.0), (
            f"Expected 4 observations, got {ale.observations}"
        )
        assert ale.Gamma_size == 3

        # With 3 leaves, all bipartitions are trivial (size 1 or 2),
        # so Bip_counts should be non-empty
        assert len(ale.Bip_counts) > 0, "Bip_counts should not be empty"


class TestSaveLoadRoundtrip:

    def test_save_load_roundtrip(self):
        """Save ALE to file, load it back, compare all fields."""
        trees = [
            "((A:0.5,B:0.3):0.2,(C:0.4,D:0.6):0.1);",
            "((A:0.5,C:0.3):0.2,(B:0.4,D:0.6):0.1);",
            "((A:0.5,B:0.3):0.2,(C:0.4,D:0.6):0.1);",
            "(((A:0.5,B:0.3):0.1,C:0.4):0.2,D:0.6);",
        ]
        original = observe_ALE_from_strings(trees)

        with tempfile.NamedTemporaryFile(suffix=".ale", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            save_ALE_to_file(original, tmp_path)
            loaded = load_ALE_from_file(tmp_path)

            # Compare observations
            assert loaded.observations == pytest.approx(original.observations), (
                f"observations mismatch: {loaded.observations} vs {original.observations}"
            )

            # Compare Gamma_size
            assert loaded.Gamma_size == original.Gamma_size

            # Compare leaf_ids
            assert loaded.leaf_ids == original.leaf_ids

            # Compare Bip_counts
            assert set(loaded.Bip_counts.keys()) == set(original.Bip_counts.keys())
            for key in original.Bip_counts:
                assert loaded.Bip_counts[key] == pytest.approx(
                    original.Bip_counts[key], rel=1e-6
                ), f"Bip_counts[{key}] mismatch"

            # Compare non-empty Dip_counts entries (empty trailing entries
            # may differ because decompose pushes empty dicts that are not
            # serialised, so load_state pads to last_leafset_id+1 only).
            orig_nonempty = {i: d for i, d in enumerate(original.Dip_counts) if d}
            load_nonempty = {i: d for i, d in enumerate(loaded.Dip_counts) if d}
            assert set(orig_nonempty.keys()) == set(load_nonempty.keys()), (
                f"Dip_counts non-empty index mismatch"
            )
            for idx in orig_nonempty:
                for key in orig_nonempty[idx]:
                    assert load_nonempty[idx][key] == pytest.approx(
                        orig_nonempty[idx][key], rel=1e-6
                    )

            # Compare set_ids mapping
            for bitset, sid in original.set_ids.items():
                assert bitset in loaded.set_ids, (
                    f"Bitset {bitset} missing from loaded set_ids"
                )
                assert loaded.set_ids[bitset] == sid

        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)


class TestPBip:

    def test_p_bip_simple(self):
        """Test bipartition probability calculation on a small example."""
        # With fewer than 4 leaves, p_bip should always return 1.0
        trees_3 = ["((A:1,B:1):1,C:1);"] * 5
        ale_3 = observe_ALE_from_strings(trees_3)

        for g_id in ale_3.Bip_counts:
            assert ale_3.p_bip(g_id) == pytest.approx(1.0), (
                f"p_bip({g_id}) should be 1.0 for 3 leaves"
            )

        # With 4+ leaves, probabilities depend on observation frequencies
        trees_4 = [
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,C:1):1,(B:1,D:1):1);",
        ]
        ale_4 = observe_ALE_from_strings(trees_4)

        # Probabilities should be between 0 and 1
        for g_id in ale_4.Bip_counts:
            p = ale_4.p_bip(g_id)
            assert 0.0 <= p <= 1.0, f"p_bip({g_id})={p} out of range"

        # Leaf-level bipartitions (size 1) should have probability 1.0
        # because they appear in every tree
        for g_id, bitset in ale_4.id_sets.items():
            size = ale_4.set_sizes.get(g_id, 0)
            if size == 1:
                assert ale_4.p_bip(g_id) == pytest.approx(1.0), (
                    f"Leaf bipartition p_bip({g_id}) should be 1.0"
                )


class TestMppTree:

    def test_mpp_tree_simple(self):
        """Test mpp tree on a small example with a clear majority topology."""
        # Give one topology a strong majority
        trees = (
            ["((A:1,B:1):1,(C:1,D:1):1);"] * 8
            + ["((A:1,C:1):1,(B:1,D:1):1);"] * 2
        )
        ale = observe_ALE_from_strings(trees)
        mpp_str, mpp_pp = ale.mpp_tree()

        assert len(mpp_str) > 0, "mpp_tree returned empty string"
        assert mpp_pp > 0.0, f"mpp posterior probability should be > 0"

        # Parse and verify leaf count
        root = parse_newick(mpp_str)
        leaves = get_leaves(root)
        leaf_names = sorted(n.name for n in leaves)

        assert len(leaf_names) == 4
        assert set(leaf_names) == {"A", "B", "C", "D"}


class TestCountTrees:

    def test_count_trees(self):
        """Test count_trees function returns reasonable values."""
        # With 3 leaves, there are exactly 3 unrooted topologies,
        # but amalgamated count depends on observed bipartitions
        trees_3 = [
            "((A:1,B:1):1,C:1);",
            "((A:1,C:1):1,B:1);",
            "((B:1,C:1):1,A:1);",
        ]
        ale_3 = observe_ALE_from_strings(trees_3)
        count_3 = ale_3.count_trees()

        # count_trees should return a positive number
        assert count_3 > 0, f"count_trees returned {count_3}, expected > 0"

        # For 4 leaves with all topologies observed
        trees_4 = [
            "((A:1,B:1):1,(C:1,D:1):1);",
            "((A:1,C:1):1,(B:1,D:1):1);",
            "((A:1,D:1):1,(B:1,C:1):1);",
        ]
        ale_4 = observe_ALE_from_strings(trees_4)
        count_4 = ale_4.count_trees()

        assert count_4 > 0, f"count_trees returned {count_4}, expected > 0"
        # With 4 leaves, count_trees counts amalgamable rooted trees (not
        # unrooted topologies).  Each of the 3 unrooted topologies has 5
        # rootings, giving 15 rooted trees total.
        assert count_4 == pytest.approx(15.0), (
            f"count_trees for 4 leaves (all topologies observed) should be 15, got {count_4}"
        )
