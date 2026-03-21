"""Tests for ALEobserve -- comparing Python output with C++ output."""

import os
import shutil
import subprocess

import pytest

from ALE_python.ale import ApproxPosterior
from ALE_python.ale_util import observe_ALE_from_file, load_ALE_from_file
from ALE_python.newick import parse_newick, get_leaves


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_cpp_observe(cpp_bin, gene_trees_path, tmp_dir, burnin=1000):
    """Run the C++ ALEobserve and return the path to the .ale file."""
    gene_trees_copy = os.path.join(tmp_dir, os.path.basename(gene_trees_path))
    shutil.copy(gene_trees_path, gene_trees_copy)
    exe = os.path.join(cpp_bin, "ALEobserve")
    result = subprocess.run(
        [exe, gene_trees_copy, f"burnin={burnin}"],
        capture_output=True,
        text=True,
        cwd=tmp_dir,
        timeout=300,
    )
    assert result.returncode == 0, f"C++ ALEobserve failed:\n{result.stderr}"
    ale_path = gene_trees_copy + ".ale"
    assert os.path.isfile(ale_path), f"Expected .ale file not created: {ale_path}"
    return ale_path


def _load_cpp_ale(cpp_bin, gene_trees_path, tmp_dir, burnin=1000):
    """Run C++ ALEobserve and load the resulting .ale file."""
    ale_path = _run_cpp_observe(cpp_bin, gene_trees_path, tmp_dir, burnin=burnin)
    return load_ALE_from_file(ale_path)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestObserve:

    def test_observe_creates_ale_file(self, gene_trees, tmp_dir):
        """Run Python ALEobserve on example data with burnin=1000 and verify
        that an .ale file is created."""
        gene_trees_copy = os.path.join(tmp_dir, os.path.basename(gene_trees))
        shutil.copy(gene_trees, gene_trees_copy)

        ale = observe_ALE_from_file(gene_trees_copy, burnin=1000)
        out_path = gene_trees_copy + ".ale"
        ale.save_state(out_path)

        assert os.path.isfile(out_path), "Python ALEobserve did not create .ale file"
        assert os.path.getsize(out_path) > 0, ".ale file is empty"

    def test_observe_observation_count(self, gene_trees):
        """Load the .ale and verify observations == 8504."""
        ale = observe_ALE_from_file(gene_trees, burnin=1000)
        assert ale.observations == pytest.approx(8504.0), (
            f"Expected 8504 observations, got {ale.observations}"
        )

    def test_observe_matches_cpp_leaf_count(self, gene_trees, cpp_bin, tmp_dir):
        """Compare number of leaves (36) in both Python and C++ .ale files."""
        py_ale = observe_ALE_from_file(gene_trees, burnin=1000)
        cpp_ale = _load_cpp_ale(cpp_bin, gene_trees, tmp_dir, burnin=1000)

        assert py_ale.Gamma_size == 36, (
            f"Python Gamma_size={py_ale.Gamma_size}, expected 36"
        )
        assert cpp_ale.Gamma_size == 36, (
            f"C++ Gamma_size={cpp_ale.Gamma_size}, expected 36"
        )
        assert py_ale.Gamma_size == cpp_ale.Gamma_size

    def test_observe_matches_cpp_bip_counts(self, gene_trees, cpp_bin, tmp_dir):
        """Load both .ale files and compare that Bip_counts are similar.

        Python and C++ may assign different integer IDs to the same
        bipartitions, so we compare by bitset content rather than by ID.

        Leaf bipartitions (size 1 or Gamma_size-1) are excluded because
        the C++ code does not store their counts in Bip_counts during
        decompose; instead it special-cases them in p_bip().
        """
        py_ale = observe_ALE_from_file(gene_trees, burnin=1000)
        cpp_ale = _load_cpp_ale(cpp_bin, gene_trees, tmp_dir, burnin=1000)

        def _is_leaf_size(size, gamma_size):
            return size <= 1 or size >= gamma_size - 1

        def _bitset_to_leafset(ale, bitset):
            """Convert a bitset integer to a frozenset of leaf names."""
            names = set()
            for name, lid in ale.leaf_ids.items():
                if (bitset >> lid) & 1:
                    names.add(name)
            return frozenset(names)

        # Build leaf-name-set -> count maps for non-leaf bipartitions
        py_by_leaves = {}
        for kid, count in py_ale.Bip_counts.items():
            size = py_ale.set_sizes.get(kid, 0)
            if _is_leaf_size(size, py_ale.Gamma_size):
                continue
            leafset = _bitset_to_leafset(py_ale, py_ale.id_sets[kid])
            py_by_leaves[leafset] = count

        cpp_by_leaves = {}
        for kid, count in cpp_ale.Bip_counts.items():
            size = cpp_ale.set_sizes.get(kid, 0)
            if _is_leaf_size(size, cpp_ale.Gamma_size):
                continue
            leafset = _bitset_to_leafset(cpp_ale, cpp_ale.id_sets[kid])
            cpp_by_leaves[leafset] = count

        assert len(py_by_leaves) == len(cpp_by_leaves), (
            f"Non-leaf bipartition count mismatch: {len(py_by_leaves)} Python vs {len(cpp_by_leaves)} C++"
        )

        py_only = set(py_by_leaves.keys()) - set(cpp_by_leaves.keys())
        cpp_only = set(cpp_by_leaves.keys()) - set(py_by_leaves.keys())
        assert not py_only and not cpp_only, (
            f"Bipartition mismatch: {len(py_only)} only in Python, {len(cpp_only)} only in C++"
        )

        for leafset in py_by_leaves:
            py_val = py_by_leaves[leafset]
            cpp_val = cpp_by_leaves[leafset]
            assert py_val == pytest.approx(cpp_val, rel=1e-6), (
                f"Bip_counts for {leafset}: Python={py_val}, C++={cpp_val}"
            )

    def test_observe_mpp_tree_has_all_leaves(self, gene_trees):
        """Check the mpp tree string contains all 36 leaf names."""
        ale = observe_ALE_from_file(gene_trees, burnin=1000)
        mpp_str, mpp_pp = ale.mpp_tree()

        assert len(mpp_str) > 0, "mpp_tree returned empty string"
        assert mpp_pp > 0.0, f"mpp posterior probability should be > 0, got {mpp_pp}"

        # Parse the mpp tree and check leaf count
        root = parse_newick(mpp_str)
        leaves = get_leaves(root)
        leaf_names = sorted(n.name for n in leaves)

        assert len(leaf_names) == 36, (
            f"Expected 36 leaves in mpp tree, got {len(leaf_names)}"
        )

        # Verify all expected leaf names are present
        expected_leaves = sorted(ale.leaf_ids.keys())
        assert leaf_names == expected_leaves, (
            f"Leaf name mismatch between mpp tree and ale object"
        )

    def test_observe_singleton_tree_file(self, tmp_dir):
        """Singleton input should produce a stable one-leaf ALE."""
        source_dir = os.path.join(tmp_dir, "src")
        os.makedirs(source_dir, exist_ok=True)
        singleton_path = os.path.join(source_dir, "singleton.trees")
        with open(singleton_path, "w") as fh:
            fh.write("A;\n")

        py_ale = observe_ALE_from_file(singleton_path)

        assert py_ale.observations == pytest.approx(1.0)
        assert py_ale.Gamma_size == 1
        assert py_ale.get_leaf_names() == ["A"]
        assert py_ale.Bip_bls[1] == pytest.approx(1.0)

        roundtrip_path = os.path.join(tmp_dir, "singleton.ale")
        py_ale.save_state(roundtrip_path)
        loaded = load_ALE_from_file(roundtrip_path)
        assert loaded.observations == pytest.approx(py_ale.observations)
        assert loaded.Gamma_size == py_ale.Gamma_size
        assert loaded.get_leaf_names() == py_ale.get_leaf_names()

    def test_observe_singleton_then_tree_returns_singleton(self, tmp_dir):
        """If a singleton appears first, observe should return that one-leaf ALE."""
        source_dir = os.path.join(tmp_dir, "src")
        os.makedirs(source_dir, exist_ok=True)
        mixed_path = os.path.join(source_dir, "singleton_then_tree.trees")
        with open(mixed_path, "w") as fh:
            fh.write("A;\n")
            fh.write("(A,B);\n")

        py_ale = observe_ALE_from_file(mixed_path)

        assert py_ale.observations == pytest.approx(1.0)
        assert py_ale.Gamma_size == 1
        assert py_ale.get_leaf_names() == ["A"]
