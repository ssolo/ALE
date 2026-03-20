"""Tests for ALEml_undated -- comparing Python output with C++ output."""

import math
import os
import re
import shutil
import subprocess

import pytest

from ALE_python.ale import ApproxPosterior
from ALE_python.ale_util import load_ALE_from_file, observe_ALE_from_file
from ALE_python.exodt import ExODTModel
from ALE_python.newick import parse_newick, get_leaves


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _prepare_inputs(species_tree, gene_trees, tmp_dir, burnin=1000):
    """Copy input files to tmp_dir, run Python observe, return paths."""
    s_copy = os.path.join(tmp_dir, os.path.basename(species_tree))
    g_copy = os.path.join(tmp_dir, os.path.basename(gene_trees))
    shutil.copy(species_tree, s_copy)
    shutil.copy(gene_trees, g_copy)

    ale = observe_ALE_from_file(g_copy, burnin=burnin)
    ale_path = g_copy + ".ale"
    ale.save_state(ale_path)
    return s_copy, ale_path


def _run_python_ml_undated_fixed(species_tree_path, ale_path, delta=0.01,
                                  tau=0.01, lambda_=0.1):
    """Run Python ALEml_undated with fixed rates (no optimization) and return
    the log-likelihood."""
    with open(species_tree_path) as f:
        Sstring = f.readline().strip()

    ale = load_ALE_from_file(ale_path)

    model = ExODTModel()
    model.set_model_parameter("BOOTSTRAP_LABELS", "yes")
    model.set_model_parameter("undatedBL", 0)
    model.set_model_parameter("reldate", 0)
    model.construct_undated(Sstring)

    model.set_model_parameter("seq_beta", 1.0)
    model.set_model_parameter("O_R", 1.0)
    model.set_model_parameter("delta", delta)
    model.set_model_parameter("tau", tau)
    model.set_model_parameter("lambda", lambda_)

    model.calculate_undatedEs()
    lk = model.pun(ale, False, False)
    if lk <= 0:
        return -1e50
    return math.log(lk)


def _run_cpp_ml_undated_fixed(cpp_bin, species_tree_path, ale_path, tmp_dir,
                               delta=0.01, tau=0.01, lambda_=0.1):
    """Run C++ ALEml_undated with fixed rates and extract the LL from output."""
    exe = os.path.join(cpp_bin, "ALEml_undated")
    result = subprocess.run(
        [exe, species_tree_path, ale_path,
         f"delta={delta}", f"tau={tau}", f"lambda={lambda_}",
         "sample=1"],
        capture_output=True,
        text=True,
        cwd=tmp_dir,
        timeout=600,
    )
    assert result.returncode == 0, f"C++ ALEml_undated failed:\n{result.stderr}"

    # Extract LL from stdout
    for line in result.stdout.splitlines():
        if line.startswith("LL="):
            return float(line.split("=")[1].strip())
    # Also check in the .uml_rec file
    rec_files = [f for f in os.listdir(tmp_dir) if f.endswith(".uml_rec")]
    for rec_file in rec_files:
        with open(os.path.join(tmp_dir, rec_file)) as fh:
            for line in fh:
                m = re.search(r">logl:\s*([-\d.eE+]+)", line)
                if m:
                    return float(m.group(1))

    pytest.fail("Could not extract LL from C++ ALEml_undated output")


def _run_python_ml_undated_optimized(species_tree_path, ale_path):
    """Run Python ALEml_undated with optimization and return the LL."""
    from ALE_python.cli import ALEml_undated as cli_ml_undated

    # The CLI function takes argv and writes output files.
    # We call it from the directory containing the ale file.
    cwd_save = os.getcwd()
    out_dir = os.path.dirname(ale_path)
    os.chdir(out_dir)
    try:
        cli_ml_undated([species_tree_path, ale_path, "sample=10", "seed=42"])
    finally:
        os.chdir(cwd_save)

    # Extract LL from output .uml_rec file
    basename_s = os.path.basename(species_tree_path)
    basename_ale = os.path.basename(ale_path)
    rec_name = basename_s + "_" + basename_ale + ".uml_rec"
    rec_path = os.path.join(out_dir, rec_name)
    assert os.path.isfile(rec_path), f"Expected output file not found: {rec_path}"

    with open(rec_path) as fh:
        for line in fh:
            m = re.search(r">logl:\s*([-\d.eE+]+)", line)
            if m:
                return float(m.group(1))

    pytest.fail("Could not extract LL from Python ALEml_undated output")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestMLUndated:

    def test_ml_undated_log_likelihood(self, species_tree, gene_trees,
                                        cpp_bin, tmp_dir):
        """Run Python ALEml_undated with fixed rates delta=0.01 tau=0.01
        lambda=0.1, verify LL is close to C++ LL."""
        s_copy, ale_path = _prepare_inputs(species_tree, gene_trees, tmp_dir)

        py_ll = _run_python_ml_undated_fixed(s_copy, ale_path)
        cpp_ll = _run_cpp_ml_undated_fixed(cpp_bin, s_copy, ale_path, tmp_dir)

        # Allow some tolerance for floating-point differences between
        # the Python and C++ implementations
        assert py_ll == pytest.approx(cpp_ll, abs=1.0), (
            f"Python LL={py_ll}, C++ LL={cpp_ll}"
        )

    def test_ml_undated_optimized_rates(self, species_tree, gene_trees, tmp_dir):
        """Run optimization, check LL is around -34.6 (within 1.0 tolerance)."""
        s_copy, ale_path = _prepare_inputs(species_tree, gene_trees, tmp_dir)
        ll = _run_python_ml_undated_optimized(s_copy, ale_path)

        assert ll == pytest.approx(-34.6, abs=1.0), (
            f"Optimized LL={ll}, expected around -34.6"
        )

    def test_ml_undated_output_files(self, species_tree, gene_trees, tmp_dir):
        """Check that .uml_rec and .uTs files are created with correct format."""
        s_copy, ale_path = _prepare_inputs(species_tree, gene_trees, tmp_dir)

        from ALE_python.cli import ALEml_undated as cli_ml_undated

        cwd_save = os.getcwd()
        os.chdir(tmp_dir)
        try:
            cli_ml_undated([s_copy, ale_path, "sample=10", "seed=42"])
        finally:
            os.chdir(cwd_save)

        basename_s = os.path.basename(s_copy)
        basename_ale = os.path.basename(ale_path)
        radical = basename_s + "_" + basename_ale

        rec_path = os.path.join(tmp_dir, radical + ".uml_rec")
        ts_path = os.path.join(tmp_dir, radical + ".uTs")

        assert os.path.isfile(rec_path), f".uml_rec file not found: {rec_path}"
        assert os.path.isfile(ts_path), f".uTs file not found: {ts_path}"

        # Verify .uml_rec has expected structure
        with open(rec_path) as fh:
            content = fh.read()
        assert "#ALEml_undated" in content, ".uml_rec missing header"
        assert "S:\t" in content, ".uml_rec missing species tree line"
        assert ">logl:" in content, ".uml_rec missing log-likelihood line"
        assert "rate of" in content, ".uml_rec missing rate header"
        assert "reconciled G-s:" in content, ".uml_rec missing reconciled trees section"

        # Verify .uTs has expected structure
        with open(ts_path) as fh:
            ts_content = fh.read()
        assert "#from\tto\tfreq." in ts_content, ".uTs missing header"

    def test_ml_undated_event_counts(self, species_tree, gene_trees, tmp_dir):
        """Run with sample=100, verify total events are reasonable.
        S should be 35 = n_leaves - 1 per sample."""
        s_copy, ale_path = _prepare_inputs(species_tree, gene_trees, tmp_dir)

        from ALE_python.cli import ALEml_undated as cli_ml_undated

        cwd_save = os.getcwd()
        os.chdir(tmp_dir)
        try:
            cli_ml_undated([s_copy, ale_path, "sample=100", "seed=42"])
        finally:
            os.chdir(cwd_save)

        basename_s = os.path.basename(s_copy)
        basename_ale = os.path.basename(ale_path)
        radical = basename_s + "_" + basename_ale
        rec_path = os.path.join(tmp_dir, radical + ".uml_rec")

        # Parse event counts from the .uml_rec file
        with open(rec_path) as fh:
            lines = fh.readlines()

        # Find the "Total" line after "# of\tDuplications\tTransfers\tLosses\tSpeciations"
        total_d, total_t, total_l, total_s = None, None, None, None
        for i, line in enumerate(lines):
            if line.startswith("Total "):
                parts = line.strip().split("\t")
                # Format: "Total \tD\tT\tL\tS"
                total_d = float(parts[1])
                total_t = float(parts[2])
                total_l = float(parts[3])
                total_s = float(parts[4])
                break

        assert total_s is not None, "Could not find Total event counts in .uml_rec"

        # Speciations should be exactly n_leaves - 1 = 35 per sample (averaged)
        assert total_s == pytest.approx(35.0, abs=1.0), (
            f"Expected ~35 speciations per sample, got {total_s}"
        )

        # D, T, L should all be non-negative
        assert total_d >= 0, f"Negative duplication count: {total_d}"
        assert total_t >= 0, f"Negative transfer count: {total_t}"
        assert total_l >= 0, f"Negative loss count: {total_l}"

        # Total events should be reasonable (not zero, not astronomical)
        total_events = total_d + total_t + total_l + total_s
        assert total_events > 35.0, (
            f"Total events ({total_events}) suspiciously low"
        )
        assert total_events < 1000.0, (
            f"Total events ({total_events}) suspiciously high"
        )
