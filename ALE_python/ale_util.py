"""I/O utilities for ALE (Amalgamated Likelihood Estimation).

Ported from C++ ALE_util.cpp. Provides functions for observing, loading,
and saving ApproxPosterior objects from/to files and strings.
"""

from __future__ import annotations

import os
from typing import Union

from .ale import ApproxPosterior


def observe_ALE_from_file(
    fname_or_fnames: Union[str, list[str]],
    burnin: int = 0,
    every: int = 1,
    until: int = -1,
) -> ApproxPosterior:
    """Read tree(s) from one or more files and build an ApproxPosterior.

    For each file, lines containing '(' are treated as tree strings.
    The first *burnin* trees are skipped, then every *every*-th tree is
    kept. If *until* > 0 the collection is truncated to that many trees.

    Parameters
    ----------
    fname_or_fnames:
        A single filename or a list of filenames.
    burnin:
        Number of leading trees to discard per file.
    every:
        Sub-sampling stride (keep every *every*-th tree after burn-in).
    until:
        Maximum number of trees to use (all if <= 0).

    Returns
    -------
    ApproxPosterior
        The observed approximate posterior object.
    """
    if isinstance(fname_or_fnames, str):
        fnames = [fname_or_fnames]
    else:
        fnames = list(fname_or_fnames)

    trees: list[str] = []

    for fname in fnames:
        if not os.path.isfile(fname):
            raise FileNotFoundError(
                f"Error, file {fname} does not seem accessible."
            )

        with open(fname) as fh:
            tree_i = 0
            for line in fh:
                line = line.rstrip("\n").strip()
                if not line:
                    continue
                if "(" in line:
                    tree_i += 1
                    if tree_i > burnin and tree_i % every == 0:
                        trees.append(line)
                elif line.endswith(";"):
                    # Singleton tree like "A;" — tokenize and take
                    # the first field, matching C++ boost::split on
                    # ",;: " delimiters.
                    name = line.rstrip(";").split(",")[0].split(":")[0].split()[0]
                    if name:
                        tree_i += 1
                        if tree_i > burnin and tree_i % every == 0:
                            trees.append(name)

    if not trees:
        raise ValueError("No trees found in the provided file(s).")

    if until > 0:
        trees = trees[:until]

    ale = ApproxPosterior(trees[0])
    ale.observation(trees)
    return ale


def observe_ALE_from_strings(trees: list[str]) -> ApproxPosterior:
    """Build an ApproxPosterior from a list of tree strings.

    Parameters
    ----------
    trees:
        Newick tree strings.

    Returns
    -------
    ApproxPosterior
    """
    if not trees:
        raise ValueError("Tree list must not be empty.")

    ale = ApproxPosterior(trees[0])
    ale.observation(trees)
    return ale


def observe_ALE_from_string(tree: str) -> ApproxPosterior:
    """Build an ApproxPosterior from a single tree string.

    Parameters
    ----------
    tree:
        A Newick tree string.

    Returns
    -------
    ApproxPosterior
    """
    return observe_ALE_from_strings([tree])


def load_ALE_from_file(fname: str) -> ApproxPosterior:
    """Load a previously saved ApproxPosterior from a ``.ale`` file.

    Parameters
    ----------
    fname:
        Path to the saved state file.

    Returns
    -------
    ApproxPosterior
    """
    ale = ApproxPosterior()
    ale.load_state(fname)
    return ale


def save_ALE_to_file(ale: ApproxPosterior, fname: str) -> str:
    """Save an ApproxPosterior to a file.

    Parameters
    ----------
    ale:
        The ApproxPosterior instance to save.
    fname:
        Destination file path.

    Returns
    -------
    str
        The filename that was written.
    """
    ale.save_state(fname)
    return fname
