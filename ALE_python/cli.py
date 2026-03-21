"""Main command-line interface that dispatches to individual ALE programs.

Usage: python -m ALE_python <program_name> [args...]

Programs: ALEobserve, ALEml_undated, ALEmcmc_undated, ALEcount,
          ls_leaves, CCPscore, ALEadd, ALEevaluate_undated
"""

import math
import os
import random
import sys

from .ale import approx_posterior as ApproxPosterior
from .ale_util import (
    load_ALE_from_file,
    observe_ALE_from_file,
    observe_ALE_from_string,
    save_ALE_to_file,
)
from .exodt import ExODTModel
from .newick import get_leaves, parse_newick, to_newick

ALE_VERSION = "1.0"


# ---------------------------------------------------------------------------
# Argument parsing helpers
# ---------------------------------------------------------------------------

def _print_subcommand_help(name):
    """Print the argparse help for a specific subcommand and exit."""
    parser = _build_argparse()
    parser.parse_args([name, "--help"])


def _normalize_argv(argv):
    """Convert ``--key value`` long options to ``key=value`` form.

    This allows callers to use either the original C++ style
    (``delta=0.01``, ``rate_multiplier:tau_to:43:0.0``) or standard
    long options (``--delta 0.01``, ``--rate-multiplier tau_to:43:0.0``).
    Positional arguments (not starting with ``--``) are passed through
    unchanged.  Bare flags (``--MLOR``, ``--reldate``) become ``MLOR``
    and ``reldate``.
    """
    result = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg.startswith("--"):
            key = arg[2:].replace("-", "_")
            # Peek: is the next token a value (doesn't start with --)?
            if i + 1 < len(argv) and not argv[i + 1].startswith("--"):
                result.append(f"{key}={argv[i + 1]}")
                i += 2
            else:
                # Boolean flag like --MLOR or --reldate
                result.append(key)
                i += 1
        else:
            result.append(arg)
            i += 1
    return result


def _parse_kwargs(argv):
    """Parse key=value and key:value arguments from an argv list.

    Returns (positional_args, keyword_dict).
    """
    positional = []
    kwargs = {}
    for arg in argv:
        if "=" in arg:
            key, _, val = arg.partition("=")
            kwargs[key] = val
        elif ":" in arg:
            # Some arguments use colon separators (e.g. S_branch_lengths:0.2)
            key, _, val = arg.partition(":")
            kwargs[key] = val
        else:
            positional.append(arg)
    return positional, kwargs


# ---------------------------------------------------------------------------
# ALEobserve
# ---------------------------------------------------------------------------

def ALEobserve(argv):
    """Observe gene trees and build an ALE file.

    argv[0] = gene tree file, additional positional args are more files.
    burnin=N keyword arg (default 0).
    """
    print(f"ALEobserve using ALE v{ALE_VERSION}")

    if not argv:
        _print_subcommand_help("ALEobserve")
        return 1

    positional, kwargs = _parse_kwargs(argv)
    burnin = int(kwargs.get("burnin", 0))

    first_file = positional[0].strip()
    ale_files = [f.strip() for f in positional]
    ale_name = first_file + ".ale"

    ale = observe_ALE_from_file(ale_files, burnin=burnin)

    print(f"# observe {ale.observations} tree(s) from: {' '.join(ale_files)}")
    print(f"{burnin} burn in per file discarded.")
    ale.save_state(ale_name)
    print(f"# saved in {ale_name}")
    print("# mpp tree from sample: ")
    mpp_tree_str, _mpp_pp = ale.mpp_tree()
    print(mpp_tree_str)
    return 0


# ---------------------------------------------------------------------------
# ALEml_undated
# ---------------------------------------------------------------------------

def ALEml_undated(argv):
    """Maximum-likelihood reconciliation under the undated DTL model."""
    import numpy as np

    print(f"ALEml_undated using ALE v{ALE_VERSION}")

    if len(argv) < 2:
        _print_subcommand_help("ALEml_undated")
        return 1

    S_treefile = argv[0]
    ale_file_arg = argv[1]

    if not os.path.isfile(S_treefile):
        print(f"Error, file {S_treefile} does not seem accessible.")
        sys.exit(1)

    with open(S_treefile) as f:
        Sstring = f.readline().strip()
    print(f"Read species tree from: {S_treefile}..")

    ale = load_ALE_from_file(ale_file_arg)
    print(f"Read summary of tree sample for {ale.observations} trees from: {ale_file_arg}..")

    # Output file radical: basename of S_treefile + "_" + basename of ale_file
    ale_file = os.path.basename(S_treefile) + "_" + os.path.basename(ale_file_arg)

    model = ExODTModel()

    samples = 100
    O_R = 1.0
    beta = 1.0
    delta_fixed = False
    tau_fixed = False
    lambda_fixed = False
    DT_fixed = False
    MLOR = False
    no_T = False

    delta = 1e-2
    tau = 1e-2
    lambda_ = 1e-1
    DT_ratio = 0.05
    fraction_missing_file = ""
    output_species_tree = ""
    rate_multipliers = {}
    ml_branch_multipliers = []  # [(branch_id, rate_name)] for per-branch optimization

    model.set_model_parameter("undatedBL", 0)
    model.set_model_parameter("reldate", 0)

    # Parse optional arguments
    for arg in argv[2:]:
        if "=" in arg or ":" in arg:
            sep = "=" if "=" in arg else ":"
            tokens = arg.replace("=", ":").split(":")
            key = tokens[0]
        else:
            continue

        if key == "sample":
            samples = int(tokens[1])
        elif key == "separators":
            model.set_model_parameter("gene_name_separators", tokens[1])
        elif key == "delta":
            delta = float(tokens[1])
            delta_fixed = True
            print(f"# delta fixed to {delta}")
        elif key == "tau":
            tau = float(tokens[1])
            tau_fixed = True
            if tau < 1e-10:
                no_T = True
                print("# tau fixed to no transfer!")
                tau = 1e-19
            else:
                print(f"# tau fixed to {tau}")
        elif key == "lambda":
            lambda_ = float(tokens[1])
            lambda_fixed = True
            print(f"# lambda fixed to {lambda_}")
        elif key == "DT":
            DT_ratio = float(tokens[1])
            DT_fixed = True
            model.set_model_parameter("DT_ratio", DT_ratio)
            print(f"# D/T ratio fixed to {model.scalar_parameter['DT_ratio']}")
        elif key == "O_R":
            O_R = float(tokens[1])
            print(f"# O_R set to {O_R}")
        elif key == "beta":
            beta = float(tokens[1])
            print(f"# beta set to {beta}")
        elif key == "fraction_missing":
            fraction_missing_file = tokens[1]
            print(f"# File containing fractions of missing genes set to {fraction_missing_file}")
        elif key == "S_branch_lengths":
            model.set_model_parameter("undatedBL", 1)
            if len(tokens) == 1 or tokens[1] == "":
                model.set_model_parameter("root_BL", 1)
                print("# unsing branch lengths of input S tree as rate multipliers with 1 at root! ")
            else:
                root_rm = float(tokens[1])
                model.set_model_parameter("root_BL", root_rm)
                print(f"# unsing branch lengths of input S tree as rate multipliers with {root_rm} at root! ")
        elif key == "reldate":
            print("Respecting realtive ages from input S tree, please make sure input S tree is ultrametric!")
            model.set_model_parameter("reldate", 1)
        elif key == "MLOR":
            print("Optimizing root origination multiplier.")
            MLOR = True
        elif key == "rate_multiplier":
            rate_name = tokens[1]
            e = int(tokens[2])
            rm = float(tokens[3])
            if rm >= -1:
                print(f"# rate multiplier for rate {rate_name} on branch with ID {e} set to {rm}")
                rate_multipliers.setdefault("rate_multiplier_" + rate_name, {})[e] = rm
            else:
                print(f"# rate multiplier for rate {rate_name} on branch with ID {e} to be optimized")
                ml_branch_multipliers.append((e, "rate_multiplier_" + rate_name))
        elif key == "output_species_tree":
            val = tokens[1].lower()
            if val in ("y", "ye", "yes"):
                output_species_tree = ale_file + ".spTree"
                print(f"# outputting the annotated species tree to {output_species_tree}")
        elif key == "seed":
            seed_val = int(tokens[1])
            print(f"Set random seed to {seed_val}")
            random.seed(seed_val)
            np.random.seed(seed_val)

    model.set_model_parameter("BOOTSTRAP_LABELS", "yes")
    model.construct_undated(Sstring, fraction_missing_file)

    # Apply rate multipliers
    for rm_name, rm_dict in rate_multipliers.items():
        for e, rm_val in rm_dict.items():
            model.vector_parameter[rm_name][e] = rm_val

    model.set_model_parameter("seq_beta", beta)
    model.set_model_parameter("O_R", O_R)
    model.set_model_parameter("delta", delta)
    model.set_model_parameter("tau", tau)
    model.set_model_parameter("lambda", lambda_)

    model.calculate_undatedEs()
    print("Reconciliation model initialised, starting DTL rate optimisation..")

    # Optimize if not all rates are fixed, or if branch multipliers need optimization
    if not (delta_fixed and tau_fixed and lambda_fixed and not MLOR and not ml_branch_multipliers):
        print("#optimizing rates")

        from scipy.optimize import minimize

        def neg_log_lk(params):
            idx = 0
            d_val = delta
            t_val = tau
            l_val = lambda_
            o_val = O_R

            if not delta_fixed and not DT_fixed:
                d_val = params[idx]
                idx += 1
            if not tau_fixed and not DT_fixed:
                t_val = params[idx]
                idx += 1
            if not lambda_fixed:
                l_val = params[idx]
                idx += 1
            if DT_fixed:
                t_val = params[idx]
                idx += 1
                d_val = t_val * model.scalar_parameter.get("DT_ratio", DT_ratio)
            if MLOR:
                o_val = params[idx]
                idx += 1

            for branch_e, rm_name in ml_branch_multipliers:
                multiplier = params[idx]
                idx += 1
                model.vector_parameter[rm_name][branch_e] = max(multiplier, 1e-7)

            model.set_model_parameter("delta", max(d_val, 1e-10))
            model.set_model_parameter("tau", max(t_val, 1e-10))
            model.set_model_parameter("lambda", max(l_val, 1e-10))
            if MLOR:
                model.set_model_parameter("O_R", max(o_val, 1e-10))
            model.calculate_undatedEs()
            lk = model.pun(ale, False, no_T)
            if lk <= 0:
                return 1e50
            return -math.log(lk)

        # Build initial parameter vector
        x0 = []
        param_names = []
        if not delta_fixed and not DT_fixed:
            x0.append(delta)
            param_names.append("delta")
            print("#optimizing delta rate")
        if not tau_fixed and not DT_fixed:
            x0.append(tau)
            param_names.append("tau")
            print("#optimizing tau rate")
        if not lambda_fixed:
            x0.append(lambda_)
            param_names.append("lambda")
            print("#optimizing lambda rate")
        if DT_fixed:
            x0.append(tau)
            param_names.append("tau_DT")
            print("#optimizing delta and tau rates with fixed D/T ratio")
        if MLOR:
            x0.append(1.0)
            param_names.append("O_R")
            print("#optimizing O_R")

        for branch_e, rm_name in ml_branch_multipliers:
            x0.append(1.0)
            param_names.append(f"rm_{rm_name}_{branch_e}")
            print(f"#optimizing for branch {branch_e} ratemultiplier {rm_name}")

        x0 = np.array(x0)

        # Build per-parameter bounds matching C++ constraints
        rate_lo = 1e-6 if no_T else 1e-10
        bounds = []
        for pname in param_names:
            if pname == "O_R":
                bounds.append((1e-10, 1000.0))
            elif pname.startswith("rm_"):
                bounds.append((1e-7, 10000.0))
            else:
                bounds.append((rate_lo, 100.0))
        result = minimize(
            neg_log_lk,
            x0,
            method="Nelder-Mead",
            bounds=bounds,
            options={"maxiter": 10000, "xatol": 1e-6, "fatol": 1e-6, "adaptive": True},
        )

        # Extract optimized parameters
        idx = 0
        if not delta_fixed and not DT_fixed:
            delta = result.x[idx]
            idx += 1
        if not tau_fixed and not DT_fixed:
            tau = result.x[idx]
            idx += 1
        if not lambda_fixed:
            lambda_ = result.x[idx]
            idx += 1
        if DT_fixed:
            tau = result.x[idx]
            idx += 1
            delta = tau * model.scalar_parameter.get("DT_ratio", DT_ratio)
        if MLOR:
            O_R = result.x[idx]
            idx += 1

        ml_rm_strings = []
        for branch_e, rm_name in ml_branch_multipliers:
            multiplier = result.x[idx]
            idx += 1
            model.vector_parameter[rm_name][branch_e] = multiplier
            ml_rm_strings.append(f"{rm_name}\t{branch_e}\t{multiplier};")

        mlll = -result.fun
    else:
        mlll = math.log(model.pun(ale, False, no_T))

    print()
    print(f"ML rates:  delta={delta}; tau={tau}; lambda={lambda_}; O_R={O_R}.")
    if ml_branch_multipliers:
        print("ML rate multipliers:")
        for s in ml_rm_strings:
            print(s)
    print(f"LL={mlll}")

    # Set final rates for sampling
    model.set_model_parameter("delta", delta)
    model.set_model_parameter("tau", tau)
    model.set_model_parameter("lambda", lambda_)
    model.set_model_parameter("O_R", O_R)
    model.calculate_undatedEs()
    model.pun(ale, False, no_T)

    # Sample reconciled gene trees
    print("Sampling reconciled gene trees..")
    sample_strings = []
    total_events = {"D": 0.0, "T": 0.0, "L": 0.0, "S": 0.0}
    for i in range(int(samples)):
        # Reset event counts per sample
        model.MLRec_events.clear()
        model.Ttokens = []
        sample_tree = model.sample_undated(no_T)
        sample_strings.append(sample_tree)
        for key in total_events:
            total_events[key] += model.MLRec_events.get(key, 0.0)

    # Write .uml_rec output
    outname = ale_file + ".uml_rec"
    with open(outname, "w") as fout:
        fout.write(f"#ALEml_undated using ALE v{ALE_VERSION} by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;\n\n")
        s_tree_str = model.string_parameter.get("S_with_ranks", model.string_parameter.get("S_un", Sstring))
        fout.write(f"S:\t{s_tree_str}\n")
        fout.write("\n")
        fout.write(f"Input ale from:\t{ale_file}\n")
        fout.write(f">logl: {mlll}\n")
        fout.write("rate of\t Duplications\tTransfers\tLosses\n")
        fout.write(f"ML \t{delta}\t{tau}\t{lambda_}\n")
        fout.write(f"{int(samples)} reconciled G-s:\n\n")
        for s in sample_strings:
            fout.write(s + "\n")
        fout.write("# of\t Duplications\tTransfers\tLosses\tSpeciations\n")
        div = samples if samples > 0 else 1
        fout.write(
            f"Total \t{total_events['D'] / div}\t"
            f"{total_events['T'] / div}\t"
            f"{total_events['L'] / div}\t"
            f"{total_events['S'] / div}\n"
        )
        fout.write("\n")
        fout.write("# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\textinction_prob\tpresence\tLL\n")
        fout.write(model.counts_string_undated(samples))

    # Output species tree
    if output_species_tree:
        with open(output_species_tree, "w") as fout:
            s_tree_str = model.string_parameter.get("S_with_ranks", model.string_parameter.get("S_un", Sstring))
            fout.write(s_tree_str + "\n")

    print(f"Results in: {outname}")

    # Write transfer file
    t_name = ale_file + ".uTs"
    with open(t_name, "w") as tout:
        tout.write("#from\tto\tfreq.\n")
        for e in range(model.last_branch):
            for f in range(model.last_branch):
                if model.T_to_from[e][f] > 0:
                    if e < model.last_leaf:
                        e_name = f"{model._node_name[model._id_nodes[e].id]}({e})"
                    else:
                        e_name = str(e)
                    if f < model.last_leaf:
                        f_name = f"{model._node_name[model._id_nodes[f].id]}({f})"
                    else:
                        f_name = str(f)
                    tout.write(f"\t{e_name}\t{f_name}\t{model.T_to_from[e][f] / div}\n")
    print(f"Transfers in: {t_name}")
    return 0


# ---------------------------------------------------------------------------
# ALEmcmc_undated
# ---------------------------------------------------------------------------

def _scale_double_constrained(value, maxi, lam):
    """Scaling move on a positive real.  Returns (new_value, hastings_ratio)."""
    u = random.random()
    scaling_factor = math.exp(lam * (u - 0.5))
    new_value = value * scaling_factor
    if new_value < 0.00001:
        new_value = 0.00001
    if new_value > maxi:
        new_value = maxi
    return new_value, scaling_factor


def _compute_exponential_log_probability(param, value):
    """Log probability under Exponential(param)."""
    if param <= 0:
        return 0.0
    return math.log(param) - param * value


def _compute_log_prior(o, d, t, l, prior_o, prior_d, prior_t, prior_l):
    pp = 0.0
    pp += _compute_exponential_log_probability(prior_o, o)
    pp += _compute_exponential_log_probability(prior_d, d)
    pp += _compute_exponential_log_probability(prior_t, t)
    pp += _compute_exponential_log_probability(prior_l, l)
    return pp


def _compute_log_lk(model, ale, o, d, t, l):
    model.set_model_parameter("O_R", o)
    model.set_model_parameter("delta", d)
    model.set_model_parameter("tau", t)
    model.set_model_parameter("lambda", l)
    model.calculate_undatedEs()
    lk = model.pun(ale)
    if lk <= 0:
        return -1e100
    return math.log(lk)


def ALEmcmc_undated(argv):
    """MCMC sampling of reconciliation under the undated DTL model."""
    print(f"ALEmcmc using ALE v{ALE_VERSION}")

    if len(argv) < 2:
        _print_subcommand_help("ALEmcmc_undated")
        return 1

    S_treefile = argv[0]
    ale_file_arg = argv[1]

    if not os.path.isfile(S_treefile):
        print(f"Error, file {S_treefile} does not seem accessible.")
        sys.exit(1)

    with open(S_treefile) as f:
        Sstring = f.readline().strip()
    print(f"Read species tree from: {S_treefile}..")

    ale = load_ALE_from_file(ale_file_arg)
    print(f"Read summary of tree sample for {ale.observations} trees from: {ale_file_arg}..")

    ale_file = os.path.basename(ale_file_arg)

    model = ExODTModel()

    samples = 100
    prior_origination = 1.0
    prior_delta = 0.01
    prior_tau = 0.01
    prior_lambda = 0.1
    sampling_rate = 1
    beta = 1.0
    fraction_missing_file = ""
    output_species_tree = ""
    rate_multipliers = {}

    model.set_model_parameter("undatedBL", 0)
    model.set_model_parameter("reldate", 0)

    for arg in argv[2:]:
        if "=" in arg or ":" in arg:
            tokens = arg.replace("=", ":").split(":")
            key = tokens[0]
        else:
            continue

        if key == "sample":
            samples = int(tokens[1])
        elif key == "separators":
            model.set_model_parameter("gene_name_separators", tokens[1])
        elif key == "delta":
            prior_delta = float(tokens[1])
            print(f"# priorDelta fixed to {prior_delta}")
        elif key == "tau":
            prior_tau = float(tokens[1])
            print(f"# priorTau fixed to {prior_tau}")
        elif key == "lambda":
            prior_lambda = float(tokens[1])
            print(f"# priorLambda fixed to {prior_lambda}")
        elif key == "O_R":
            prior_origination = float(tokens[1])
            print(f"# priorOrigination set to {prior_origination}")
        elif key == "beta":
            beta = float(tokens[1])
            print(f"# beta set to {beta}")
        elif key == "sampling_rate":
            sampling_rate = int(tokens[1])
            print(f"# sampling_rate set to {sampling_rate}")
        elif key == "fraction_missing":
            fraction_missing_file = tokens[1]
            print(f"# File containing fractions of missing genes set to {fraction_missing_file}")
        elif key == "output_species_tree":
            val = tokens[1].lower()
            if val in ("y", "ye", "yes"):
                output_species_tree = ale_file + ".spTree"
                print(f"# outputting the annotated species tree to {output_species_tree}")
        elif key == "S_branch_lengths":
            model.set_model_parameter("undatedBL", 1)
            if len(tokens) == 1 or tokens[1] == "":
                model.set_model_parameter("root_BL", 1)
                print("# unsing branch lengths of input S tree as rate multipliers with 1 at root! ")
            else:
                root_rm = float(tokens[1])
                model.set_model_parameter("root_BL", root_rm)
                print(f"# unsing branch lengths of input S tree as rate multipliers with {root_rm} at root! ")
        elif key == "rate_multiplier":
            rate_name = tokens[1]
            e = int(tokens[2])
            rm = float(tokens[3])
            print(f"# rate multiplier for rate {rate_name} on branch with ID {e} set to {rm}")
            rate_multipliers.setdefault("rate_multiplier_" + rate_name, {})[e] = rm
        elif key == "reldate":
            print("Respecting realtive ages from input S tree, please make sure input S tree is ultrametric!")
            model.set_model_parameter("reldate", 1)

    model.set_model_parameter("BOOTSTRAP_LABELS", "yes")
    model.set_model_parameter("seq_beta", beta)

    model.construct_undated(Sstring, fraction_missing_file)

    # Apply rate multipliers after construct_undated so vectors are initialized
    for rm_name, rm_dict in rate_multipliers.items():
        for e, rm_val in rm_dict.items():
            model.vector_parameter[rm_name][e] = rm_val

    # Draw initial values from exponential priors
    current_origination = random.expovariate(prior_origination) if prior_origination > 0 else 1.0
    current_delta = random.expovariate(prior_delta) if prior_delta > 0 else 0.01
    current_tau = random.expovariate(prior_tau) if prior_tau > 0 else 0.01
    current_lambda = random.expovariate(prior_lambda) if prior_lambda > 0 else 0.1

    new_origination = current_origination
    new_delta = current_delta
    new_tau = current_tau
    new_lambda = current_lambda

    current_log_lk = _compute_log_lk(model, ale, current_origination, current_delta, current_tau, current_lambda)
    current_log_prior = _compute_log_prior(
        current_origination, current_delta, current_tau, current_lambda,
        prior_origination, prior_delta, prior_tau, prior_lambda,
    )

    print(f"Initial logLK: {current_log_lk} and logPrior: {current_log_prior}")
    print("Reconciliation model initialised, starting DTL rate sampling..")

    # Move setup
    ORIGINATION_ID = 0
    DELTA_ID = 1
    LAMBDA_ID = 2
    TAU_ID = 3
    move_weights = [1.0, 1.0, 1.0, 1.0]
    max_sum_dtl = 10.0
    max_origination = 1000000.0
    scale_move_params = [0.1, 1.0, 10.0]

    # CSV trace file
    mcmc_outname = ale_file + "_umcmc.csv"
    mcmc_fh = open(mcmc_outname, "w")
    mcmc_fh.write("Iteration\tLogLk\tLogPrior\tOrigination\tDelta\tTau\tLambda\n")
    print("Iteration\tLogLk\tLogPrior\tOrigination\tDelta\tTau\tLambda")

    def _pick_weighted(weights):
        total = sum(weights)
        r = random.random() * total
        cumsum = 0.0
        for i, w in enumerate(weights):
            cumsum += w
            if r < cumsum:
                return i
        return len(weights) - 1

    def _do_mcmc_step():
        nonlocal current_origination, current_delta, current_tau, current_lambda
        nonlocal new_origination, new_delta, new_tau, new_lambda
        nonlocal current_log_lk, current_log_prior

        move = _pick_weighted(move_weights)
        scale = scale_move_params[_pick_weighted([1, 1, 1])]
        hastings_ratio = 1.0

        if move == ORIGINATION_ID:
            new_origination, hastings_ratio = _scale_double_constrained(current_origination, max_origination, scale)
        elif move == DELTA_ID:
            new_delta, hastings_ratio = _scale_double_constrained(current_delta, max_sum_dtl - current_lambda - current_tau, scale)
        elif move == LAMBDA_ID:
            new_lambda, hastings_ratio = _scale_double_constrained(current_lambda, max_sum_dtl - current_delta - current_tau, scale)
        elif move == TAU_ID:
            new_tau, hastings_ratio = _scale_double_constrained(current_tau, max_sum_dtl - current_lambda - current_delta, scale)

        new_log_lk = _compute_log_lk(model, ale, new_origination, new_delta, new_tau, new_lambda)
        new_log_prior = _compute_log_prior(
            new_origination, new_delta, new_tau, new_lambda,
            prior_origination, prior_delta, prior_tau, prior_lambda,
        )

        log_ratio = (new_log_lk + new_log_prior) - (current_log_lk + current_log_prior)
        # Clamp exponent to prevent OverflowError; large positive values
        # mean the proposal is overwhelmingly better, so auto-accept.
        acceptance_prob = math.exp(min(log_ratio, 700.0)) * hastings_ratio

        if random.random() < acceptance_prob:
            current_origination = new_origination
            current_delta = new_delta
            current_tau = new_tau
            current_lambda = new_lambda
            current_log_lk = new_log_lk
            current_log_prior = new_log_prior
        else:
            new_origination = current_origination
            new_delta = current_delta
            new_tau = current_tau
            new_lambda = current_lambda

    # BURNIN
    burnin_length = 100
    print(f"BURNIN during {burnin_length} iterations.")
    print("LogLk\tLogPrior\tOrigination\tDelta\tTau\tLambda")
    for i in range(burnin_length):
        _do_mcmc_step()
        print(f"{i}\t{current_log_lk}\t{current_log_prior}\t{current_origination}\t{current_delta}\t{current_tau}\t{current_lambda}")

    # MCMC
    total_iterations = int(samples * sampling_rate)
    print(f"MCMC during {total_iterations} iterations.")
    print("LogLk\tLogPrior\tOrigination\tDelta\tTau\tLambda")

    sample_strings = []
    num_speciations = 0.0
    num_duplications = 0.0
    num_transfers = 0.0
    num_losses = 0.0
    t_to_from_accum = {}  # {("from_name", "to_name"): count}

    for i in range(total_iterations):
        _do_mcmc_step()

        if i % sampling_rate == 0:
            model.MLRec_events.clear()
            model.reset_T_to_from()
            model.Ttokens = []
            # Recompute with current params for sampling
            model.set_model_parameter("O_R", current_origination)
            model.set_model_parameter("delta", current_delta)
            model.set_model_parameter("tau", current_tau)
            model.set_model_parameter("lambda", current_lambda)
            model.calculate_undatedEs()
            model.pun(ale)
            sample_tree = model.sample_undated()
            sample_strings.append(sample_tree)
            num_speciations += model.MLRec_events.get("S", 0)
            num_duplications += model.MLRec_events.get("D", 0)
            num_transfers += model.MLRec_events.get("T", 0)
            num_losses += model.MLRec_events.get("L", 0)
            # Accumulate T_to_from into running aggregate
            for e in range(model.last_branch):
                for f in range(model.last_branch):
                    if model.T_to_from[e][f] > 0:
                        e_name = model._node_name[model._id_nodes[e].id] if e < model.last_leaf else str(e)
                        f_name = model._node_name[model._id_nodes[f].id] if f < model.last_leaf else str(f)
                        key = (e_name, f_name)
                        t_to_from_accum[key] = t_to_from_accum.get(key, 0.0) + model.T_to_from[e][f]
            print(f"{i}\t{current_log_lk}\t{current_log_prior}\t{current_origination}\t{current_delta}\t{current_tau}\t{current_lambda}")

        mcmc_fh.write(f"{i}\t{current_log_lk}\t{current_log_prior}\t{current_origination}\t{current_delta}\t{current_tau}\t{current_lambda}\n")

    mcmc_fh.close()

    # Write .umcmc_rec output
    outname = ale_file + ".umcmc_rec"
    with open(outname, "w") as fout:
        fout.write(f"#ALEmcmc_undated using ALE v{ALE_VERSION} by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;\n\n")
        s_tree_str = model.string_parameter.get("S_with_ranks", model.string_parameter.get("S_un", Sstring))
        fout.write(f"S:\t{s_tree_str}\n")
        fout.write("\n")
        fout.write(f"Input ale from:\t{ale_file}\n")
        fout.write("\n")
        fout.write(f"{int(samples)} reconciled G-s:\n\n")
        for s in sample_strings:
            fout.write(s + "\n")
        fout.write("# of\t Duplications\tTransfers\tLosses\tSpeciations\n")
        div = samples if samples > 0 else 1
        fout.write(
            f"Total \t{num_duplications / div}\t"
            f"{num_transfers / div}\t"
            f"{num_losses / div}\t"
            f"{num_speciations / div}\n"
        )
        fout.write("\n")
        fout.write("# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\textinction_prob\tpresence\tLL\n")
        fout.write(model.counts_string_undated(samples))

    if output_species_tree:
        with open(output_species_tree, "w") as fout:
            s_tree_str = model.string_parameter.get("S_with_ranks", model.string_parameter.get("S_un", Sstring))
            fout.write(s_tree_str + "\n")

    print(f"Results in: {outname}")

    # Transfer file (from accumulated counts across all samples)
    t_name = ale_file + "_mcmc.uTs"
    with open(t_name, "w") as tout:
        tout.write("#from\tto\tfreq.\n")
        div = samples if samples > 0 else 1
        for (e_name, f_name), count in sorted(t_to_from_accum.items()):
            tout.write(f"\t{e_name}\t{f_name}\t{count / div}\n")
    print(f"Transfers in: {t_name}")
    return 0


# ---------------------------------------------------------------------------
# ALEcount
# ---------------------------------------------------------------------------

def ALEcount(argv):
    """Load an ALE file and print the number of amalgamated trees."""
    if not argv:
        _print_subcommand_help("ALEcount")
        return 1

    ale_file = argv[0]
    ale = load_ALE_from_file(ale_file)
    print(ale.count_trees())
    return 0


# ---------------------------------------------------------------------------
# ls_leaves
# ---------------------------------------------------------------------------

def ls_leaves(argv):
    """For each file: parse tree, list leaves with counts."""
    if not argv:
        _print_subcommand_help("ls_leaves")
        return 1

    names = {}
    for fname in argv:
        with open(fname) as f:
            tree_str = f.readline().strip()
        root = parse_newick(tree_str)
        leaves = get_leaves(root)
        for leaf in leaves:
            name = leaf.name
            names[name] = names.get(name, 0) + 1

    for name in sorted(names.keys()):
        print(f"{name} {names[name]}")
    return 0


# ---------------------------------------------------------------------------
# CCPscore
# ---------------------------------------------------------------------------

def CCPscore(argv):
    """Load an ALE file, read a tree, and print log(p(tree))."""
    if len(argv) < 2:
        _print_subcommand_help("CCPscore")
        return 1

    ale_file = argv[0]
    tree_file = argv[1]

    ale = load_ALE_from_file(ale_file)
    with open(tree_file) as f:
        tree_str = f.readline().strip()
    p_val = ale.p(tree_str)
    if p_val > 0:
        print(math.log(p_val))
    else:
        print("-inf")
    return 0


# ---------------------------------------------------------------------------
# ALEadd
# ---------------------------------------------------------------------------

def ALEadd(argv):
    """Load an existing ALE file and add new trees."""
    print(f"ALEadd using ALE v{ALE_VERSION}")

    if len(argv) < 2:
        _print_subcommand_help("ALEadd")
        return 1

    ale_file = argv[0].strip()
    trees_file = argv[1].strip()
    ale_name = ale_file

    burnin = 0
    every = 1
    until = -1
    weight = 1.0

    for arg in argv[2:]:
        if "=" in arg:
            tokens = arg.split("=")
            key = tokens[0]
            val = tokens[1]
            if key == "burnin":
                burnin = int(val)
            elif key == "every":
                every = int(val)
            elif key == "until":
                until = int(val)
            elif key == "weight":
                weight = float(val)
            elif key == "outfile":
                ale_name = val

    ale = load_ALE_from_file(ale_file)
    print(".")

    # Read trees from file
    trees = []
    tree_i = 0
    with open(trees_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if "(" in line:
                if tree_i >= burnin and tree_i % every == 0:
                    trees.append(line)
                tree_i += 1

    print("..")

    observe_trees = trees
    if until > 0:
        observe_trees = trees[:until]

    ale.observation(observe_trees, weight=weight)

    print(f"# {len(observe_trees)} new tree(s) observed with weight {weight} from: {trees_file}")
    print(f"; {burnin} trees burnin discarded.")
    print(f"# .ale with {ale.observations} tree(s) from: {ale_file} and {trees_file}")
    ale.save_state(ale_name)
    print(f"# saved in {ale_name}")
    return 0


# ---------------------------------------------------------------------------
# ALEevaluate_undated
# ---------------------------------------------------------------------------

def ALEevaluate_undated(argv):
    """Evaluate a single gene tree under the undated DTL model."""
    print(f"ALEestimate using ALE v{ALE_VERSION}")

    if len(argv) < 2:
        _print_subcommand_help("ALEevaluate_undated")
        return 1

    species_tree_file = argv[0].strip()
    gene_tree_file = argv[1].strip()

    # Read species tree
    species_tree_str = ""
    with open(species_tree_file) as f:
        for line in f:
            line = line.strip()
            if "(" in line:
                species_tree_str = line
    print(f"\n\tRead species tree from: {species_tree_file}")

    # Read gene tree
    gene_tree_str = ""
    with open(gene_tree_file) as f:
        for line in f:
            line = line.strip()
            if "(" in line:
                gene_tree_str = line

    ale = ApproxPosterior(gene_tree_str)
    ale.observation([gene_tree_str])
    print(f"\n\tObserved {ale.observations} gene tree(s) from: {gene_tree_file}")

    model = ExODTModel()

    samples = 100
    O_R = 1.0
    beta = 1.0
    delta = 0.01
    tau = 0.01
    lambda_ = 0.1
    fraction_missing_file = ""
    output_files = False

    for arg in argv[2:]:
        if "=" in arg:
            tokens = arg.split("=")
            key = tokens[0]
            val = tokens[1]
            if key == "sample":
                samples = int(val)
            elif key == "separators":
                model.set_model_parameter("gene_name_separators", val)
            elif key == "delta":
                delta = float(val)
                print(f"\n\tDelta fixed to {delta}")
            elif key == "tau":
                tau = float(val)
                print(f"\n\tTau fixed to {tau}")
            elif key == "lambda":
                lambda_ = float(val)
                print(f"Lambda fixed to {lambda_}")
            elif key == "O_R":
                O_R = float(val)
                print(f"\n\tO_R set to {O_R}")
            elif key == "beta":
                beta = float(val)
                print(f"\n\tBeta set to {beta}")
            elif key == "fraction_missing":
                fraction_missing_file = val
                print(f"\n\tFile containing fractions of missing genes set to {fraction_missing_file}")
            elif key == "outputFiles":
                if val.lower() in ("y", "yes"):
                    output_files = True

    model.set_model_parameter("BOOTSTRAP_LABELS", "yes")
    model.construct_undated(species_tree_str, fraction_missing_file)

    model.set_model_parameter("seq_beta", beta)
    model.set_model_parameter("O_R", O_R)
    model.set_model_parameter("delta", delta)
    model.set_model_parameter("tau", tau)
    model.set_model_parameter("lambda", lambda_)

    model.calculate_undatedEs()
    loglk = math.log(max(model.pun(ale, True), sys.float_info.min))
    print(f"\n\tReconciliation model likelihood computed, logLk: {loglk}")

    if output_files:
        print("\n\tSampling reconciled gene trees..")
        sample_strings = []
        total_events = {"D": 0.0, "T": 0.0, "L": 0.0, "S": 0.0}
        for i in range(int(samples)):
            model.MLRec_events.clear()
            model.Ttokens = []
            sample_tree = model.sample_undated()
            sample_strings.append(sample_tree)
            for key in total_events:
                total_events[key] += model.MLRec_events.get(key, 0.0)

        ale_name = os.path.basename(gene_tree_file)
        outname = ale_name + ".uml_rec"
        with open(outname, "w") as fout:
            fout.write(f"#ALEevaluate using ALE v{ALE_VERSION}; CC BY-SA 3.0;\n\n")
            s_tree_str = model.string_parameter.get("S_with_ranks", model.string_parameter.get("S_un", species_tree_str))
            fout.write(f"S:\t{s_tree_str}\n")
            fout.write("\n")
            fout.write(f"Gene tree from:\t{gene_tree_file}\n")
            fout.write(f">logl: {loglk}\n")
            fout.write("rate of\t Duplications\tTransfers\tLosses\n")
            fout.write(f"\t{delta}\t{tau}\t{lambda_}\n")
            fout.write("\n")
            fout.write(f"{int(samples)} reconciled G-s:\n\n")
            for s in sample_strings:
                fout.write(s + "\n")
            fout.write("# of\t Duplications\tTransfers\tLosses\tSpeciations\n")
            div = samples if samples > 0 else 1
            fout.write(
                f"Total \t{total_events['D'] / div}\t"
                f"{total_events['T'] / div}\t"
                f"{total_events['L'] / div}\t"
                f"{total_events['S'] / div}\n"
            )
            fout.write("\n")
            fout.write("# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies\n")
            fout.write(model.counts_string_undated(samples))

        print(f"Results in: {outname}")

        # Transfer file
        t_name = ale_name + ".uTs"
        with open(t_name, "w") as tout:
            tout.write("#from\tto\tfreq.\n")
            for e in range(model.last_branch):
                for f in range(model.last_branch):
                    if model.T_to_from[e][f] > 0:
                        if e < model.last_leaf:
                            e_name = model._node_name[model._id_nodes[e].id]
                        else:
                            e_name = str(e)
                        if f < model.last_leaf:
                            f_name = model._node_name[model._id_nodes[f].id]
                        else:
                            f_name = str(f)
                        tout.write(f"\t{e_name}\t{f_name}\t{model.T_to_from[e][f] / div}\n")
        print(f"Transfers in: {t_name}")

    return 0


# ---------------------------------------------------------------------------
# Dispatch table and main entry point
# ---------------------------------------------------------------------------

PROGRAMS = {
    "ALEobserve": ALEobserve,
    "ALEml_undated": ALEml_undated,
    "ALEmcmc_undated": ALEmcmc_undated,
    "ALEcount": ALEcount,
    "ls_leaves": ls_leaves,
    "CCPscore": CCPscore,
    "ALEadd": ALEadd,
    "ALEevaluate_undated": ALEevaluate_undated,
}


def _build_argparse():
    """Build an argparse parser with subcommands for help/documentation.

    The actual argument parsing is done by each subcommand's bespoke
    parser for backward compatibility with the C++ key=value syntax.
    argparse is used only for --help and subcommand dispatch.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m ALE_python",
        description=(
            "ALE (Amalgamated Likelihood Estimation) — "
            "Python port of the phylogenetic reconciliation toolkit. "
            "Options accept both C++ style (delta=0.01) and long-option "
            "style (--delta 0.01). "
            "Use 'python -m ALE_python <command> --help' for per-command help."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"ALE {ALE_VERSION}",
    )
    subparsers = parser.add_subparsers(dest="command", title="commands")

    # --- ALEobserve ---
    p = subparsers.add_parser(
        "ALEobserve",
        help="Build an .ale file from a sample of gene trees",
        description=(
            "Read a sample of gene trees in Newick format and construct an "
            "approximate posterior (ALE) file summarising clade frequencies. "
            "The output .ale file is used as input by ALEml_undated, "
            "ALEmcmc_undated, and other tools."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Both key=value and --key value styles are accepted for options.\n\n"
            "examples:\n"
            "  python -m ALE_python ALEobserve gene_trees.newicks\n"
            "  python -m ALE_python ALEobserve gene_trees.newicks burnin=1000\n"
            "  python -m ALE_python ALEobserve gene_trees.newicks --burnin 1000"
        ),
    )
    p.add_argument("gene_trees", nargs="+", help="one or more Newick gene tree files")
    p.add_argument("burnin", nargs="?", default="burnin=0", help="burnin=N — discard the first N trees per file (default: 0)")

    # --- ALEml_undated ---
    p = subparsers.add_parser(
        "ALEml_undated",
        help="Maximum-likelihood reconciliation under the undated DTL model",
        description=(
            "Optimise duplication, transfer, and loss (DTL) rates by maximum "
            "likelihood under the undated reconciliation model, then sample "
            "reconciled gene trees. Reads a species tree and an .ale file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Both key=value and --key value styles are accepted for options.\n\n"
            "options:\n"
            "  sample=N                  number of reconciled trees to sample (default: 100)\n"
            "  seed=INT                  random seed for reproducibility\n"
            "  separators=STR            gene-name separator characters\n"
            "  delta=FLOAT               fix duplication rate (skip optimisation)\n"
            "  tau=FLOAT                 fix transfer rate (tau<1e-10 disables transfers)\n"
            "  lambda=FLOAT              fix loss rate\n"
            "  O_R=FLOAT                 origination-at-root multiplier (default: 1.0)\n"
            "  MLOR                      optimise the root origination multiplier\n"
            "  DT=FLOAT                  fix duplication/transfer ratio\n"
            "  beta=FLOAT                weight of sequence evidence (default: 1.0)\n"
            "  fraction_missing=FILE     file with per-species missing-gene fractions\n"
            "  S_branch_lengths:ROOT_BL  use species-tree branch lengths as rate multipliers\n"
            "  reldate                   respect relative dates from an ultrametric species tree\n"
            "  output_species_tree=y     write annotated species tree to a .spTree file\n"
            "  rate_multiplier:RATE:BRANCH:VALUE\n"
            "                            set or optimise a per-branch rate multiplier\n"
            "                            (value >= -1 → fixed; value < -1 → optimised)\n\n"
            "examples:\n"
            "  python -m ALE_python ALEml_undated species.nwk genes.ale\n"
            "  python -m ALE_python ALEml_undated species.nwk genes.ale sample=1000 seed=42\n"
            "  python -m ALE_python ALEml_undated species.nwk genes.ale --delta 0.01 --tau 0.01 --lambda 0.1"
        ),
    )
    p.add_argument("species_tree", help="species tree in Newick format")
    p.add_argument("ale_file", help=".ale file from ALEobserve")
    p.add_argument("options", nargs="*", help="key=value options (see below)")

    # --- ALEmcmc_undated ---
    p = subparsers.add_parser(
        "ALEmcmc_undated",
        help="MCMC sampling of reconciliation under the undated DTL model",
        description=(
            "Run a Metropolis-Hastings MCMC chain to sample DTL rates and "
            "reconciled gene trees under the undated model. Priors are "
            "exponential distributions parametrised by the given rate values."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Both key=value and --key value styles are accepted for options.\n\n"
            "options:\n"
            "  sample=N                  total MCMC iterations (default: 100)\n"
            "  sampling_rate=N           record a sample every N iterations (default: 1)\n"
            "  separators=STR            gene-name separator characters\n"
            "  delta=FLOAT               prior rate for duplications (default: 0.01)\n"
            "  tau=FLOAT                 prior rate for transfers (default: 0.01)\n"
            "  lambda=FLOAT              prior rate for losses (default: 0.1)\n"
            "  O_R=FLOAT                 prior for origination at root (default: 1.0)\n"
            "  beta=FLOAT                weight of sequence evidence (default: 1.0)\n"
            "  fraction_missing=FILE     file with per-species missing-gene fractions\n"
            "  S_branch_lengths:ROOT_BL  use species-tree branch lengths as rate multipliers\n"
            "  rate_multiplier:RATE:BRANCH:VALUE\n"
            "                            fix a per-branch rate multiplier\n"
            "  output_species_tree=y     write annotated species tree\n"
            "  reldate                   respect relative dates from an ultrametric species tree\n\n"
            "examples:\n"
            "  python -m ALE_python ALEmcmc_undated species.nwk genes.ale sample=10000\n"
            "  python -m ALE_python ALEmcmc_undated species.nwk genes.ale --sample 5000 --sampling-rate 10"
        ),
    )
    p.add_argument("species_tree", help="species tree in Newick format")
    p.add_argument("ale_file", help=".ale file from ALEobserve")
    p.add_argument("options", nargs="*", help="key=value options (see below)")

    # --- ALEcount ---
    p = subparsers.add_parser(
        "ALEcount",
        help="Print the number of amalgamated trees in an .ale file",
        description="Load an .ale file and print the number of distinct amalgamated gene tree topologies.",
    )
    p.add_argument("ale_file", help=".ale file to count trees from")

    # --- ls_leaves ---
    p = subparsers.add_parser(
        "ls_leaves",
        help="List leaf names from one or more tree files",
        description="Parse each Newick tree file and print sorted leaf names with occurrence counts.",
    )
    p.add_argument("tree_files", nargs="+", help="one or more Newick tree files")

    # --- CCPscore ---
    p = subparsers.add_parser(
        "CCPscore",
        help="Score a tree under the conditional clade probability model",
        description="Load an .ale file and a tree, then print log(P(tree)) under the CCP model.",
    )
    p.add_argument("ale_file", help=".ale file with clade posteriors")
    p.add_argument("tree_file", help="Newick tree file to score")

    # --- ALEadd ---
    p = subparsers.add_parser(
        "ALEadd",
        help="Add new tree observations to an existing .ale file",
        description=(
            "Load an existing .ale file and observe additional gene trees from "
            "a Newick file, updating the clade frequencies."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Both key=value and --key value styles are accepted for options.\n\n"
            "options:\n"
            "  weight=FLOAT   observation weight for new trees (default: 1.0)\n"
            "  burnin=N       discard first N trees (default: 0)\n"
            "  every=N        keep every Nth tree (default: 1)\n"
            "  until=N        use at most N trees (default: all)\n"
            "  outfile=PATH   output .ale file name (default: overwrite input)\n\n"
            "examples:\n"
            "  python -m ALE_python ALEadd existing.ale new_trees.newicks\n"
            "  python -m ALE_python ALEadd existing.ale new_trees.newicks --burnin 100 --weight 0.5"
        ),
    )
    p.add_argument("ale_file", help="existing .ale file")
    p.add_argument("gene_trees", help="Newick gene tree file to add")
    p.add_argument("options", nargs="*", help="key=value options (see below)")

    # --- ALEevaluate_undated ---
    p = subparsers.add_parser(
        "ALEevaluate_undated",
        help="Evaluate a single gene tree under the undated DTL model",
        description=(
            "Compute the reconciliation likelihood of a single gene tree "
            "against a species tree under the undated DTL model with fixed rates."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Both key=value and --key value styles are accepted for options.\n\n"
            "options:\n"
            "  sample=N           number of reconciled trees to sample (default: 100)\n"
            "  separators=STR     gene-name separator characters\n"
            "  delta=FLOAT        duplication rate (default: 0.01)\n"
            "  tau=FLOAT          transfer rate (default: 0.01)\n"
            "  lambda=FLOAT       loss rate (default: 0.1)\n"
            "  O_R=FLOAT          origination-at-root multiplier (default: 1.0)\n"
            "  beta=FLOAT         weight of sequence evidence (default: 1.0)\n"
            "  fraction_missing=FILE  per-species missing-gene fractions\n"
            "  outputFiles=y      write .uml_rec and .uTs output files\n\n"
            "examples:\n"
            "  python -m ALE_python ALEevaluate_undated species.nwk gene.nwk\n"
            "  python -m ALE_python ALEevaluate_undated species.nwk gene.nwk --delta 0.02 --outputFiles y"
        ),
    )
    p.add_argument("species_tree", help="species tree in Newick format")
    p.add_argument("gene_tree", help="gene tree in Newick format")
    p.add_argument("options", nargs="*", help="key=value options (see below)")

    return parser


def main():
    # If invoked with no args, print help; for --help/--version let argparse handle it
    if len(sys.argv) < 2:
        parser = _build_argparse()
        parser.print_help()
        return
    if sys.argv[1] in ("-h", "--help", "--version"):
        parser = _build_argparse()
        parser.parse_args()
        return

    program_name = sys.argv[1]

    # If the subcommand itself is followed by --help, use argparse
    if "--help" in sys.argv[2:] or "-h" in sys.argv[2:]:
        parser = _build_argparse()
        parser.parse_args()
        return

    if program_name not in PROGRAMS:
        print(f"Error: unknown program '{program_name}'")
        print(f"Available programs: {', '.join(sorted(PROGRAMS.keys()))}")
        sys.exit(1)

    # Normalize --key value args to key=value, then dispatch to the
    # original bespoke parser for backward compatibility
    remaining_args = _normalize_argv(sys.argv[2:])
    ret = PROGRAMS[program_name](remaining_args)
    sys.exit(ret or 0)
