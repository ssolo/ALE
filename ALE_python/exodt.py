"""exODT model -- core reconciliation model for ALE (undated).

Ported from the C++ files: exODT.h, exODT.cpp, undated.cpp.
Focuses on the UNDATED model (most commonly used).

All code by Szollosi GJ et al.; ssolo@elte.hu; GNU GPL 3.0;
Python port.
"""

import math
import re
import sys
from collections import defaultdict

from .newick import parse_newick, get_leaves, get_all_nodes, is_leaf
from .fraction_missing import read_fraction_missing_file

# Smallest positive float -- used to avoid log(0) and division by zero
EPSILON = sys.float_info.min


class ExODTModel:
    """Species-tree / gene-tree reconciliation model (undated variant).

    This class contains the description of a species tree with parameters of
    duplication, transfer and loss.  Given an ``ApproxPosterior`` object it can
    compute the probability of the approximate posterior under the DTL model
    and sample reconciled gene trees.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self):
        # Parameters --------------------------------------------------
        self.scalar_parameter: dict[str, float] = {}
        self.vector_parameter: dict[str, list[float]] = {}
        self.string_parameter: dict[str, str] = {}

        # Species tree topology ---------------------------------------
        self.father: dict[int, int] = {}          # node_id -> parent_id
        self.daughter: dict[int, int] = {}         # node_id -> left child
        self.son: dict[int, int] = {}              # node_id -> right child
        self.extant_species: dict[int, str] = {}   # leaf_id -> species name

        self.last_branch: int = 0
        self.last_leaf: int = 0
        self.last_rank: int = 0
        self.root_i: int = -1

        # Node objects keyed by name / branch id ----------------------
        self._name_node: dict[str, object] = {}
        self._node_name: dict[int, str] = {}       # node.id -> name
        self._node_ids: dict[int, int] = {}         # newick-node.id -> branch id
        self._id_nodes: dict[int, object] = {}      # branch id -> newick Node

        # Undated model arrays ----------------------------------------
        self.fm: list[float] = []
        self.uE: list[float] = []
        self.mPTE: float = 0.0
        self.mPTE_ancestral_correction: list[float] = []
        self.uq: list[list[float]] = []
        self.mPTuq: list[float] = []
        self.mPTuq_ancestral_correction: list[list[float]] = []

        self.PD: list[float] = []
        self.wT: list[float] = []
        self.PL: list[float] = []
        self.PS: list[float] = []
        self.rmD: list[float] = []
        self.rmT: list[float] = []
        self.rmL: list[float] = []
        self.tau_norm: list[float] = []

        self.ancestors: list[list[int]] = []
        self.below: dict[int, dict[int, int]] = {}
        self.ancestral: dict[int, dict[int, int]] = {}

        # ALE pointer and gene-clade bookkeeping ---------------------
        self.ale_pointer = None
        self.g_ids: list[int] = []
        self.g_id_sizes: list[int] = []
        self.g_id2i: dict[int, int] = {}
        self.gid_sps: dict[int, str] = {}

        # Event tracking ----------------------------------------------
        self.MLRec_events: dict[str, float] = defaultdict(float)
        self.branch_counts: dict[str, list[float]] = {}
        self.T_to_from: list[list[float]] = []
        self.Ttokens: list[str] = []

        # Rank-to-label map (bootstrap values) ------------------------
        self.rank2label: dict[int, int] = {}

        # Set default parameters (mirrors C++ constructor) ------------
        self.string_parameter["BOOTSTRAP_LABELS"] = "no"
        self.string_parameter["gene_name_separators"] = "_@"
        self.scalar_parameter["species_field"] = 0
        self.scalar_parameter["event_node"] = 0
        self.scalar_parameter["min_bip_count"] = -1
        self.scalar_parameter["min_branch_lenghts"] = 0
        self.scalar_parameter["stem_length"] = 1
        self.scalar_parameter["D"] = 3
        self.scalar_parameter["grid_delta_t"] = 0.005
        self.scalar_parameter["min_D"] = 3
        self.scalar_parameter["DD"] = 10
        self.scalar_parameter["O_R"] = 1.0
        self.scalar_parameter["seq_beta"] = 1.0
        self.scalar_parameter["undatedBL"] = 0
        self.scalar_parameter["reldate"] = 0
        self.scalar_parameter["root_BL"] = 0

    # ------------------------------------------------------------------
    # Tree height helper (average of the two children, like the C++)
    # ------------------------------------------------------------------

    def _height(self, node) -> float:
        """Compute height of a node: 0 for leaves, otherwise average of
        (child branch length + child height) over both children."""
        if is_leaf(node):
            return 0.0
        sons = node.children
        h0 = (sons[0].branch_length or 0.0) + self._height(sons[0])
        h1 = (sons[1].branch_length or 0.0) + self._height(sons[1])
        return 0.5 * (h0 + h1)

    # ------------------------------------------------------------------
    # construct_undated
    # ------------------------------------------------------------------

    def construct_undated(self, s_string: str, fraction_missing_file: str = ""):
        """Construct the undated reconciliation model from a Newick species
        tree string and (optionally) a fraction-missing file.

        Mirrors ``exODT_model::construct_undated`` in *undated.cpp*.
        """
        # Reset maps
        self.daughter.clear()
        self.son.clear()
        self._name_node.clear()
        self._node_name.clear()
        self._node_ids.clear()
        self._id_nodes.clear()

        self.string_parameter["S_un"] = s_string

        # Parse the species tree
        root = parse_newick(s_string)
        all_nodes = get_all_nodes(root)  # post-order

        # Build name_node / node_name maps
        for nd in all_nodes:
            if is_leaf(nd):
                self._name_node[nd.name] = nd
                self._node_name[nd.id] = nd.name
            else:
                leaf_names = sorted(lf.name for lf in get_leaves(nd))
                name = ".".join(leaf_names) + "."
                self._name_node[name] = nd
                self._node_name[nd.id] = name

        # Register leaves (alphabetical order, as C++ iterates std::map)
        self.last_branch = 0
        self.last_leaf = 0

        self.vector_parameter["BL_rate_multiplier"] = []
        self.vector_parameter["rate_multiplier_tau_to"] = []
        self.vector_parameter["rate_multiplier_tau_from"] = []
        self.vector_parameter["rate_multiplier_delta"] = []
        self.vector_parameter["rate_multiplier_lambda"] = []
        self.vector_parameter["rate_multiplier_O"] = []
        self.wT = []
        self.rmD = []
        self.rmT = []
        self.rmL = []

        saw = set()  # set of newick-node ids already registered

        # Sort leaf names alphabetically (std::map ordering)
        leaf_names_sorted = sorted(
            name for name, nd in self._name_node.items() if is_leaf(nd)
        )

        for lname in leaf_names_sorted:
            node = self._name_node[lname]
            self.extant_species[self.last_branch] = node.name
            self._node_ids[node.id] = self.last_branch
            self._id_nodes[self.last_branch] = node
            self.last_branch += 1
            self.last_leaf += 1
            saw.add(node.id)
            # Leaves have no children
            self.daughter[self.last_branch] = -1
            self.son[self.last_branch] = -1
            self.vector_parameter["BL_rate_multiplier"].append(
                node.branch_length if node.branch_length is not None else 0.0
            )
            self.vector_parameter["rate_multiplier_tau_to"].append(1.0)
            self.vector_parameter["rate_multiplier_tau_from"].append(1.0)
            self.wT.append(1.0)
            self.rmD.append(1.0)
            self.rmT.append(1.0)
            self.rmL.append(1.0)
            self.vector_parameter["rate_multiplier_delta"].append(1.0)
            self.vector_parameter["rate_multiplier_lambda"].append(1.0)
            self.vector_parameter["rate_multiplier_O"].append(1.0)

        # Ad-hoc post-order: propagate from leaves upward
        next_generation = [
            self._name_node[lname] for lname in leaf_names_sorted
        ]

        while next_generation:
            new_generation = []
            for node in next_generation:
                if node.parent is not None:
                    parent = node.parent
                    sons = parent.children
                    sister = sons[1] if sons[0] is node else sons[0]
                    if parent.id not in self._node_ids and sister.id in saw:
                        self._node_ids[parent.id] = self.last_branch
                        self._id_nodes[self.last_branch] = parent
                        # Note: C++ stores node->getDistanceToFather() which
                        # is the branch length of the *child* that triggered
                        # this registration.  Replicating that behaviour.
                        self.vector_parameter["BL_rate_multiplier"].append(
                            node.branch_length
                            if node.branch_length is not None
                            else 0.0
                        )
                        self.vector_parameter["rate_multiplier_tau_to"].append(1.0)
                        self.vector_parameter["rate_multiplier_tau_from"].append(1.0)
                        self.wT.append(1.0)
                        self.rmD.append(1.0)
                        self.rmT.append(1.0)
                        self.rmL.append(1.0)
                        self.vector_parameter["rate_multiplier_delta"].append(1.0)
                        self.vector_parameter["rate_multiplier_lambda"].append(1.0)
                        self.vector_parameter["rate_multiplier_O"].append(1.0)
                        self.last_branch += 1
                        saw.add(parent.id)
                        new_generation.append(parent)
            next_generation = new_generation

        # Build ``below`` matrix:  below[e][f] = 1 iff
        #   height(father(e)) < height(f)
        self.below.clear()
        for e in range(self.last_branch - 1):
            self.below[e] = {}
            node_e = self._id_nodes[e]
            h_father_e = self._height(node_e.parent) if node_e.parent is not None else 0.0
            for f in range(self.last_branch - 1):
                node_f = self._id_nodes[f]
                h_f = self._height(node_f)
                self.below[e][f] = 1 if h_father_e < h_f else 0

        # Extra rate-multiplier entries for the stem above root
        # (index == last_branch, which is one past the last registered)
        # In C++ this writes to vector_parameter["BL_rate_multiplier"][last_branch]
        # which extends the vector by one.
        bl_vec = self.vector_parameter["BL_rate_multiplier"]
        if len(bl_vec) <= self.last_branch:
            bl_vec.append(self.scalar_parameter.get("root_BL", 0.0))
        else:
            bl_vec[self.last_branch] = self.scalar_parameter.get("root_BL", 0.0)
        self.vector_parameter["rate_multiplier_tau_to"].append(1.0)
        self.vector_parameter["rate_multiplier_tau_from"].append(1.0)
        self.wT.append(1.0)
        self.vector_parameter["rate_multiplier_delta"].append(1.0)
        self.vector_parameter["rate_multiplier_lambda"].append(1.0)
        self.vector_parameter["rate_multiplier_O"].append(1.0)

        # Build ``ancestral`` and ``ancestors`` ----------------------------
        self.ancestors = [[] for _ in range(self.last_branch)]
        self.ancestral = {}
        for e in range(self.last_branch):
            self.ancestral[e] = {}
            for f in range(self.last_branch):
                self.ancestral[e][f] = 0

        for nd in all_nodes:
            if nd.id not in self._node_ids:
                continue
            e = self._node_ids[nd.id]
            walker = nd
            while True:
                f = self._node_ids[walker.id]
                if not self.ancestral[e][f]:
                    self.ancestors[e].append(f)
                self.ancestral[e][f] = 1
                if walker.parent is None:
                    break
                walker = walker.parent

        # If reldate mode, also mark ``below`` entries as ancestral
        if self.scalar_parameter.get("reldate", 0):
            for e in range(self.last_branch):
                for f in range(self.last_branch):
                    if self.below.get(e, {}).get(f, 0) == 1:
                        if not self.ancestral[e][f]:
                            self.ancestors[e].append(f)
                        self.ancestral[e][f] = 1

        # Set daughter / son for internal nodes
        for name, nd in self._name_node.items():
            if not is_leaf(nd):
                sons = nd.children
                self.daughter[self._node_ids[nd.id]] = self._node_ids[sons[0].id]
                self.son[self._node_ids[nd.id]] = self._node_ids[sons[1].id]

        # Generate S_with_ranks: species tree with branch IDs as labels
        # Mirrors C++: node->setBranchProperty("ID", rank) then treeToParenthesis
        def _to_newick_with_ranks(nd):
            if is_leaf(nd):
                return nd.name
            child_strs = [_to_newick_with_ranks(c) for c in nd.children]
            rank = self._node_ids[nd.id]
            return "(" + ",".join(child_strs) + ")" + str(rank)
        self.string_parameter["S_with_ranks"] = _to_newick_with_ranks(root) + ";"

        # Initialise branch_counts
        count_keys = [
            "Os", "Ds", "Ts", "Tfroms", "Ls", "count",
            "presence", "saw", "O_LL", "copies", "singleton",
        ]
        for key in count_keys:
            self.branch_counts[key] = [0.0] * self.last_branch

        # Initialise T_to_from matrix
        self.T_to_from = [[0.0] * self.last_branch for _ in range(self.last_branch)]

        # last_rank mirrors C++
        self.last_rank = self.last_branch

        # Set N=1 (undated does not use population size)
        self.set_model_parameter("N", 1.0)

        # Fraction missing
        self.vector_parameter["fraction_missing"] = [0.0] * self.last_leaf
        if fraction_missing_file:
            frac_map = read_fraction_missing_file(fraction_missing_file)
            idx = 0
            for lname in leaf_names_sorted:
                node = self._name_node[lname]
                species = node.name
                if species in frac_map:
                    self.vector_parameter["fraction_missing"][idx] = frac_map[species]
                idx += 1

    # ------------------------------------------------------------------
    # set_model_parameter (three overloads merged into one)
    # ------------------------------------------------------------------

    def set_model_parameter(self, name: str, value):
        """Set a model parameter.

        *value* can be a ``str``, a ``float``/``int``, or a ``list[float]``.
        Mirrors the three C++ overloads.
        """
        if isinstance(value, str):
            self.string_parameter[name] = value
            return

        if isinstance(value, (list, tuple)):
            self._set_model_parameter_vector(name, list(value))
            return

        # Scalar
        value = float(value)
        if name in ("delta", "tau", "lambda"):
            N = self.vector_parameter.get("N", [1.0])[0]
            self.vector_parameter[name] = []
            for _branch in range(self.last_branch):
                if name == "tau":
                    self.vector_parameter[name].append(value / N)
                else:
                    self.vector_parameter[name].append(value)
            if name == "tau":
                self.scalar_parameter[name + "_avg"] = value / N
            else:
                self.scalar_parameter[name + "_avg"] = value
        elif name in ("N", "Delta_bar", "Lambda_bar"):
            self.vector_parameter[name] = [value] * self.last_rank
        else:
            self.scalar_parameter[name] = value

    def _set_model_parameter_vector(self, name: str, value_vector: list[float]):
        if name in ("delta", "tau", "lambda"):
            N = self.vector_parameter.get("N", [1.0])[0]
            self.vector_parameter[name] = []
            avg = 0.0
            for branch in range(self.last_branch):
                if name == "tau":
                    v = value_vector[branch] / N
                else:
                    v = value_vector[branch]
                self.vector_parameter[name].append(v)
                avg += v
            self.scalar_parameter[name + "_avg"] = avg / max(self.last_branch, 1)
        else:
            self.vector_parameter[name] = []
            for rank in range(self.last_rank):
                self.vector_parameter[name].append(value_vector[rank])

    # ------------------------------------------------------------------
    # calculate_undatedEs
    # ------------------------------------------------------------------

    def calculate_undatedEs(self):
        """Compute extinction probabilities for each branch (undated model).

        Mirrors ``exODT_model::calculate_undatedEs`` in *undated.cpp*.
        """
        self.uE = []
        self.fm = []
        self.mPTE_ancestral_correction = []
        self.PD = []
        self.tau_norm = []
        self.PL = []
        self.PS = []

        use_bl = bool(self.scalar_parameter.get("undatedBL", 0))

        # Compute raw rates per branch
        for e in range(self.last_branch):
            if use_bl:
                self.wT[e] = (
                    self.vector_parameter["rate_multiplier_tau_to"][e]
                    * self.vector_parameter["BL_rate_multiplier"][e]
                )
                self.rmT[e] = (
                    self.vector_parameter["tau"][e]
                    * self.vector_parameter["rate_multiplier_tau_from"][e]
                    * self.vector_parameter["BL_rate_multiplier"][e]
                )
                self.rmD[e] = (
                    self.vector_parameter["delta"][e]
                    * self.vector_parameter["rate_multiplier_delta"][e]
                    * self.vector_parameter["BL_rate_multiplier"][e]
                )
                self.rmL[e] = (
                    self.vector_parameter["lambda"][e]
                    * self.vector_parameter["rate_multiplier_lambda"][e]
                    * self.vector_parameter["BL_rate_multiplier"][e]
                )
            else:
                self.wT[e] = self.vector_parameter["rate_multiplier_tau_to"][e]
                self.rmT[e] = (
                    self.vector_parameter["tau"][e]
                    * self.vector_parameter["rate_multiplier_tau_from"][e]
                )
                self.rmD[e] = (
                    self.vector_parameter["delta"][e]
                    * self.vector_parameter["rate_multiplier_delta"][e]
                )
                self.rmL[e] = (
                    self.vector_parameter["lambda"][e]
                    * self.vector_parameter["rate_multiplier_lambda"][e]
                )

        tau_sum = sum(self.wT[f] for f in range(self.last_branch))

        for e in range(self.last_branch):
            P_D = self.rmD[e]
            P_T = self.rmT[e]
            P_L = self.rmL[e]
            P_S = 1.0

            total = P_D + P_T + P_L + P_S
            P_D /= total
            P_T /= total
            P_L /= total
            P_S /= total

            self.PD.append(P_D)

            # tau_norm[e] = (tau_sum - sum_wT_ancestors) / P_T
            tn = tau_sum
            for f in self.ancestors[e]:
                tn -= self.wT[f]
            if P_T > 0:
                tn /= P_T
            else:
                tn = 1.0  # avoid division by zero
            self.tau_norm.append(tn)

            self.PL.append(P_L)
            self.PS.append(P_S)
            self.uE.append(0.0)

            if e < self.last_leaf:
                self.fm.append(self.vector_parameter["fraction_missing"][e])
            else:
                self.fm.append(0.0)

            self.mPTE_ancestral_correction.append(0.0)

        # Iterative computation of extinction probabilities (4 iterations)
        self.mPTE = 0.0
        for iteration in range(4):
            new_mPTE = 0.0

            if iteration > 0:
                for e in range(self.last_branch):
                    self.mPTE_ancestral_correction[e] = 0.0
                    for f in self.ancestors[e]:
                        self.mPTE_ancestral_correction[e] += self.wT[f] * self.uE[f]

            for e in range(self.last_branch):
                if e < self.last_leaf:
                    # Leaf: no speciation
                    self.uE[e] = (
                        self.PL[e]
                        + self.PD[e] * self.uE[e] * self.uE[e]
                        + self.uE[e]
                        * (self.mPTE - self.mPTE_ancestral_correction[e])
                        / self.tau_norm[e]
                    )
                else:
                    f = self.daughter[e]
                    g = self.son[e]
                    self.uE[e] = (
                        self.PL[e]
                        + self.PS[e] * self.uE[f] * self.uE[g]
                        + self.PD[e] * self.uE[e] * self.uE[e]
                        + self.uE[e]
                        * (self.mPTE - self.mPTE_ancestral_correction[e])
                        / self.tau_norm[e]
                    )
                new_mPTE += self.wT[e] * self.uE[e]
            self.mPTE = new_mPTE

        # One more update incorporating fraction missing
        new_mPTE = 0.0
        for e in range(self.last_branch):
            if e < self.last_leaf:
                self.uE[e] = (1.0 - self.fm[e]) * self.uE[e] + self.fm[e]
            else:
                f = self.daughter[e]
                g = self.son[e]
                self.uE[e] = (
                    self.PL[e]
                    + self.PS[e] * self.uE[f] * self.uE[g]
                    + self.PD[e] * self.uE[e] * self.uE[e]
                    + self.uE[e]
                    * (self.mPTE - self.mPTE_ancestral_correction[e])
                    / self.tau_norm[e]
                )
            new_mPTE += self.wT[e] * self.uE[e]
        self.mPTE = new_mPTE

    # ------------------------------------------------------------------
    # Gene-name to species-name mapping helper
    # ------------------------------------------------------------------

    def _gene_to_species(self, gene_name: str) -> str:
        """Extract species name from a gene name using the configured
        separators and species field."""
        seps = self.string_parameter.get("gene_name_separators", "_@")
        # Build regex pattern that splits on any separator character
        pattern = "[" + re.escape(seps) + "]+"
        tokens = re.split(pattern, gene_name)
        field = int(self.scalar_parameter.get("species_field", 0))
        if field == -1:
            return "_".join(tokens[1:])
        return tokens[field]

    # ------------------------------------------------------------------
    # pun  --  compute log-likelihood under undated model
    # ------------------------------------------------------------------

    def pun(self, ale, verbose: bool = False, no_T: bool = False) -> float:
        """Compute the (non-log) likelihood of *ale* under the undated model.

        Returns the likelihood value (not log-transformed).
        Mirrors ``exODT_model::pun`` in *undated.cpp*.

        Parameters
        ----------
        ale : ApproxPosterior
            The approximate posterior object.
        verbose : bool
            If True, print gene-to-species mapping.
        no_T : bool
            If True, disable transfers.
        """
        survive = 0.0
        root_sum = 0.0
        O_norm = 0.0
        self.mPTuq_ancestral_correction = []
        self.uq = []
        self.mPTuq = []
        self.ale_pointer = ale

        # Build ordered clade list (small to large)
        self.g_ids = []
        self.g_id_sizes = []
        for size in sorted(ale.size_ordered_bips.keys()):
            for g_id in ale.size_ordered_bips[size]:
                self.g_ids.append(g_id)
                self.g_id_sizes.append(size)
        # Root bipartition (-1) handled separately
        self.g_ids.append(-1)
        self.g_id_sizes.append(ale.Gamma_size)
        self.root_i = len(self.g_ids) - 1

        # Gene <-> species mapping
        self.gid_sps.clear()
        species_set = set(self.extant_species.values())
        if verbose:
            print("\nGene\t:\tSpecies")
        for idx in range(len(self.g_ids)):
            g_id = self.g_ids[idx]
            if self.g_id_sizes[idx] == 1:
                # Find the single leaf id in the bitset
                leaf_id = None
                for bit_i in range(ale.Gamma_size + 1):
                    if (ale.id_sets[g_id] >> bit_i) & 1:
                        leaf_id = bit_i
                        break
                gene_name = ale.id_leaves[leaf_id]
                species_name = self._gene_to_species(gene_name)
                self.gid_sps[g_id] = species_name
                if species_name not in species_set:
                    print(
                        f"Error: gene name {gene_name} is associated to "
                        f"species name {species_name} that cannot be found "
                        f"in the species tree.",
                        file=sys.stderr,
                    )
                    sys.exit(-1)
                if verbose:
                    print(f"{gene_name}\t:\t{species_name}")

        # Build g_id2i and initialise uq / mPTuq arrays
        self.g_id2i = {}
        for i in range(len(self.g_ids)):
            g_id = self.g_ids[i]
            self.g_id2i[g_id] = i

            if i >= len(self.uq):
                self.uq.append([0.0] * self.last_branch)
                self.mPTuq_ancestral_correction.append([0.0] * self.last_branch)
                self.mPTuq.append(0.0)
            else:
                self.mPTuq[i] = 0.0
                for e in range(self.last_branch):
                    self.uq[i][e] = 0.0
                    self.mPTuq_ancestral_correction[i][e] = 0.0

        # Main iteration (4 rounds for convergence)
        for _iteration in range(4):
            for i in range(len(self.g_ids)):
                new_mPTuq = 0.0

                g_id = self.g_ids[i]
                is_a_leaf = self.g_id_sizes[i] == 1

                # Build partition lists
                gp_is = []
                gpp_is = []
                p_part = []

                if g_id != -1:
                    # Normal clade: iterate over Dip_counts
                    for (gp_id, gpp_id), count in ale.Dip_counts[g_id].items():
                        gp_is.append(self.g_id2i[gp_id])
                        gpp_is.append(self.g_id2i[gpp_id])
                        if ale.Bip_counts.get(g_id, 0) <= self.scalar_parameter["min_bip_count"]:
                            p_part.append(0.0)
                        else:
                            p_part.append(
                                ale.p_dip(g_id, gp_id, gpp_id)
                                ** self.scalar_parameter["seq_beta"]
                            )
                else:
                    # Root bipartition: enumerate all bipartitions
                    bip_parts = {}
                    for gp_id in ale.Bip_counts:
                        gamma = ale.id_sets[gp_id]
                        not_gamma = gamma ^ ale.Gamma
                        if not_gamma not in ale.set_ids:
                            continue
                        gpp_id = ale.set_ids[not_gamma]
                        key = frozenset([gp_id, gpp_id])
                        bip_parts[key] = 1
                    for key in bip_parts:
                        parts = sorted(key)
                        gp_id = parts[0]
                        gpp_id = parts[1] if len(parts) > 1 else parts[0]
                        gp_is.append(self.g_id2i[parts[0]])
                        gpp_is.append(self.g_id2i[parts[1]] if len(parts) > 1 else self.g_id2i[parts[0]])
                        bip_count = ale.Bip_counts.get(gp_id, 0)
                        if bip_count <= self.scalar_parameter.get("min_bip_count", -1) and not ale.Gamma_size < 4:
                            p_part.append(0.0)
                        else:
                            p_part.append(
                                ale.p_bip(gp_id)
                                ** self.scalar_parameter["seq_beta"]
                            )

                # Inner loop over branches
                for e in range(self.last_branch):
                    uq_sum = 0.0

                    # S-leaf and G-leaf match
                    if (
                        e < self.last_leaf
                        and is_a_leaf
                        and self.extant_species[e] == self.gid_sps.get(g_id, "")
                    ):
                        uq_sum += self.PS[e] * 1.0

                    # G internal: enumerate partitions
                    if not is_a_leaf:
                        n_parts = len(gp_is)
                        for pi in range(n_parts):
                            gp_i = gp_is[pi]
                            gpp_i = gpp_is[pi]
                            pp = p_part[pi]

                            if not (e < self.last_leaf):
                                f = self.daughter[e]
                                g = self.son[e]
                                # Speciation event
                                uq_sum += self.PS[e] * (
                                    self.uq[gp_i][f] * self.uq[gpp_i][g]
                                    + self.uq[gp_i][g] * self.uq[gpp_i][f]
                                ) * pp

                            # Duplication event
                            uq_sum += self.PD[e] * (
                                self.uq[gp_i][e] * self.uq[gpp_i][e]
                            ) * pp

                            # Transfer event
                            if not no_T:
                                uq_sum += (
                                    self.uq[gp_i][e]
                                    * (
                                        self.mPTuq[gpp_i]
                                        - self.mPTuq_ancestral_correction[gpp_i][e]
                                    )
                                    / self.tau_norm[e]
                                    + self.uq[gpp_i][e]
                                    * (
                                        self.mPTuq[gp_i]
                                        - self.mPTuq_ancestral_correction[gp_i][e]
                                    )
                                    / self.tau_norm[e]
                                ) * pp

                    # SL (speciation-loss) event
                    if not (e < self.last_leaf):
                        f = self.daughter[e]
                        g = self.son[e]
                        uq_sum += self.PS[e] * (
                            self.uq[i][f] * self.uE[g]
                            + self.uq[i][g] * self.uE[f]
                        )

                    # DL (duplication-loss) event
                    uq_sum += self.PD[e] * (self.uq[i][e] * self.uE[e] * 2)

                    # TL (transfer-loss) event
                    if not no_T:
                        uq_sum += (
                            (
                                self.mPTuq[i]
                                - self.mPTuq_ancestral_correction[i][e]
                            )
                            / self.tau_norm[e]
                            * self.uE[e]
                            + self.uq[i][e]
                            * (self.mPTE - self.mPTE_ancestral_correction[e])
                            / self.tau_norm[e]
                        )

                    if uq_sum < EPSILON:
                        uq_sum = EPSILON
                    self.uq[i][e] = uq_sum
                    new_mPTuq += self.wT[e] * uq_sum

                    # Update ancestral correction for this clade
                    self.mPTuq_ancestral_correction[i][e] = 0.0
                    for f in self.ancestors[e]:
                        self.mPTuq_ancestral_correction[i][e] += self.wT[f] * uq_sum

                self.mPTuq[i] = new_mPTuq

            # Root summation
            survive = 0.0
            root_sum = 0.0
            O_norm = 0.0

            # Check for single-origination mode
            single_O = any(
                self.vector_parameter["rate_multiplier_O"][e] < 0
                for e in range(self.last_branch)
            )
            if single_O:
                for e in range(self.last_branch):
                    if self.vector_parameter["rate_multiplier_O"][e] > 0:
                        self.vector_parameter["rate_multiplier_O"][e] = 0.0
                    else:
                        self.vector_parameter["rate_multiplier_O"][e] = 1.0

            for e in range(self.last_branch):
                O_p = self.vector_parameter["rate_multiplier_O"][e]
                if e == (self.last_branch - 1) and O_p == 1.0:
                    O_p = self.scalar_parameter["O_R"]
                O_norm += O_p
                root_sum += self.uq[self.root_i][e] * O_p
                survive += 1.0 - self.uE[e]

            for e in range(self.last_branch):
                O_p = self.vector_parameter["rate_multiplier_O"][e]
                if e == (self.last_branch - 1) and O_p == 1.0:
                    O_p = self.scalar_parameter["O_R"]
                if O_p > 0 and O_norm > 0:
                    self.branch_counts["O_LL"][e] = (
                        math.log(max(self.uq[self.root_i][e], EPSILON))
                        + math.log(max(O_p, EPSILON))
                        - math.log(O_norm)
                    )

        # Return non-log likelihood (C++ returns this directly)
        if survive <= 0 or O_norm <= 0:
            return EPSILON
        return root_sum / survive / O_norm * self.last_branch

    # ------------------------------------------------------------------
    # sample_undated  --  sample a reconciled gene tree (entry point)
    # ------------------------------------------------------------------

    def sample_undated(self, no_T: bool = False) -> str:
        """Sample a reconciled gene tree under the undated model.

        Returns a Newick-like string describing the reconciled tree.
        Mirrors ``exODT_model::sample_undated()`` (no args) in *undated.cpp*.
        """
        import random

        r = random.random()

        root_sum = 0.0
        O_norm = 0.0
        for e in range(self.last_branch):
            self.branch_counts["saw"][e] = 0
            O_p = self.vector_parameter["rate_multiplier_O"][e]
            if e == (self.last_branch - 1) and O_p == 1.0:
                O_p = self.scalar_parameter["O_R"]
            O_norm += O_p
            root_sum += self.uq[len(self.g_ids) - 1][e] * O_p + EPSILON

        root_resum = 0.0
        for e in range(self.last_branch):
            O_p = self.vector_parameter["rate_multiplier_O"][e]
            if e == (self.last_branch - 1) and O_p == 1.0:
                O_p = self.scalar_parameter["O_R"]
            root_resum += self.uq[self.root_i][e] * O_p + EPSILON
            if r * root_sum < root_resum:
                self.register_O(e)
                return self._sample_undated_recursive(e, self.root_i, "O", "", no_T) + ";"

        return "-!=-"

    # ------------------------------------------------------------------
    # _sample_undated_recursive  --  recursive sampling
    # ------------------------------------------------------------------

    def _sample_undated_recursive(
        self,
        e: int,
        i: int,
        last_event: str,
        branch_string: str = "",
        no_T: bool = False,
    ) -> str:
        """Recursively sample events for a reconciled gene tree.

        Mirrors ``exODT_model::sample_undated(int e, int i, ...)`` in
        *undated.cpp*.
        """
        import random

        r = random.random()
        ale = self.ale_pointer

        is_a_leaf = False
        g_id = self.g_ids[i]
        if self.g_id_sizes[i] == 1:
            is_a_leaf = True

        # Branch length
        if g_id in ale.Bip_counts and ale.Bip_counts[g_id] > 0:
            bl = max(
                ale.Bip_bls[g_id] / ale.Bip_counts[g_id],
                self.scalar_parameter["min_branch_lenghts"],
            )
        else:
            bl = max(
                ale.Bip_bls.get(g_id, 0.0) / ale.observations,
                self.scalar_parameter["min_branch_lenghts"],
            )
        branch_length = str(bl)

        # Build partition lists (same logic as pun)
        gp_is = []
        gpp_is = []
        p_part = []

        if g_id != -1:
            for (gp_id, gpp_id), count in ale.Dip_counts[g_id].items():
                gp_is.append(self.g_id2i[gp_id])
                gpp_is.append(self.g_id2i[gpp_id])
                if ale.Bip_counts.get(g_id, 0) <= self.scalar_parameter["min_bip_count"]:
                    p_part.append(0.0)
                else:
                    p_part.append(
                        ale.p_dip(g_id, gp_id, gpp_id)
                        ** self.scalar_parameter["seq_beta"]
                    )
        else:
            bip_parts = {}
            for gp_id in ale.Bip_counts:
                gamma = ale.id_sets[gp_id]
                not_gamma = gamma ^ ale.Gamma
                if not_gamma not in ale.set_ids:
                    continue
                gpp_id = ale.set_ids[not_gamma]
                key = frozenset([gp_id, gpp_id])
                bip_parts[key] = 1
            for key in bip_parts:
                parts = sorted(key)
                gp_is.append(self.g_id2i[parts[0]])
                gpp_is.append(self.g_id2i[parts[1]] if len(parts) > 1 else self.g_id2i[parts[0]])
                gp_id = parts[0]
                bip_count = ale.Bip_counts.get(gp_id, 0)
                if bip_count <= self.scalar_parameter.get("min_bip_count", -1) and not ale.Gamma_size < 4:
                    p_part.append(0.0)
                else:
                    p_part.append(
                        ale.p_bip(gp_id) ** self.scalar_parameter["seq_beta"]
                    )

        # Compute total uq_sum for sampling
        uq_sum = 0.0

        # S-leaf and G-leaf
        if e < self.last_leaf and is_a_leaf and self.extant_species[e] == self.gid_sps.get(g_id, ""):
            uq_sum += self.PS[e] * 1.0 + EPSILON

        # G internal
        if not is_a_leaf:
            n_parts = len(gp_is)
            for pi in range(n_parts):
                gp_i = gp_is[pi]
                gpp_i = gpp_is[pi]
                pp = p_part[pi]
                if not (e < self.last_leaf):
                    f = self.daughter[e]
                    g = self.son[e]
                    uq_sum += self.PS[e] * self.uq[gp_i][f] * self.uq[gpp_i][g] * pp + EPSILON
                    uq_sum += self.PS[e] * self.uq[gp_i][g] * self.uq[gpp_i][f] * pp + EPSILON
                # D event
                uq_sum += self.PD[e] * (self.uq[gp_i][e] * self.uq[gpp_i][e] * 2) * pp + EPSILON
                # T event
                for f in range(self.last_branch):
                    if not self.ancestral[e][f] and not no_T:
                        uq_sum += self.uq[gp_i][e] * (self.wT[f] / self.tau_norm[e]) * self.uq[gpp_i][f] * pp + EPSILON
                        uq_sum += self.uq[gpp_i][e] * (self.wT[f] / self.tau_norm[e]) * self.uq[gp_i][f] * pp + EPSILON

        if not (e < self.last_leaf):
            f = self.daughter[e]
            g = self.son[e]
            uq_sum += self.PS[e] * self.uq[i][f] * self.uE[g] + EPSILON
            uq_sum += self.PS[e] * self.uq[i][g] * self.uE[f] + EPSILON

        uq_sum += self.PD[e] * (self.uq[i][e] * self.uE[e] * 2) + EPSILON

        for f in range(self.last_branch):
            if not self.ancestral[e][f] and not no_T:
                uq_sum += (self.wT[f] / self.tau_norm[e]) * self.uq[i][f] * self.uE[e] + EPSILON
                uq_sum += (self.wT[f] / self.tau_norm[e]) * self.uE[f] * self.uq[i][e] + EPSILON

        # Branch label
        if not (e < self.last_leaf):
            estr = str(e)
        else:
            estr = self.extant_species[e]

        # Now resample: walk through events in same order, accumulating
        uq_resum = 0.0

        # S-leaf and G-leaf
        if e < self.last_leaf and is_a_leaf and self.extant_species[e] == self.gid_sps.get(g_id, ""):
            uq_resum += self.PS[e] * 1.0 + EPSILON
            if r * uq_sum < uq_resum:
                self.register_leafu(e, last_event)
                return ale.set2name(ale.id_sets[g_id]) + branch_string + ":" + branch_length

        # G internal
        if not is_a_leaf:
            n_parts = len(gp_is)
            for pi in range(n_parts):
                gp_i = gp_is[pi]
                gpp_i = gpp_is[pi]
                pp = p_part[pi]

                if not (e < self.last_leaf):
                    f = self.daughter[e]
                    g = self.son[e]

                    # S event (gp->f, gpp->g)
                    uq_resum += self.PS[e] * self.uq[gp_i][f] * self.uq[gpp_i][g] * pp + EPSILON
                    if r * uq_sum < uq_resum:
                        self.register_Su(e, last_event)
                        return (
                            "("
                            + self._sample_undated_recursive(f, gp_i, "S", "", no_T)
                            + ","
                            + self._sample_undated_recursive(g, gpp_i, "S", "", no_T)
                            + ")."
                            + estr
                            + branch_string
                            + ":"
                            + branch_length
                        )

                    # S event (gp->g, gpp->f)
                    uq_resum += self.PS[e] * self.uq[gp_i][g] * self.uq[gpp_i][f] * pp + EPSILON
                    if r * uq_sum < uq_resum:
                        self.register_Su(e, last_event)
                        return (
                            "("
                            + self._sample_undated_recursive(g, gp_i, "S", "", no_T)
                            + ","
                            + self._sample_undated_recursive(f, gpp_i, "S", "", no_T)
                            + ")."
                            + estr
                            + branch_string
                            + ":"
                            + branch_length
                        )

                # D event
                uq_resum += self.PD[e] * (self.uq[gp_i][e] * self.uq[gpp_i][e] * 2) * pp + EPSILON
                if r * uq_sum < uq_resum:
                    self.register_D(e)
                    return (
                        "("
                        + self._sample_undated_recursive(e, gp_i, "D", "", no_T)
                        + ","
                        + self._sample_undated_recursive(e, gpp_i, "D", "", no_T)
                        + ").D@"
                        + estr
                        + branch_string
                        + ":"
                        + branch_length
                    )

                # T event
                for f in range(self.last_branch):
                    if not self.ancestral[e][f] and not no_T:
                        fstr = self.extant_species[f] if f < self.last_leaf else str(f)

                        uq_resum += self.uq[gp_i][e] * (self.wT[f] / self.tau_norm[e]) * self.uq[gpp_i][f] * pp + EPSILON
                        if r * uq_sum < uq_resum:
                            self.register_Tfrom(e)
                            self.register_Tto(f)
                            self.register_T_to_from(e, f)
                            t_token = f"{estr}>{fstr}|{ale.set2name(ale.id_sets[self.g_ids[gpp_i]])}"
                            self.Ttokens.append(t_token)
                            return (
                                "("
                                + self._sample_undated_recursive(e, gp_i, "S", "", no_T)
                                + ","
                                + self._sample_undated_recursive(f, gpp_i, "T", "", no_T)
                                + ").T@"
                                + estr
                                + "->"
                                + fstr
                                + branch_string
                                + ":"
                                + branch_length
                            )

                        uq_resum += self.uq[gpp_i][e] * (self.wT[f] / self.tau_norm[e]) * self.uq[gp_i][f] * pp + EPSILON
                        if r * uq_sum < uq_resum:
                            self.register_Tfrom(e)
                            self.register_Tto(f)
                            self.register_T_to_from(e, f)
                            t_token = f"{estr}>{fstr}|{ale.set2name(ale.id_sets[self.g_ids[gp_i]])}"
                            self.Ttokens.append(t_token)
                            return (
                                "("
                                + self._sample_undated_recursive(e, gpp_i, "S", "", no_T)
                                + ","
                                + self._sample_undated_recursive(f, gp_i, "T", "", no_T)
                                + ").T@"
                                + estr
                                + "->"
                                + fstr
                                + branch_string
                                + ":"
                                + branch_length
                            )

        # SL event
        if not (e < self.last_leaf):
            f = self.daughter[e]
            g = self.son[e]
            uq_resum += self.PS[e] * self.uq[i][f] * self.uE[g] + EPSILON
            if r * uq_sum < uq_resum:
                self.register_Su(e, last_event)
                self.register_L(g)
                return self._sample_undated_recursive(f, i, "S", "." + estr + branch_string, no_T)

            uq_resum += self.PS[e] * self.uq[i][g] * self.uE[f] + EPSILON
            if r * uq_sum < uq_resum:
                self.register_Su(e, last_event)
                self.register_L(f)
                return self._sample_undated_recursive(g, i, "S", "." + estr + branch_string, no_T)

        # DL event
        uq_resum += self.PD[e] * (self.uq[i][e] * self.uE[e] * 2) + EPSILON
        if r * uq_sum < uq_resum:
            return self._sample_undated_recursive(e, i, "S", branch_string, no_T)

        # TL event
        for f in range(self.last_branch):
            if not self.ancestral[e][f] and not no_T:
                fstr = self.extant_species[f] if f < self.last_leaf else str(f)

                uq_resum += (self.wT[f] / self.tau_norm[e]) * self.uq[i][f] * self.uE[e] + EPSILON
                if r * uq_sum < uq_resum:
                    self.register_Tfrom(e)
                    self.register_Tto(f)
                    self.register_T_to_from(e, f)
                    self.register_L(e)
                    return self._sample_undated_recursive(
                        f, i, "T", ".T@" + estr + "->" + fstr + branch_string, no_T
                    )

                uq_resum += (self.wT[f] / self.tau_norm[e]) * self.uE[f] * self.uq[i][e] + EPSILON
                if r * uq_sum < uq_resum:
                    return self._sample_undated_recursive(e, i, "S", "", no_T)

        print("sum error!", file=sys.stderr)
        return "-!=-"

    # ------------------------------------------------------------------
    # Event registration helpers
    # ------------------------------------------------------------------

    def register_O(self, e: int):
        """Register an origination event on branch *e*."""
        if e > -1:
            self.branch_counts["count"][e] += 1
            self.branch_counts["Os"][e] += 1

    def register_D(self, e: int):
        """Register a duplication event on branch *e*."""
        self.MLRec_events["D"] += 1
        if e > -1:
            self.branch_counts["Ds"][e] += 1

    def register_Tto(self, e: int):
        """Register a transfer-to event on branch *e*."""
        self.MLRec_events["T"] += 1
        if e > -1:
            self.branch_counts["Ts"][e] += 1

    def register_Tfrom(self, e: int):
        """Register a transfer-from event on branch *e*."""
        if e > -1:
            self.branch_counts["Tfroms"][e] += 1

    def register_L(self, e: int):
        """Register a loss event on branch *e*."""
        self.MLRec_events["L"] += 1
        if e > -1:
            self.branch_counts["Ls"][e] += 1

    def register_S(self, e: int):
        """Register a speciation event on branch *e* (dated model)."""
        self.MLRec_events["S"] += 1
        if e > -1:
            f = self.daughter[e]
            g = self.son[e]
            self.branch_counts["copies"][e] += 1
            self.branch_counts["count"][f] += 1
            self.branch_counts["count"][g] += 1

    def register_Su(self, e: int, last_event: str):
        """Register a speciation event on branch *e* (undated model)."""
        self.MLRec_events["S"] += 1
        if e > -1:
            f = self.daughter[e]
            g = self.son[e]
            if last_event in ("S", "O"):
                self.branch_counts["singleton"][e] += 1
            self.branch_counts["copies"][e] += 1
            if self.branch_counts["saw"][e] == 0:
                self.branch_counts["presence"][e] += 1
            self.branch_counts["saw"][e] = 1
            self.branch_counts["count"][f] += 1
            self.branch_counts["count"][g] += 1

    def register_leafu(self, e: int, last_event: str):
        """Register reaching a leaf (undated model)."""
        if e > -1:
            self.branch_counts["copies"][e] += 1
            if self.branch_counts["saw"][e] == 0:
                self.branch_counts["presence"][e] += 1
            self.branch_counts["saw"][e] = 1
            if last_event in ("S", "O"):
                self.branch_counts["singleton"][e] += 1

    def register_leaf(self, e: int):
        """Register reaching a leaf (dated model)."""
        if e > -1:
            self.branch_counts["copies"][e] += 1

    def register_T_to_from(self, e: int, f: int):
        """Record a transfer from branch *e* to branch *f*."""
        self.T_to_from[e][f] += 1

    def reset_T_to_from(self):
        """Reset the transfer-direction matrix to zeros."""
        for e in range(self.last_branch):
            for f in range(self.last_branch):
                self.T_to_from[e][f] = 0

    def register_Ttoken(self, token: str):
        """Append a transfer token string."""
        self.Ttokens.append(token)

    # ------------------------------------------------------------------
    # counts_string_undated
    # ------------------------------------------------------------------

    def counts_string_undated(self, samples: float = 1.0) -> str:
        """Format event counts as a tab-separated string.

        Mirrors ``exODT_model::counts_string_undated`` in *undated.cpp*.
        """
        lines = []
        for e in range(self.last_branch):
            is_leaf_branch = e < self.last_leaf
            if is_leaf_branch:
                named = f"{self.extant_species[e]}({e})"
                prefix = "S_terminal_branch"
            else:
                named = str(e)
                prefix = "S_internal_branch"

            line = (
                f"{prefix}\t{named}\t"
                f"{self.branch_counts['Ds'][e] / samples}\t"
                f"{self.branch_counts['Ts'][e] / samples}\t"
                f"{self.branch_counts['Ls'][e] / samples}\t"
                f"{self.branch_counts['Os'][e] / samples}\t"
                f"{self.branch_counts['copies'][e] / samples}\t"
                f"{self.branch_counts['singleton'][e] / samples}\t"
                f"{self.uE[e]}\t"
                f"{self.branch_counts['presence'][e] / samples}\t"
                f"{self.branch_counts['O_LL'][e]}"
            )
            lines.append(line)
        return "\n".join(lines) + "\n"
