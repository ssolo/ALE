"""
Core approx_posterior class for the ALE Python port.

Direct port of ALE.h / ALE.cpp. Uses Python integers as bitsets
(bit i set means species i is in the set). Leaf IDs are 1-based
(bit 0 is unused).
"""

import math
import random

from .newick import parse_newick, get_leaves, get_all_nodes, is_leaf, to_newick


def double_factorial(n):
    """Compute n!! = n * (n-2) * (n-4) * ... * 1 (or 2).

    Returns 1 for n <= 0.
    """
    if n <= 0:
        return 1.0
    result = 1.0
    k = n
    while k > 0:
        result *= k
        k -= 2
    return result


def _popcount(x):
    """Count set bits in integer *x*."""
    return bin(x).count("1")


def _bits_to_list(x, size):
    """Return sorted list of set bit positions in *x* (up to bit *size*)."""
    out = []
    for i in range(1, size + 1):
        if x & (1 << i):
            out.append(i)
    return out


# ---------------------------------------------------------------------------
# Unrooting helper
# ---------------------------------------------------------------------------

def _build_unrooted_adjacency(root):
    """Build an unrooted adjacency representation from a parsed tree.

    The C++ code calls ``G->unroot()`` which, for a rooted binary tree,
    removes the root node and connects both subtrees into a trifurcation
    at one of the root's children.

    Returns
    -------
    neighbor : dict[int, list[int]]
        node-id -> list of adjacent node-ids
    node_map : dict[int, Node]
        node-id -> Node object
    """
    all_nodes = get_all_nodes(root)
    node_map = {n.id: n for n in all_nodes}
    neighbor = {}

    # Determine if the tree is rooted (root has exactly 2 children)
    if len(root.children) == 2:
        # Unroot: merge the root's two children.  The root is removed and
        # the two children become neighbours of each other in addition to
        # their existing children.
        left, right = root.children
        # Build normal adjacency for every node except root
        for n in all_nodes:
            if n is root:
                continue
            adj = []
            for c in n.children:
                adj.append(c.id)
            if n.parent is not None and n.parent is not root:
                adj.append(n.parent.id)
            neighbor[n.id] = adj

        # Connect left <-> right through the removed root
        neighbor[left.id].append(right.id)
        neighbor[right.id].append(left.id)
    else:
        # Already unrooted (or multifurcation at root)
        for n in all_nodes:
            adj = []
            for c in n.children:
                adj.append(c.id)
            if n.parent is not None:
                adj.append(n.parent.id)
            neighbor[n.id] = adj

    return neighbor, node_map


# ---------------------------------------------------------------------------
# approx_posterior
# ---------------------------------------------------------------------------

class approx_posterior:
    """Approximate posterior on phylogenetic trees via clade and conditional
    clade probabilities.  Direct port of the C++ ``approx_posterior`` class.
    """

    def __init__(self, tree_string=None):
        self.observations = 0.0
        self.constructor_string = ""
        self.alpha = 0.0
        self.beta = 0.0
        self.leaf_ids = {}       # name -> id (1-based)
        self.id_leaves = {}      # id -> name
        self.Bip_counts = {}     # bipartition_id -> count
        self.Bip_bls = {}        # bipartition_id -> sum of branch lengths
        self.Dip_counts = [{}]   # indexed by g_id; keys are (gp_id, gpp_id)
        self.set_ids = {}        # bitset -> bipartition_id
        self.id_sets = {}        # bipartition_id -> bitset
        self.Gamma = 0           # bitset with all species
        self.Gamma_size = 0
        self.K_Gamma = 0.0
        self.N_Gamma = 0.0
        self.last_leafset_id = 0
        self.set_sizes = {}      # bipartition_id -> set size
        self.size_ordered_bips = {}  # size -> [bipartition_ids]
        self.name_separator = "+"

        if tree_string is not None:
            self.constructor_string = tree_string
            self.construct(tree_string)

    # ------------------------------------------------------------------
    # construct
    # ------------------------------------------------------------------

    def construct(self, tree_string):
        """Parse a Newick tree to initialise leaf mappings and Gamma."""
        self.last_leafset_id = 0
        self.observations = 0.0

        tree_string = tree_string.strip()

        # Get leaf names -- from a tree or comma-separated list
        if tree_string.startswith("("):
            root = parse_newick(tree_string)
            leaves = get_leaves(root)
            leaf_names = sorted(n.name for n in leaves)
        else:
            leaf_names = sorted(
                s.strip() for s in tree_string.split(",") if s.strip()
            )

        # Assign 1-based IDs in alphabetical order (matching C++ std::map)
        self.leaf_ids = {}
        self.id_leaves = {}
        for idx, name in enumerate(leaf_names, start=1):
            self.leaf_ids[name] = idx
            self.id_leaves[idx] = name

        self.alpha = 0.0
        self.beta = 0.0
        self.Gamma_size = len(leaf_names)

        # Gamma bitset: bits 1..Gamma_size set, bit 0 unused
        self.Gamma = 0
        for i in range(1, self.Gamma_size + 1):
            self.Gamma |= (1 << i)

        # Number of bipartitions of Gamma
        self.K_Gamma = 2.0 ** (self.Gamma_size - 1) - 1

        # Number of unrooted trees on Gamma_size leaves
        if self.Gamma_size < 3:
            self.N_Gamma = 1.0
        else:
            self.N_Gamma = double_factorial(2 * self.Gamma_size - 5)

        # Pre-allocate Dip_counts entries (matching C++ push_back(temp))
        self.Dip_counts = [{}]
        while len(leaf_names) + 1 > len(self.Dip_counts):
            self.Dip_counts.append({})

    # ------------------------------------------------------------------
    # set2id
    # ------------------------------------------------------------------

    def set2id(self, leaf_set):
        """Return existing ID or create a new one for *leaf_set* (int bitset)."""
        if leaf_set in self.set_ids:
            sid = self.set_ids[leaf_set]
            if sid != 0:
                return sid

        self.last_leafset_id += 1
        self.set_ids[leaf_set] = self.last_leafset_id
        self.id_sets[self.last_leafset_id] = leaf_set
        # Ensure Dip_counts has enough entries
        while self.last_leafset_id >= len(self.Dip_counts):
            self.Dip_counts.append({})
        self.Bip_bls[self.last_leafset_id] = 0.0
        return self.last_leafset_id

    # ------------------------------------------------------------------
    # decompose
    # ------------------------------------------------------------------

    def decompose(self, G_string, weight=1.0):
        """Parse a tree and extract bipartitions, updating counts."""
        # Push an empty dict at the start (matching C++ behaviour)
        self.Dip_counts.append({})

        root = parse_newick(G_string)
        neighbor, node_map = _build_unrooted_adjacency(root)

        # Build directed edges (from_id, to_id) -> count
        dedges = {}
        for from_id, adj in neighbor.items():
            for to_id in adj:
                dedges[(from_id, to_id)] = 0

        flat_names = {}  # (from_id, to_id) -> bitset

        # Name all leaf edges
        for dedge in list(dedges.keys()):
            from_id, to_id = dedge
            node = node_map[from_id]
            if is_leaf(node):
                leaf_bit = 1 << self.leaf_ids[node.name]
                flat_names[dedge] = leaf_bit

                # Register leaf set and accumulate branch length
                g_id = self.set2id(leaf_bit)
                bl = node.branch_length if node.branch_length is not None else 0.0
                self.Bip_bls[g_id] += bl

                # Mark done
                dedges[dedge] = -1

                # Increment counts for outgoing edges from 'to'
                for adj_id in neighbor[to_id]:
                    if adj_id != from_id:
                        dedge_out = (to_id, adj_id)
                        if dedge_out in dedges:
                            dedges[dedge_out] += 1

        # Special case: only 2 leaves
        leaf_count = sum(1 for n in node_map.values() if is_leaf(n))
        if leaf_count == 2:
            self.Bip_counts[1] = self.Bip_counts.get(1, 0.0) + weight
            return

        # Iteratively resolve internal edges
        edges_left = any(v != -1 for v in dedges.values())
        while edges_left:
            for dedge in list(dedges.keys()):
                if dedges[dedge] != 2:
                    continue

                from_id, to_id = dedge
                # Find the two incoming edges
                dedges_in = []
                for adj_id in neighbor[from_id]:
                    if adj_id != to_id:
                        dedges_in.append((adj_id, from_id))

                leaf_set_in_1 = flat_names[dedges_in[0]]
                leaf_set_in_2 = flat_names[dedges_in[1]]
                combined = leaf_set_in_1 | leaf_set_in_2
                flat_names[dedge] = combined

                g_id = self.set2id(combined)

                tmp_id1 = self.set2id(leaf_set_in_1)
                tmp_id2 = self.set2id(leaf_set_in_2)
                parts = (min(tmp_id1, tmp_id2), max(tmp_id1, tmp_id2))

                # Branch length
                from_node = node_map[from_id]
                to_node = node_map[to_id]
                if (from_node.parent is not None
                        and from_node.parent.id == to_id):
                    bl = (from_node.branch_length
                          if from_node.branch_length is not None else 0.0)
                elif (to_node.parent is not None
                      and to_node.parent.id == from_id):
                    bl = (to_node.branch_length
                          if to_node.branch_length is not None else 0.0)
                else:
                    # After unrooting, the edge between the two former
                    # children of the root may not have a direct
                    # parent-child relationship.  Use the sum of their
                    # branch lengths to the (removed) root.
                    bl_from = (from_node.branch_length
                               if from_node.branch_length is not None else 0.0)
                    bl_to = (to_node.branch_length
                             if to_node.branch_length is not None else 0.0)
                    bl = bl_from + bl_to

                self.Bip_bls[g_id] = self.Bip_bls.get(g_id, 0.0) + bl

                # Update counts
                if g_id >= len(self.Dip_counts):
                    while g_id >= len(self.Dip_counts):
                        self.Dip_counts.append({})
                self.Dip_counts[g_id][parts] = (
                    self.Dip_counts[g_id].get(parts, 0.0) + weight
                )
                self.Bip_counts[g_id] = (
                    self.Bip_counts.get(g_id, 0.0) + weight
                )

                # Mark done
                dedges[dedge] = -1
                for adj_id in neighbor[to_id]:
                    if adj_id != from_id:
                        dedge_out = (to_id, adj_id)
                        if dedge_out in dedges:
                            dedges[dedge_out] += 1

            edges_left = any(v != -1 for v in dedges.values())

    # ------------------------------------------------------------------
    # observation
    # ------------------------------------------------------------------

    def observation(self, trees, weight=1.0):
        """Observe a list of tree strings, updating all counts."""
        for tree_str in trees:
            self.decompose(tree_str, weight=weight)
            self.observations += weight

        # Rebuild set_sizes and size_ordered_bips
        self.set_sizes = {}
        self.size_ordered_bips = {}
        for sid, bitset in self.id_sets.items():
            size = _popcount(bitset)
            self.set_sizes[sid] = size
            self.size_ordered_bips.setdefault(size, []).append(sid)

    # ------------------------------------------------------------------
    # Probability functions
    # ------------------------------------------------------------------

    def Bi(self, n2):
        """Number of tree topologies for a bipartition of sizes n1, n2."""
        n1 = self.Gamma_size - n2
        if n1 == 1 or n2 == 1:
            return double_factorial(2 * self.Gamma_size - 5)
        n1 = max(2, n1)
        n2 = max(2, n2)
        return (double_factorial(2 * n1 - 3)
                * double_factorial(2 * n2 - 3))

    def Tri(self, n2, n3):
        """Number of tree topologies for a trifurcation of sizes n1, n2, n3."""
        n1 = self.Gamma_size - n2 - n3
        n1 = max(2, n1)
        n2 = max(2, n2)
        n3 = max(2, n3)
        return (double_factorial(2 * n1 - 3)
                * double_factorial(2 * n2 - 3)
                * double_factorial(2 * n3 - 3))

    def p_bip(self, g_id):
        """Bipartition probability by ID."""
        if self.Gamma_size < 4:
            return 1.0

        Bip_count = 0.0
        if g_id not in self.Bip_counts or g_id == -10 or g_id == 0:
            Bip_count = 0.0
        else:
            Bip_count = self.Bip_counts[g_id]

        if g_id not in self.set_sizes or g_id == -10:
            Bip_count = 0.0
        elif (self.set_sizes[g_id] == 1
              or self.set_sizes[g_id] == self.Gamma_size - 1):
            Bip_count = self.observations

        if self.alpha > 0:
            size = self.set_sizes.get(g_id, 0)
            return (Bip_count + self.alpha / self.N_Gamma
                    * self.Bi(size)) / (
                        self.observations + self.alpha)
        else:
            return Bip_count / self.observations

    def p_bip_by_set(self, gamma):
        """Bipartition probability by bitset."""
        if self.Gamma_size < 4:
            return 1.0
        if gamma in self.set_ids:
            g_id = self.set_ids[gamma]
        else:
            g_id = -10
        return self.p_bip(g_id)

    def p_dip(self, g_id, gp_id, gpp_id):
        """Conditional clade probability (by IDs, with gp_id < gpp_id)."""
        if self.Gamma_size < 4:
            return 1.0

        beta_switch = 1.0
        Dip_count = 0.0
        Bip_count = 0.0

        if g_id not in self.Bip_counts or g_id == -10 or g_id == 0:
            beta_switch = 0.0
            Bip_count = 0.0
            Dip_count = 0.0
        else:
            parts = (min(gp_id, gpp_id), max(gp_id, gpp_id))
            Bip_count = self.Bip_counts[g_id]

            if (gp_id == -10 or gpp_id == -10
                    or gp_id == 0 or gpp_id == 0
                    or g_id >= len(self.Dip_counts)
                    or parts not in self.Dip_counts[g_id]):
                Dip_count = 0.0
            else:
                Dip_count = self.Dip_counts[g_id][parts]

        if (g_id not in self.set_sizes
                or self.set_sizes[g_id] == 1
                or self.set_sizes[g_id] == self.Gamma_size - 1):
            Bip_count = self.observations

        if self.alpha > 0 or self.beta > 0:
            g_size = self.set_sizes.get(g_id, 1)
            gp_size = self.set_sizes.get(gp_id, 1)
            gpp_size = self.set_sizes.get(gpp_id, 1)
            numerator = (Dip_count
                         + self.alpha / self.N_Gamma
                         * self.Tri(gp_size, gpp_size)
                         + beta_switch * self.beta
                         / (2.0 ** (g_size - 1) - 1))
            denominator = (Bip_count
                           + self.alpha / self.N_Gamma
                           * self.Bi(g_size)
                           + beta_switch * self.beta)
            return numerator / denominator
        else:
            if Bip_count == 0:
                return 0.0
            return Dip_count / Bip_count

    def p_dip_by_set(self, gamma, gammap, gammapp):
        """Conditional clade probability by bitsets."""
        if self.Gamma_size < 4:
            return 1.0

        g_id = self.set_ids[gamma] if gamma in self.set_ids else -10
        gp_id = self.set_ids[gammap] if gammap in self.set_ids else -10
        gpp_id = self.set_ids[gammapp] if gammapp in self.set_ids else -10

        if gpp_id > gp_id:
            return self.p_dip(g_id, gp_id, gpp_id)
        else:
            return self.p_dip(g_id, gpp_id, gp_id)

    # ------------------------------------------------------------------
    # recompose
    # ------------------------------------------------------------------

    def recompose(self, G_string):
        """Compute conditional clade probabilities for a given tree.

        Returns a dict mapping bitset -> probability.
        """
        return_map = {}

        root = parse_newick(G_string)
        neighbor, node_map = _build_unrooted_adjacency(root)

        dedges = {}
        for from_id, adj in neighbor.items():
            for to_id in adj:
                dedges[(from_id, to_id)] = 0

        flat_names = {}
        q = {}

        # Name leaves
        for dedge in list(dedges.keys()):
            from_id, to_id = dedge
            node = node_map[from_id]
            if is_leaf(node):
                leaf_bit = 1 << self.leaf_ids[node.name]
                flat_names[dedge] = leaf_bit
                q[dedge] = 1.0
                return_map[leaf_bit] = 1.0
                dedges[dedge] = -1
                for adj_id in neighbor[to_id]:
                    if adj_id != from_id:
                        dedge_out = (to_id, adj_id)
                        if dedge_out in dedges:
                            dedges[dedge_out] += 1

        edges_left = any(v != -1 for v in dedges.values())
        while edges_left:
            for dedge in list(dedges.keys()):
                if dedges[dedge] != 2:
                    continue

                from_id, to_id = dedge
                dedges_in = []
                for adj_id in neighbor[from_id]:
                    if adj_id != to_id:
                        dedges_in.append((adj_id, from_id))

                leaf_set_in_1 = flat_names[dedges_in[0]]
                leaf_set_in_2 = flat_names[dedges_in[1]]
                combined = leaf_set_in_1 | leaf_set_in_2
                flat_names[dedge] = combined

                q[dedge] = (q[dedges_in[0]] * q[dedges_in[1]]
                            * self.p_dip_by_set(combined,
                                                leaf_set_in_1,
                                                leaf_set_in_2))
                return_map[combined] = q[dedge]

                dedges[dedge] = -1
                for adj_id in neighbor[to_id]:
                    if adj_id != from_id:
                        dedge_out = (to_id, adj_id)
                        if dedge_out in dedges:
                            dedges[dedge_out] += 1

            edges_left = any(v != -1 for v in dedges.values())

        return return_map

    # ------------------------------------------------------------------
    # p  (tree probability)
    # ------------------------------------------------------------------

    def p(self, tree_string):
        """Probability of a tree given the approximate posterior."""
        rec_map = self.recompose(tree_string)
        p_val = 0.0
        for gamma, val in rec_map.items():
            p_val = val
            # not_gamma: flip all bits except bit 0
            not_gamma = gamma ^ self.Gamma
            p_val *= rec_map.get(not_gamma, 0.0) * self.p_bip_by_set(gamma)
            if math.isnan(p_val):
                p_val = 0.0
            break
        return p_val

    # ------------------------------------------------------------------
    # mpp_tree  (maximum posterior probability tree)
    # ------------------------------------------------------------------

    def mpp_tree(self):
        """Return (newick_string, max_posterior_probability)."""
        qmpp = {}

        # DP from small to large
        for size in sorted(self.size_ordered_bips.keys()):
            for g_id in self.size_ordered_bips[size]:
                if size == 1:
                    qmpp[g_id] = 1.0
                else:
                    max_cp = 0.0
                    if g_id < len(self.Dip_counts):
                        for (gp_id, gpp_id), _ in self.Dip_counts[g_id].items():
                            cp = (self.p_dip(g_id, gp_id, gpp_id)
                                  * qmpp.get(gp_id, 0.0)
                                  * qmpp.get(gpp_id, 0.0))
                            if cp > max_cp:
                                max_cp = cp
                    qmpp[g_id] = max_cp

        # Find best root bipartition
        max_pp = 0.0
        sum_pp = 0.0
        max_bip = -1
        max_not_bip = -1

        for g_id in self.Bip_counts:
            gamma = self.id_sets[g_id]
            not_gamma = gamma ^ self.Gamma
            if not_gamma not in self.set_ids:
                continue
            not_g_id = self.set_ids[not_gamma]
            pp = (qmpp.get(g_id, 0.0)
                  * qmpp.get(not_g_id, 0.0)
                  * self.p_bip(g_id))
            sum_pp += pp
            if pp > max_pp:
                max_pp = pp
                max_bip = g_id
                max_not_bip = not_g_id

        if max_bip == -1 or sum_pp == 0:
            return ("", 0.0)

        # Root support and branch length
        support = max_pp / sum_pp * 2  # we looked at everything twice
        bl = min(
            self.Bip_bls.get(max_bip, 0.0)
            / max(self.Bip_counts.get(max_bip, 1.0), 1e-100),
            0.99,
        )
        bs_str = f"{support}:{bl}"

        left = self.mpp_backtrack(max_bip, qmpp)
        right = self.mpp_backtrack(max_not_bip, qmpp)
        max_tree = f"({left},{right}){bs_str};\n"

        return (max_tree, max_pp)

    def mpp_backtrack(self, g_id, qmpp):
        """Recursive backtrack for mpp_tree."""
        # Leaf
        if self.set_sizes.get(g_id, 0) == 1:
            bl = self.Bip_bls.get(g_id, 0.0) / max(self.observations, 1e-100)
            bitset = self.id_sets[g_id]
            leaf_id = 0
            for i in range(1, self.Gamma_size + 1):
                if bitset & (1 << i):
                    leaf_id = i
                    break
            return f"{self.id_leaves[leaf_id]}:{bl}"

        # Internal node: find best split
        max_cp = 0.0
        sum_cp = 0.0
        max_gp_id = -1
        max_gpp_id = -1

        if g_id < len(self.Dip_counts):
            for (gp_id, gpp_id), _ in self.Dip_counts[g_id].items():
                cp = (self.p_dip(g_id, gp_id, gpp_id)
                      * qmpp.get(gp_id, 0.0)
                      * qmpp.get(gpp_id, 0.0))
                sum_cp += cp
                if cp > max_cp:
                    max_cp = cp
                    max_gp_id = gp_id
                    max_gpp_id = gpp_id

        if sum_cp == 0:
            support_str = "0"
        else:
            support_str = str(max_cp / sum_cp)
        bl = (self.Bip_bls.get(g_id, 0.0)
              / max(self.Bip_counts.get(g_id, 1.0), 1e-100))

        left = self.mpp_backtrack(max_gp_id, qmpp)
        right = self.mpp_backtrack(max_gpp_id, qmpp)
        return f"({left},{right}){support_str}:{bl}"

    # ------------------------------------------------------------------
    # count_trees
    # ------------------------------------------------------------------

    def count_trees(self):
        """Count amalgamated trees with the complete leaf set."""
        g_id_count = {}

        for size in sorted(self.size_ordered_bips.keys()):
            for g_id in self.size_ordered_bips[size]:
                if size == 1:
                    g_id_count[g_id] = 1.0
                else:
                    g_id_count[g_id] = 0.0
                    if g_id < len(self.Dip_counts):
                        for (gp_id, gpp_id), _ in self.Dip_counts[g_id].items():
                            g_id_count[g_id] += (
                                g_id_count.get(gp_id, 0.0)
                                * g_id_count.get(gpp_id, 0.0)
                            )

        count = 0.0
        for g_id in self.Bip_counts:
            gamma = self.id_sets[g_id]
            not_gamma = gamma ^ self.Gamma
            if not_gamma not in self.set_ids:
                continue
            gamma_size = _popcount(gamma)
            not_gamma_size = _popcount(not_gamma)
            val = (g_id_count.get(self.set_ids[gamma], 0.0)
                   * g_id_count.get(self.set_ids[not_gamma], 0.0))
            if gamma_size > not_gamma_size:
                count += val
            elif gamma_size == not_gamma_size:
                count += val / 2.0

        return count

    # ------------------------------------------------------------------
    # random_tree / random_split
    # ------------------------------------------------------------------

    def random_tree(self):
        """Sample a random tree from the approximate posterior."""
        total = sum(self.Bip_counts.values())
        if total == 0:
            return ""
        rnd = random.random()
        cumsum = 0.0
        g_id = None
        for gid, cnt in self.Bip_counts.items():
            cumsum += cnt
            g_id = gid
            if cumsum > total * rnd:
                break

        gamma = self.id_sets[g_id]
        not_gamma = gamma ^ self.Gamma
        return (f"({self.random_split(gamma)}:1,"
                f"{self.random_split(not_gamma)}:1);\n")

    def random_split(self, gamma):
        """Recursive random split of a clade."""
        gamma_v = _bits_to_list(gamma, self.Gamma_size)
        gamma_size = len(gamma_v)

        if gamma_size == 1:
            return self.id_leaves[gamma_v[0]]

        rnd = random.random()
        p_sum = 0.0
        g_id = self.set_ids.get(gamma, 0)
        beta_switch = 1.0
        Bip_count = 0.0

        if not g_id:
            beta_switch = 0.0
            Bip_count = 0.0

        gammap = 0
        gammapp = 0

        for gp_size in range(1, gamma_size // 2 + 1):
            saw = 0
            found = False

            if g_id and g_id < len(self.Dip_counts):
                for (gp_id_k, gpp_id_k), _ in self.Dip_counts[g_id].items():
                    gp_bits = self.id_sets.get(gp_id_k, 0)
                    gpp_bits = self.id_sets.get(gpp_id_k, 0)
                    gp_id_size = _popcount(gp_bits)
                    gpp_id_size = _popcount(gpp_bits)
                    this_size = min(gp_id_size, gpp_id_size)

                    if this_size == gp_size:
                        p_sum += self.p_dip_by_set(
                            gamma,
                            self.id_sets[gp_id_k],
                            self.id_sets[gpp_id_k],
                        )
                        if rnd < p_sum:
                            p_sum = -1
                        saw += 1

                    if p_sum < 0:
                        gammap = self.id_sets[gp_id_k]
                        gammapp = self.id_sets[gpp_id_k]
                        found = True
                        break

            if found:
                break

            # Unobserved partitions
            if g_id:
                Bip_count = self.Bip_counts.get(g_id, 0.0)
            if gamma_size == 1 or gamma_size == self.Gamma_size - 1:
                Bip_count = self.observations

            nbip = math.comb(gamma_size, gp_size)
            if gamma_size - gp_size == gp_size:
                nbip //= 2

            denom = (Bip_count
                     + self.alpha / self.N_Gamma * self.Bi(gamma_size)
                     + beta_switch * self.beta)
            if denom != 0:
                p_sum += ((0
                           + self.alpha / self.N_Gamma
                           * self.Tri(gp_size, gamma_size - gp_size)
                           + beta_switch * self.beta
                           / (2.0 ** (gamma_size - 1) - 1))
                          / denom * (nbip - saw))

            if rnd < p_sum:
                p_sum = -1

            if p_sum < 0:
                # Pick random unsampled partition
                while True:
                    chosen = random.sample(gamma_v, gp_size)
                    gammap = 0
                    for v in chosen:
                        gammap |= (1 << v)
                    gammapp = gamma ^ gammap

                    gp_id_t = self.set_ids.get(gammap, 0)
                    gpp_id_t = self.set_ids.get(gammapp, 0)
                    parts = (min(gp_id_t, gpp_id_t), max(gp_id_t, gpp_id_t))
                    if (g_id < 0 or g_id >= len(self.Dip_counts)
                            or parts not in self.Dip_counts[g_id]
                            or self.Dip_counts[g_id][parts] == 0):
                        break
                break

        return (f"({self.random_split(gammap)}:1,"
                f"{self.random_split(gammapp)}:1)")

    # ------------------------------------------------------------------
    # save_state / load_state
    # ------------------------------------------------------------------

    def save_state(self, fname):
        """Write ALE file in the exact C++ format."""
        with open(fname, "w") as fout:
            fout.write("#constructor_string\n")
            fout.write(self.constructor_string.strip() + "\n")

            fout.write("#observations\n")
            fout.write(f"{self.observations}\n")

            fout.write("#Bip_counts\n")
            for gid in sorted(self.Bip_counts.keys()):
                fout.write(f"{gid}\t{self.Bip_counts[gid]}\n")

            fout.write("#Bip_bls\n")
            for gid in sorted(self.Bip_bls.keys()):
                fout.write(f"{gid}\t{self.Bip_bls[gid]}\n")

            fout.write("#Dip_counts\n")
            for index in range(len(self.Dip_counts)):
                for (gp_id, gpp_id), count in self.Dip_counts[index].items():
                    fout.write(f"{index}\t{gp_id}\t{gpp_id}\t{count}\n")

            fout.write("#last_leafset_id\n")
            fout.write(f"{self.last_leafset_id}\n")

            fout.write("#leaf-id\n")
            for name in sorted(self.leaf_ids.keys()):
                fout.write(f"{name}\t{self.leaf_ids[name]}\n")

            fout.write("#set-id\n")
            for bitset in sorted(self.set_ids.keys()):
                sid = self.set_ids[bitset]
                fout.write(f"{sid}\t:")
                for i in range(self.Gamma_size + 1):
                    if bitset & (1 << i):
                        fout.write(f"\t{i}")
                fout.write("\n")

            fout.write("#END\n")

    def load_state(self, fname):
        """Read ALE file in the exact C++ format."""
        reading = "#nothing"

        with open(fname, "r") as fin:
            for line in fin:
                line = line.rstrip("\n")
                if "#" in line:
                    reading = line.strip()
                elif reading == "#constructor_string":
                    tree_string = line.strip()
                    self.constructor_string = tree_string
                    self.construct(tree_string)
                    reading = "#nothing"
                elif reading == "#observations":
                    self.observations = float(line.strip())
                elif reading == "#Bip_counts":
                    tokens = line.strip().split()
                    if len(tokens) >= 2:
                        self.Bip_counts[int(tokens[0])] = float(tokens[1])
                elif reading == "#Bip_bls":
                    tokens = line.strip().split()
                    if len(tokens) >= 2:
                        self.Bip_bls[int(tokens[0])] = float(tokens[1])
                elif reading == "#Dip_counts":
                    tokens = line.strip().split()
                    if len(tokens) >= 4:
                        idx = int(tokens[0])
                        parts = (int(tokens[1]), int(tokens[2]))
                        count = float(tokens[3])
                        while idx >= len(self.Dip_counts):
                            self.Dip_counts.append({})
                        self.Dip_counts[idx][parts] = count
                elif reading == "#last_leafset_id":
                    self.last_leafset_id = int(line.strip())
                elif reading == "#leaf-id":
                    tokens = line.strip().split()
                    if len(tokens) >= 2:
                        name = tokens[0]
                        lid = int(tokens[1])
                        self.leaf_ids[name] = lid
                        self.id_leaves[lid] = name
                elif reading == "#set-id":
                    fields = line.strip().split(":")
                    if len(fields) >= 2:
                        set_id = int(fields[0].strip())
                        tokens = fields[1].strip().split()
                        bitset = 0
                        for t in tokens:
                            bitset |= (1 << int(t))
                        self.set_ids[bitset] = set_id
                        self.id_sets[set_id] = bitset

        # Pad Dip_counts to match last_leafset_id + 1 so that indices
        # corresponding to empty dicts (which save_state does not write)
        # are present after loading.
        while len(self.Dip_counts) < self.last_leafset_id + 1:
            self.Dip_counts.append({})

        # Root bipartition (bits 1..Gamma_size set)
        root_bits = 0
        for i in range(1, self.Gamma_size + 1):
            root_bits |= (1 << i)
        self.id_sets[-1] = root_bits
        self.set_ids[root_bits] = -1

        # Rebuild set_sizes and size_ordered_bips
        self.set_sizes = {}
        self.size_ordered_bips = {}
        for sid, bitset in self.id_sets.items():
            size = _popcount(bitset)
            self.set_sizes[sid] = size
            self.size_ordered_bips.setdefault(size, []).append(sid)

    # ------------------------------------------------------------------
    # Utility methods
    # ------------------------------------------------------------------

    def get_leaf_names(self):
        """Return all leaf names."""
        return list(self.leaf_ids.keys())

    def set_alpha(self, a):
        """Set the alpha correction parameter."""
        self.alpha = a

    def set_beta(self, b):
        """Set the beta correction parameter."""
        self.beta = b

    def set2name(self, leaf_set):
        """Return human-readable name for a bitset."""
        parts = []
        for i in range(1, self.Gamma_size + 1):
            if leaf_set & (1 << i):
                parts.append(self.id_leaves[i])
        return self.name_separator.join(parts)

    def compute_ordered_vector_of_clades(self):
        """Return (ids, sizes) ordered by size, with root (-1) appended."""
        ids = []
        id_sizes = []
        for size in sorted(self.size_ordered_bips.keys()):
            for gid in self.size_ordered_bips[size]:
                ids.append(gid)
                id_sizes.append(size)
        ids.append(-1)
        id_sizes.append(self.Gamma_size)
        return ids, id_sizes


# Alias used by ale_util and cli modules
ApproxPosterior = approx_posterior
