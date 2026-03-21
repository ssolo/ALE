class Node:
    _next_id = 0

    def __init__(self, name="", branch_length=None, parent=None):
        self.name = name
        self.branch_length = branch_length
        self.children = []
        self.parent = parent
        self.id = Node._next_id
        Node._next_id += 1

    def add_child(self, child):
        child.parent = self
        self.children.append(child)
        return child

    def __repr__(self):
        if self.name:
            return f"Node({self.name!r})"
        return f"Node(id={self.id}, children={len(self.children)})"


def is_leaf(node):
    return len(node.children) == 0


def get_leaves(node):
    leaves = []
    _collect_leaves(node, leaves)
    return leaves


def _collect_leaves(node, leaves):
    if is_leaf(node):
        leaves.append(node)
    else:
        for child in node.children:
            _collect_leaves(child, leaves)


def get_all_nodes(node):
    nodes = []
    _collect_postorder(node, nodes)
    return nodes


def _collect_postorder(node, nodes):
    for child in node.children:
        _collect_postorder(child, nodes)
    nodes.append(node)


def copy_tree(node, parent=None):
    new_node = Node(
        name=node.name,
        branch_length=node.branch_length,
        parent=parent,
    )
    new_node.id = node.id
    for child in node.children:
        new_child = copy_tree(child, parent=new_node)
        new_node.children.append(new_child)
    return new_node


def parse_newick(string):
    string = string.strip()
    if string.endswith(";"):
        string = string[:-1]
    tokens = _tokenize(string)
    root, _ = _parse_tokens(tokens, 0)
    _assign_ids(root)
    return root


def _tokenize(s):
    tokens = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c in "(),":
            tokens.append(c)
            i += 1
        elif c == ":":
            tokens.append(":")
            i += 1
            start = i
            while i < n and s[i] not in "(),;:":
                i += 1
            tokens.append(s[start:i])
        else:
            start = i
            while i < n and s[i] not in "(),;:":
                i += 1
            tokens.append(s[start:i])
    return tokens


def _parse_tokens(tokens, pos):
    node = Node()

    if pos < len(tokens) and tokens[pos] == "(":
        pos += 1  # consume '('
        child, pos = _parse_tokens(tokens, pos)
        node.add_child(child)

        while pos < len(tokens) and tokens[pos] == ",":
            pos += 1  # consume ','
            child, pos = _parse_tokens(tokens, pos)
            node.add_child(child)

        if pos < len(tokens) and tokens[pos] == ")":
            pos += 1  # consume ')'

    # Parse label (name or bootstrap) after ')' or as a leaf name
    if pos < len(tokens) and tokens[pos] not in "(),:":
        node.name = tokens[pos]
        pos += 1

    # Parse branch length
    if pos < len(tokens) and tokens[pos] == ":":
        pos += 1  # consume ':'
        if pos < len(tokens):
            try:
                node.branch_length = float(tokens[pos])
            except ValueError:
                node.branch_length = None
            pos += 1

    return node, pos


def _assign_ids(node):
    counter = [0]
    _assign_ids_postorder(node, counter)


def _assign_ids_postorder(node, counter):
    for child in node.children:
        _assign_ids_postorder(child, counter)
    node.id = counter[0]
    counter[0] += 1


def to_newick(node):
    return _to_newick_recursive(node) + ";"


def _to_newick_recursive(node):
    result = ""
    if node.children:
        child_strs = [_to_newick_recursive(c) for c in node.children]
        result = "(" + ",".join(child_strs) + ")"
    if node.name:
        result += node.name
    if node.branch_length is not None:
        result += ":" + str(node.branch_length)
    return result
