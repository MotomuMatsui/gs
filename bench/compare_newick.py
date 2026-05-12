#!/usr/bin/env python3
# Compare two gs2 Newick outputs by topology and EP-support tolerance.
#
# gs2 Newick grammar (no branch lengths, leaves are integer IDs, internal nodes
# carry a single floating-point EP support value after the close paren):
#     tree   := node ';'
#     node   := leaf | '(' node (',' node)+ ')' support?
#     leaf   := INT
#     support:= FLOAT
#
# A "clade" is the frozenset of leaf IDs under an internal node.  Two trees are
# topology-equal iff their clade sets are equal.  When topologies match, EP
# support is compared per-clade with abs tolerance.

import sys

def tokenize(s):
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c in '(),;':
            yield (c, c)
            i += 1
        elif c.isspace():
            i += 1
        else:
            j = i
            while j < n and s[j] not in '(),;':
                j += 1
            yield ('TOK', s[i:j])
            i = j

class Node:
    __slots__ = ('children', 'leaf_id', 'support')
    def __init__(self):
        self.children = []
        self.leaf_id = None   # int, set on leaves only
        self.support = None   # float, set on internal nodes if a token followed ')'

def parse_newick(s):
    # Recursive descent over the token stream.  No error recovery — we assume
    # well-formed gs2 output (parens are balanced, every internal node either
    # has a trailing FLOAT token or nothing).
    #
    # Defensive prefix/suffix strip: some build/runtime configurations leak a
    # `which mmseqs` path to stdout before the Newick (dash mishandles the
    # `&> /dev/null` redirection in main.cpp's PATH probe).  We anchor on the
    # first '(' and the last ';' to ignore any surrounding noise.
    lo = s.find('(')
    hi = s.rfind(';')
    if lo < 0 or hi < lo:
        raise ValueError("no Newick tree found (no '(' and ';' bracket)")
    toks = list(tokenize(s[lo:hi]))
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else (None, None)

    def consume():
        t = toks[pos[0]]
        pos[0] += 1
        return t

    def parse_node():
        kind, val = peek()
        node = Node()
        if kind == '(':
            consume()
            node.children.append(parse_node())
            while peek()[0] == ',':
                consume()
                node.children.append(parse_node())
            assert peek()[0] == ')', f"expected ')' at pos {pos[0]}"
            consume()
            # Optional trailing label/support
            if peek()[0] == 'TOK':
                node.support = float(consume()[1])
        else:
            # leaf
            assert kind == 'TOK', f"expected leaf at pos {pos[0]}: {kind!r}"
            consume()
            node.leaf_id = int(val)
        return node

    root = parse_node()
    return root

def collect_clades(root):
    """Walk the tree once, return:
       leaves: frozenset of all leaf IDs
       clades: dict[frozenset[int], float | None]
               keyed by clade leaf set, value = EP support (None if absent).
       Trivial clades (singleton leaves, full leaf set) are excluded — they
       carry no useful information and the full-leaf clade has no support.
    """
    clades = {}
    all_leaves = set()

    def visit(node):
        if node.leaf_id is not None:
            all_leaves.add(node.leaf_id)
            return frozenset([node.leaf_id])
        leaves = set()
        for c in node.children:
            leaves |= visit(c)
        fs = frozenset(leaves)
        if len(fs) > 1:
            clades[fs] = node.support
        return fs

    visit(root)
    # Drop the root clade (== all leaves) — uninformative
    full = frozenset(all_leaves)
    clades.pop(full, None)
    return full, clades

def compare(path_a, path_b, eps):
    with open(path_a) as f: tree_a = parse_newick(f.read())
    with open(path_b) as f: tree_b = parse_newick(f.read())

    leaves_a, clades_a = collect_clades(tree_a)
    leaves_b, clades_b = collect_clades(tree_b)

    errors = []
    if leaves_a != leaves_b:
        only_a = sorted(leaves_a - leaves_b)[:10]
        only_b = sorted(leaves_b - leaves_a)[:10]
        errors.append(f"leaf-set mismatch: a={len(leaves_a)} b={len(leaves_b)} "
                      f"only-a≈{only_a} only-b≈{only_b}")
        return errors

    set_a = set(clades_a.keys())
    set_b = set(clades_b.keys())

    if set_a != set_b:
        only_a = len(set_a - set_b)
        only_b = len(set_b - set_a)
        errors.append(f"topology mismatch: {only_a} clades only in {path_a}, "
                      f"{only_b} only in {path_b}")
        # Show a couple of examples
        for c in list(set_a - set_b)[:3]:
            errors.append(f"  only in a: {sorted(c)[:6]}…")
        for c in list(set_b - set_a)[:3]:
            errors.append(f"  only in b: {sorted(c)[:6]}…")
        return errors

    # Topology matches; check EP within eps.
    worst = 0.0
    worst_clade = None
    diffs = []
    for c in set_a:
        sa, sb = clades_a[c], clades_b[c]
        if sa is None and sb is None:
            continue
        if sa is None or sb is None:
            errors.append(f"support presence mismatch on clade {sorted(c)[:6]}…: "
                          f"a={sa} b={sb}")
            continue
        d = abs(sa - sb)
        diffs.append(d)
        if d > worst:
            worst = d
            worst_clade = c

    if errors:
        return errors

    over = sum(1 for d in diffs if d > eps)
    if over > 0:
        errors.append(f"EP support diff > {eps}: {over}/{len(diffs)} clades. "
                      f"worst={worst:.4f} on clade size={len(worst_clade)}")
    return errors, worst, len(diffs)

def main():
    import argparse
    p = argparse.ArgumentParser(description="Compare two gs2 Newick files.")
    p.add_argument('a')
    p.add_argument('b')
    p.add_argument('--eps', type=float, default=0.10,
                   help="Max allowed |EP_a - EP_b| per clade (default 0.10).")
    p.add_argument('--strict', action='store_true',
                   help="Topology must match AND no EP diff > eps.")
    args = p.parse_args()

    result = compare(args.a, args.b, args.eps)
    if isinstance(result, list):
        # Topology mismatch — fatal
        print("FAIL (topology):")
        for e in result:
            print("  " + e)
        sys.exit(2)
    errors, worst, n = result
    if errors:
        print("FAIL (EP tolerance):")
        for e in errors:
            print("  " + e)
        sys.exit(1 if not args.strict else 2)
    print(f"OK: topology matches ({n} non-trivial clades), max EP diff = {worst:.4f} (eps={args.eps})")
    sys.exit(0)

if __name__ == '__main__':
    main()
