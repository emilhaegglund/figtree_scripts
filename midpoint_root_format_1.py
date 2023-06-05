from ete3 import Tree
import sys


tree = Tree(sys.argv[1], format=1)
## QUICK-FIX to resolve issue with moving support values
# on midpoint root.

# Store clades and support values from orig phylo
node_support = {}
for node in tree.traverse():
    if not node.is_leaf() and not node.is_root():
        node_leafs = frozenset(node.get_leaves())
        hash_node_leafs = hash(node_leafs)
        if hash_node_leafs not in node_support.keys():
            node_support[hash_node_leafs] = node.name

# Midpoint root
R = tree.get_midpoint_outgroup()
tree.set_outgroup(R)

# Find support values that has been wrongly placed
wrongly_placed_supports = []
error_nodes = []
for node in tree.traverse():
    if not node.is_leaf() and not node.is_root():
        node_leafs = frozenset(node.get_leaves())
        hash_node_leafs = hash(node_leafs)
        if hash_node_leafs not in node_support.keys():
            wrongly_placed_supports.append(node.name)
            error_nodes.append(node)
            print(node.name)

# Find the node with missing support
for node in tree.traverse():
    if not node.is_leaf() and not node.is_root():
        if node.name == '':
            missing_node = node
            break

wrongly_placed_supports = [s for s in wrongly_placed_supports if s != '']
# Starting from node with missing value,
# add values from the list and continue up
# until the list is empty
while len(wrongly_placed_supports) != 0:
    missing_node.name = wrongly_placed_supports[-1]
    print(missing_node.name)
    wrongly_placed_supports.pop(-1)
    missing_node = missing_node.up

root = tree.get_tree_root()
root_children =  root.get_children()
for node in root_children:
    if node not in error_nodes:
        tmp_support = node.name
for node in root_children:
    if node in error_nodes:
        node.name = tmp_support

tree.write(outfile=sys.argv[2], format=1)
