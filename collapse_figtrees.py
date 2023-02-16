"""
Script to midpoint-root phylogenies and collapse clades for
visualization in FigTree
"""
import argparse
import pandas as pd
from ete3 import Tree


def loc_write_newick(rootnode, format_root_node=True, is_leaf_fn=None):
    """
    This function is from the ete-toolkit. It has been rewritten so that
    it writes the correct format for reading with FigTree.
    Iteratively export a tree structure and returns its NHX
    representation."""
    newick = []
    leaf = is_leaf_fn if is_leaf_fn else lambda n: not bool(n.children)
    for postorder, node in rootnode.iter_prepostorder(is_leaf_fn=is_leaf_fn):
        if postorder:
            newick.append(")")
            if node.up is not None or format_root_node:
                if node.collapse:
                    newick.append(
                        node.name
                        + "[&label="
                        + str(node.support)
                        + ',!collapse={"collapsed",'
                        + str(node.length)
                        + '},gtdb_taxonomy="'
                        + node.label
                        + '",!name="'
                        + node.label
                        + '"]:'
                        + str(node.dist)
                    )
                else:
                    newick.append(
                        node.name
                        + "[&label="
                        + str(node.support)
                        + "]:"
                        + str(node.dist)
                    )
        else:
            if node is not rootnode and node != node.up.children[0]:
                newick.append(",")
            if leaf(node):
                newick.append("'" + node.name + "':" + str(node.dist))
            else:
                newick.append("(")

    newick.append(";")
    return "".join(newick)


def write_tree(nexus_tree, figtree_settings, output):
    with open(output, "w") as f:
        f.write("#NEXUS\n")
        f.write("begin trees;\n")
        f.write("tree tree_1 = [&R] " + nexus_tree + "\n")
        f.write("end;\n")
        # Write FigTree Display settings block
        for line in figtree_settings:
            f.write(line)

        return True


def parse_command_line():

    args = argparse.ArgumentParser(
        prog="collapse_figtrees",
        description="Midpoint-root phylogeny and collapse clades for visualization in FigTree",
    )
    args.add_argument(
        "--tree",
        "-t",
        required=True,
        help="Path to phylogeny to midpoint-root and collapse",
    ),
    args.add_argument(
        "--annotation",
        "-a",
        required=True,
        help="""
            TSV-file with annotation for the taxa, first column must contain
            the accessions of taxa in phylogeny, remaining columns can be
            taxonomic information or other data.""",
    )
    args.add_argument(
        "--column",
        "-c",
        help="Name of the column in the TSV to use for determine monophyly, default will use second column",
    )
    args.add_argument(
        "--taxa",
        "-n",
        type=int,
        default=3,
        help="Minimal number of monophyletic taxa to collapse a node",
    )
    args.add_argument("--figtree-settings", "-s", required=True)
    args.add_argument("--output", "-o", required=True, help="Path to output file")

    return args.parse_args()


args = parse_command_line()
tree = Tree(args.tree)
annotation = pd.read_csv(args.annotation, sep="\t")


if args.column:
    column = args.column
else:
    column = annotation.columns[1]
    print(column)

output = args.output
min_taxa = args.taxa
figtree_settings = args.figtree_settings

# Midpoint-root
R = tree.get_midpoint_outgroup()
tree.set_outgroup(R)

# This is reference dist for triangles
dist = tree.get_farthest_leaf()

# Loop over the tree to find what nodes to collapse
for node in tree.traverse():
    # Chek if ancestor was collapsed
    ancestors = node.get_ancestors()
    ancestor_collapsed = False
    for a_node in ancestors:
        if a_node.collapse:
            ancestor_collapsed = True
    leafs = node.get_leaves()
    leaf_names = [l.name for l in leafs]

    # Count the number of taxa in the clade
    node_annotation = annotation[annotation["taxa"].isin(leaf_names)]
    n_taxa = node_annotation.shape[0]

    # Check if the node is a monophyletic clade
    taxonomical_units = set(node_annotation[column].to_list())
    if len(taxonomical_units) == 1:
        monophyletic = True
    else:
        monophyletic = False

    # Collapse node
    if monophyletic and n_taxa >= min_taxa and not ancestor_collapsed:
        node.collapse = True
        root_dist = node.get_distance(tree.get_tree_root())
        far_leaf = node.get_farthest_leaf()
        length = dist[1] - (root_dist + far_leaf[1])
        node.length = length
        node.label = (
            node_annotation[column].iat[0] + " (" + str(node_annotation.shape[0]) + ")"
        )
        print("Collapsed " + node.label)
    else:
        node.collapse = False

# Read FigTree-settings
figtree_settings_block = []
with open(figtree_settings, "r") as f:
    for line in f:
        figtree_settings_block.append(line)

# Convert to nexus-format and write to file
nexus_tree = loc_write_newick(tree)
write_tree(nexus_tree, figtree_settings_block, output)
