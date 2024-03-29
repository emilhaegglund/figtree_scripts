"""
Script to midpoint-root phylogenies and collapse clades for
visualization in FigTree
"""
import argparse
import sys
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
            if node.is_leaf():
                print("Node is leag")
            if node.up is not None or format_root_node:
                if "/" in node.name:
                    node.support = float(node.name.split("/")[1])
                if node.collapse:
                    newick.append(
                        #node.name
                        "[&label="
                        + '"' + str(node.name) + '"'
                        + ',!collapse={"collapsed",'
                        + str(node.length)
                        + '},!name="'
                        + node.label
                        + '"]:'
                        + str(node.dist)
                    )
                else:
                    newick.append(
                        #node.name
                        '[&label='
                        + '"' + str(node.name) + '"'
                        + ']:'
                        + str(node.dist)
                    )
        else:
            if node is not rootnode and node != node.up.children[0]:
                newick.append(",")
            if leaf(node):
                newick.append("'" + node.label + "':" + str(node.dist))
            else:
                newick.append("(")

    newick.append(";")
    return "".join(newick)

def create_label(node_annotation, n_taxa, column, leaf):
    """
    Create annotation string
    """
    label = []
    if leaf:
        annotation_columns = node_annotation.columns
    else:
        annotation_columns = node_annotation.columns[1:]
    for c in annotation_columns:
        label.append(node_annotation[c].iat[0])
        if c == column and not leaf:
            break
    label_str = "_".join(label[1:])

    if not leaf:
        label_str = label_str + "_(" + str(n_taxa) + ")"

    return label_str

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
            the accession and the following columns should be
            taxonomic information on the format domain, phylum, class, order, family, genus, and species""",
    )
    args.add_argument(
        "--column",
        "-c",
        required=True,
        help="Name of the column in the TSV to use for determine monophyly, default will use second column",
    )
    args.add_argument(
        "--headers",
        help="Comma-separated list of header to use for collapsing and annotation"
    )
    args.add_argument(
        "--taxa",
        "-n",
        type=int,
        default=3,
        help="Minimal number of monophyletic taxa to collapse a node",
    )
    args.add_argument("--figtree-settings", "-s", required=True, help="""FigTree Display settings, the file in the example directory of the repo can be used.""")
    args.add_argument("--midpoint-root", "-m", action='store_true', help="Midpoint root the phylogeny before collapsing nodes.")
    args.add_argument("--output", "-o", required=True, help="Path to output file")

    return args.parse_args()


def filter_tsv(df, taxa_columns):
    """
    Filter the dataframe so it only contains columns
    to use for collapse and annotation.

    Returns the dataframe with the columns in the same
    order as the items in the taxa_columns.
    """
    accession_column = df.columns[0]
    df = df[[accession_column] + taxa_columns]
    return df


def check_tsv(df, taxa_columns):
    """
    Check that the selected columns to use for collapsing and
    annotation is present in the dataframe.
    """
    if set(taxa_columns).issubset(df.columns):
        return True
    else:
        return False

args = parse_command_line()
tree = Tree(args.tree, format=1)
annotation = pd.read_csv(args.annotation, sep="\t")

if args.column not in annotation.columns:
    sys.exit(f"Column '{args.column}' is not present in the annotation file")

column = args.column
if args.headers:
    taxa_columns = args.headers.split(",")
else:
    taxa_columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

if check_tsv(annotation, taxa_columns):
    annotation = filter_tsv(annotation, taxa_columns)
else:
    sys.exit("TSV-file does not include required columns")


if annotation.columns[0] != 'taxa':
    annotation.rename(columns={annotation.columns[0]:'taxa'}, inplace=True)
output = args.output
min_taxa = args.taxa
figtree_settings = args.figtree_settings

if args.midpoint_root:
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
    er_nodes = []
    for node in tree.traverse():
        if not node.is_leaf() and not node.is_root():
            node_leafs = frozenset(node.get_leaves())
            hash_node_leafs = hash(node_leafs)
            if hash_node_leafs not in node_support.keys():
                wrongly_placed_supports.append(node.name)
                er_nodes.append(node)
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
        if node not in er_nodes:
            tmp_support = node.name
    for node in root_children:
        if node in er_nodes:
            node.name = tmp_support

##### End fix #####

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
        node.label = (create_label(node_annotation, n_taxa=n_taxa, column=column, leaf=False))
        print("Collapsed " + node.label)
    else:
        node.collapse = False
        node.label = (create_label(node_annotation, n_taxa=1, column=column, leaf=True))
        print("Not collapsed " + node.label)
        if not monophyletic:
            print("Not monophyletic")
        if min_taxa > n_taxa:
            print('Not enough taxa')


# Read FigTree-settings
figtree_settings_block = []
with open(figtree_settings, "r") as f:
    for line in f:
        figtree_settings_block.append(line)

# Convert to nexus-format and write to file
nexus_tree = loc_write_newick(tree)
write_tree(nexus_tree, figtree_settings_block, output)
