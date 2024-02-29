from ete3 import Tree
import random
import os,sys

species_tree_file = "s_tree_10.txt"
num_gene_trees = 1000
duprate = 2
numerator = 1
denominator = 0
numerator_fix = 0
denominator_fix = 1

#species_tree_file = sys.argv[1]
#num_gene_trees = int(sys.argv[2])
#duprate = int(sys.argv[3])
#numerator = int(sys.argv[4])
#denominator = int(sys.argv[5])
#numerator_fix = int(sys.argv[6])
#denominator_fix = int(sys.argv[7])

def change_dup_rates(species_tree, rate):
    for node in species_tree.traverse():
        if node.dist == 0 and not node.is_root():
            node.dist = rate
    #print(species_tree.write())

def duplicate_node(subtree_root):
    new_subtree1 = subtree_root.copy()
    new_subtree2 = subtree_root.copy()
    subtree_root.children = []
    subtree_root.name = subtree_root.name + "_Dup"
    subtree_root.add_child(new_subtree1)
    subtree_root.add_child(new_subtree2)
    return new_subtree1


def dodup(subtree_root, pr, duprate):
    if subtree_root.is_leaf():
        return
    random_pr = random.random()
    # if subtree_root.name == "D":
    # print(pr)
    # print(random_pr)
    if random_pr <= pr:
        # print("dup")
        duplicate_node(subtree_root)
        if duprate != 0:
            subtree_root.get_children()[0].dist = pr / duprate
            subtree_root.get_children()[1].dist = pr / duprate
            dodup(subtree_root.get_children()[0], pr / duprate, duprate)
            dodup(subtree_root.get_children()[1], pr / duprate, duprate)
        else:
            subtree_root.get_children()[0].dist = 0
            subtree_root.get_children()[1].dist = 0
            dodup(subtree_root.get_children()[0], 0, duprate)
            dodup(subtree_root.get_children()[1], 0, duprate)
    else:
        dodup(subtree_root.get_children()[0], subtree_root.get_children()[0].dist, duprate)
        dodup(subtree_root.get_children()[1], subtree_root.get_children()[1].dist, duprate)


def count_leaves(tree):
    leaf_count = {}
    for leaf in tree.iter_leaves():
        leaf_name = leaf.name
        leaf_count[leaf_name] = leaf_count.get(leaf_name, 0) + 1
    return leaf_count


def do_loss(gene_tree, leaf_name, pr, numerator, numerator_fix, denominator, denominator_fix):
    current_pr = pr - numerator
    leaves = list(gene_tree.iter_leaves())
    random.shuffle(leaves)
    for leaf in leaves:
        if leaf_name == leaf.name:
            random_pr = random.random()
            if random_pr <= current_pr / pr:
                leaf.detach()
                if current_pr != 0 and numerator_fix == 0:
                    current_pr = current_pr - numerator
                if pr != 0 and denominator_fix == 0:
                    pr = pr - denominator


def rename_internal_nodes(tree):
    counter = 1  
    for node in tree.traverse("preorder"):
        if not node.is_leaf():
            node_name_parts = node.name.split('_')
            if len(node_name_parts) == 2:
                # Rename nodes with two words
                node.name = f"{counter}_{node_name_parts[0]}_{node_name_parts[1]}"
                counter += 1
            elif len(node_name_parts) == 1:
                # Rename nodes with one word
                node.name = f"{counter}_{node_name_parts[0]}_Sp"
                counter += 1
        elif node.is_leaf():
            if "'" in node.name:
                node.name = node.name.replace("'", "")
                node.name = f"{node.name}_0_0"
            else:
                parent = node.up
                parent1 = node.up
                node.detach()
                while parent1.is_leaf():
                    parent1 = parent.up
                    parent.detach()
                    parent = parent1


def remove_single_child_internal_nodes(tree):
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and len(node.children) == 1:
            child = node.children[0]
            parent = node.up
            # Remove the internal node
            node.detach()
            # Link the parent to the child
            parent.add_child(child)


def generate_gene_trees(species_tree, num_gene_trees):
    gene_trees = []
    for _ in range(num_gene_trees):
        gene_tree = species_tree.copy()
        root = gene_tree.get_tree_root()
        dodup(root, root.dist, duprate)
        leaf_counts = count_leaves(gene_tree)
        for leaf_name, count in leaf_counts.items():
            do_loss(gene_tree, leaf_name, count, numerator, numerator_fix, denominator, denominator_fix)
            # print(f"{leaf_name}: {count}")
        rename_internal_nodes(gene_tree)
        remove_single_child_internal_nodes(gene_tree)
        gene_trees.append(gene_tree)
    return gene_trees


def main():
    species_tree = Tree(species_tree_file, format=1)
    num_nodes_sp = len(species_tree.get_tree_root().get_descendants()) + 1
    rate = 1/(num_nodes_sp*num_gene_trees)
    change_dup_rates(species_tree, rate)
    #print(species_tree.write())
    gene_trees = generate_gene_trees(species_tree, num_gene_trees)

    with open("rated_" + species_tree_file, "w") as output_file:
            output_file.write(species_tree.write(format=1, format_root_node=lambda node: f"[&name={node.name}]") + "\n")

    with open("gene_trees_" + species_tree_file, "w") as output_file:
        for i, gene_tree in enumerate(gene_trees, start=1):
            # output_file.write(f"> GeneTree_{i}\n")
            output_file.write(gene_tree.write(format=1, format_root_node=lambda node: f"[&name={node.name}]") + "\n")


if __name__ == "__main__":
    main()
