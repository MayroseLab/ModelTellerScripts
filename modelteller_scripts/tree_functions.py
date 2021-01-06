import itertools
import ete3
from msa_functions import get_msa_from_file, count_substitutions
from utils import compute_entropy, lists_diff
from definitions import *
from ete3 import Tree
import ete3.coretype.tree
import numpy as np


def get_newick_tree(tree):
	"""
	:param tree: newick tree string or txt file containing one tree
	:return:	tree: a string of the tree in ete3.Tree format
	"""
	if type(tree) == str:
		if os.path.exists(tree):
			with open(tree, 'r') as tree_fpr:
				tree = tree_fpr.read().strip()
		tree = Tree(tree, format=1)
	return tree


def rescale_tree_branch_lengths(tree, factor):
	"""
	:param tree: newick tree string or txt file containing one tree
	:param factor: the factor by which to multiply all branch lengths in tree
	:return:	reformatted_tree: a string of the scaled tree in Newick format
	"""
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	for node in tree_root.iter_descendants(): # the root dist is 1.0, we don't want it
		node.dist = node.dist * factor
	return tree.write(format=1, dist_formatter="%.10f")


def reroot_tree(tree, outgroup_name):
	tree = get_newick_tree(tree)
	tree.set_outgroup(tree & outgroup_name)
	return tree


def scale_tree_to_length(tree, target_dist, outgroup_name=None, scaling_method="tbl"):
	"""
	:param tree:  newick tree string or txt file containing one tree OR ete3.Tree object
	:param target_dist: numeric, the desired total tree distance
	:param outgroup_name: the name of a tree node if tree needs to be rooted,
	                      otherwise would calculate the distance from the inferred root (acoorsding to the newick order)
	:param scaling_method: "height" for longest distance from root to leaf, "tbl" for total branch lengths
	:return: a newick string of the rescaled tree
	"""
	t = get_newick_tree(tree)

	if outgroup_name: # re-root tree
		t = reroot_tree(t, outgroup_name)
	root = t.get_tree_root()

	if scaling_method.lower() == "tbl":
		dist = get_total_branch_lengths(root)
	elif scaling_method.lower() == "height":
		dist = get_tree_height(root)

	scaling_factor = target_dist / dist
	rescaled_tree = rescale_tree_branch_lengths(tree, scaling_factor)

	return rescaled_tree


def get_tree_height(tree_root):
	"""
	:param tree_root: ete3 node; because we traverse only its descendants
	(its dist is 1.0, we don't want it)
	:return: longest distance from root to leaf
	"""
	# the two configurations are the same - I compared! (Shiran)
	# current_length = 0
	# for leaf in tree_root:
	# 	current_length = max(current_length, tree_root.get_distance(leaf))
	return tree_root.get_farthest_leaf()[1]


def get_frac_of_cherries(tree):
	"""
	McKenzie, Andy, and Mike Steel. "Distributions of cherries for two models of trees."
	 Mathematical biosciences 164.1 (2000): 81-92.
	:param tree:
	:return:
	"""
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	leaves = list(tree_root.iter_leaves())
	cherries_cnt = 0
	for leaf1, leaf2 in itertools.combinations(leaves, 2):
		if leaf1.up is leaf2.up:
			cherries_cnt += 1
	return 2*cherries_cnt/len(leaves)


def get_leaves_branches(tree):
	"""
	:param tree:
	:return: a list of pendant edges lengths, i.e., the brnches that lead to leaves
	"""
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	leaves_bl = []
	for node in tree_root.iter_leaves():
		leaves_bl.append(node.dist)
	return leaves_bl


def get_sackin_index(tree):
	"""
	:param tree: ete3 tree; not rerooted!
	:return: sackin_index
	formula Sackin: https://academic.oup.com/sysbio/article/21/2/225/1713481
	"""
	node_depth_dict = {}
	sackin_index_lst = []
	for node in tree.traverse(strategy="levelorder"):
		if node.is_root():
			node_depth_dict[node] = 1
		else:
			node_depth_dict[node] = node_depth_dict[node.up] + 1

		if node.is_leaf():
			sackin_index_lst.append(node_depth_dict[node])

	return np.mean(sackin_index_lst)


def get_colless_index(tree):
	"""
	:param tree: ete3 tree; not rerooted!
	:return: colless_index
	formula Colless: https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1558-5646.1992.tb01171.x
	"""
	ndes_dict = {}
	colless_index = 0
	ntaxa = 0
	for node in tree.traverse(strategy="postorder"):
		if node.is_leaf():
			ndes_dict[node] = 1
			ntaxa += 1
		else:
			ndes_dict[node] = ndes_dict[node.children[0]] + ndes_dict[node.children[1]]
			colless_index += abs(ndes_dict[node.children[0]] - ndes_dict[node.children[1]])

	return colless_index/float((ntaxa-1)*(ntaxa-2)/2)


def get_stemminess_indexes(tree):
	"""
		:param tree: ete3 tree; not rerooting!
		:return: cumulative stemminess index Fiala and Sokal 1985
				 noncumulative stemminess index Rohlf 1990
		formula cumulative stemminess: https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1558-5646.1985.tb00398.x
		formula noncumulative stemminess: https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1558-5646.1990.tb03855.x
		"""
	subtree_blsum_dict = {}
	nodes_height_dict = {}
	stem85_index_lst = []
	stem90_index_lst = []
	for node in tree.traverse(strategy="postorder"):
		if node.is_leaf():
			subtree_blsum_dict[node] = 0
			nodes_height_dict[node] = 0
		elif node.is_root():
			continue
		else:
			subtree_blsum_dict[node] = subtree_blsum_dict[node.children[0]] + subtree_blsum_dict[node.children[1]] + \
			                           node.children[0].dist + node.children[1].dist
			nodes_height_dict[node] = max(nodes_height_dict[node.children[0]] + node.children[0].dist,
			                              nodes_height_dict[node.children[1]] + node.children[1].dist)
			stem85_index_lst.append(node.dist/(subtree_blsum_dict[node] + node.dist))
			stem90_index_lst.append(node.dist/(nodes_height_dict[node]) + node.dist)

	return np.mean(stem85_index_lst), np.mean(stem90_index_lst)


def pendant_to_internal_bl_ratio(tree):
	"""
	:param tree:
	:return: sort of stemminess index, i.e., how long are the inner branches compared to the external ones
	"""
	tree = get_newick_tree(tree)
	pendant_bl_sum = np.sum(get_leaves_branches(tree))
	internal_bl_sum = get_total_branch_lengths(tree) - pendant_bl_sum

	if pendant_bl_sum == 0:
		return 0
	if internal_bl_sum == 0:
		return 100000
	return pendant_bl_sum/float(internal_bl_sum)


def get_branch_lengths(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return: total branch lengths
	"""
	# TBL
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	branches = []
	for node in tree_root.iter_descendants(): # the root dist is 1.0, we don't want it
		branches.append(node.dist)
	return branches


def get_branch_lengths_estimates(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return:
	"""
	# TBL
	branches = get_branch_lengths(tree)
	entropy = compute_entropy(branches)

	return max(branches), min(branches), np.mean(branches), np.std(branches), entropy


def get_diameters_estimates(tree_filepath, actual_bl=True):
	"""
	if not actual_bl - function changes the tree! send only filepath
	:param tree_filepath: tree file or newick tree string;
	:param actual_bl: True to sum actual dists, False for num of branches
	:return: min, max, mean, and std of tree diameters
	"""
	# tree = copy.deepcopy(get_newick_tree(tree)) # do not deepcopy! when trees are large it exceeds recursion depth
	if not actual_bl:
		assert isinstance(tree_filepath, str)
	tree = get_newick_tree(tree_filepath)
	tree_root = tree.get_tree_root()
	if not actual_bl:
		for node in tree_root.iter_descendants():
			node.dist = 1.0
	tree_diams = []
	leaves = list(tree_root.iter_leaves())
	for leaf1, leaf2 in itertools.combinations(leaves, 2):
		tree_diams.append(leaf1.get_distance(leaf2))
	entropy = compute_entropy(tree_diams)

	return max(tree_diams), min(tree_diams), np.mean(tree_diams), np.std(tree_diams), entropy


def get_distance_from_tips(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return: mean and std of terminal branches
	"""
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	terminal_branches = []
	for leaf in tree_root.iter_leaves():
		terminal_branches.append(leaf.dist)
	return np.mean(terminal_branches), np.std(terminal_branches)


def get_total_branch_lengths(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return: total branch lengths
	"""
	tree = get_newick_tree(tree) #if it's a file
	branches = get_branch_lengths(tree)
	return sum(branches)


def rename_ids(tree_str, conversion_dict):
	"""
	:param tree_str:  Newick format tree string
	:param conversion_dict: {current_id: new_id}
	:return: an updated tree string
	"""
	for sp in conversion_dict:
		tree_str = re.sub(sp + ":", conversion_dict[sp] + ":", tree_str)
	return tree_str


def calc_branch_length(tree):
	branch_lengths = get_branch_lengths(tree)
	N = len(branch_lengths)
	total_BL = sum(branch_lengths)
	mean_BL = total_BL/N

	return total_BL, mean_BL


def assign_simulated_seqs_to_nodes(tree, msa_file, ancestral_seqs_file):
	"""
	:param tree: tree file or ete3 tree
	:param msa_file: sequences names should match tree leaves
	:param ancestral_seqs_file: sequences names should match tree inner nodes + root
	:return: ete3 tree, with feature "seq" inside each node
	To extract this value, use: node.seq
	"""
	tree = get_newick_tree(tree)
	msa = get_msa_from_file(msa_file)

	for rec in msa:
		rec_name = rec.id
		(tree & rec_name).add_feature("seq", rec._seq)
	for line in open(ancestral_seqs_file):
		rec = line.strip().split("\t")
		(tree & rec[0]).add_feature("seq", rec[1])

	return tree


def assign_simulated_seqs_to_indelible_tree(sim_dir, simulated_msa_filename):
	"""
	:param sim_dir: the directory in which indelible was simulated
	:param simulated_msa_filename: the name of the simulated msa file (to extract the simulated msa + ancestral seqs)
	:return: ete3 tree, with feature "seq" inside each node
	To extract this value, use: node.seq
	"""
	sim_msa_file = sim_dir + simulated_msa_filename
	trees_file = sim_dir + "trees.txt"
	ancestral_seqs_file = sim_dir + simulated_msa_filename[:-4] + "_ANCESTRAL.phy"
	tree = indelible_executions.get_indelible_tree_from_treesfile(trees_file)

	return assign_simulated_seqs_to_nodes(tree, sim_msa_file, ancestral_seqs_file)


def copy_nodes_names_from_similar_topology(ref_topology, tree):
	tree = get_newick_tree(tree)
	ref_topology = get_newick_tree(ref_topology)
	for node in tree.traverse("postorder"):
		if node.is_leaf():
			continue
		children = node.children
		node_in_ref = (ref_topology & children[0].name).up

		# assert that nodes in ref and in tree have the same children
		ref_children = node_in_ref.children
		assert ".".join(sorted([child.name for child in children])) == ".".join(sorted([child.name for child in ref_children]))

		node.name = node_in_ref.name
	return tree


def get_cumulative_num_of_substitutions(tree_with_seqs):
	"""
	:param tree_with_seqs: ete3 tree in which every node is assigned with a 'seq' feature
	:return: cumulative number of sequences along the tree (sum over average number of substitutions per branch)
	 over nucleotides and codons substitutions
	"""
	cum_differences_nuc = 0
	cum_differences_codon = 0

	for node in tree_with_seqs.traverse("postorder"):
		if node.name != "ROOT":
			seq1 = node.seq
			seq2 = node.up.seq
			cum_differences_nuc += count_substitutions(seq1, seq2, is_codons=False)
			try:
				cum_differences_codon += count_substitutions(seq1, seq2, is_codons=True)
			except AssertionError:
				cum_differences_codon = -1
		# bl_codons = node.dist

	return cum_differences_nuc, cum_differences_codon


def convert_hexadec_branches_to_decimal(tree_str):
	reformatted_number = re.compile(r'\d+[\.\d+]*e\-\d+') # reformat small numbers
	for number in reformatted_number.findall(tree_str):
		reformatted_num = str("{0:.10f}".format(float(number)))
		tree_str = tree_str.replace(number, reformatted_num, 1)

	return tree_str


def get_internal_and_external_leaves_relative_to_subroot(tree_root, subroot):
	all_leaves = tree_root.get_leaves()
	subtree_leaves = subroot.get_leaves()
	other_leaves = lists_diff(all_leaves, subtree_leaves)
	return subtree_leaves, other_leaves


def prune_tree_at_node(tree, node_as_subroot):
	assert(isinstance(node_as_subroot, ete3.TreeNode))
	subtree_leaves, other_leaves = get_internal_and_external_leaves_relative_to_subroot(tree, node_as_subroot)
	t1 = tree.prune(subtree_leaves, preserve_branch_length=True)
	t2 = tree.prune(other_leaves, preserve_branch_length=True)
	return t1, t2


def get_largest_branch(tree):
	tree = get_newick_tree(tree)
	tree_nodes = list(tree.traverse("levelorder"))
	max_bl_node = max(tree_nodes, key=lambda node: node.dist)
	return max_bl_node
