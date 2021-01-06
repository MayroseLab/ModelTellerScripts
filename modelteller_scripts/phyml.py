import os,sys

import itertools

sys.path.append("/groups/itay_mayrose/shiranabad/MsMl/code/")
__author__ = 'Shiran'

from utils import is_file_empty#, untargz_files_by_pattern, targz_files_by_pattern
# from getters import *
# from data_processing.traverse_data_dirs import traverse_data_dirs
# from utilities import createJobFile
from defs import *


############################### additional parameters ###############################
PHYML_PINV_TAGS = {True: "-v e",
				   False: ""}
PHYML_GAMMA_TAGS = {True: "-a e -c 4",
					False: "-c 1"}
PHYML_OPT_TAGS = {"ml": "-o tlr -s NNI",
				  "bionj": "-o lr",
				  "fixed": "-o lr",
				  "rates": "-o r",
				  "noopt": "-o n"}
################################# model parameters #################################
PHYML_BASE_FREQS_TAGS = {True: "-f m",
						 False: "-f 0.25,0.25,0.25,0.25"}
PHYML_SUBS_RATES_TAGS = {1: "-m 000000",
						 2: "-m 010010",
							 3: "-m 012345"}
PHYML_MODEL_TAGS = {"JC": [PHYML_SUBS_RATES_TAGS[1], PHYML_BASE_FREQS_TAGS[False]],
					"F81": [PHYML_SUBS_RATES_TAGS[1], PHYML_BASE_FREQS_TAGS[True]],
					"K80": [PHYML_SUBS_RATES_TAGS[2], PHYML_BASE_FREQS_TAGS[False]],
					"HKY": [PHYML_SUBS_RATES_TAGS[2], PHYML_BASE_FREQS_TAGS[True]],
					"SYM": [PHYML_SUBS_RATES_TAGS[3], PHYML_BASE_FREQS_TAGS[False]],
					"GTR": [PHYML_SUBS_RATES_TAGS[3], PHYML_BASE_FREQS_TAGS[True]]}

###################################### general ######################################
PHYML_GENERAL_TAGS = "-d nt -n 1 -b 0 --no_memory_check"
PHYML_SCRIPT = "/groups/itay_mayrose/shiranabad/applications/jmodeltest-2.1.7/exe/phyml/PhyML_3.0_linux64"


def create_phyml_exec_line(msa_file_full_path, base_model, pinv, gamma, topology="ml", tree_file=None, run_id=None):
	run_id = (base_model + ("+I" if pinv else "") + ("+G" if gamma else "")) if run_id is None else run_id
	execution_tags = " ".join(PHYML_MODEL_TAGS[base_model] +
							  [PHYML_PINV_TAGS[pinv], PHYML_GAMMA_TAGS[gamma],
							   PHYML_OPT_TAGS[topology], PHYML_GENERAL_TAGS, "--run_id " + run_id])
	if tree_file:
		execution_tags += " -u " + tree_file

	return " ".join([PHYML_SCRIPT, "-i", msa_file_full_path, execution_tags])


def create_phyml_exec_line_full_model(msa_file_full_path, full_model, topology="ml", tree_file=None, run_id=None):
	run_id = full_model if run_id is None else run_id
	pinv = "+I" in full_model
	gamma = "+G" in full_model
	base_model = re.sub(r'\+.*', '', full_model)

	return create_phyml_exec_line(msa_file_full_path, base_model, pinv, gamma, topology, tree_file, run_id)


def run_phyml(msa_filepath, full_model, topology="ml", tree_file=None, run_id=None, output_dir=None):
	"""
	:param msa_filepath:
	:param full_model: e.g., JC, HKY+I+G
	:param topology: "ml" or "bionj"
	:param tree_file:
	:param run_id: if not given, take the full model
	:return: the stats/tree output filepath (i.e., msa_filepath+"_phyml_stats/tree_"+run_id+".txt"
	"""
	run_id = full_model if run_id is None else run_id
	phyml_exec_line = create_phyml_exec_line_full_model(msa_filepath, full_model, topology, tree_file, run_id)
	output_filename = msa_filepath + "_phyml_{}_" + run_id + ".txt"

	stats_file, tree_file = output_filename.format("stats"), output_filename.format("tree")
	final_stats_file, final_tree_file = stats_file, tree_file
	if output_dir:
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)
		final_stats_file = output_dir + stats_file.split("/")[-1]
		final_tree_file = output_dir + tree_file.split("/")[-1]
	if is_file_empty(stats_file) and is_file_empty(final_stats_file):
		os.system(phyml_exec_line)

	if is_file_empty(final_stats_file) and output_dir:
		os.rename(stats_file, final_stats_file)
		os.rename(tree_file, final_tree_file)

	return final_stats_file, final_tree_file


def parse_phyml_stats_file(phyml_stats_filepath):
	"""
	:param dirpath: where phylip and phyml stats outputs are located
	:param model: could be: JC, F81, K80, HKY, SYM, GTR ( and then '+I+G' or each of them)
	:return: dictionary with the attributes - string typed. if parameter was not estimated, empty string
	"""

	def get_frequency(nuc1, nuc2, gtr_file_reader):
		freq = re.search("(?<=  " + nuc1 + " <-> " + nuc2 + ") {2,4}[0-9\.]*", gtr_file_reader).group(0)
		freq = re.sub(" ", "", freq)
		return float(freq)

	def read_matrix(gtr_file_reader):
		mat = re.search("\. Instantaneous rate matrix :\s+\[A-+C-+G-+T-+\](\s+.*[0-9]\.[0-9]{5} *){4}$",
						gtr_file_reader, re.M).group(0)
		matrix_reader = re.findall("-?[0-9]\.[0-9]{5}", mat, re.M)
		matrix = [[], [], [], []]
		for i in range(0, 16):
			element = matrix_reader[i]
			matrix[i // 4].append(float(element))
		return matrix

	def search_for_value(lookup_string, file_reader):
		value = re.search("(?<= " + lookup_string + ")[0-9\.\-]+", file_reader)
		if value:
			return value.group(0)
		else:
			return ""

	res_dict = dict.fromkeys(["parsimony", "tree_size",
							  "fA", "fC", "fG", "fT",
							  "pInv", "gamma", "Tstv",
							  "mu_rate"], "")

	#open file if exists
	with open(phyml_stats_filepath, 'r') as fpr:
		phyml_content = fpr.read()

	# general measures
	res_dict["parsimony"] = search_for_value("Parsimony: \t{4}", phyml_content)
	res_dict["tree_size"] = search_for_value("Tree size: \t{4}", phyml_content)
	res_dict["gamma"] = search_for_value("Gamma shape parameter: \t{2}", phyml_content)
	res_dict["pInv"] = search_for_value("Proportion of invariant: \t{2}", phyml_content)
	res_dict["Tstv"] = search_for_value("Transition/transversion ratio: \t{1}", phyml_content)
	res_dict["logL"] = search_for_value("Log-likelihood: \t{3}", phyml_content)

	# nucleotide frequencies
	res_dict["fA"] = float(re.search("(?<=f\(A\)\= )[\.0-9]+",phyml_content).group(0))
	res_dict["fC"] = float(re.search("(?<=f\(C\)\= )[\.0-9]+",phyml_content).group(0))
	res_dict["fG"] = float(re.search("(?<=f\(G\)\= )[\.0-9]+",phyml_content).group(0))
	res_dict["fT"] = float(re.search("(?<=f\(T\)\= )[\.0-9]+",phyml_content).group(0))
	nuc_freq = [res_dict["fA"], res_dict["fC"], res_dict["fG"], res_dict["fT"]]

	# relative rate parameters
	freq_AC = get_frequency("A", "C", phyml_content)
	freq_AG = get_frequency("A", "G", phyml_content)
	freq_AT = get_frequency("A", "T", phyml_content)
	freq_CG = get_frequency("C", "G", phyml_content)
	freq_CT = get_frequency("C", "T", phyml_content)
	freq_GT = get_frequency("G", "T", phyml_content)

	res_dict["rel_subAC"], res_dict["rel_subAG"], res_dict["rel_subAT"], res_dict["rel_subCG"], res_dict["rel_subCT"], res_dict["rel_subGT"] =\
		freq_AC, freq_AG, freq_AT, freq_CG, freq_CT, freq_GT

	subs = [[None, freq_AC, freq_AG, freq_AT],
			[None, None, freq_CG, freq_CT],
			[None, None, None, freq_GT]]
	matrix = read_matrix(phyml_content)

	i = 0
	j = 1
	while (matrix[i][j] == 0) and (i < 4): # get to the first rate that != 0
		if j == 3:
			i += 1
			j = i+1
		else:
			j += 1

	# relevant mainly to GTR
	res_dict["mu_rate"] = matrix[i][j]/float(nuc_freq[j]*subs[i][j])
	sub_rates = zip(["sub" + x + y for x, y in itertools.product('ACGT', repeat=2)], matrix[0] + matrix[1] + matrix[2] + matrix[3])
	sub_rates = list(sub_rates)
	res_dict.update(sub_rates)

	return res_dict


def get_phyml_running_time(phyml_stats_filepath):
	try:
		with open(phyml_stats_filepath, 'r') as fpr:
			phyml_content = fpr.read()
		re_matches = re.search("\((\d+) seconds\)", phyml_content)
		return int(re_matches.group(1))
	except:
		return


def run_phyml_for_all_models(msa_filepath, dest_path, model=None, fixed_treefile=None):
	"""
	:param msa_filepath: full path to msa file
	:param dest_path: destination to save phyml output files
	:param model: a model, a list of models, or all models if None
	:param fixed_treefile: if provided, fixed topology
	:return:
	"""
	_, msa_filename = os.path.split(msa_filepath)
	if model in ALL_PHYML_MODELS:
		run_on_models = [model]
	elif isinstance(model, list):
		run_on_models = model
	else:
		run_on_models = ALL_PHYML_MODELS
	final_output_pattern = SEP.join([dest_path, msa_filename + "_phyml_tree_{}.txt"])
	for full_model in run_on_models:
		final_output = final_output_pattern.format(full_model)
		if is_file_empty(final_output):
			stats_filename, tree_filename = run_phyml(msa_filepath, full_model,
													  topology="ml" if (fixed_treefile is None) else "fixed",
													  tree_file=fixed_treefile)

			os.rename(tree_filename, final_output)
			os.rename(stats_filename, re.sub("tree", "stats", final_output))
		else:
			print("phyml exists", final_output)
	return re.sub("tree", "stats", final_output_pattern), final_output_pattern