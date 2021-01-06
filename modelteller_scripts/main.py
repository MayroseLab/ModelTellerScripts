import argparse
import shutil
import traceback
from time import sleep

from Bio import AlignIO
import Bio.Align.AlignInfo
import pandas as pd

import phyml
import tree_functions
from msa_functions import make_msa_generic
import numpy as np
from web_utils import *
import pickle
import sklearn

sys.path.append('/bioseq/modelteller/auxiliaries/')
from html_editor import edit_results_html, notify_by_email, post_html_editing, show_final_tree_in_html, show_features_contribution_in_html
import compute_features
import h2o
from definitions import  *

MAPPING_FEATURE_NAMES_FILE = "/bioseq/modelteller/modelteller_scripts/features_labels_mapping.csv" 
CODE_PATH = "/bioseq/modelteller/modelteller_scripts/"
global RESULT_MSG

def get_ml_features_to_include():
	with open(ML_FEATURES_FILE) as fpr:
		features = (fpr.read().strip()).split(",")
	# features = features.tolist() + ["model_matrix", "model_F", "model_I", "model_G"]
	return features


def read_alignment_file(msa_file, user_tree_file):
	"""
	:param msa_file: the path to an MSA file, one of biopython's formats
	:return: msa_objs_list - an iterator of the alignments (AlignIO objects) or None if the file is corrupted or the format is invalid
	"""
	msa_objs_list = None
	# identify format and retrieve all MSAs
	for aln_format in ["clustal", "emboss", "fasta", "fasta-m10", "ig", "maf", "mauve", "nexus", "phylip-relaxed", "phylip-sequential", "stockholm"]:
		try:
			msa_objs_list = list(AlignIO.parse(msa_file, format=aln_format))
			if len(msa_objs_list) == 0:
				raise Exception()
			logger.info("The MSA file is format: " + aln_format)
			break
		except Exception:
			msa_objs_list = None
	if msa_objs_list is None:
		set_results_message("File format isn't supported or sequences are not in the same length. Please upload a valid MSA format.")
		raise Exception()

	# validate all MSAs
	for i, msa in enumerate(msa_objs_list):
		msa_info = Bio.Align.AlignInfo.SummaryInfo(msa)
		aln_letters = msa_info._get_all_letters()
		for let in aln_letters:
			if not (let.lower() in "acgt-"):
				set_results_message("Warning! Alignment number " + str(i + 1) + " consists of characters that are not A,C,G,T,-\n")
				break
		if len(aln_letters) > 10: # say, acgt, gap, RYN
			print(aln_letters)
			set_results_message("There are too many non-nucleotide characters in your alignment. Please validate your input.")
			raise Exception()

	tree_list = None
	if user_tree_file:
		with open(user_tree_file) as fpr:
			tree_list = fpr.read().strip().splitlines(keepends=False)
		# First assert that the trees are valid
		new_tree_list = []
		for i, tree in enumerate(tree_list): #Here, the line numbering matches those in the input file
			if re.search("\S", tree):
				try:
					new_tree_list.append(tree_functions.get_newick_tree(tree.strip()))
				except:
					set_results_message("In line {} there is no valid Newick tree".format(i))
		tree_list = new_tree_list
		# assert that the number of trees matches the number of alignments
		if len(msa_objs_list)!= len(tree_list):
			set_results_message("Number of trees does not match number of alignments")
			raise Exception()
		# assert that each tree matches the corresponding MSA
		for i, tree_obj in enumerate(tree_list):
			leaves = sorted([node.name for node in tree_obj.get_leaves()])
			seq_names = sorted([rec.id for rec in msa_objs_list[i]])
			if len(leaves) != len(seq_names) or (not all(x == y for x,y  in zip(seq_names,leaves))):
				set_results_message("The tips of the tree and the MSA sequences names do not match" + ((" for MSA number " + str(i+1)) if (i>0 and len(msa_objs_list)>1) else ""))
				raise Exception()

	return msa_objs_list, tree_list


def build_execution_dirtree(msa_objs_list, user_trees_list, results_path):
	"""
	:param msa_objs_list: an iterator of the alignments (AlioIO objects)
	:param results_path:
	saves each msa in iterator in its own directory (dirname is numeric, chronologically) by the name "msa.phy"
	"""
	exec_dir = results_path + "data/"
	if not os.path.exists(exec_dir):
		os.mkdir(exec_dir)
	msas_list = []
	tree_files_list = []

	for i, msa in enumerate(msa_objs_list):
		if not os.path.exists(exec_dir + str(i)):
			os.mkdir(exec_dir + str(i))

		## write msa to file
		make_msa_generic(msa, dest_filename=exec_dir + str(i) + SEP + "msa.phy")
		msas_list.append((msa, exec_dir + str(i) + SEP + "msa.phy"))

		## write tree to file
		if user_trees_list:
			with open(exec_dir + str(i) + SEP + "tree.txt", "w") as fpw:
				fpw.write(user_trees_list[i].write(format=1, dist_formatter="%.10f"))
			tree_files_list.append(exec_dir + str(i) + SEP + "tree.txt")

	if len(tree_files_list) == 0:
		tree_files_list = None

	return msas_list, tree_files_list


def init_execution(msa_file, results_path, user_tree_file):
	"""
	:param msa_file: A file with an alignment or several ones of the same format (unknow format)
	:param results_path: job results directory
	"""
	msa_objs_list, user_tree_list = read_alignment_file(msa_file, user_tree_file) # list of alignIO objs
	# write every msa/tree to a different directory for additional analysis
	msas_list, tree_list = build_execution_dirtree(msa_objs_list, user_tree_list, results_path)
	return msas_list, tree_list


def report_results(status_ok=True, results_path="", job_name="", msg=''):

	i = 0
	msas_results = []
	# aggregate all results paths to a list
	while os.path.exists(f"{results_path}/msa_{i}.csv"):
		msas_results.append(f"{results_path}/msa_{i}.csv")
		i += 1

	if not msas_results: # no results
		if msg == "":
			msg = "Unknown error occurred. Please contact us for help. We apologize for the inconvenience."
		else:
			msg = "Error occurred: " + msg +". Please fix the input and submit again."
		edit_results_html(False, msas_results, f"{results_path}/output.html", run_number=job_name, msg=msg)
        
        # Josef: send failed email
		if os.path.exists(results_path + "user_email.txt"):
			online_res_path = f"https://modelteller.tau.ac.il/results/{job_name}/output.html"                 
			with open(results_path + "user_email.txt") as fpr:
				addressee = fpr.read().strip()
			notify_by_email(addressee, status_ok, online_res_path, msg=msg + "\nCheck out faliure information at: " + online_res_path)
                            
	else: # at least one result was detected
		edit_results_html(status_ok, msas_results, f"{results_path}/output.html", run_number=job_name, msg=msg)

		if os.path.exists(results_path + "user_email.txt"):
			online_res_path = f"https://modelteller.tau.ac.il/results/{job_name}/output.html"                 
			with open(results_path + "user_email.txt") as fpr:
				addressee = fpr.read().strip()
			notify_by_email(addressee, status_ok, online_res_path, msg=msg + "\nCheck out your results at: " + online_res_path
			                + ".\nModelTeller will now reconstruct a maximum-likelihood tree using the chosen model. This might take some time.")


def predict_h2o(ext_df, features_to_include):
	# in case another instance of h2o is currently running
	flag_file = CODE_PATH + "H2O_currently_running.txt"
	while os.path.exists(flag_file):
		sleep(15)
	fpr = open(flag_file, "w")
	fpr.flush()
	fpr.close()

	drf_model_path = MODEL_PATH

	try:
		h2o.init(max_mem_size="8G")
		covtype_df = h2o.H2OFrame(ext_df[features_to_include])
		drf_model = h2o.load_model(drf_model_path)

		ext_df["pred_Bs"] = drf_model.predict(covtype_df).as_data_frame()
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
	finally:
		h2o.cluster().shutdown(prompt=False)
		os.remove(flag_file)
	return


def predict_sklearn(ext_df, features_to_include):
	with open(MODEL_PATH, "rb") as pklr:
		clf = pickle.load(pklr)

	try:
		ext_df["pred_Bs"] = clf.predict(ext_df[features_to_include])
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)

	return


def predict_sklearn_wTreeInterpreter(ext_df, features_to_include):
	with open(MODEL_PATH, "rb") as pklr:
		clf = pickle.load(pklr)

	try:
		from treeinterpreter import treeinterpreter
		ext_df["pred_Bs"], bias, contribution = treeinterpreter.predict(clf, ext_df[features_to_include])

		return contribution
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)

	return


def main(results_path, msas_list, run_mode, user_tree_list, prediction_function):
	"""
	:param msas_list:
	:param run_mode: 0 for regular ModelTeller, 1 for a fixed GTR+I+G tree, 2 for a user tree
	:param user_tree_list: if run_mode=2
	:return:
	"""
	ext_df = compute_features.prepare_features_df(msas_list, run_mode, user_tree_list)
	ext_df.to_csv(results_path + "ext_df_before.csv")
	# predict_h2o(ext_df, get_ml_features_to_include())
	features_list = get_ml_features_to_include()
	feature_contribution_matrix = prediction_function(ext_df, features_list)

	probs_df = ext_df.pivot(index='index', columns='model', values='pred_Bs')[ALL_PHYML_MODELS]
	ranked_df = pd.DataFrame.rank(probs_df, axis=1, method="min")
	ext_df["model_rank"] = ranked_df.stack().values

	mapping_dict = pd.read_csv(MAPPING_FEATURE_NAMES_FILE, index_col=0).to_dict()['label']
	mapping_dict = {x.replace("GTR+I+G_", ""): mapping_dict[x] for x in mapping_dict}
	ext_df.rename(mapper=mapping_dict, axis="columns", inplace=True)

	ext_df.to_csv(results_path + "ext_df_aftermapping.csv")
	for i in range(ext_df["index"].max()+1):
		df_res = ext_df.loc[ext_df["index"]==i, ["model_rank", "model"]]
		df_res.sort_values(by="model_rank", inplace=True)
		df_res["model"].reset_index(drop=True).to_csv(results_path + "msa_" + str(i) + ".csv", header=False, index=False)

	try:
		features_contribution_df = pd.DataFrame(feature_contribution_matrix, ext_df["model"], columns=features_list)
		features_contribution_df.rename(mapper=mapping_dict, axis="columns", inplace=True)
		features_contribution_df.abs().to_csv(results_path + "feature_contribution.csv")
	except:
		pass # feature_contribution_matrix is None :)


def reconstruct_ml_trees(results_path):
	i = 0
	cur_dir = f"{results_path}data/{i}/"
	while os.path.exists(cur_dir):
		user_tree_file = cur_dir + "tree.txt"
		tree_topology_mode = "fixed"
		if not os.path.exists(user_tree_file):
			user_tree_file = None
			tree_topology_mode = "ml"
		with open(f"{results_path}msa_{i}.csv") as fpr:
			selected_model = fpr.readline().strip()
		opt_phyml_stats_filepath, opt_phyml_tree_filepath = phyml.run_phyml(cur_dir + "msa.phy", selected_model, topology=tree_topology_mode,
		                                                                    tree_file=user_tree_file)
		shutil.copyfile(opt_phyml_tree_filepath, f"{results_path}output_tree_{i}.txt")
		show_final_tree_in_html(f"{results_path}output_tree_{i}.txt", f"{results_path}output.html", i=i)

		i += 1
		cur_dir = f"{results_path}data/{i}/"


if __name__ == '__main__':
	logger = logging.getLogger('ModelTeller main script')
	init_commandline_logger(logger)

	parser = argparse.ArgumentParser(description='run modelteller online')
	parser.add_argument('--msa_filepath', '-m', default=None, help='A file with an alignment or several ones of the same format')
	parser.add_argument('--job_name', '-j', default=None)
	parser.add_argument('--run_mode', '-p') # 0 for regular ModelTeller, 1 for a fixed GTR+I+G tree, 2 for a user tree
	parser.add_argument('--feature_contribution', '-f') # 0 for no, 1 for yes
	parser.add_argument('--user_tree_file', '-u', default=None) # if p=2
	args = parser.parse_args()

	global JOB_NAME
	global RESULTS_PATH
	global MODEL_PATH
	global ML_FEATURES_FILE
	JOB_NAME = args.job_name
	RESULTS_PATH = RESULTS_HEAD_DIR + JOB_NAME + SEP

	run_mode = args.run_mode
	compute_feature_contribution = bool(int(args.feature_contribution))
	user_tree_file = args.user_tree_file
	assert run_mode in "012"
	if run_mode == "2":
		assert os.path.exists(user_tree_file)

	if False: # for h2o
		if run_mode != "1":
			MODEL_DIR_PATH = "/groups/itay_mayrose/shiranabad/MsMl/summary_files/ML_REG/GTRagain_nodecimals_trees50_simulated0_10fold_ALLFEATURES_h2oFilFeat__BsrankNone_RFNoneNone_ml/"
		else:
			MODEL_DIR_PATH = "/groups/itay_mayrose/shiranabad/MsMl/summary_files/ML_REG/GTRagain_nodecimals_trees50_simulated0_10fold_ALLFEATURES_h2oFilFeat_ml_gtrig_BsrankNone_RFNoneNone_fixed-GTRIG/"
		MODEL_PATH = MODEL_DIR_PATH + "Bs_model/Bs_model"
		ML_FEATURES_FILE = MODEL_DIR_PATH + "Bs_features_order.txt"
		prediction_func = predict_h2o
	else: #sklearn
		if run_mode != "1":
			MODEL_DIR_PATH = "/groups/itay_mayrose/shiranabad/MsMl/summary_files/ML_REG/web_verGTRagain_nodecimals_trees50_simulated1_10fold_ImpAnalysisFilFeat100_rates_gtrig_BsrankNone_RFNoneNone_ml/"
		else:
			MODEL_DIR_PATH = "/groups/itay_mayrose/shiranabad/MsMl/summary_files/ML_REG/web_verGTRagain_nodecimals_trees50_simulated1_10fold_ImpAnalysisFilFeat100_rates_gtrig_BsrankNone_RFNoneNone_fixed-GTRIG/"
		MODEL_PATH = MODEL_DIR_PATH + "Bs_predictor.pkl"
		ML_FEATURES_FILE = MODEL_DIR_PATH + "Bs_features_order.txt"
		if not compute_feature_contribution:
			prediction_func = predict_sklearn
		else:
			prediction_func = predict_sklearn_wTreeInterpreter

	try:
		msas_list, trees_list = init_execution(args.msa_filepath, RESULTS_PATH, user_tree_file)
		main(RESULTS_PATH, msas_list, run_mode, trees_list, prediction_func)
		report_results(status_ok=True, results_path=RESULTS_PATH, job_name=JOB_NAME, msg=get_results_message() + '\nModelTeller run finished successfully!')
		reconstruct_ml_trees(RESULTS_PATH)
		if compute_feature_contribution:
			show_features_contribution_in_html(RESULTS_PATH + "feature_contribution.csv", f"{RESULTS_PATH}output.html")
		post_html_editing(f"{RESULTS_PATH}/output.html")
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
		report_results(status_ok=False, results_path=RESULTS_PATH, job_name=JOB_NAME,
					   msg=get_results_message())
		# with the following error:\n' + get_error_message())
		# exc_value should be the relevant one to present to user
