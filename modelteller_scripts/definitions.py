#imports
import sys, logging, os, re


SEP = "/"
RESULTS_HEAD_DIR = "/bioseq/data/results/modelteller/"
PICKLE_PATH = "/bioseq/modelteller/Bs_predictor.pkl"

global RESULT_MSG
RESULT_MSG = ""

PHYLIP_FORMAT = "phylip-relaxed"
BASE_MODELS = ["JC", "F81", "K80", "HKY", "SYM", "GTR"]
MODELS_TAGS = ["", "+I", "+G", "+I+G"]
ALL_PHYML_MODELS = [base_model + tag for base_model in BASE_MODELS for tag in MODELS_TAGS]

BASE_TREE_TYPE_BIONJ = "BIONJ"
BASE_TREE_TYPE_ML = "ML"

PHYML_STATS_SUFFIX = "_phyml_stats{0}.txt"
PHYML_TREE_SUFFIX = "_phyml_tree{0}.txt"
PHYML_LK_SUFFIX = "_phyml_lk.txt"


def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)