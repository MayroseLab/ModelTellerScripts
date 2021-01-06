#imports
import numpy as np
import re
import matplotlib.pyplot as plt
#import seaborn as sns
import pandas as pd
import scipy, os, sys, math, random, time, copy, argparse, platform, subprocess, logging, socket
from shutil import copyfile
import shutil


sys.path.append("/groups/itay_mayrose/shiranabad/MsMl/code/")
__author__ = 'Shiran'


if platform.system() == 'Linux':
	#linux paths
	JMT_APP = "/groups/itay_mayrose/shiranabad/applications/jmodeltest-2.1.7"
	RAW_DATABASES_PATH = "/groups/itay_mayrose/shiranabad/MsMl/databases/"
	DATA_PATH = "/groups/itay_mayrose/shiranabad/MsMl/data/"
	SEP = "/"
	PHYML_APP = "/groups/itay_mayrose/shiranabad/applications/jmodeltest-2.1.7/exe/phyml/PhyML_3.0_linux64"
	CODE_PATH = "/groups/itay_mayrose/shiranabad/MsMl/code/"
	SUMMARY_FILES_DIR = "/groups/itay_mayrose/shiranabad/MsMl/summary_files/"
	FIGURES_DIR = SUMMARY_FILES_DIR + "/figures/"
else:
	#windows paths
	dropbox_path = "D:\\Dropbox\\lab\\"
	JMT_APP = dropbox_path + "\\project\\jmodeltest-2.1.5"
	# DATABASES_PATH = dropbox_path + "\\project\\ppln_test\\"
	DATA_PATH = dropbox_path + "/MsMl/databases/"
	SUMMARY_FILES_DIR = dropbox_path + "/MsMl/summary_files/"
	FIGURES_DIR = dropbox_path + "/MsMl/figures/"
	SEP = "\\"
	PHYML_APP = dropbox_path + "\\project\\code\\modelSelectionCode\\PhyML-3.1_win32.exe"

PLOIDB_DBNAME = "ploiDB"
PROTDB_DBNAME = "protDBs"
SELECTOME_DBNAME = "selectome"
PANDIT_DBNAME = "PANDIT"
ORTHOMAM_DBNAME = "orthoMam"
RFAM_DBNAME = "Rfam"
TREEBASE_DBNAME = "TreeBASE"

PLOIDB_PATH = DATA_PATH + PLOIDB_DBNAME + SEP
PROTDB_PATH = DATA_PATH + PROTDB_DBNAME + SEP
SELECTOME_PATH = DATA_PATH + SELECTOME_DBNAME + SEP
PANDIT_PATH = DATA_PATH + PANDIT_DBNAME + SEP
ORTHOMAM_PATH = DATA_PATH + ORTHOMAM_DBNAME + SEP
RFAM_PATH = DATA_PATH + RFAM_DBNAME + SEP
TREEBASE_PATH = DATA_PATH + TREEBASE_DBNAME + SEP

ORIG_MSA_PHYLIP_FILENAME = "seqs-organism-concat_phy.phy"
MSA_PHYLIP_FILENAME = "real_msa.phy"
MSA_PHYLIP_SEQUENTIAL_FILENAME = "real_msa_sequential.phy"
MSA_NEXUS_FILENAME = "real_msa.nex"
FASTA_FORMAT = "fasta"
PHYLIP_FORMAT = "phylip-relaxed"
NEXUS_FORMAT = "nexus"
BASE_MODELS = ["JC", "F81", "K80", "HKY", "SYM", "GTR"]
MODELS_TAGS = ["", "+I", "+G", "+I+G"]
ALL_PHYML_MODELS = [base_model + tag for base_model in BASE_MODELS for tag in MODELS_TAGS]


############# simulated data
N_SIMULATION = 1
SIMULATED_MSA_FILENAME_FASTA = "simulated_msa.fas" # add .format(re.sub("\+", "-", model))
SIMULATED_MSA_FILENAME_PHYLIP = "simulated_msa.phy" # add .format(re.sub("\+", "-", model))
SIMULATED_MSA_FILENAME_NEXUS = "simulated_msa.nex" # add .format(re.sub("\+", "-", model))
SIMULATED_MSA_FILENAME_PAML = "simulated_msa_paml.phy"

ML_RESULTS_DIR = SUMMARY_FILES_DIR + "ML/"
MLREG_RESULTS_DIR = SUMMARY_FILES_DIR + "ML_REG/"
MLREG_split_RESULTS_DIR = SUMMARY_FILES_DIR + "ML_REG_split/"

############ jModelTest
JMODELTEST_OUTPUT_FILENAME = "jmt.txt"
JMODELTEST_FIXEDTREE_OUTPUT_FILENAME = "jmt_fixed.txt"
JMODELTEST_SUMMARY_FILENAME = "jmt_results.csv"
JMODELTEST_FIXEDTREE_SUMMARY_FILENAME = "jmt_fixedTree_results.csv"
WEIGHTS_FILENAME = "weights_{0}_{1}.csv"
JMODELTEST_FIXEDTREE_WEIGHTS_FILENAME = "jmt_fixedTree_weights_{0}.csv"

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