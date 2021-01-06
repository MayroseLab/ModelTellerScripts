import sys, logging, os, math,re
from Bio import AlignIO, Alphabet
from definitions import *
import numpy as np


__author__ = 'Shiran'


def change_path_permissions_to_777(path):
	os.chmod(path, 0o777)
	for root, dirs, files in os.walk(path):
		for dir in dirs:
			try:
				os.chmod(os.path.join(root, dir), 0o777)
			except:
				pass
		for file in files:
			try:
				os.chmod(os.path.join(root, file), 0o777)
			except:
				pass


def convert_phylipInterleaved_to_sequential_relaxed(msa_file, output_file):
	with open(msa_file, "rU") as input_handle:
		alignments = AlignIO.read(input_handle, PHYLIP_FORMAT)

	with open(output_file, "w") as output_handle:
		out_handle = AlignIO.PhylipIO.SequentialPhylipWriter(output_handle)
		out_handle.write_alignment(alignments, id_width=30)


def get_avg(l):
	avg = sum(l)/len(l)
	return avg


def get_var(l):
	avg = get_avg(l)
	dist_from_avg = list(map(lambda x: (x - avg)**2, l))
	var = sum(dist_from_avg)/len(dist_from_avg)
	return var


def get_std(l):
	var = get_var(l)
	std = math.sqrt(var)
	return std


def sum_of_squares(l):
	return sum([c**2 for c in l])


def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False


def format_float(num):
	n = str(num)
	if "e-" in n:
		k = n.split("e-")[-1]
		n = ("{0:." + k + "f}").format(num)
	return n


def is_file_empty(filepath):
	"""
	:param filepath:
	:return: True if filepath doesn't exist or is empty
	"""
	if os.path.exists(filepath):
		with open(filepath) as fpr:
			if re.search("\S", fpr.read()):
				return False
	return True


def compute_entropy(lst, epsilon=0.000001):
	if np.sum(lst) != 0:
		lst_norm = np.array(lst)/np.sum(lst)
	else:
		lst_norm = np.array(lst) + epsilon
	entropy = -1*sum(np.log2(lst_norm)*lst_norm)
	if np.isnan(entropy):
		lst_norm += epsilon
		entropy = -1*sum(np.log2(lst_norm)*lst_norm)
	return entropy


def lists_diff(li1, li2):
	"""
	:param li1, li2: two lists
	:return: diff between lists
	"""
	return (list(set(li1) - set(li2)))