import copy
from Bio.Align.AlignInfo import SummaryInfo
from utils import *
from Bio import AlignIO, Seq
from collections import Counter
__author__ = 'Shiran'


def make_msa_generic(msa, dest_filename=None):
	"""
	:param msa: an AlignIO object
	:param back_transcribe:
	:param dest_filename:
	:return:
	"""
	if dest_filename:
		with open(dest_filename, "w") as fpw:
			AlignIO.write(msa, fpw, format=PHYLIP_FORMAT)
	# return msa


def get_chars_content(sequence, sel_chars, all_chars):
	"""returns the fraction of sel_chars within all_chars in sequence"""
	#chars = a string of characters concatenated within square parenthesis to look for in sequence
	count = len(re.findall(sel_chars, sequence, re.I))
	total_base_count = len(re.findall(all_chars, sequence, re.I))
	fraction = count/total_base_count
	return fraction


def get_chars_statistic_per_alignment(msa_seqs, sel_chars, all_chars):
	seqs_chars_content = []

	for rec in msa_seqs:
		seq = rec._seq._data
		chars_content = get_chars_content(seq, sel_chars, all_chars)
		seqs_chars_content.append(chars_content)

	avg = get_avg(seqs_chars_content)
	std = get_std(seqs_chars_content)

	return avg, std


def remove_gaps_from_sequence(seq):
	return re.sub("[^agctAGCT]+", "", seq, re.I)


def is_column_invariant(col_string, thres=100):
	# return true if column is invariant above a threshold (percentage)
	return get_column_var_percent(col_string) >= thres


def get_column_var_percent(col_string):
	col_gapless = remove_gaps_from_sequence(col_string)
	#if something that is not gaps but not AGCT appears in the column
	if len(col_gapless) == 0:
		return 0
	freqs = (col_gapless.upper().count(x) for x in "ACGT")
	return round((max(freqs)/len(col_gapless)) * 100)


def get_pInv_above_thres(msa, thres=100):
	"""
	:param msa:
	:param thresholds: a list of percentages
	:return:  a list - for every percentage, how many sites above this conservation thresholds
	"""
	msa_length = msa.get_alignment_length()
	invariant_sites = 0
	for col_i in range(0, msa_length):
		invariant_sites += is_column_invariant(msa[:,col_i], thres)

	return invariant_sites/float(msa_length)


def count_fully_conserved_fraction(msa):
	"""
	:param msa:
	:param thresholds: a list of percentages
	:return:  a list - for every percentage, how many sites above this conservation thresholds
	"""
	msa_length = msa.get_alignment_length()
	invariant_sites = 0
	for col_i in range(0, msa_length):
		col_gapless = remove_gaps_from_sequence(msa[:,col_i])
		# if something that is not gaps but not AGCT appears in the column
		if len(col_gapless) == 0:
			continue
		if col_gapless.count(col_gapless[0]) == len(col_gapless):
			invariant_sites += 1

	return invariant_sites/msa_length


def calculate_column_entropy(col_string):
	# column_entropy = - sum(for every nucleotide x) {count(x)*log2(Prob(nuc x in col i))}
	col_gapless = remove_gaps_from_sequence(col_string).upper()
	col_entropy = 0
	for x in ['A', 'G', 'C', 'T']:
		count_x = str.count(col_gapless, x)
		if count_x == 0:
			entropy_x = 0
		else:
			prob_x = count_x/len(col_gapless)
			entropy_x = count_x*math.log2(prob_x)
		col_entropy += entropy_x

	return -col_entropy


def get_msa_avg_entropy(msa):
	msa_length = msa.get_alignment_length()
	sum_entropy = 0
	for col_i in range(0, msa_length):
		sum_entropy += calculate_column_entropy(msa[:,col_i])

	return sum_entropy/msa_length


def count_distinct_seqs(msa, only_middle=False):
	start_idx = 0
	end_idx = msa.get_alignment_length()
	if only_middle:
		start_idx = end_idx//3
		end_idx = start_idx*2
	seqs_list = []
	for i in range(len(msa)):
		rec = msa[i]
		seqs_list.append(remove_gaps_from_sequence(str(rec._seq[start_idx:end_idx])))
	return len(set(seqs_list))


def get_columns_gc_variance(msa):
	msa_length = msa.get_alignment_length()
	gc_cnts = []
	for col_i in range(0, msa_length):
		gc_cnts.append(get_chars_content(msa[:,col_i],"[GC]", "[ATGC]"))

	return get_avg(gc_cnts),get_var(gc_cnts)


def calculate_bollback_multinomial(msa):
	msa_length = msa.get_alignment_length()
	cols = []
	for col_i in range(0, msa_length):
		cols.append(msa[:,col_i])
	counts = Counter(cols)

	multinomial = 0
	for k in counts:
		c = counts[k]
		multinomial += c*math.log(c)
	multinomial -= msa_length*math.log(msa_length)
	return multinomial, len(counts), len(counts)/msa_length


def count_patterns(msa):
	msa_length = msa.get_alignment_length()
	col_patterns = []
	for col_i in range(0, msa_length):
		x = 0
		col = msa[:,col_i]
		uniques = set(list(col))
		for item in col:
			if item in uniques:
				col = re.sub(item, str(x), col)
				x += 1
				uniques -= set(item)
			if len(uniques) == 0:
				break
		col_patterns.append(col)
	n_unique = len(set(col_patterns))
	return n_unique, n_unique/msa_length


def infer_pairwise_substitution_matrix(seq1, seq2):
	substitution_count_dictionary = {"AC": 0, "AG": 0, "AT": 0, "CG": 0, "CT": 0, "GT": 0, "AA": 0, "GG": 0, "CC": 0,
	                                 "TT": 0, "1s": 0, "2s": 0} #1s for one space vs nucleotide, 2s for 2 spaces
	pa_length = 0
	for i in range(0, len(seq1)):
		ch1 = min(seq1[i].upper(), seq2[i].upper())
		ch2 = max(seq1[i].upper(), seq2[i].upper())
		if ch1 not in ["A", "G", "C", "T"] and ch2 not in ["A", "G", "C", "T"]:
			substitution_count_dictionary["2s"] +=1
		elif ch1 in ["A", "G", "C", "T"] and ch2 in ["A", "G", "C", "T"]: #both nucleotides
			substitution_count_dictionary[ch1+ch2] += 1
			pa_length +=1
		else: #1 space
			substitution_count_dictionary["1s"] += 1

	return substitution_count_dictionary, pa_length


def compute_pairwise_substitution_rates(seq1, seq2):
	MATCH_SCORE = 1
	MISMATCH_SCORE = -1
	GAP_SCORE = -1

	transition_rate, transversion_rate, match, mismatch, gap = 0, 0, 0, 0, 0

	substitution_count_dictionary, pa_length = infer_pairwise_substitution_matrix(seq1, seq2)
	if pa_length != 0:
		transition_rate = float(substitution_count_dictionary["AG"] + substitution_count_dictionary["CT"]) / pa_length
		transversion_rate = float(substitution_count_dictionary["AC"] + substitution_count_dictionary["AT"] +
								  substitution_count_dictionary["CG"] + substitution_count_dictionary["GT"]) / pa_length
		match = (substitution_count_dictionary["AA"]+substitution_count_dictionary["CC"]+substitution_count_dictionary["GG"]+substitution_count_dictionary["TT"])*MATCH_SCORE
		mismatch = (substitution_count_dictionary["AC"]+substitution_count_dictionary["AG"]+substitution_count_dictionary["AT"] +
				   substitution_count_dictionary["CG"]+substitution_count_dictionary["CT"]+substitution_count_dictionary["GT"]+substitution_count_dictionary["1s"])*MISMATCH_SCORE
	gap = substitution_count_dictionary["1s"]*GAP_SCORE

	return transition_rate, transversion_rate, match+mismatch+gap,\
			substitution_count_dictionary["AC"], substitution_count_dictionary["AG"], substitution_count_dictionary["AT"], \
			substitution_count_dictionary["CG"], substitution_count_dictionary["CT"], \
			substitution_count_dictionary["GT"]


def compute_base_frequencies(msa):
	freqs = []
	seqs_msa = list(msa)
	allchars = "".join([rec._seq._data for rec in seqs_msa])

	for nuc in "ACGT":
		freqs.append(len(re.findall(nuc, allchars, re.I)))
	freqs = {"freq_" + nuc : freq / sum(freqs) for nuc, freq in zip("ACGT", freqs)}
	return freqs


def calculate_substitution_rates(msa):
	seqs_msa = list(msa)
	transition_rates = []
	transversion_rates = []
	sop_score = 0
	ac_cnt = ag_cnt = at_cnt = cg_cnt = ct_cnt = gt_cnt = 0
	for i in range(0, len(seqs_msa)-1):
		for j in range(i+1,len(seqs_msa)):
			seq1 = seqs_msa[i]._seq._data
			seq2 = seqs_msa[j]._seq._data
			transition_rate, transversion_rate, pair_sop, ac, ag, at, cg, ct, gt = compute_pairwise_substitution_rates(seq1, seq2)
			transition_rates.append(transition_rate)
			transversion_rates.append(transversion_rate)
			sop_score += pair_sop
			ac_cnt, ag_cnt, at_cnt, cg_cnt, ct_cnt, gt_cnt = \
				ac_cnt + ac, ag_cnt + ag, at_cnt + at, cg_cnt + cg, ct_cnt + ct, gt_cnt + gt

	subs_sum = sum([ac_cnt, ag_cnt, at_cnt, cg_cnt, ct_cnt, gt_cnt])
	return {"transition_avg": get_avg(transition_rates),
			"transversion_avg": get_avg(transversion_rates),
			"sop_score": sop_score},\
		   {c: x/subs_sum if subs_sum != 0 else 0 for c,x in
			zip(["ac_subs", "ag_subs", "at_subs", 'cg_subs', 'ct_subs', 'gt_subs'],
				[ac_cnt, ag_cnt, at_cnt, cg_cnt, ct_cnt, gt_cnt])}


def get_msa_from_file(msa_file_path):
	#open file if exists
	if not os.path.exists(msa_file_path):
		return None
	try:
		msa = AlignIO.read(msa_file_path, PHYLIP_FORMAT)
	except:
		return None
	return msa


def get_msa_properties(msa):
	"""
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	if isinstance(msa, str):
		msa = get_msa_from_file(msa)
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()

	return ntaxa, nchars


def remove_nonACTGU_sites(msa, dest_filename=None):
	"""
	removes sites that contain non ACGT characters (i.e., N's, gaps etc.)
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	if isinstance(msa, str):
		msa = get_msa_from_file(msa)
	nchars = msa.get_alignment_length()
	new_msa = copy.deepcopy(msa)
	for col_i in range(nchars - 1, -1, -1):
		col = msa[:, col_i]
		if not re.search("[ACGTUacgtu]", col, re.IGNORECASE):
			new_msa = new_msa[:, :col_i] + new_msa[:, col_i + 1:]
	if dest_filename:
		with open(dest_filename, "w") as fpw:
			AlignIO.write(new_msa, fpw, format=PHYLIP_FORMAT)
	return new_msa


def remove_gapped_sites(msa, thres):
	"""
	removes sites that contain more than thres gaps in a column
	:param msa: a Bio.Align.MultipleSeqAlignment object
	:param thres: a floating number in [0,1], the maximal fraction of
	sites that may hold gaps within every alignment column
	:return: a new msa without sites that have more than thres gaps
	"""
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()
	new_msa = copy.deepcopy(msa)
	for col_i in range(nchars - 1, -1, -1):
		if msa[:, col_i].count("-") / ntaxa > thres:
			new_msa = new_msa[:, :col_i] + new_msa[:, col_i + 1:]
	return new_msa


def validate_msa_characters(msa, legal_chars="-acgtn"):
	msa_info = SummaryInfo(msa)
	msa_letters = msa_info._get_all_letters()
	for let in msa_letters:
		if let not in legal_chars:
			return False
	return True


def count_substitutions(seq1, seq2, is_codons):
	"""
	:param seq1:
	:param seq2:
	:param is_codons: if True, count every triplet substitution as a single one, else - every individual position
	:return: average number of substitutions
	"""
	m = len(seq1)
	jump = 3 if is_codons else 1
	assert (len(seq1) == len(seq2)) and (m%jump==0)
	cnt = 0
	for i in range(0, m, jump):
		cnt += int(seq1[i: i+jump] != seq2[i: i+jump])
	return cnt/(m/jump)


def reduce_msa_to_seqs_by_name(msa, keep_names_lst):
	new_msa = []
	all_names = [rec.id for rec in list(msa)]
	for name in keep_names_lst:
		new_msa.append(msa[all_names.index(name), :])
	#remove positions that are just gaps after removal of sequences
	new_msa = remove_nonACTGU_sites(AlignIO.MultipleSeqAlignment(new_msa))
	return new_msa

if __name__ == '__main__':
	msa_file = "/groups/itay_mayrose/shiranabad/MsMl/data/protDBs/balibase/RV30_BBS30021/real_msa.phy"
	msa = get_msa_from_file(msa_file)
	count_patterns(msa)