#!/bin/python

import sys
import numpy as np
import copy
import pickle
from tqdm import tqdm
from TranscriptClass import *


def ReadSAD_region(filename, p_threshold = 0.01):
	predictions = {}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if float(strs[12]) >= p_threshold:
			continue
		predictions[strs[0]] = [(int(strs[3]), int(strs[4])), (int(strs[8]), int(strs[9]))] # under-expression region, over-expression region
	fp.close()
	return predictions


def ReadLongReadNucmer(filename, dist_threshold = 20000):
	# treat each alignment as a transcript
	longread_trans = []
	tmp_trans = Transcript_t("", "", "", "", True, -1, -1)
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount < 6:
			continue
		strs = [x for x in line.strip().split() if x != ""]
		# if a read has changed, starting make a new transcript
		# this_strand = (int(strs[3]) < int(strs[4]))
		# if strs[-1] != tmp_trans.TransID or tmp_trans.Chr != strs[-2] or tmp_trans.Strand != this_strand or np.abs(tmp_trans.Exons[-1][1] - int(strs[0])) > dist_threshold:
		if strs[-1] != tmp_trans.TransID or tmp_trans.Chr != strs[-2] or np.abs(tmp_trans.Exons[-1][1] - int(strs[0])) > dist_threshold:
			if tmp_trans.TransID != "":
				tmp_trans.EndPos = np.max([e[1] for e in tmp_trans.Exons])
				longread_trans.append( copy.deepcopy(tmp_trans) )
			# tmp_trans = Transcript_t(strs[-1], "", "", strs[-2], this_strand, int(strs[0]), -1)
			tmp_trans = Transcript_t(strs[-1], "", "", strs[-2], True, int(strs[0]), -1)
		# adding the aligned regions as exons to the transcript
		tmp_trans.Exons.append( (int(strs[0]), int(strs[1])) )
	# add the last transcript group
	if tmp_trans.TransID != "":
		tmp_trans.EndPos = np.max([e[1] for e in tmp_trans.Exons])
		longread_trans.append( copy.deepcopy(tmp_trans) )
	fp.close()
	# make each long read transcript have non-overlapping exons
	for i in range(len(longread_trans)):
		sorted_exons = sorted(longread_trans[i].Exons)
		nonoverlap = []
		for e in sorted_exons:
			if len(nonoverlap) != 0 and e[0] < nonoverlap[-1][1]:
				assert( nonoverlap[-1][0] <= e[0] )
				nonoverlap[-1][1] = max(nonoverlap[-1][1], e[1])
			else:
				nonoverlap.append( [e[0], e[1]] )
		nonoverlap = [(e[0], e[1]) for e in nonoverlap]
		# if longread_trans[i].Strand:
		# 	longread_trans[i].Exons = nonoverlap
		# else:
		# 	longread_trans[i].Exons = nonoverlap[::-1]
		longread_trans[i].Exons = nonoverlap
	return longread_trans


def RemoveExactMatch(longread_trans, transcripts, fuzzy = 20):
	# note that longread_trans is a list, but transcripts is a map
	# definition of exact match: intron chain match + sboth start and ending positions are within 200bp range of the matched reference transcript
	new_longread_trans = []
	with tqdm(total=len(longread_trans)) as pbar:
		for tl in longread_trans:
			is_exact_match = False
			related_ref_trans = [v for t,v in transcripts.items() if v.Chr == tl.Chr and max(v.StartPos, tl.StartPos) < min(v.EndPos, tl.EndPos)]
			for tr in related_ref_trans:
				if len(tl.Exons) != len(tr.Exons):
					continue
				# check intron chain
				sorted_tl = sorted(tl.Exons)
				sorted_tr = sorted(tr.Exons)
				has_violation = False
				for i in range(len(sorted_tr)):
					if i == 0:
						if np.abs(sorted_tl[i][1] - sorted_tr[i][1]) > fuzzy:
							has_violation = True
							break
					elif i + 1 == len(sorted_tr):
						if np.abs(sorted_tl[i][0] - sorted_tr[i][0]) > fuzzy:
							has_violation = True
							break
					else:
						if np.abs(sorted_tl[i][0] - sorted_tr[i][0]) > fuzzy or np.abs(sorted_tl[i][1] - sorted_tr[i][1]):
							has_violation = True
							break
				# check start / ending position
				if not has_violation:
					if np.abs(tl.StartPos - tr.StartPos) < 200 and np.abs(tl.EndPos - tr.EndPos) < 200:
						is_exact_match = True
						break
			if tl.TransID == "PB2015.919.2|chr1:145096393-145117731(+)|i2a_c42094/f2p200/2505":
				print(is_exact_match)
			if not is_exact_match:
				new_longread_trans.append( tl )
			pbar.update(1)
	return new_longread_trans


def ValidateSAD(predictions, transcripts, longread_trans, pro_threshold = 0.75):
	# definition of being validated: cover > 75% over-expression region, exclude > 75% under-expression region
	status = {tname:False for tname in predictions.keys()}
	correspondence = {tname:"." for tname in predictions.keys()}
	with tqdm(total=len(predictions)) as pbar:
		for tname, p in predictions.items():
			t = transcripts[tname]
			# convert the over and under-expression from transcript coordinate to genome coordination
			under = []
			over = []
			l = 0
			for i in range(len(t.Exons)):
				# whether overlap with under-expression region
				if max(l, p[0][0]) < min(l + t.Exons[i][1] - t.Exons[i][0], p[0][1]):
					if t.Strand:
						under.append( (t.Exons[i][0] + max(l, p[0][0]) - l, t.Exons[i][0] + min(l + t.Exons[i][1] - t.Exons[i][0], p[0][1]) - l) )
					else:
						under.append( (t.Exons[i][1] - min(l + t.Exons[i][1] - t.Exons[i][0], p[0][1]) + l, t.Exons[i][1] - max(l, p[0][0]) + l) )
				# whether overlap with over-expression region
				if max(l, p[1][0]) < min(l + t.Exons[i][1] - t.Exons[i][0], p[1][1]):
					if t.Strand:
						over.append( (t.Exons[i][0] + max(l, p[1][0]) - l, t.Exons[i][0] + min(l + t.Exons[i][1] - t.Exons[i][0], p[1][1]) - l) )
					else:
						over.append( (t.Exons[i][1] - min(l + t.Exons[i][1] - t.Exons[i][0], p[1][1]) + l, t.Exons[i][1] - max(l, p[1][0]) + l) )
				l += t.Exons[i][1] - t.Exons[i][0]
			# compare the over and under-expression region with the long read transcripts
			for tl in longread_trans:
				if tl.Chr != t.Chr or max(tl.StartPos, t.StartPos) >= min(tl.EndPos, t.EndPos):
					continue
				# calculate overlapping length with the over and under-expression region
				overlap_over = 0
				overlap_under = 0
				for e in tl.Exons:
					for e1 in over:
						if max(e[0], e1[0]) < min(e[1], e1[1]):
							# print([e[0], e[1], e1[0], e1[1], min(e[1], e1[1]) - max(e[0], e1[0])])
							overlap_over += min(e[1], e1[1]) - max(e[0], e1[0])
					for e2 in under:
						if max(e[0], e2[0]) < min(e[1], e2[1]):
							overlap_under += min(e[1], e2[1]) - max(e[0], e2[0])
				if overlap_over > pro_threshold * np.sum([e[1] - e[0] for e in over]) and overlap_under < (1 - pro_threshold) * np.sum([e[1] - e[0] for e in under]):
					status[tname] = True
					correspondence[tname] = tl.TransID
					break
			pbar.set_postfix(status = status[tname])
			pbar.update(1)
	return status, correspondence


def WriteStatus(output_file, status, correspondence):
	fp = open(output_file, 'w')
	fp.write("# Name\tStatus\tCorresponding_long_read\n")
	for k,v in status.items():
		fp.write("{}\t{}\t{}\n".format(k, v, correspondence[k]))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python LongreadsNucmer.py <gtf_file> <pvalue_file> <nucmer_file> <output_file>")
	else:
		gtf_file = sys.argv[1]
		pvalue_file = sys.argv[2]
		nucmer_file = sys.argv[3]
		output_file = sys.argv[4]

		transcripts = ReadGTF(gtf_file)
		predictions = ReadSAD_region(pvalue_file)
		longread_trans = ReadLongReadNucmer(nucmer_file)

		longread_trans = RemoveExactMatch(longread_trans, transcripts)
		print(".".join(nucmer_file.split(".")[:-1]) + "_noexactmatch.pickle")
		pickle.dump( longread_trans, open(".".join(nucmer_file.split(".")[:-1]) + "_noexactmatch.pickle", 'wb') )
		status, correspondence = ValidateSAD(predictions, transcripts, longread_trans)
		WriteStatus(output_file, status, correspondence)
