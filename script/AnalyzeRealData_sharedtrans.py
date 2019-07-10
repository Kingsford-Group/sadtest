#!/bin/python

# analyze real data SAD novel isoform prediction
# questions to answer:
# 	- the shared transcripts of unadjustable anomalies, what is the transcript lengths? compared to the whole gene transcript distribution?
# 	- the under-expression region (or deleted region in novel isoform), how many exons does it span? Or hown many exons are left out in the novel one?
# 	- what is the start / end proportion of over-expressed region in transcript?
# 	- what is the exon length distribution of the omitting-exon of potential novel isoform?

import sys
import copy
import gzip
import struct
import pickle
import numpy as np
from TranscriptClass import *
from pathlib import Path


def ReadSignificantTranscripts(filename, p_threshold = 0.01):
	trans_list = []
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if float(strs[12]) < p_threshold:
			trans_list.append(strs[0])
	fp.close()
	return trans_list


def GetSharedTranscripts(InFolders, pvaluefilename):
	shared_trans_set = set()
	counter = 0
	for folder in InFolders:
		p = Path(folder)
		files = [str(f) for f in p.glob("*/"+pvaluefilename)]
		for f in files:
			trans_list = ReadSignificantTranscripts(f)
			if counter == 0:
				shared_trans_set = set(trans_list)
			else:
				shared_trans_set = shared_trans_set & set(trans_list)
			counter += 1
	print("There are {} shared transcripts under the unadjustable anomaly category.".format(len(shared_trans_set)))
	return list(shared_trans_set)


def GetUnionTranscripts(InFolders, pvaluefilename):
	union_trans_set = set()
	for folder in InFolders:
		p = Path(folder)
		files = [str(f) for f in p.glob("*/"+pvaluefilename)]
		for f in files:
			trans_list = ReadSignificantTranscripts(f)
			union_trans_set = union_trans_set | set(trans_list)
	print("There are {} union transcripts under the unadjustable anomaly category.".format(len(union_trans_set)))
	return list(union_trans_set)


def CheckNumSpanningExons(t, start, end):
	assert(isinstance(t, Transcript_t))
	exoncount = 0
	coverlen = 0
	for e in t.Exons:
		if max(coverlen, start) < min(coverlen+e[1]-e[0], end):
			exoncount += 1
		coverlen += e[1] - e[0]
	return exoncount


def CheckNumSpanningExons_Length(t, start, end):
	assert(isinstance(t, Transcript_t))
	lengths = []
	coverlen = 0
	for e in t.Exons:
		if max(coverlen, start) < min(coverlen+e[1]-e[0], end):
			lengths.append( e[1]-e[0] )
		coverlen += e[1] - e[0]
	assert(len(lengths) >= 1)
	return lengths


def CheckNumContainingExons_Length(t, start, end):
	assert(isinstance(t, Transcript_t))
	lengths = []
	coverlen = 0
	for e in t.Exons:
		if start <= coverlen and end >= coverlen+e[1]-e[0]:
			lengths.append( e[1]-e[0] )
		coverlen += e[1] - e[0]
	return lengths


def RetrieveAffectIsoform(pvaluefile, SharedTransList):
	SharedTransSet = set(SharedTransList)
	IsoformInfo = {}
	fp = open(pvaluefile, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		if strs[0] in SharedTransSet:
			# if not (int(strs[3]) >= 0 and int(strs[4]) >= 0 and int(strs[6]) >= 0 and int(strs[7]) >= 0):
			if not (int(strs[3]) >= 0 and int(strs[4]) >= 0 and int(strs[8]) >= 0 and int(strs[9]) >= 0):
				print(pvaluefile)
				print(line)
			# assert(int(strs[3]) >= 0 and int(strs[4]) >= 0 and int(strs[6]) >= 0 and int(strs[7]) >= 0)
			assert(int(strs[3]) >= 0 and int(strs[4]) >= 0 and int(strs[8]) >= 0 and int(strs[9]) >= 0)
			# IsoformInfo[strs[-1]] = [strs[0], (int(strs[3]), int(strs[4])), (int(strs[6]), int(strs[7]))] # transcript name, under-expression region, over-expression region
			IsoformInfo[strs[0]] = [(int(strs[3]), int(strs[4])), (int(strs[8]), int(strs[9]))] # under-expression region, over-expression region
	fp.close()
	assert( len(IsoformInfo) == len(SharedTransList) )
	return IsoformInfo


def WriteTransLengthTable(outputfile, TransLength, SharedTransList, UnionTransList):
	fp = open(outputfile, 'w')
	fp.write("# TransName\tTransLength\tType\n")
	SharedTransSet = set(SharedTransList)
	UnionTransSet = set(UnionTransList)
	for t,v in TransLength.items():
		if t in SharedTransSet:
			fp.write("{}\t{}\t{}\n".format(t, v, "SharedSig"))
		elif t in UnionTransSet:
			fp.write("{}\t{}\t{}\n".format(t, v, "NonSharedSig"))
		else:
			fp.write("{}\t{}\t{}\n".format(t, v, "NotSig"))
	fp.close()


def WriteNumSpanningExons(outputfile, InFolders, pvaluefilename, Transcripts, SharedTransList):
	fp = open(outputfile, 'w')
	fp.write("# TransName\tDataSet\tSampleID\tNumSpanningExons\n")
	for folder in InFolders:
		p = Path(folder)
		files = [str(f) for f in p.glob("*/"+pvaluefilename)]
		dataset = files[0].split("/")[-3]
		print("under-expression region spanning exons: "+dataset)
		for f in files:
			print("\t"+f)
			sampleID = f.split("/")[-2].split("_")[-1]
			IsoformInfo = RetrieveAffectIsoform(f, SharedTransList)
			# count spanning exon for under-expression region
			ExonCount = {}
			for t,v in IsoformInfo.items():
				ExonCount[t] = CheckNumSpanningExons(Transcripts[t], v[0][0], v[0][1])
			# write to file
			for t in SharedTransList:
				fp.write("{}\t{}\t{}\t{}\n".format(t, dataset, sampleID, ExonCount[t]))
	fp.close()


def WriteOverRegionProportion(outputfile, InFolders, pvaluefilename, TransLength, SharedTransList):
	fp = open(outputfile, 'w')
	fp.write("# TransName\tDataSet\tSampleID\tStartProportion\tEndProportion\n")
	for folder in InFolders:
		p = Path(folder)
		files = [str(f) for f in p.glob("*/"+pvaluefilename)]
		dataset = files[0].split("/")[-3]
		print("Over-expression region proportion: "+dataset)
		for f in files:
			print("\t"+f)
			sampleID = f.split("/")[-2].split("_")[-1]
			IsoformInfo = RetrieveAffectIsoform(f, SharedTransList)
			# calculate the proportion and write to file
			for t in SharedTransList:
				overregion = IsoformInfo[t][1]
				start_prop = overregion[0] / TransLength[t]
				end_prop = overregion[1] / TransLength[t]
				fp.write("{}\t{}\t{}\t{}\t{}\n".format(t, dataset, sampleID, start_prop, end_prop))
	fp.close()


def WriteUnderRegionProportion(outputfile, InFolders, pvaluefilename, TransLength, SharedTransList):
	fp = open(outputfile, 'w')
	fp.write("# TransName\tDataSet\tSampleID\tStartProportion\tEndProportion\n")
	for folder in InFolders:
		p = Path(folder)
		files = [str(f) for f in p.glob("*/"+pvaluefilename)]
		dataset = files[0].split("/")[-3]
		print("Under-expression region proportion: "+dataset)
		for f in files:
			print("\t"+f)
			sampleID = f.split("/")[-2].split("_")[-1]
			IsoformInfo = RetrieveAffectIsoform(f, SharedTransList)
			# calculate the proportion and write to file
			for t in SharedTransList:
				underregion = IsoformInfo[t][0]
				start_prop = underregion[0] / TransLength[t]
				end_prop = underregion[1] / TransLength[t]
				fp.write("{}\t{}\t{}\t{}\t{}\n".format(t, dataset, sampleID, start_prop, end_prop))
	fp.close()


def WriteLengthSpanningExons(outputfile, InFolders, pvaluefilename, Transcripts, GeneTransMap, SharedTransList):
	fp = open(outputfile, 'w')
	fp.write("# TransName\tExonLength\tLabel\n")
	TransExonLength_omit = {}
	TransExonLength_contain = {}
	for folder in InFolders:
		p = Path(folder)
		files = [str(f) for f in p.glob("*/"+pvaluefilename)]
		dataset = files[0].split("/")[-3]
		print("under-expression region exon lengths: "+dataset)
		for f in files:
			print("\t"+f)
			sampleID = f.split("/")[-2].split("_")[-1]
			IsoformInfo = RetrieveAffectIsoform(f, SharedTransList)
			# count spanning exon for under-expression region
			for t,v in IsoformInfo.items():
				# save the omitted exon length
				lengths = CheckNumSpanningExons_Length(Transcripts[t], v[0][0], v[0][1])
				if t in TransExonLength_omit:
					if len(lengths) == 1:
						TransExonLength_omit[t] += [(elength, "Single") for elength in lengths]
					else:
						TransExonLength_omit[t] += [(elength, "Multiple") for elength in lengths]
				else:
					if len(lengths) == 1:
						TransExonLength_omit[t] = [(elength, "Single") for elength in lengths]
					else:
						TransExonLength_omit[t] = [(elength, "Multiple") for elength in lengths]
				# save the containing exon length
				lengths = CheckNumContainingExons_Length(Transcripts[t], v[1][0], v[1][1])
				if len(lengths) > 0:
					if t in TransExonLength_contain:
						TransExonLength_contain[t] += [(elength, "Contain") for elength in lengths]
					else:
						TransExonLength_contain[t] = [(elength, "Contain") for elength in lengths]
	# write novel isoform omitting exon lengths to file
	for t,v in TransExonLength_omit.items():
		TransExonLength_omit[t] = list(set(v))
	for t,v in TransExonLength_omit.items():
		for elength in v:
			fp.write("{}\t{}\t{}\n".format(t, elength[0], elength[1]))
	# write novel isoform containing exon lengths to file
	for t,v in TransExonLength_contain.items():
		TransExonLength_contain[t] = list(set(v))
	for t,v in TransExonLength_contain.items():
		for elength in v:
			fp.write("{}\t{}\t{}\n".format(t, elength[0], elength[1]))
	# compare to all the exons
	for g,v in GeneTransMap.items():
		alllengths = sum([Transcripts[t].Exons for t in v], [])
		alllengths = list(set(alllengths))
		alllengths = [e[1] - e[0] for e in alllengths]
		for elength in alllengths:
			fp.write("{}\t{}\tAll\n".format(g, elength))
	fp.close()


if __name__=="__main__":
	Transcripts = ReadGTF("/home/cong/Documents/SADnewresult/gencode.v26.annotation.gtf")
	[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)
	InFolders = ["/home/cong/Documents/SADnewresult/GEUVADIS/", "/home/cong/Documents/SADnewresult/HumanBodyMap/"]

	pvaluefilename = "test_pvalue_overall_sorted"
	SharedTransList = GetSharedTranscripts(InFolders, pvaluefilename)
	UnionTransList = GetUnionTranscripts(InFolders, pvaluefilename)

	# shared gene length vs all gene length distribution
	if not Path("/home/cong/Documents/SADnewresult/SharedTrans_Length.txt").exists():
		TransLength = GetTransLength(Transcripts)
		WriteTransLengthTable("/home/cong/Documents/SADnewresult/SharedTrans_Length.txt", TransLength, SharedTransList, UnionTransList)

	# number of spanning exons of the under-expression region for each corresponding isoform for each sample
	if not Path("/home/cong/Documents/SADnewresult/SharedTrans_NumSpanningExons.txt").exists():
		# pvaluefilename = "test_pvalue_overall_sorted_uniqgene"
		pvaluefilename = "test_pvalue_overall_sorted"
		WriteNumSpanningExons("/home/cong/Documents/SADnewresult/SharedTrans_NumSpanningExons.txt", InFolders, pvaluefilename, Transcripts, SharedTransList)

	# the start / end proportion in transcript length
	if not Path("/home/cong/Documents/SADnewresult/SharedTrans_OverExpProportion.txt").exists():
		TransLength = GetTransLength(Transcripts)
		# pvaluefilename = "test_pvalue_overall_sorted_uniqgene"
		pvaluefilename = "test_pvalue_overall_sorted"
		WriteOverRegionProportion("/home/cong/Documents/SADnewresult/SharedTrans_OverExpProportion.txt", InFolders, pvaluefilename, TransLength, SharedTransList)
	
	if not Path("/home/cong/Documents/SADnewresult/SharedTrans_UnderExpProportion.txt").exists():
		TransLength = GetTransLength(Transcripts)
		# pvaluefilename = "test_pvalue_overall_sorted_uniqgene"
		pvaluefilename = "test_pvalue_overall_sorted"
		WriteUnderRegionProportion("/home/cong/Documents/SADnewresult/SharedTrans_UnderExpProportion.txt", InFolders, pvaluefilename, TransLength, SharedTransList)

	# what is the exon length distribution of the omitting-exon of novel isoform
	if not Path("/home/cong/Documents/SADnewresult/SharedTrans_ExonLengths.txt").exists():
		# pvaluefilename = "test_pvalue_overall_sorted_uniqgene"
		pvaluefilename = "test_pvalue_overall_sorted"
		WriteLengthSpanningExons("/home/cong/Documents/SADnewresult/SharedTrans_ExonLengths.txt", InFolders, pvaluefilename, Transcripts, GeneTransMap, SharedTransList)
