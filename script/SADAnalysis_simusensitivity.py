#!/bin/python

import sys
import numpy as np
from pathlib import Path


def CountDetection(SimuDir, SADFolderName, PredFile, EvaFile, _type, _numevents, _id):
	if (not Path(SimuDir+"/salmon/"+_type+"_"+_id+"_"+str(_numevents)+"/"+SADFolderName+"/"+PredFile).exists()) or (not Path(SimuDir+"/salmon/"+_type+"_"+_id+"_"+str(_numevents)+"/"+SADFolderName+"/"+EvaFile).exists()):
		return -1
	# read prediction file, count the number of significant predictions
	fp = open(SimuDir+"/salmon/"+_type+"_"+_id+"_"+str(_numevents)+"/"+SADFolderName+"/"+PredFile, 'r')
	countsignificant = 0
	ind_adj = -1
	for line in fp:
		strs = line.strip().split("\t")
		if line[0] == '#':
			#ind_adj = strs.index("AdjustedPvalue")
			ind_adj = strs.index("MinAdjPValue")
			assert(ind_adj > -1)
		else:
			if float(strs[ind_adj]) < 0.01:
				countsignificant += 1
	fp.close()
	# read evaluation file, count unique events corresponding to the significant predictions
	fp = open(SimuDir+"/salmon/"+_type+"_"+_id+"_"+str(_numevents)+"/"+SADFolderName+"/"+EvaFile, 'r')
	linecount = 0
	Events = []
	for line in fp:
		linecount += 1
		if linecount < countsignificant + 1:
			strs = line.strip().split("\t")
			if strs[4] != ".":
				Events += strs[4].split(",")
	fp.close()
	Events = set(Events)
	return len(Events)


def CountSimulation(SimuDir, Suffix, _type, _numevents, _id):
	# count fusion
	fp = open(SimuDir+"/groundtruth_"+_type+"_"+_id+"_"+str(_numevents)+"/fusion_"+Suffix+".txt", 'r')
	lines = fp.readlines()
	numfusion = len(lines)
	# count remove
	fp = open(SimuDir+"/groundtruth_"+_type+"_"+_id+"_"+str(_numevents)+"/removelist_"+Suffix+".txt", 'r')
	lines = fp.readlines()
	numremove = len(lines)
	return (numfusion+numremove)


def WriteSensitivity(OutputFile, Result):
	fp = open(OutputFile, 'w')
	fp.write("# Type\tNumEvents\tID\tNumSimulated\tNumDetected\tSensitivity\n")
	for re in Result:
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(re[0], re[1], re[2], re[3], re[4], 1.0*re[4]/re[3]))
	fp.close()


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python3 AnalyzeRealData_sharedtrans.py <prepare data directory>")
	else:
		prepdir = sys.argv[1]

		Type = ["PC", "Full"]
		NumEvents = [200, 500, 1000, 1500]
		RefQuantID = ["GEU", "GM12878", "K562"]
		SimuDir = prepdir + "/Simulation"
		SADFolderName = "sad"
		OutputFile = prepdir + "/Simulation/Sensitivity.txt"

		Result = []
		for _type in Type:
			for _numevents in NumEvents:
				for _id in RefQuantID:
					numsimulated = CountSimulation(SimuDir, "existjunc", _type, _numevents, _id)
					numdetected = CountDetection(SimuDir, SADFolderName, "sad_unadjustable_pvalue_sorted.tsv", "evaluation_overall_exist", _type, _numevents, _id)
					if numdetected > 0:
						Result.append([_type, _numevents, _id, numsimulated, numdetected])
		WriteSensitivity(OutputFile, Result)
