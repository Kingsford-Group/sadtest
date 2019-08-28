#!/bin/python

import sys
import numpy as np
from pathlib import Path


def MergeHits(infiles, outputfile):
	status = {}
	correspondence = {}
	for fname in infiles:
		fp = open(fname, 'r')
		for line in fp:
			if line[0] == '#':
				continue
			strs = line.strip().split("\t")
			if strs[0] in status:
				if (not status[strs[0]]) and (strs[1] == "True"):
					status[strs[0]] = True
					correspondence[strs[0]] = strs[2]
			else:
				status[strs[0]] = (strs[1] == "True")
				correspondence[strs[0]] = strs[2]
		fp.close()
		print("{}\t{}".format(fname, len([k for k,v in status.items() if v]) ))
		assert( len(status) == len(correspondence) )
	# write output
	fp = open(outputfile, 'w')
	for k,v in status.items():
		fp.write("{}\t{}\t{}\n".format(k, v, correspondence[k]))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python MergeLongreadNucmer.py <hitprefix> <outputfile>")
	else:
		hitprefix = sys.argv[1]
		outputfile = sys.argv[2]

		folder = "/".join(hitprefix.split("/")[:-1])
		prefix = hitprefix.split("/")[-1]
		infiles = [str(x) for x in list(Path(folder).glob(prefix + "*"))]
		MergeHits(infiles, outputfile)
