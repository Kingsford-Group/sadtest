#!/bin/python

import sys
import numpy as np
from PlotDisttributionTrans import *
from TranscriptClass import *
from pathlib import Path


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python3 plot_examples.py <prepare data directory> <figure output directory>")
	else:
		prepdir = sys.argv[1]
		outdir = sys.argv[2]

		# ADJUSTABLE ANOMALY EXAMPLE 1
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/correction.dat").exists())
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/startpos.dat").exists())
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/sad/sad_expectedbinnorm.dat").exists())
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/sad/sad_covariance.dat").exists())

		Expected = ReadRawCorrection(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/correction.dat")
		Observed = ReadRawStartPos(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/startpos.dat")
		AdjExpected = ReadAjustedTheoDistribution(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/sad/sad_expectedbinnorm.dat", Expected)
		AdjObserved = ReadLPAdjustment(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/sad/sad_adjusted")
		[Covariances, LenClass] = ReadCovariance(prepdir + "/HumanBodyMap/salmon_Full_ERR030885/sad/sad_covariance.dat")
		Transcripts = ReadGTF(prepdir + "/gencode.v26.annotation.gtf")

		t1="ENST00000545682.5"
		t2="ENST00000537601.5"
		sharedexons = GetSharedExon(Transcripts, [t1, t2])
		exp_t1 = ConvertCoordinate(AdjExpected[t1], Transcripts, sharedexons, t1)
		obs_t1 = ConvertCoordinate(Observed[t1], Transcripts, sharedexons, t1)
		exp_t2 = ConvertCoordinate(AdjExpected[t2], Transcripts, sharedexons, t2)
		obs_t2 = ConvertCoordinate(Observed[t2], Transcripts, sharedexons, t2)
		obs_adj_t1 = ConvertCoordinate(AdjObserved[t1], Transcripts, sharedexons, t1)
		obs_adj_t2 = ConvertCoordinate(AdjObserved[t2], Transcripts, sharedexons, t2)

		fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, sharex='col', figsize=(16, 6.5))
		ax1 = PlotDist_ribbon_v2(ax1, exp_t1, obs_t1, Covariances[LenClass[t1]], "gene $\it{"+Transcripts[t1].GeneName+"}$, transcript "+t1, binsize=25, legend_bottom=None, coordinate="gene", ytop_lim=0.15)
		ax2 = PlotDist_ribbon_v2(ax2, exp_t2, obs_t2, Covariances[LenClass[t2]], "gene $\it{"+Transcripts[t2].GeneName+"}$, transcript "+t2, binsize=25, legend_bottom=True, coordinate="gene", ytop_lim=None)
		ax3 = PlotDist_ribbon_v2(ax3, exp_t1, obs_adj_t1, Covariances[LenClass[t1]], "reassigned: gene $\it{"+Transcripts[t1].GeneName+"}$, transcript "+t1, binsize=25, legend_bottom=None, coordinate="gene", ytop_lim=0.15)
		ax4 = PlotDist_ribbon_v2(ax4, exp_t2, obs_adj_t2, Covariances[LenClass[t2]], "reassigned: gene $\it{"+Transcripts[t2].GeneName+"}$, transcript "+t2, binsize=25, legend_bottom=True, coordinate="gene", ytop_lim=None)
		fig.subplots_adjust(left=0.05, right=0.975, hspace=0.3, bottom = 0.15)
		fig.savefig(outdir + "/misquant_hubmap_kidney_ENSG00000172663.pdf", transparent=True)

		
		# ADJUSTABLE ANOMALY EXAMPLE 2
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188088/correction.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188088/startpos.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188088/sad/sad_expectedbinnorm.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188088/sad/sad_covariance.dat").exists())

		Expected = ReadRawCorrection(prepdir + "/GEUVADIS/salmon_Full_ERR188088/correction.dat")
		Observed = ReadRawStartPos(prepdir + + "/GEUVADIS/salmon_Full_ERR188088/startpos.dat")
		AdjExpected = ReadAjustedTheoDistribution(prepdir + "/GEUVADIS/salmon_Full_ERR188088/sad/sad_expectedbinnorm.dat", Expected)
		AdjObserved = ReadLPAdjustment(prepdir +"/GEUVADIS/salmon_Full_ERR188088/sad/sad_adjusted")
		[Covariances, LenClass] = ReadCovariance(prepdir +"/GEUVADIS/salmon_Full_ERR188088/sad/sad_covariance.dat")
		Transcripts = ReadGTF(prepdir + "/gencode.v26.annotation.gtf")

		t1 = "ENST00000532808.5"
		t2 = "ENST00000263464.7"
		sharedexons = GetSharedExon(Transcripts, [t1, t2])
		exp_t1 = ConvertCoordinate(AdjExpected[t1], Transcripts, sharedexons, t1)
		obs_t1 = ConvertCoordinate(Observed[t1], Transcripts, sharedexons, t1)
		exp_t2 = ConvertCoordinate(AdjExpected[t2], Transcripts, sharedexons, t2)
		obs_t2 = ConvertCoordinate(Observed[t2], Transcripts, sharedexons, t2)
		obs_adj_t1 = ConvertCoordinate(AdjObserved[t1], Transcripts, sharedexons, t1)
		obs_adj_t2 = ConvertCoordinate(AdjObserved[t2], Transcripts, sharedexons, t2)

		fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, sharex='col', figsize=(16, 6.5))
		ax1 = PlotDist_ribbon_v2(ax1, exp_t1, obs_t1, Covariances[LenClass[t1]], "gene $\it{"+Transcripts[t1].GeneName+"}$, transcript "+t1, binsize=50, legend_bottom=None, coordinate="gene", ytop_lim=0.055)
		ax2 = PlotDist_ribbon_v2(ax2, exp_t2, obs_t2, Covariances[LenClass[t2]], "gene $\it{"+Transcripts[t2].GeneName+"}$, transcript "+t2, binsize=50, legend_bottom=True, coordinate="gene", ytop_lim=None)
		ax3 = PlotDist_ribbon_v2(ax3, exp_t1, obs_adj_t1, Covariances[LenClass[t1]], "reassigned: gene $\it{"+Transcripts[t1].GeneName+"}$, transcript "+t1, binsize=50, legend_bottom=None, coordinate="gene", ytop_lim=0.055)
		ax4 = PlotDist_ribbon_v2(ax4, exp_t2, obs_adj_t2, Covariances[LenClass[t2]], "reassigned: gene $\it{"+Transcripts[t2].GeneName+"}$, transcript "+t2, binsize=50, legend_bottom=True, coordinate="gene", ytop_lim=None)
		fig.subplots_adjust(left=0.05, right=0.975, hspace=0.3, bottom = 0.15)
		fig.savefig(outdir + "/misquant_geuvadis_ERR188088_ENSG00000023445.pdf", transparent=True)


		# UNADJUSTABLE ANOMALY EXAMPLES
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/correction.dat").exists())
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/startpos.dat").exists())
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/sad/sad_expectedbinnorm.dat").exists())
		assert(Path(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/sad/sad_covariance.dat").exists())
		
		Expected = ReadRawCorrection(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/correction.dat")
		Observed = ReadRawStartPos(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/startpos.dat")
		AdjExpected = ReadAjustedTheoDistribution(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/sad/sad_expectedbinnorm.dat", Expected)
		AdjObserved = ReadLPAdjustment(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/sad/sad_adjusted")
		[Covariances, LenClass] = ReadCovariance(prepdir + "/HumanBodyMap/salmon_Full_ERR030886/sad/sad_covariance.dat")
		Transcripts = ReadGTF(prepdir + "/gencode.v26.annotation.gtf")

		t1 = "ENST00000292211.4"
		t2 = "ENST00000273317.4"
		fig, ax = plt.subplots(2, figsize=(8, 6.5))
		ax[0] = PlotDist_ribbon_v2(ax[0], AdjExpected[t1], Observed[t1], Covariances[LenClass[t1]], "gene $\it{"+Transcripts[t1].GeneName+"}$, transcript "+t1, legend_bottom=None, coordinate="transcript")
		ax[1] = PlotDist_ribbon_v2(ax[1], AdjExpected[t2], Observed[t2], Covariances[LenClass[t2]], "gene $\it{"+Transcripts[t2].GeneName+"}$, transcript "+t2, legend_bottom=True, coordinate="transcript")
		fig.subplots_adjust(left=0.1, right=0.975, hspace=0.5, bottom = 0.15)
		ax[0].text(-500, 0.12, "A", fontsize=14, fontweight='bold')
		ax[1].text(-1800, 0.037, "B", fontsize=14, fontweight='bold')
		fig.savefig(outdir + "/novelseq_hubmap_heart.pdf", transparent=True)

		
		# UNIQUE DETECTIONS IN RSEM OR IN SALMON
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188265/correction.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188265/startpos.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188265/sad/sad_expectedbinnorm.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/salmon_Full_ERR188265/sad/sad_covariance.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/rsem_Full_ERR188265/correction.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/rsem_Full_ERR188265/startpos.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/rsem_Full_ERR188265/sad/sad_expectedbinnorm.dat").exists())
		assert(Path(prepdir + "/GEUVADIS/rsem_Full_ERR188265/sad/sad_covariance.dat").exists())

		Expected1 = ReadRawCorrection(prepdir + "/GEUVADIS/salmon_Full_ERR188265/correction.dat")
		Observed1 = ReadRawStartPos(prepdir + "/GEUVADIS/salmon_Full_ERR188265/startpos.dat")
		AdjExpected1 = ReadAjustedTheoDistribution(prepdir + "/GEUVADIS/salmon_Full_ERR188265/test_expectedbinnorm.dat", Expected1)
		[Covariances1, LenClass1] = ReadCovariance(prepdir + "/GEUVADIS/salmon_Full_ERR188265/test_covariance.dat")
		Transcripts = ReadGTF(prepdir + "/home/cong/Documents/SADresult/gencode.v26.annotation.gtf")

		Expected2 = ReadRawCorrection(prepdir + "/GEUVADIS/rsem_Full_ERR188265/correction.dat")
		Observed2 = ReadRawStartPos(prepdir + "/GEUVADIS/rsem_Full_ERR188265/startpos.dat")
		AdjExpected2 = ReadAjustedTheoDistribution(prepdir + "/GEUVADIS/rsem_Full_ERR188265/test_expectedbinnorm.dat", Expected1)
		[Covariances2, LenClass2] = ReadCovariance(prepdir + "/GEUVADIS/rsem_Full_ERR188265/test_covariance.dat")

		# ENST00000425389.2: unique in RSEM because salmon expected distribution is more refined and sequence specific
		t = "ENST00000425389.2"
		fig, ax = plt.subplots(2, figsize=(8, 6.5))
		ax[0] = PlotDist_ribbon_v2(ax[0], AdjExpected1[t], Observed1[t], Covariances1[LenClass1[t]], "Salmon: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=None, coordinate="transcript")
		ax[1] = PlotDist_ribbon_v2(ax[1], AdjExpected2[t], Observed2[t], Covariances2[LenClass2[t]], "RSEM: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=True, coordinate="transcript")
		fig.subplots_adjust(left=0.1, right=0.975, hspace=0.5, bottom = 0.15)
		ax[0].text(-500, 0.11, "A", fontsize=14, fontweight='bold')
		ax[1].text(-500, 0.11, "B", fontsize=14, fontweight='bold')
		fig.savefig(outdir + "/diff_salmon_rsem_1.pdf", transparent=True)

		# ENST00000380381.3: unique in RSEM because different observed distribution
		t = "ENST00000380381.3"
		fig, ax = plt.subplots(2, figsize=(8, 6.5))
		ax[0] = PlotDist_ribbon_v2(ax[0], AdjExpected1[t], Observed1[t], Covariances1[LenClass1[t]], "Salmon: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=None, coordinate="transcript")
		ax[1] = PlotDist_ribbon_v2(ax[1], AdjExpected2[t], Observed2[t], Covariances2[LenClass2[t]], "RSEM: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=True, coordinate="transcript")
		fig.subplots_adjust(left=0.1, right=0.975, hspace=0.5, bottom = 0.15)
		ax[0].text(-135, 0.4, "A", fontsize=14, fontweight='bold')
		ax[1].text(-135, 0.75, "B", fontsize=14, fontweight='bold')
		fig.savefig(outdir + "/diff_salmon_rsem_2.pdf", transparent=True)

		# ENST00000339647.5: unique in Salmon because the different and unfixable read assignment
		t = "ENST00000339647.5"
		fig, ax = plt.subplots(2, figsize=(8, 6.5))
		ax[0] = PlotDist_ribbon_v2(ax[0], AdjExpected1[t], Observed1[t], Covariances1[LenClass1[t]], "Salmon: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=None, coordinate="transcript")
		ax[1] = PlotDist_ribbon_v2(ax[1], AdjExpected2[t], Observed2[t], Covariances2[LenClass2[t]], "RSEM: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=True, coordinate="transcript")
		fig.subplots_adjust(left=0.1, right=0.975, hspace=0.5, bottom = 0.15)
		ax[0].text(-450, 0.125, "A", fontsize=14, fontweight='bold')
		ax[1].text(-450, 0.07, "B", fontsize=14, fontweight='bold')
		fig.savefig(outdir + "/diff_salmon_rsem_3.pdf", transparent=True)

		# ENST00000527673.1: unique in salmon, similar observed distribution but RSEM has a larger covariance
		t = "ENST00000527673.1"
		fig, ax = plt.subplots(2, figsize=(8, 6.5))
		ax[0] = PlotDist_ribbon_v2(ax[0], AdjExpected1[t], Observed1[t], Covariances1[LenClass1[t]], "Salmon: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=None, coordinate="transcript")
		ax[1] = PlotDist_ribbon_v2(ax[1], AdjExpected2[t], Observed2[t], Covariances2[LenClass2[t]], "RSEM: gene $\it{"+Transcripts[t].GeneName+"}$, transcript "+t, legend_bottom=True, coordinate="transcript")
		fig.subplots_adjust(left=0.1, right=0.975, hspace=0.5, bottom = 0.15)
		ax[0].text(-125, 0.5, "A", fontsize=14, fontweight='bold')
		ax[1].text(-125, 0.5, "B", fontsize=14, fontweight='bold')
		fig.savefig(outdir + "/diff_salmon_rsem_4.pdf", transparent=True)

