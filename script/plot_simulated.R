library(ggplot2)
library(cowplot)
library(assertthat)

options(stringsAsFactors = FALSE)

ReadSADEvaluationFile <- function(evafilename, pvaluefilename, binsize, pthresh = 0.01) {
	# read p value file, determine the maxmum number of records to read
	t1 <- read.table(pvaluefilename, header=F, sep="\t")
	colnames(t1) <- c("Name", "Coverage", "DeletionScorePos", "DeletioScoreNeg", "RawPvalue", "AdjustedPvalue", "Choice")
	maxsize <- sum(t1$AdjustedPvalue < pthresh)
	assert_that(maxsize > binsize)
	# read evaluation file
	t2 <- read.table(evafilename, header=T, sep="\t")
	# store the accuracy for each bin
	res <- NULL
	for (i in seq(binsize, maxsize, binsize)) {
		tp <- sum(t2$Hit[1:i])
		accuracy <- tp/i
		tmp <- data.frame(bin=c(i), TP=c(tp), Accuracy=c(accuracy), DiffPredictions=c(0))
		if (is.null(res)) {
			res <- tmp
		} else {
			res <- rbind(res, tmp)
		}
	}
	# store the accuracy for the actual number of predictions
	if (maxsize %% binsize == 0) {
		res[nrow(res), "DiffPredictions"] <- maxsize
	} else {
		tp <- sum(t2$Hit[1:maxsize])
		accuracy <- tp/maxsize
		tmp <- data.frame(bin=c(floor(maxsize/binsize+1)*binsize), TP=c(tp), Accuracy=c(accuracy), DiffPredictions=c(maxsize))
		res <- rbind(res, tmp)
	}
	# return result
	return(res)
}


ReadOtherEvaluationFile <- function(evafilename, maxsize, binsize) {
	# read evaluation file
	t2 <- read.table(evafilename, header=T, sep="\t")
	# store the accuracy for each bin
	res <- NULL
	for (i in seq(binsize, maxsize, binsize)) {
		tp <- sum(t2$Hit[1:i])
		accuracy <- tp/i
		tmp <- data.frame(bin=c(i), TP=c(tp), Accuracy=c(accuracy), DiffPredictions=c(0))
		if (is.null(res)) {
			res <- tmp
		} else {
			res <- rbind(res, tmp)
		}
	}
	# store the accuracy for the actual number of predictions
	if (maxsize %% binsize == 0) {
		res[nrow(res), "DiffPredictions"] <- maxsize
	} else {
		tp <- sum(t2$Hit[1:maxsize])
		accuracy <- tp/maxsize
		tmp <- data.frame(bin=c(floor(maxsize/binsize+1)*binsize), TP=c(tp), Accuracy=c(accuracy), DiffPredictions=c(maxsize))
		res <- rbind(res, tmp)
	}
	# return result
	return(res)
}

RestrictedSADQuantification <- function(exprdifffile, salmonquant, SADquant) {
	# read expression_difference, salmonquant, SADquant files
	t0 <- read.table(exprdifffile, header=F, sep="\t")
	colnames(t0) <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads", "nsimulated", "diffTPM")
	t1 <- read.table(salmonquant, header=T, sep="\t")
	t2 <- read.table(SADquant, header=F, sep="\t")
	colnames(t2) <- c("Name", "Length", "NumReads")
	# only keep the overlapping transcripts with SADquant
	t0 <- t0[match(t2$Name, t0$Name), ]
	assert_that(!any(is.na(t0)))
	t1 <- t1[match(t2$Name, t1$Name), ]
	assert_that(!any(is.na(t1)))
	# calculate correlations of salmon
	salmon_pearson <- cor(t0$nsimulated, t1$NumReads, method="pearson")
	salmon_spearman <- cor(t0$nsimulated, t1$NumReads, method="spearman")
	salmon_ard <- abs(t0$nsimulated-t1$NumReads) / (t0$nsimulated+t1$NumReads)
	salmon_ard[which(is.na(salmon_ard))] <- 0
	salmon_ardmean <- mean(salmon_ard)
	salmon_ardmedian <- median(salmon_ard)
	# calculate correlations of SAD
	sad_pearson <- cor(t0$nsimulated, t2$NumReads, method="pearson")
	sad_spearman <- cor(t0$nsimulated, t2$NumReads, method="spearman")
	sad_ard <- abs(t0$nsimulated-t2$NumReads) / (t0$nsimulated+t2$NumReads)
	sad_ard[which(is.na(sad_ard))] <- 0
	sad_ardmean <- mean(sad_ard)
	sad_ardmedian <- median(sad_ard)
	# return result
	res <- data.frame(salmon_pearson=salmon_pearson, salmon_spearman=salmon_spearman, salmon_ardmean=salmon_ardmean, salmon_ardmedian=salmon_ardmedian, 
		sad_pearson=sad_pearson, sad_spearman=sad_spearman, sad_ardmean=sad_ardmean, sad_ardmedian=sad_ardmedian)
	return(res)
}


RestrictedSADQuantification_isorange <- function(exprdifffile, salmonquant, SADquant, Tx2Gene, lb, ub) {
	# read expression_difference, salmonquant, SADquant files
	t0 <- read.table(exprdifffile, header=F, sep="\t")
	colnames(t0) <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads", "nsimulated", "diffTPM")
	t1 <- read.table(salmonquant, header=T, sep="\t")
	t2 <- read.table(SADquant, header=F, sep="\t")
	colnames(t2) <- c("Name", "Length", "NumReads")
	# add the corresponding gene name
	t2$GeneName <- Tx2Gene[match(t2$Name, Tx2Gene$TransName), "GeneName"]
	gene_count <- data.frame(table(t2$GeneName))
	genelist <- gene_count[gene_count$Freq >= lb & gene_count$Freq < ub, "Var1"]
	t2 <- t2[t2$GeneName %in% genelist, ]
	# only keep the overlapping transcripts with SADquant
	t0 <- t0[match(t2$Name, t0$Name), ]
	assert_that(!any(is.na(t0)))
	t1 <- t1[match(t2$Name, t1$Name), ]
	assert_that(!any(is.na(t1)))
	# calculate correlations of salmon
	salmon_pearson <- cor(t0$nsimulated, t1$NumReads, method="pearson")
	salmon_spearman <- cor(t0$nsimulated, t1$NumReads, method="spearman")
	salmon_ard <- abs(t0$nsimulated-t1$NumReads) / (t0$nsimulated+t1$NumReads)
	salmon_ard[which(is.na(salmon_ard))] <- 0
	salmon_ardmean <- mean(salmon_ard)
	salmon_ardmedian <- median(salmon_ard)
	# calculate correlations of SAD
	sad_pearson <- cor(t0$nsimulated, t2$NumReads, method="pearson")
	sad_spearman <- cor(t0$nsimulated, t2$NumReads, method="spearman")
	sad_ard <- abs(t0$nsimulated-t2$NumReads) / (t0$nsimulated+t2$NumReads)
	sad_ard[which(is.na(sad_ard))] <- 0
	sad_ardmean <- mean(sad_ard)
	sad_ardmedian <- median(sad_ard)
	# return result
	res <- data.frame(salmon_pearson=salmon_pearson, salmon_spearman=salmon_spearman, salmon_ardmean=salmon_ardmean, salmon_ardmedian=salmon_ardmedian, 
		sad_pearson=sad_pearson, sad_spearman=sad_spearman, sad_ardmean=sad_ardmean, sad_ardmedian=sad_ardmedian)
	return(res)
}


FullSADQuantification <- function(exprdifffile, salmonquant, SADquant) {
	# read expression_difference, salmonquant, SADquant files
	t0 <- read.table(exprdifffile, header=F, sep="\t")
	colnames(t0) <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads", "nsimulated", "diffTPM")
	t1 <- read.table(salmonquant, header=T, sep="\t")
	t2 <- read.table(SADquant, header=F, sep="\t")
	colnames(t2) <- c("Name", "Length", "NumReads")
	# adding salmon quant of the rest of transcripts back to SAD
	t3 <- t1
	t3[match(t2$Name, t3$Name), "NumReads"] <- t2$NumReads
	t2 <- t3
	# calculate correlations of salmon
	salmon_pearson <- cor(t0$nsimulated, t1$NumReads, method="pearson")
	salmon_spearman <- cor(t0$nsimulated, t1$NumReads, method="spearman")
	salmon_ard <- abs(t0$nsimulated-t1$NumReads) / (t0$nsimulated+t1$NumReads)
	salmon_ard[which(is.na(salmon_ard))] <- 0
	salmon_ardmean <- mean(salmon_ard)
	salmon_ardmedian <- median(salmon_ard)
	# calculate correlations of SAD
	sad_pearson <- cor(t0$nsimulated, t2$NumReads, method="pearson")
	sad_spearman <- cor(t0$nsimulated, t2$NumReads, method="spearman")
	sad_ard <- abs(t0$nsimulated-t2$NumReads) / (t0$nsimulated+t2$NumReads)
	sad_ard[which(is.na(sad_ard))] <- 0
	sad_ardmean <- mean(sad_ard)
	sad_ardmedian <- median(sad_ard)
	# return result
	res <- data.frame(salmon_pearson=salmon_pearson, salmon_spearman=salmon_spearman, salmon_ardmean=salmon_ardmean, salmon_ardmedian=salmon_ardmedian, 
		sad_pearson=sad_pearson, sad_spearman=sad_spearman, sad_ardmean=sad_ardmean, sad_ardmedian=sad_ardmedian)
	return(res)
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	print("Rscript plot_simulated.R <prepare data directory> <output figure directory>")
} else {
	indir <- args[1]
	outdir <- args[2]

	Type <- c("PC", "Full")
	NumEvents <- c(200, 500, 1000, 1500)
	RefQuantID <- c("GEU", "GM12878", "K562")
	LocalDir <- paste0(indir, "/Simulation", sep="")

	#################### unadjustable anomaly ####################
	binsize <- 50
	ResultType <- c()
	ResultNumEvents <- c()
	ResultRefQuantID <- c()
	ResultBin <- c()
	ResultAccuracy_SAD <- c()
	ResultAccuracy_stringtie <- c()
	ResultAccuracy_scallop <- c()
	ResultIsMaxPred <- c()

	for (t in Type) {
		for (n in NumEvents) {
			for (ID in RefQuantID) {
				SADevaFile <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/sad/evaluation_genelevel_overall_exist", sep="")
				SADpvalueFile <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/sad/sad_unadjustable_pvalue_sorted_uniqgene.tsv", sep="")
				StringtieFile <- paste0(LocalDir, "/Star/", t, "_", ID, "_", n, "/evaluation_stringtie_existing_novelisoforms_existjunc", sep="")
				ScallopFile <- paste0(LocalDir, "/Star/", t, "_", ID, "_", n, "/evaluation_scallop_existing_novelisoforms_existjunc", sep="")

				# check whether SAD run successfully and whether evaluation file exists
				if (!file.exists(SADevaFile) || !file.exists(SADpvalueFile)) {
					next
				}

				# read evaluation result
				res_SAD <- ReadSADEvaluationFile(SADevaFile, SADpvalueFile, binsize)
				maxsize <- res_SAD$DiffPredictions[nrow(res_SAD)]
				assert_that(maxsize != 0)
				res_stringtie <- ReadOtherEvaluationFile(StringtieFile, maxsize, binsize)
				res_scallop <- ReadOtherEvaluationFile(ScallopFile, maxsize, binsize)
				assert_that(nrow(res_SAD) == nrow(res_stringtie))
				assert_that(nrow(res_SAD) == nrow(res_scallop))

				ResultType <- c( ResultType, rep(t, nrow(res_SAD)) )
				ResultNumEvents <- c( ResultNumEvents, rep(n, nrow(res_SAD)) )
				ResultRefQuantID <- c( ResultRefQuantID, rep(ID, nrow(res_SAD))  )
				ResultBin <- c( ResultBin, res_SAD$bin )
				ResultAccuracy_SAD <- c( ResultAccuracy_SAD, res_SAD$Accuracy )
				ResultAccuracy_stringtie <- c( ResultAccuracy_stringtie, res_stringtie$Accuracy )
				ResultAccuracy_scallop <- c( ResultAccuracy_scallop, res_scallop$Accuracy )
				ResultIsMaxPred <- c(ResultIsMaxPred, res_SAD$DiffPredictions )
				print(c(t, n, ID, nrow(ResultAccuracy_SAD), nrow(ResultAccuracy_stringtie), nrow(ResultAccuracy_scallop)))
			}
		}
	}
	t <- data.frame(ResultType, ResultNumEvents, ResultRefQuantID, ResultBin, ResultAccuracy_SAD, ResultAccuracy_stringtie, ResultAccuracy_scallop, ResultIsMaxPred)
	colnames(t) <- c("Type", "NumEvents", "ID", "Bin", "Accuracy_SAD", "Accuracy_Stringtie", "Accuracy_SCALLOP", "IsMaxPred")
	t$Type <- as.character(t$Type)
	t[t$Type=="PC", "Type"] <- "Protein-coding"

	p_stringtie <- ggplot(t[t$IsMaxPred != 0, ]) + geom_point(aes(x=Accuracy_Stringtie, y=Accuracy_SAD, color=Type, shape=as.factor(NumEvents)), size=2.5) + 
		scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + geom_abline(slope=1, intercept=0) + background_grid() + 
		labs(shape="Number events", color="Annotation type", title="Unannotated isoform detection precision against StringTie", x="StringTie precision", y="SAD precision")
	p_scallop <- ggplot(t[t$IsMaxPred != 0, ]) + geom_point(aes(x=Accuracy_SCALLOP, y=Accuracy_SAD, color=Type, shape=as.factor(NumEvents)), size=2.5) + 
		scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) + geom_abline(slope=1, intercept=0) + background_grid() + 
		labs(shape="Number events", color="Annotation type", title="Unannotated isoform detection precision against Scallop", x="Scallop precision", y="SAD precision")

	t <- read.table(paste0(LocalDir, "/Sensitivity.txt", sep=""), header=F, sep="\t")
	colnames(t) <- c("Type", "NumEvents", "ID", "NumSimulated", "NumDetected", "Sensitivity")
	t$Type <- as.character(t$Type)
	t[t$Type=="PC", "Type"] <- "Protein-coding"
	p9 <- ggplot(t) + geom_boxplot(aes(x=as.factor(NumEvents), y=Sensitivity, color=Type), width=0.3) + labs(x="Number events", title="Sensitivity of SAD", color="Annotation type") + background_grid()


	#################### adjustable anomaly ####################
	ResultType <- c()
	ResultNumEvents <- c()
	ResultRefQuantID <- c()
	ResultPearson_Salmon <- c()
	ResultPearson_SAD <- c()
	ResultSpearman_Salmon <- c()
	ResultSpearman_SAD <- c()
	ResultARDMean_Salmon <- c()
	ResultARDMean_SAD <- c()
	ResultARDMedian_Salmon <- c()
	ResultARDMedian_SAD <- c()

	for (t in Type) {
		for (n in NumEvents) {
			for (ID in RefQuantID) {
				ExprDiffFile <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/sad/expression_difference.txt", sep="")
				SalmonQuant <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/sad/quant.sf", sep="")
				SADQuant <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/sad/sad_adjusted_quantification.tsv", sep="")
	
				# check whether SAD run successfully and whether evaluation file exists
				if (!file.exists(SADQuant)) {
					next
				}

				# read evaluation result
				res <- RestrictedSADQuantification(ExprDiffFile, SalmonQuant, SADQuant)

				ResultType <- c(ResultType, t)
				ResultNumEvents <- c(ResultNumEvents, n)
				ResultRefQuantID <- c(ResultRefQuantID, ID)
				ResultPearson_Salmon <- c(ResultPearson_Salmon, res$salmon_pearson)
				ResultPearson_SAD <- c(ResultPearson_SAD, res$sad_pearson)
				ResultSpearman_Salmon <- c(ResultSpearman_Salmon, res$salmon_spearman)
				ResultSpearman_SAD <- c(ResultSpearman_SAD, res$sad_spearman)
				ResultARDMean_Salmon <- c(ResultARDMean_Salmon, res$salmon_ardmean)
				ResultARDMean_SAD <- c(ResultARDMean_SAD, res$sad_ardmean)
				ResultARDMedian_Salmon <- c(ResultARDMedian_Salmon, res$salmon_ardmedian)
				ResultARDMedian_SAD <- c(ResultARDMedian_SAD, res$sad_ardmedian)
			}
		}
	}

	t <- data.frame(ResultType, ResultNumEvents, ResultRefQuantID, ResultPearson_Salmon, ResultPearson_SAD, ResultSpearman_Salmon, ResultSpearman_SAD, 
		ResultARDMean_Salmon, ResultARDMean_SAD, ResultARDMedian_Salmon, ResultARDMedian_SAD)
	colnames(t) <- c("Type", "NumEvents", "ID", "PearsonR_Salmon", "PearsonR_SAD", "SpearmanR_Salmon", "SpearmanR_SAD", 
		"ARDMean_Salmon", "ARDMean_SAD", "ARDMedian_Salmon", "ARDMedian_SAD")
	t$Type <- as.character(t$Type)
	t[t$Type=="PC", "Type"] <- "Protein-coding"
	df <- data.frame(Type=rep("Correlation", nrow(t)), Metric=rep("Pearson correlation (SAD - Salmon)", nrow(t)), Value=t$PearsonR_SAD-t$PearsonR_Salmon)
	df <- rbind( df, data.frame(Type=rep("Correlation", nrow(t)), Metric=rep("Spearman correlation (SAD - Salmon)", nrow(t)), Value=t$SpearmanR_SAD-t$SpearmanR_Salmon) )
	df <- rbind( df, data.frame(Type=rep("Distance", nrow(t)), Metric=rep("Mean ARD (Salmon - ARD)", nrow(t)), Value=t$ARDMean_Salmon-t$ARDMean_SAD) )
	df <- rbind( df, data.frame(Type=rep("Distance", nrow(t)), Metric=rep("Median ARD (Salmon - ARD)", nrow(t)), Value=t$ARDMedian_Salmon-t$ARDMedian_SAD) )
	p_violin <- ggplot(df) + geom_violin(aes(x="", y=Value, color=Metric), size=0.8) + background_grid() + scale_y_continuous(limits=c(-0.1, 0.5), breaks=seq(-0.05, 0.5, 0.05)) +
		labs(x="", y="Accuracy difference", title="Accuracy difference between Salmon and SAD on the adjusted transcripts") + theme(axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.1))

	Tx2Gene <- read.table(paste0(indir, "/gencode.v26.annotation.tx2gene.txt", sep=""), header=T, sep="\t")	
	ResultType <- c()
	ResultNumEvents <- c()
	ResultRefQuantID <- c()
	ResultNumIsoforms <- c()
	ResultPearson_Salmon <- c()
	ResultPearson_SAD <- c()
	ResultSpearman_Salmon <- c()
	ResultSpearman_SAD <- c()
	ResultARDMean_Salmon <- c()
	ResultARDMean_SAD <- c()
	ResultARDMedian_Salmon <- c()
	ResultARDMedian_SAD <- c()

	for (t in Type) {
		for (n in NumEvents) {
			for (ID in RefQuantID) {
				ExprDiffFile <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/expression_difference.txt", sep="")
				SalmonQuant <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/quant.sf", sep="")
				SADQuant <- paste0(LocalDir, "/salmon/", t, "_", ID, "_", n, "/sad/sad_adjusted_quantification.tsv", sep="")

				# check whether SAD run successfully and whether evaluation file exists
				if (!file.exists(SADQuant)) {
					next
				}

				# read evaluation result
				for (numiso in seq(2, 5)) {
					if (numiso != 5) {
						res <- RestrictedSADQuantification_isorange(ExprDiffFile, SalmonQuant, SADQuant, Tx2Gene, numiso, numiso+1)
					} else {
						res <- RestrictedSADQuantification_isorange(ExprDiffFile, SalmonQuant, SADQuant, Tx2Gene, numiso, 100000)
					}

					if (!any(is.na(res))) {
						ResultType <- c(ResultType, t)
						ResultNumEvents <- c(ResultNumEvents, n)
						ResultRefQuantID <- c(ResultRefQuantID, ID)
						if (numiso != 5) {
							ResultNumIsoforms <- c(ResultNumIsoforms, as.character(numiso))
						} else {
							ResultNumIsoforms <- c(ResultNumIsoforms, ">=5")
						}
						ResultPearson_Salmon <- c(ResultPearson_Salmon, res$salmon_pearson)
						ResultPearson_SAD <- c(ResultPearson_SAD, res$sad_pearson)
						ResultSpearman_Salmon <- c(ResultSpearman_Salmon, res$salmon_spearman)
						ResultSpearman_SAD <- c(ResultSpearman_SAD, res$sad_spearman)
						ResultARDMean_Salmon <- c(ResultARDMean_Salmon, res$salmon_ardmean)
						ResultARDMean_SAD <- c(ResultARDMean_SAD, res$sad_ardmean)
						ResultARDMedian_Salmon <- c(ResultARDMedian_Salmon, res$salmon_ardmedian)
						ResultARDMedian_SAD <- c(ResultARDMedian_SAD, res$sad_ardmedian)
					}
				}

				# evaluation result for all adjusted
				res <- RestrictedSADQuantification_isorange(ExprDiffFile, SalmonQuant, SADQuant, Tx2Gene, 0, 100000)
				if (!any(is.na(res))) {
					ResultType <- c(ResultType, t)
					ResultNumEvents <- c(ResultNumEvents, n)
					ResultRefQuantID <- c(ResultRefQuantID, ID)
					ResultNumIsoforms <- c(ResultNumIsoforms, "overall")
					ResultPearson_Salmon <- c(ResultPearson_Salmon, res$salmon_pearson)
					ResultPearson_SAD <- c(ResultPearson_SAD, res$sad_pearson)
					ResultSpearman_Salmon <- c(ResultSpearman_Salmon, res$salmon_spearman)
					ResultSpearman_SAD <- c(ResultSpearman_SAD, res$sad_spearman)
					ResultARDMean_Salmon <- c(ResultARDMean_Salmon, res$salmon_ardmean)
					ResultARDMean_SAD <- c(ResultARDMean_SAD, res$sad_ardmean)
					ResultARDMedian_Salmon <- c(ResultARDMedian_Salmon, res$salmon_ardmedian)
					ResultARDMedian_SAD <- c(ResultARDMedian_SAD, res$sad_ardmedian)
				}
				print(c(t, n, ID))
			}
		}
	}
	df <- data.frame(ResultType, ResultNumEvents, ResultRefQuantID, ResultNumIsoforms, ResultPearson_Salmon, ResultPearson_SAD, ResultSpearman_Salmon, ResultSpearman_SAD, 
		ResultARDMean_Salmon, ResultARDMean_SAD, ResultARDMedian_Salmon, ResultARDMedian_SAD)
	colnames(df) <- c("Type", "NumEvents", "ID", "NumIsoforms", "PearsonR_Salmon", "PearsonR_SAD", "SpearmanR_Salmon", "SpearmanR_SAD", 
		"ARDMean_Salmon", "ARDMean_SAD", "ARDMedian_Salmon", "ARDMedian_SAD")
	df$NumIsoforms <- factor(df$NumIsoforms, levels=c("overall", "2", "3", "4", ">=5"))
	p_meanard <- ggplot(df) + geom_violin(aes(x = NumIsoforms, y = ARDMean_Salmon-ARDMean_SAD)) + background_grid() +
		labs(x="Number of involved isoforms within each gene", y="Mean ARD (Salmon - ARD)", title="Accuracy comparison of expression quantification \non the transcripts involved in re-assignment")

	
	#################### merge plots ####################
	ptmp1 <- plot_grid(p_stringtie + theme(legend.position="none"), p_scallop + theme(legend.position="none"), nrow=1, labels=c("A", "B"))
	ptmp1 <- plot_grid(ptmp1, get_legend(p_stringtie + theme(legend.position="bottom")), nrow=2, rel_heights=c(1, 0.1))
	ptmp2 <- plot_grid(p9 + theme(legend.position="bottom"), p_meanard, labels=c("C", "D"))
	p <- ggdraw() + 
		draw_plot(ptmp1, 0, 0.66, 1, 0.34) + 
		draw_plot(ptmp2, 0, 0.32, 1, 0.34) + 
		draw_plot(p_violin, 0.15, 0, 0.69, 0.32) + 
		draw_plot_label(c("E"), c(0.145), c(0.32))
	save_plot(paste0(outdir, "/combine_simulate.pdf", sep=""), p, base_height = 14, base_aspect_ratio = 0.9)
}
