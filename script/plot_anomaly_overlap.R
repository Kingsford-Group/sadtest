library(ggplot2)
library(cowplot)
library(assertthat)

options(stringsAsFactors = FALSE)


ReadMeta <- function(metafile) {
	Centers  <- c()
	IDs <- c()
	con = file(metafile, "r")
	line <- readLines(con, n = 1)
	while (length(line) > 0) {
		strs <- strsplit(line, "\t")[[1]]
		Centers <- c(Centers, strs[2])
		IDs <- c(IDs, strs[length(strs)])
		line <- readLines(con, n = 1)
	}
	close(con)
	t_meta <- data.frame(Center=Centers, ID=IDs)
	return(t_meta)
}

ReadNovelPrediction <- function(folder, ID) {
	filename <- paste0(folder, "/salmon_Full_", ID, "/sad/sad_unadjustable_pvalue_sorted.tsv", sep="")
	t <- read.table(filename, header=F, sep="\t")
	# colnames(t) <- c("Name", "Coverage", "DeletionScorePos", "DeletionScoreNeg", "RawPvalue", "AdjustedPvalue", "Choice")
	colnames(t) <- c("Name", "Coverage", "AnomalyScorePos", "RegionStartPos", "RegionEndPos", "PValue_Pos", "AdjPValue_Pos",
		"AnomalyScoreNeg", "RegionStartNeg", "RegionEndNeg", "PValue_Neg", "AdjPValue_Neg", "AdjustedPvalue", "Choice")
	t <- t[t$AdjustedPvalue < 0.01, ]
	return(t$Name)
}

ReadAssemblyCompare <- function(folder, ID) {
	scallopcompfile <- paste0(folder, "/salmon_Full_", ID, "/sad/sad_unadjustable_scallopcomp", sep="")
	stringtiecompfile <- paste0(folder, "/salmon_Full_", ID, "/sad/sad_unadjustable_stringtiecomp", sep="")
	t_res <- NULL
	# read scallop comparison and add to data frame
	t_scallop <- read.table(scallopcompfile, header=F, sep="\t")
	colnames(t_scallop) <- c("TransID", "GeneID", "Score", "RegionStart", "RegionEnd", "IsOverExpression", "CorrespondingAssembly", "OnEnds")
	t_res <- data.frame(ID=c(ID), Ratio=c(nrow(t_scallop[t_scallop$CorrespondingAssembly!="None", ]) / nrow(t_scallop)), Assembler=c("Scallop"), AssemblyContain=c(TRUE))
	t_res <- rbind(t_res, data.frame(ID=c(ID), Ratio=c(nrow(t_scallop[t_scallop$CorrespondingAssembly=="None", ]) / nrow(t_scallop)), Assembler=c("Scallop"), AssemblyContain=c(FALSE)) )
	# read stringtie comparison and add to data frame
	t_stringtie <- read.table(stringtiecompfile, header=F, sep="\t")
	colnames(t_stringtie) <- c("TransID", "GeneID", "Score", "RegionStart", "RegionEnd", "IsOverExpression", "CorrespondingAssembly", "OnEnds")
	t_res <- rbind(t_res, data.frame(ID=c(ID), Ratio=c(nrow(t_stringtie[t_stringtie$CorrespondingAssembly!="None", ]) / nrow(t_stringtie)), Assembler=c("StringTie"), AssemblyContain=c(TRUE)) )
	t_res <- rbind(t_res, data.frame(ID=c(ID), Ratio=c(nrow(t_stringtie[t_stringtie$CorrespondingAssembly=="None", ]) / nrow(t_stringtie)), Assembler=c("StringTie"), AssemblyContain=c(FALSE)) )
	return(t_res)
}

CompareSalmonRSEM <- function(folder, ID) {
	# check file exists
	if (!file.exists(paste0(folder, "/rsem_Full_", ID, "/sad/sad_unadjustable_pvalue_sorted.tsv", sep="")) || !file.exists(paste0(folder, "/rsem_Full_", ID, "/sad/sad_adjusted_quantification.tsv", sep=""))) {
		return(NULL)
	}
	# unadjustable 
	t_rsem <- read.table( paste0(folder, "/rsem_Full_", ID, "/sad/sad_unadjustable_pvalue_sorted.tsv", sep=""), header=F, sep="\t" )
	colnames(t_rsem) <- c("Name","Coverage","AnomalyScorePos","RegionStartPos","RegionEndPos","PValue_Pos","AdjPValue_Pos","AnomalyScoreNeg","RegionStartNeg","RegionEndNeg","PValue_Neg","AdjPValue_Neg","MinAdjPValue","Choice")
	t_salmon <- read.table( paste0(folder, "/salmon_Full_", ID, "/sad/sad_unadjustable_pvalue_sorted.tsv", sep=""), header=F, sep="\t" )
	colnames(t_salmon) <- c("Name","Coverage","AnomalyScorePos","RegionStartPos","RegionEndPos","PValue_Pos","AdjPValue_Pos","AnomalyScoreNeg","RegionStartNeg","RegionEndNeg","PValue_Neg","AdjPValue_Neg","MinAdjPValue","Choice")
	# calculate overlap and unique
	u_rsem <- setdiff(t_rsem[t_rsem$MinAdjPValue < 0.01, "Name"], t_salmon[t_salmon$MinAdjPValue < 0.01, "Name"])
	u_salmon <- setdiff(t_salmon[t_salmon$MinAdjPValue < 0.01, "Name"], t_rsem[t_rsem$MinAdjPValue < 0.01, "Name"])
	overlap <- intersect(t_rsem[t_rsem$MinAdjPValue < 0.01, "Name"], t_salmon[t_salmon$MinAdjPValue < 0.01, "Name"])
	sum_length <- length(u_rsem) + length(u_salmon) + length(overlap)
	# collect result of unadjustable anomaly overlap
	t_res <- data.frame(ID=ID, Type="unadjustable", Label="rsem unique", Number=length(u_rsem), Value=length(u_rsem)/sum_length)
	t_res <- rbind(t_res, data.frame(ID=ID, Type="unadjustable", Label="salmon unique", Number=length(u_salmon), Value=length(u_salmon)/sum_length))
	t_res <- rbind(t_res, data.frame(ID=ID, Type="unadjustable", Label="overlap", Number=length(overlap), Value=length(overlap)/sum_length))

	# adjustable anomaly
	t_rsem <- read.table( paste0(folder, "/rsem_Full_", ID, "/sad/sad_adjusted_quantification.tsv", sep=""), header=F, sep="\t" )
	colnames(t_rsem) <- c("Name", "Length", "NumReads")
	t_salmon <- read.table( paste0(folder, "/salmon_Full_", ID, "/sad/sad_adjusted_quantification.tsv", sep=""), header=F, sep="\t" )
	colnames(t_salmon) <- c("Name", "Length", "NumReads")
	# calculate overlap and unique
	u_rsem <- setdiff(t_rsem$Name, t_salmon$Name)
	u_salmon <- setdiff(t_salmon$Name, t_rsem$Name)
	overlap <- intersect(t_rsem$Name, t_salmon$Name)
	sum_length <- length(u_rsem) + length(u_salmon) + length(overlap)
	# collect the result of adjustable anomaly overlap
	t_res <- rbind(t_res, data.frame(ID=ID, Type="adjustable", Label="rsem unique", Number=length(u_rsem), Value=length(u_rsem)/sum_length) )
	t_res <- rbind(t_res, data.frame(ID=ID, Type="adjustable", Label="salmon unique", Number=length(u_salmon), Value=length(u_salmon)/sum_length) )
	t_res <- rbind(t_res, data.frame(ID=ID, Type="adjustable", Label="overlap", Number=length(overlap), Value=length(overlap)/sum_length) )

	return(t_res)
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	print("Rscript plot_anomaly_overlap.R <prepare data directory> <output figure directory>")
} else {
	indir <- args[1]
	outdir <- args[2]

	GFolder <- paste0(indir, "/GEUVADIS", sep="")
	HFolder <- paste0(indir, "/HumanBodyMap", sep="")
	t_meta <- ReadMeta( paste0(indir, "/GEUVADIS/Metadata.txt", sep=""))
	t_h_meta <- read.table(paste0(indir, "/HumanBodyMap/HumanBodyMap_SRR_Acc_List.txt", sep=""), header=F, sep="\t")
	colnames(t_h_meta) <- c("ID")

	#################### overlap with assembled transcript ####################
	t_g_compassembly <- NULL
	for (id in t_meta$ID) {
		t <- ReadAssemblyCompare(GFolder, id)
		if (is.null(t_g_compassembly)) {
			t_g_compassembly <- t
		} else {
			t_g_compassembly <- rbind(t_g_compassembly, t )
		}
	}
	t_g_compassembly[t_g_compassembly$AssemblyContain==TRUE, "AssemblyContain"] <- "Contained in assembly"
	t_g_compassembly[t_g_compassembly$AssemblyContain==FALSE, "AssemblyContain"] <- "Not contained in assembly"
	p1 <- ggplot(t_g_compassembly) + geom_boxplot(aes(x=Assembler, y=Ratio, color=AssemblyContain)) + labs(title="GEUVADIS") + theme(legend.title=element_blank()) + background_grid()

	t_h_compassembly <- NULL
	for (id in t_h_meta$ID) {
		t <- ReadAssemblyCompare(HFolder, id)
		if (is.null(t_h_compassembly)) {
			t_h_compassembly <- t
		} else {
			t_h_compassembly <- rbind(t_h_compassembly, t )
		}
	}
	t_h_compassembly[t_h_compassembly$AssemblyContain==TRUE, "AssemblyContain"] <- "Contained in assembly"
	t_h_compassembly[t_h_compassembly$AssemblyContain==FALSE, "AssemblyContain"] <- "Not contained in assembly"
	p2 <- ggplot(t_h_compassembly) + geom_boxplot(aes(x=Assembler, y=Ratio, color=AssemblyContain)) + labs(title="HumanBodyMap") + theme(legend.title=element_blank()) + background_grid()

	ptmp <- plot_grid(p1+theme(legend.position="none"), p2+theme(legend.position="none"), nrow=1, labels="AUTO")
	p <- plot_grid(ptmp, get_legend(p1), nrow=1, rel_widths=c(2, 0.5))
	save_plot(paste0(outdir, "/CompareAssembly.pdf", sep=""), p, base_aspect_ratio=2.2, base_height=5)


	#################### overlap with RSEM anomalies ####################	
	df <- NULL
	for (i in 1:length(t_meta$ID)) {
		id <- t_meta$ID[i]
		tmp <- CompareSalmonRSEM(GFolder, id)
		if (is.null(df)) {
			df <- tmp
		} else {
			df <- rbind(df, tmp)
		}
		print(id)
	}

	p1 <- ggplot(df) + geom_boxplot(aes(x=Label, y=Value, color=Type)) + labs(y="percentage among the union of RSEM and Salmon", title="anomalies in RSEM and Salmon of GEUVADIS") + 
		theme(axis.title.x=element_blank(), legend.title=element_blank()) + background_grid()

	df <- NULL
	for (i in 1:length(t_h_meta$ID)) {
		id <- t_h_meta$ID[i]
		tmp <- CompareSalmonRSEM(HFolder, id)
		if (is.null(df)) {
			df <- tmp
		} else {
			df <- rbind(df, tmp)
		}
		print(id)
	}

	p2 <- ggplot(df) + geom_boxplot(aes(x=Label, y=Value, color=Type)) + labs(y="percentage among the union of RSEM and Salmon", title="anomalies in RSEM and Salmon of HumanBodyMap") + 
		theme(axis.title.x=element_blank(), legend.title=element_blank()) + background_grid()

	p <- plot_grid(p1, p2, nrow=1, labels="AUTO")
	save_plot(paste0(outdir, "/overlap_salmon_rsem.pdf", sep=""), p, base_aspect_ratio=2.5, base_height=5)


	#################### overlap with identifiability ####################
	t <- NULL

	for (i in 1:length(t_meta$ID)) {
		id <- t_meta$ID[i]
		tmp2 <- read.table(paste0(GFolder, "/star_transcript_Full_", id, "/express_results.xprs", sep=""), header=T, sep="\t")
		# unadjustable
		tmp <- ReadNovelPrediction(GFolder, id)
		x <- tmp2[match(tmp, tmp2$target_id), "solvable"]
		if (is.null(t)) {
			t <- data.frame(ID=id, Dataset="GEUVADIS", Type="unadjustable", Identifiability=1.0*sum(x==TRUE)/length(tmp))
		} else {
			t <- rbind(t, data.frame(ID=id, Dataset="GEUVADIS", Type="unadjustable", Identifiability=1.0*sum(x==TRUE)/length(tmp)) )
		}
		# adjustable
		tmp <- read.table(paste0(GFolder, "/salmon_Full_", id, "/sad/sad_adjusted_quantification.tsv", sep=""), header=F, sep="\t")
		colnames(tmp) <- c("Name", "Length", "NumReads")
		x <- tmp2[match(tmp$Name, tmp2$target_id), "solvable"]
		t <- rbind(t, data.frame(ID=id, Dataset="GEUVADIS", Type="adjustable", Identifiability=1.0*sum(x==TRUE)/nrow(tmp)) )
		print(id)
	}

	for (i in 1:length(t_h_meta$ID)) {
		id <- t_h_meta$ID[i]
		tmp2 <- read.table(paste0(HFolder, "/star_transcript_Full_", id, "/express_results.xprs", sep=""), header=T, sep="\t")
		# unadjustable
		tmp <- ReadNovelPrediction(HFolder, id)
		x <- tmp2[match(tmp, tmp2$target_id), "solvable"]
		if (is.null(t)) {
			t <- data.frame(ID=id, Dataset="HumanBodyMap", Type="unadjustable", Identifiability=1.0*sum(x==TRUE)/length(tmp))
		} else {
			t <- rbind(t, data.frame(ID=id, Dataset="HumanBodyMap", Type="unadjustable", Identifiability=1.0*sum(x==TRUE)/length(tmp)) )
		}
		# adjustable
		tmp <- read.table(paste0(HFolder, "/salmon_Full_", id, "/sad/sad_adjusted_quantification.tsv", sep=""), header=F, sep="\t")
		colnames(tmp) <- c("Name", "Length", "NumReads")
		x <- tmp2[match(tmp$Name, tmp2$target_id), "solvable"]
		t <- rbind(t, data.frame(ID=id, Dataset="HumanBodyMap", Type="adjustable", Identifiability=1.0*sum(x==TRUE)/nrow(tmp)) )
		print(id)
	}

	p <- ggplot(t) + geom_boxplot(aes(x=Dataset, y=Identifiability, color=Type)) + background_grid() + scale_y_continuous(limits=c(0.5,1)) + 
		theme(legend.title=element_blank(), axis.title.x=element_blank()) + labs(y="percentage identifiable", title="percentage of identifiable transcripts within anomalies")
	save_plot(paste0(outdir, "/Identifiable.pdf", sep=""), p, base_aspect_ratio = 1.7)
}
