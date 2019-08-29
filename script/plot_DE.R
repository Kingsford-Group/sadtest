library(ggplot2)
library(cowplot)
library(DESeq2)
library(readr)
library(tximport)
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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	print("Rscript plot_DE.R <prepare data directory> <output figure directory>")
} else {
	indir <- args[1]
	outdir <- args[2]

	GFolder <- paste0(indir, "/GEUVADIS", sep="")
	t_meta <- ReadMeta("/home/cong/Documents/SADnewresult/GEUVADIS/Metadata.txt")
	# get the sequencing centers
	C <- unique(t_meta$Center)

	# DESeq2 to infer DE transcripts
	# salmon
	files <- paste0(GFolder, "/salmon_Full_", t_meta$ID, "/quant.sf", sep="")
	names(files) <- t_meta$ID
	txi <- tximport(files, type = "salmon", txOut = TRUE)
	dds <- DESeqDataSetFromTximport(txi, t_meta, ~Center)
	dds <- DESeq(dds)

	# LP
	txi2 <- txi
	for (id in t_meta$ID) {
		t_quant <- ReadQuantification(GFolder, id)
		txi2$counts[match(t_quant$Name, rownames(txi2$counts)), id] <- t_quant$NumReads_LP
		TPM_LP <- t_quant$NumReads_LP / t_quant$EffectiveLength
		txi2$abundance[match(t_quant$Name, rownames(txi2$abundance)), id] <- 1e6 * TPM_LP / sum(TPM_LP)
	}
	dds2 <- DESeqDataSetFromTximport(txi2, t_meta, ~Center)
	dds2 <- DESeq(dds2)

	# number of predictions
	df <- data.frame( FDR=rep(c(0.01, 0.05, 0.1), 2), Method=c(rep("Salmon", 3), rep("SAD-adjusted", 3)), NumDE=rep(0,6) )
	for (i in seq(1,nrow(df))) {
		if (df[i, "Method"] == "Salmon") {
			res <- results(dds, alpha = df[i, "FDR"])
			df[i, "NumDE"] <- sum(res$padj < df[i, "FDR"], na.rm=TRUE)
		} else {
			res2 <- results(dds2, alpha = df[i, "FDR"])
			df[i, "NumDE"] <- sum(res2$padj < df[i, "FDR"], na.rm=TRUE)
		}
	}
	# write table
	write.table(df, paste0(outdir, "num_de_count.txt", sep=""), append=F, quote=F, sep="\t", row.names=F, col.names=T)


	# reset the result to be under FDR 0.01
	res <- results(dds, alpha=0.01)
	res2 <- results(dds2, alpha=0.01)

	# get the transcripts with larger p-value
	bettertrans <- setdiff(rownames(res[!is.na(res$padj) & res$padj < 0.01, ]), rownames(res2[!is.na(res2$padj) & res2$padj < 0.01, ]))
	df <- res[rownames(res) %in% bettertrans, c("baseMean", "log2FoldChange", "lfcSE")]
	df$adj_baseMean <- res2[match(bettertrans, rownames(res2)), "baseMean"]
	df$adj_log2FoldChange <- res2[match(bettertrans, rownames(res2)), "log2FoldChange"]
	df$adj_lfcSE <- res2[match(bettertrans, rownames(res2)), "lfcSE"]

	p1 <- ggplot(as.data.frame(df)) + geom_point(aes(x=abs(log2FoldChange), y=abs(adj_log2FoldChange)), alpha=0.3) + geom_abline(intercept=0, slope=1, color="red") + scale_x_continuous(limits=c(0,12)) + 
		background_grid() + labs(title="log2 fold change of the reduced DE transcripts", x="Salmon", y = "SAD-adjusted")
	p2 <- ggplot(as.data.frame(df)) + geom_point(aes(x=lfcSE, y=adj_lfcSE), alpha=0.3) + scale_x_continuous(breaks=seq(0,3,0.4)) + scale_y_continuous(breaks=seq(0,3,0.4)) + 
		geom_abline(intercept=0, slope=1, color="red") + background_grid() + labs(title="standard error of the reduced DE transcripts", x="Salmon", y="SAD-adjusted")

	# ploting for the example DE transcripts
	# candidate ENST00000519065.5 HDAC2
	# candidate ENST00000409764.5 CYCS
	# candidate ENST00000395377.2 TFAM
	# candidate ENST00000428459.6 NDUFA13
	df <- data.frame(Name="ENST00000519065.5", ID=colnames(txi$abundance), Abundance=txi$abundance[rownames(txi$abundance) == "ENST00000519065.5"], Method="Salmon", Center=t_meta$Center)
	df <- rbind(df, data.frame(Name="ENST00000519065.5", ID=colnames(txi2$abundance), Abundance=txi2$abundance[rownames(txi2$abundance) == "ENST00000519065.5"], Method="SAD-adjusted", Center=t_meta$Center))
	df[df$Center=="UNIGE", "Center"] <- "Center 1"
	df[df$Center=="CNAG", "Center"] <- "Center 2"
	p3 <- ggplot(df) + geom_boxplot(aes(x=Method, y=Abundance, color=Center)) + background_grid() + labs(title="Abundance adjustment of ENST00000519065.5", x="", y="TPM")

	df <- data.frame(ID=colnames(txi$abundance), Abundance=txi$abundance[rownames(txi$abundance) == "ENST00000428459.6"], Method="Salmon", Center=t_meta$Center)
	df <- rbind(df, data.frame(ID=colnames(txi2$abundance), Abundance=txi2$abundance[rownames(txi2$abundance) == "ENST00000428459.6"], Method="SAD-adjusted", Center=t_meta$Center))
	df[df$Center=="UNIGE", "Center"] <- "Center 1"
	df[df$Center=="CNAG", "Center"] <- "Center 2"
	p4 <- ggplot(df) + geom_boxplot(aes(x=Method, y=Abundance, color=Center)) + background_grid() + labs(title="Abundance adjustment of ENST00000428459.6", x="", y="TPM")

	p_tmp <- plot_grid(p1, p2, p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), nrow=2, labels="AUTO")
	p <- plot_grid(p_tmp, get_legend(p3 + theme(legend.position="bottom")), nrow=2, rel_heights=c(1,0.025))
	save_plot(paste0( outdir, "/DE_improved_transcripts.pdf", sep=""), p, base_aspect_ratio=1.3, base_height=8)

	# supplementary figure to show that increase of variance and decrease of log fold change is not a bias from SAD method
	labels <- apply(txi2$counts- txi$counts, 1, function(x) any(x!=0))
	df <- res[labels, c("lfcSE", "log2FoldChange")]
	df$adj_lfcSE <- res2[labels, "lfcSE"]
	df$adj_log2FoldChange <- res2[labels, "log2FoldChange"]

	p1 <- ggplot(as.data.frame(df)) + geom_point(aes(x=abs(log2FoldChange), y=abs(adj_log2FoldChange)), alpha=0.3) + scale_x_continuous(limits=c(0,15)) + 
		scale_y_continuous(limits=c(0,15)) + geom_abline(intercept=0, slope=1, color="red") + background_grid() + labs(title="log2 fold change of the all adjusted transcripts", x="Salmon", y = "SAD-adjusted")
	p2 <- ggplot(as.data.frame(df)) + geom_point(aes(x=lfcSE, y=adj_lfcSE), alpha=0.3) + geom_abline(intercept=0, slope=1, color="red") + 
		background_grid() + labs(title="standard error of the all adjusted transcripts", x="Salmon", y="SAD-adjusted")
	p <- plot_grid(p2, p1, nrow=1, labels="AUTO")
	save_plot(paste0(outdir, "/Supp_DE.png", sep=""), p, base_aspect_ratio=2, base_height=6)

}
