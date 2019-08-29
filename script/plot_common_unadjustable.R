library(ggplot2)
library(cowplot)
library(assertthat)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	print("Rscript plot_common_unadjustable.R <output directory of AnalyzeRealData_sharedtrans.py> <output figure directory>")
} else {
	indir <- args[1]
	outdir <- args[2]

	# Analysis: the shared novel-isoform-containing genes, what is the gene lengths? compared to the whole gene length distribution?
	t <- read.table(paste0(indir, "/SharedTrans_Length.txt", sep=""), header=F, sep="\t")
	colnames(t) <- c("TransName", "Length", "Type")
	# t <- t[t$Type == "AllSig" | t$Type == "BGAllExpress", ]
	t[t$Type == "AllSig", "Type"] <- "common anomalies"
	t[t$Type == "BGAllExpress", "Type"] <- "common expressing transcripts"
	t[t$Type == "AnySig" | t$Type == "BGAnyExpress", "Type"] <- "background"
	t$Type <- factor(t$Type, levels=c("common anomalies", "common expressing transcripts", "background"))
	p1 <- ggplot(t) + geom_density(aes(x=log10(Length), color=Type, linetype=Type), size=0.8) + labs(title="Transcript length distribution", x="log10(length)", y="Density") + background_grid() + 
		theme(legend.position="bottom") + theme(legend.title=element_blank())

	# Analysis: what is the start / end proportion of over-expressed region in transcript?
	t <- read.table(paste0(indir, "/SharedTrans_OverExpProportion.txt", sep=""), header=F, sep="\t")
	colnames(t) <- c("TransName", "DataSet", "SampleID", "StartProportion", "EndProportion")
	t$Label <- rep("Over-expression region", nrow(t))
	tmp_t <- read.table(paste0(indir, "/SharedTrans_UnderExpProportion.txt", sep=""), header=F, sep="\t")
	colnames(tmp_t) <- c("TransName", "DataSet", "SampleID", "StartProportion", "EndProportion")
	tmp_t$Label <- rep("Under-expression region", nrow(tmp_t))
	t <- rbind(t, tmp_t)
	p2 <- ggplot(t) + geom_point(aes(x=StartProportion, y=EndProportion, color=Label, shape=Label), size=1.7, alpha=0.3) + scale_x_continuous(limits=c(0, 1)) + 
		scale_y_continuous(limits=c(0, 1)) + labs(title="Position of over- (under-) expression region", x="Region start relative position", y="Region end relative position") + 
		background_grid() + theme(legend.position="bottom", legend.title=element_blank(), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size=16), legend.text = element_text(size=15))
	p2 <- p2 + geom_abline(slope=-1, intercept=1, color="red") + annotate("text", x = 0.4, y = 0.1, label = "Mainly overlap first half", size=5.5) +
		annotate("text", x = 0.8, y = 0.6, label = "Mainly overlap second half", size=5.5)

	# Analysis: exon length distribution of the novel isoform omission
	t <- read.table(paste0(indir, "/SharedTrans_ExonLengths.txt", sep=""), header=F, sep="\t")
	colnames(t) <- c("Name", "Length", "Label")
	t$Label <- as.character(t$Label)
	t[t$Label=="All", "Label"] <- "background"
	t[t$Label=="Contain", "Label"] <- "over-expression region"
	t[t$Label=="Single", "Label"] <- "single-exon under-expression region"
	t[t$Label=="Multiple", "Label"] <- "multi-exons under-expression region"
	t$Label <- factor(t$Label, levels=c("single-exon under-expression region", "over-expression region", "multi-exons under-expression region", "background"))
	p3 <- ggplot(t) + geom_density(aes(x=log10(Length), color=Label, linetype=Label), size=0.8) + labs(title="Exon length distribution", x="log10(exon length)", y="Density") + background_grid() + 
		theme(legend.title=element_blank()) + theme(legend.position="bottom") + guides(color=guide_legend(nrow=2,byrow=TRUE))

	# Analysis: the under-expression region (or deleted region in novel isoform), how many exons does it span? Or hown many exons are left out in the novel one?
	t1 <- read.table(paste0(indir, "/SharedTrans_OverExpProportion.txt", sep=""), header=F, sep="\t")
	colnames(t1) <- c("TransName", "DataSet", "SampleID", "StartProportion", "EndProportion")
	t1$UnderExpRegion <- rep("Second half", nrow(t1))
	t1[0.5-t1$StartProportion < t1$EndProportion-0.5, "UnderExpRegion"] <- "First half"
	t2 <- read.table(paste0(indir, "/SharedTrans_NumSpanningExons.txt", sep=""), header=F, sep="\t")
	colnames(t2) <- c("TransName", "DataSet", "SampleID", "NumSpanningExons")
	t1$NumSpanningExons <- t2$NumSpanningExons

	t <- data.frame(NumSpanningExons=seq(1,max(t1$NumSpanningExons),1), Count=0, UnderExpRegion="First half")
	t <- rbind( t, data.frame(NumSpanningExons=seq(1,max(t1$NumSpanningExons),1), Count=0, UnderExpRegion="Second half") )
	for (i in 1:nrow(t)) {
		t[i, "Count"] <- nrow(t1[t1$NumSpanningExons == t[i, "NumSpanningExons"] & t1$UnderExpRegion == t[i, "UnderExpRegion"], ])
	}
	p4 <- ggplot(t) + geom_bar(aes(x=NumSpanningExons, y=Count, fill=UnderExpRegion), stat="identity", position="dodge") + scale_x_continuous(limits=c(1,5), breaks=seq(1,5,1)) +
		labs(title="Number exons in under-expression region across all samples", x="Number exons", fill="Under-expression region mainly overlaps") + 
		theme(legend.position="bottom", axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), plot.title = element_text(size=16), legend.title = element_text(size=15), legend.text = element_text(size=15)) + 
		background_grid() + guides(fill=guide_legend(nrow=2,byrow=TRUE))

	# save figure
	save_plot(paste0(outdir, "/supp_SharedTrans_Length.pdf", sep=""), p1, base_aspect_ratio=1.3, base_height=5)
	p <- plot_grid(p2, p4,  nrow=1, labels=c("C", "D"), label_size=15)
	save_plot(paste0(outdir, "/SharedTrans_Length_OverExpProportion.pdf", sep=""), p, base_aspect_ratio=2.4, base_height=6.2)
}
