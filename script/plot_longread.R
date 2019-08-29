library(ggplot2)
library(cowplot)
library(assertthat)

options(stringsAsFactors = FALSE)


args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	print("Rscript plot_longread.R <prepare data directory> <output figure directory>")
} else {
	indir <- args[1]
	outdir <- args[2]

	meta <- read.table(paste0(indir, "/LongRead1000G/sr_RunTable.txt", sep=""), header=T, sep="\t")

	df <- NULL
	for (id in meta$Sample) {
		t <- read.table(paste0(indir, "/LongRead1000G/", id, "/sad_unadjustable_validation_merge.txt", sep=""), header=F, sep="\t")
		colnames(t) <- c("Name", "Status", "LongReadID")
		if (is.null(df)) {
			df <- data.frame(ID=id, Hit=sum(t$Status=="True"), NumPrediction=nrow(t))
		} else {
			df <- rbind( df, data.frame(ID=id, Hit=sum(t$Status=="True"), NumPrediction=nrow(t)) )
		}
	}

	df <- unique(df)
	p <- ggplot(df) + geom_bar(aes(x=ID, y=Hit/NumPrediction), stat="identity", position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + background_grid() + 
		labs(x="sample id", y="precision", title="precision of predicting unannotated isoforms using unadjustable anomalies")
	save_plot(paste0(outdir, "/lr_validation.pdf", sep=""), p, base_aspect_ratio=1.4, base_height=5.5)
}
