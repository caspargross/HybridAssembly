#!/usr/bin/Rscript

###################################################
#  This scripts serves to analyse assembly wide   #
#  coverage  in order to find plasmids which are  #
#  overrepresented                                #
###################################################

# Load libraries
#suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(Rsamtools)) 
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomicAlignments))
#suppressPackageStartupMessages(library(GenomicFeatures))
#suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(seqinr))

# Create command line options
option_list <- list(
  make_option(c("-g", "--genome"),type="character", default=NULL,
              help="assembled genome in Fasta format' [default %default]",
              dest="genome_filename"),
  make_option(c("-b", "--bam-file"),type="character", default=NULL,
              help="alignment file in .bam format",
              dest="bam_filename"),
  make_option(c("-o", "--out"),type="character", default="analysis_",
              help="Output identifyer [default %default]",
              dest="out")
)


options(error=traceback)
parser <- OptionParser(usage = "%prog -g genome.fasta -b genome_aln.bam -o out", 
                       option_list=option_list,
                       add_help_option = TRUE)
opt = parse_args(parser)


# Testing stuff
#pr <- "/mnt/users/ahgrosc1/data/agpeters/"
#opt$input_filename <- paste0(pr, "genomes/ESBL2135a_all/ESBL2135a_unicycler_final.fasta")
#opt$bam_filename <- paste0(pr, "analysis/ESBL2135a/aligned_reads/ESBL2135a_unicycler_final_unicycler.bam")
#opt$out <- "ESBL2135a_unicycler"

# Read input files
bamFile <- BamFile(opt$bam_filename)
aln <- readGAlignments(bamFile, use.names=TRUE)
aln_dt <- data.table(as.data.frame(readGAlignments(bamFile)))
aln_dt[,readnames := names(aln),]

geno  <- readDNAStringSet(opt$input_filename)
cvg  <- coverage(aln)

# Create data.table with stats
pl  <- aln_dt[, .(nreads = .N, n_alignedbases = sum(qwidth)), by=seqnames]
pl[,coverage := mean(cvg),]
pl[,length := width(geno),]
pl[,gc := letterFrequency(geno, "GC", as.prob=TRUE),]

# Identify outliers:
bxplt <- boxplot(mean(cvg), plot=FALSE)
isOutlier <- unlist(strsplit(names(geno), " "))[3*(1:length(names(geno)))-2] %in% names(bxplt$out)

pl[,outlier := isOutlier,]

# Save outliers to file
writeXStringSet(geno[isOutlier], paste0(opt$out,'_X_coverage_outliers.fasta'))
seqs <- as.character(geno[isOutlier])
seqnames <- paste0(names(geno[isOutlier]), "Outlier = TRUE")
write.fasta(seq, seqnames, paste0(opt$out, '_coverage_outliers.fasta'))
# Save statistics to file
fwrite(pl, file=paste0(opt$out, "_aln_stats.csv"))

# Plot boxplot
pdf(paste0(opt$out, "_cov_boxplot.pdf"))
boxplot(mean(cvg), main=opt$out, ylab="mean coverage per contig")
dev.off()

# Plot length histogram
pdf(paste0(opt$out, "_length_hist.pdf"))
hist(width(geno), breaks=20, main=opt$out, ylab="Contig length frequency")
dev.off()


