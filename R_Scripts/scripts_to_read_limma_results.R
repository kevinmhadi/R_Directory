library(reshape2)
library(graphics)

WB_limma_results <- read.table(file = "Documents/Blood_OI_Limma_excl_duplicate_ID.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

WB_raw_cts <- read.table(file = "Documents/WholeBloodCounts.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

NK_limma_results_excl_outliers <- read.table(file = "Documents/NKcell_all_limmaresults.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)
NK_limma_results_incl_outliers <- read.table(file = "Documents/NKcell_all_limmaresults_incl_outliers.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)
NK_limma_results_OI_vs_CTL <- read.table(file = "Documents/NKcell_all_limmaresults_CFS_OI_vs_CTL.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)


NK_raw_cts <- read.table(file = "Documents/POYFLIB-NKCELLcounts.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

limma_results <- NK_limma_results_excl_outliers

raw_cts <- NK_raw_cts


# Whole Blood Analysis

###############

WB_limma_results <- read.table(file = "Documents/Blood_OI_Limma_excl_duplicate_ID.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

WB_raw_cts <- read.table(file = "Documents/WholeBloodCounts.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)
WB_raw_cts$element.type <- NULL
col_to_remove <- grep(pattern =  "70.CFS.A9", x = names(WB_raw_cts), value = F)
WB_raw_cts[,col_to_remove] <- NULL
CTL_col_idx <- grep("CTL", names(WB_raw_cts))
CFS_col_idx <- grep("CFS", names(WB_raw_cts))


par(mfrow = c(1,1))
# define indices of significant genes (rows) from whole blood limma results
sig_gene_idx <- which(WB_limma_results$adj.P.Val< 0.20)
# define genes above significance threshold and above average expression threshold
sig_ave_gene_idx <- which(WB_limma_results$adj.P.Val < 0.20 & WB_limma_results$AveExpr > 0)
# define significant genes from vector of indices
sel_genes <- WB_limma_results$genes[sig_ave_gene_idx]
# plot
plot(x = WB_limma_results$AveExpr, y = WB_limma_results$logFC, col = ifelse(WB_limma_results$genes %in% sel_genes, "red", "black"), pch = ifelse(WB_limma_results$genes %in% sel_genes, 19, 18), cex = ifelse(WB_limma_results$genes %in% sel_genes, 0.6, 0.2), ylab = bquote("log" [2] ~ "fold change"), xlab = bquote("log" [2] ~ "expression"), main = "Whole Blood Results")
segments(x0 = 0, y0 = -1000, x1 = 0, y1 = 1000)


mlt_WB_raw_cts <- melt(WB_raw_cts, id.vars = ("element.id"), variable.name = "sample.id", value.name = "count")
group <- data.frame(group = ifelse(grepl(pattern = "CFS", x = mlt_WB_raw_cts$sample.id),"CFS","CTL"))
mlt_WB_raw_cts <- cbind(mlt_WB_raw_cts,group)
lst_of_limma <- list(WB_limma_results, NK_limma_results_excl_outliers)



sig_gene_idx <- which(WB_limma_results$adj.P.Val<0.20)
sig_genes <- WB_limma_results[sig_gene_idx,]$genes



sig_idx_CTL <- which(mlt_WB_raw_cts$element.id %in% sig_genes & mlt_WB_raw_cts[,"group"]=="CTL")
sig_idx_CFS <- which(mlt_WB_raw_cts$element.id %in% sig_genes & mlt_WB_raw_cts[,"group"]=="CFS")


"CTL"
mean(mlt_WB_raw_cts$count[sig_idx_CTL])
"CFS"
mean(mlt_WB_raw_cts$count[sig_idx_CFS])



# by(mlt_WB_raw_cts$count, mlt_WB_raw_cts$group, FUN = mean)
dct_WB_raw_cts <- dcast(mlt_WB_raw_cts, sample.id + group ~ element.id, value.var = "count")
by(dct_WB_raw_cts[,sig_genes],dct_WB_raw_cts$group,colMeans)

# plot AveExpr vs. avg count for all genes

xmin <- 0
xmax <- 50
ymin <- -4
ymax <- 8
plot(x = wb_sig_gene_count_means, y = WB_limma_results[WB_limma_results$genes %in% sig_genes,"AveExpr"], xlim = c(0, 50), cex = 0.2, pch = 18, yaxt = "n", xlab = "Avg. Count", ylab = "Avg. Expression")
axis(2, at = ymin:ymax)

# Plot average counts of significant genes across all samples in whole blood

par(mfrow = c(1,1))
wb_sig_gene_count_means <- colMeans(dct_WB_raw_cts[,sig_genes])
wb_all_gene_count_means <- colMeans(dct_WB_raw_cts[,3:length(dct_WB_raw_cts)])
breaks <- nclass.FD(wb_sig_gene_count_means)
wb_sig_mean_counts_hist <- hist(wb_sig_gene_count_means, breaks = length(wb_sig_gene_count_means), xlim = c(0,500), xlab = "Average Counts", main = "Ave. Counts, Whole Blood, Significant Genes")
breaks <- nclass.FD(wb_all_gene_count_means)
wb_all_mean_counts_hist <- hist(wb_all_gene_count_means, breaks = breaks, xlab = "Average Counts", main = "Avg. Counts, All Genes")
max(wb_mean_counts_hist$counts)
min(wb_mean_counts_hist$counts)


breaks <- nclass.FD(wb_all_gene_count_means)
wb_all_mean_counts_hist <- hist(wb_all_gene_count_means, breaks = breaks, xlim = c(0,500), xlab = "Average Counts", main = "Avg. Counts, All Genes")



# Plot average counts of significant genes across CFS samples in whole blood
par(mfrow = c(1,1))
wb_sig_gene_count_means <- colMeans(dct_WB_raw_cts[dct_WB_raw_cts$group=="CFS",sig_genes])
wb_all_gene_count_means <- colMeans(dct_WB_raw_cts[dct_WB_raw_cts$group=="CFS",3:length(dct_WB_raw_cts)])
breaks <- nclass.FD(wb_sig_gene_count_means)
wb_sig_mean_counts_hist <- hist(wb_sig_gene_count_means, breaks = length(wb_sig_gene_count_means), xlim = c(0,500), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")

# Plot average counts of significant genes across CTL samples in whole blood

par(mfrow = c(1,1))
wb_sig_gene_count_means <- colMeans(dct_WB_raw_cts[dct_WB_raw_cts$group=="CTL",sig_genes])
wb_all_gene_count_means <- colMeans(dct_WB_raw_cts[dct_WB_raw_cts$group=="CTL",3:length(dct_WB_raw_cts)])
breaks <- nclass.FD(wb_sig_gene_count_means)
wb_sig_mean_counts_hist <- hist(wb_sig_gene_count_means, breaks = length(wb_sig_gene_count_means), xlim = c(0,500), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")


par(mfrow = c(1,1))
# define indices of significant genes (rows) from whole blood limma results to color on plot
sig_gene_idx <- which(WB_limma_results$adj.P.Val<0.20)
# define significant genes from vector of indices
sel_genes <- WB_limma_results$genes[sig_gene_idx]
# plot
plot(x = WB_limma_results$AveExpr, y = WB_limma_results$logFC, col = ifelse(WB_limma_results$genes %in% sel_genes, "red", "black"), pch = ifelse(WB_limma_results$genes %in% sel_genes, 19, 18), cex = ifelse(WB_limma_results$genes %in% sel_genes, 0.4, 0.2), ylab = bquote("log" [2] ~ "fold change"), xlab = bquote("log" [2] ~ "expression"), main = "Whole Blood Results")


#################

# NK Cell Analysis
################


# Reading in relevant files

NK_limma_results_excl_outliers <- read.table(file = "Documents/NKcell_all_limmaresults.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)
NK_limma_results_incl_outliers <- read.table(file = "Documents/NKcell_all_limmaresults_incl_outliers.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)
NK_limma_results_OI_vs_CTL <- read.table(file = "Documents/NKcell_all_limmaresults_CFS_OI_vs_CTL.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

NK_limma_results <- NK_limma_results_OI_vs_CTL

NK_raw_cts <- read.table(file = "Documents/POYFLIB-NKCELLcounts.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)
col_to_remove <- grep(pattern =  "70.CFS.A9", x = names(NK_raw_cts), value = F)
NK_raw_cts[,col_to_remove] <- NULL
NK_raw_cts$element.type <- NULL
CTL_col_idx <- grep("CTL", names(NK_raw_cts))
CFS_col_idx <- grep("CFS", names(NK_raw_cts))

# select significant genes
sig_gene_idx <- which(NK_limma_results$adj.P.Val<0.20)
sig_genes <- NK_limma_results[sig_gene_idx,]$genes


# Reshape data
mlt_NK_raw_cts <- melt(NK_raw_cts, id.vars = ("elementid"), variable.name = "sample.id", value.name = "count")
group <- data.frame(group = ifelse(grepl(pattern = "CFS", x = mlt_NK_raw_cts$sample.id),"CFS","CTL"))
mlt_NK_raw_cts <- cbind(mlt_NK_raw_cts,group)

dct_NK_raw_cts <- dcast(mlt_NK_raw_cts, sample.id + group ~ elementid, value.var = "count")
by(dct_NK_raw_cts[,sig_genes],dct_NK_raw_cts$group,colMeans)


par(mfrow = c(1,1))
# plot average counts across all samples
nk_sig_gene_count_means <- colMeans(dct_NK_raw_cts[,sig_genes])
nk_all_gene_count_means <- colMeans(dct_NK_raw_cts[,3:length(dct_NK_raw_cts)])
breaks <- nclass.FD(wb_sig_gene_count_means)
nk_sig_mean_counts_hist <- hist(nk_sig_gene_count_means, breaks = length(nk_sig_gene_count_means), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")
breaks <- nclass.FD(nk_all_gene_count_means)
nk_all_mean_counts_hist <- hist(nk_all_gene_count_means, breaks = breaks, xlim = c(0,1500), xlab = "Average Counts", main = "Avg. Counts, All Genes")
max(nk_mean_counts_hist$counts)
min(nk_mean_counts_hist$counts)


# plot AveExpr vs. avg count for all genes

xmin <- 0
xmax <- 50
ymin <- -4
ymax <- 8
plot(x = nk_sig_gene_count_means, y = NK_limma_results[NK_limma_results$genes %in% sig_genes,"AveExpr"], xlim = c(0, 50), cex = 0.2, pch = 18, xlab = "Avg. Count", ylab = "Avg. Expression")
axis(2, at = ymin:ymax)


# plot average counts in CFS

nk_sig_gene_count_means <- colMeans(dct_NK_raw_cts[dct_NK_raw_cts$group=="CFS",sig_genes])
nk_all_gene_count_means <- colMeans(dct_NK_raw_cts[dct_NK_raw_cts$group=="CFS",3:length(dct_NK_raw_cts)])
breaks <- nclass.FD(nk_sig_gene_count_means)
nk_sig_mean_counts_hist <- hist(nk_sig_gene_count_means, breaks = length(nk_sig_gene_count_means), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")

# plot average counts in CTL

nk_sig_gene_count_means <- colMeans(dct_NK_raw_cts[dct_NK_raw_cts$group=="CTL",sig_genes])
nk_all_gene_count_means <- colMeans(dct_NK_raw_cts[dct_NK_raw_cts$group=="CTL",3:length(dct_NK_raw_cts)])
breaks <- nclass.FD(nk_sig_gene_count_means)
nk_sig_mean_counts_hist <- hist(nk_sig_gene_count_means, breaks = length(nk_sig_gene_count_means), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")


par(mfrow = c(1,1))
# define indices of significant genes (rows) from NK Cell limma results to color on plot
sig_gene_idx <- which(NK_limma_results$adj.P.Val<0.20)
# define significant genes from vector of indices
sel_genes <- NK_limma_results$genes[sig_gene_idx]
# plot
plot(x = NK_limma_results$AveExpr, y = NK_limma_results$logFC, col = ifelse(NK_limma_results$genes %in% sel_genes, "red", "black"), pch = ifelse(NK_limma_results$genes %in% sel_genes, 19, 18), cex = ifelse(NK_limma_results$genes %in% sel_genes, 0.4, 0.2), ylab = bquote("log" [2] ~ "fold change"), xlab = bquote("log" [2] ~ "expression"))







# end of NK Cell Analysis

###################




# General Analysis and Plotting

##################




col_to_remove <- grep(pattern =  "70.CFS.A9", x = names(raw_cts), value = F)
raw_cts[,col_to_remove] <- NULL
raw_cts$element.type <- NULL
CTL_col_idx <- grep("CTL", names(raw_cts))
CFS_col_idx <- grep("CFS", names(raw_cts))

# select significant genes
sig_gene_idx <- which(limma_results$adj.P.Val<0.20)
sig_genes <- limma_results[sig_gene_idx,]$genes


# Reshape data
mlt_raw_cts <- melt(raw_cts, id.vars = ("element.id"), variable.name = "sample.id", value.name = "count")
group <- data.frame(group = ifelse(grepl(pattern = "CFS", x = mlt_raw_cts$sample.id),"CFS","CTL"))
mlt_raw_cts <- cbind(mlt_raw_cts,group)

dct_raw_cts <- dcast(mlt_raw_cts, sample.id + group ~ element.id, value.var = "count")
by(dct_raw_cts[,sig_genes],dct_raw_cts$group,colMeans)


par(mfrow = c(1,1))
# plot average counts across all samples
sig_gene_count_means <- colMeans(dct_raw_cts[,sig_genes])
all_gene_count_means <- colMeans(dct_raw_cts[,3:length(dct_raw_cts)])
breaks <- nclass.FD(sig_gene_count_means)
sig_mean_counts_hist <- hist(sig_gene_count_means, breaks = length(sig_gene_count_means), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")
breaks <- nclass.FD(all_gene_count_means)
all_mean_counts_hist <- hist(all_gene_count_means, breaks = breaks, xlim = c(0,1500), xlab = "Average Counts", main = "Avg. Counts, All Genes")
max(mean_counts_hist$counts)
min(mean_counts_hist$counts)


# plot AveExpr vs. avg count for all genes

xmin <- 0
xmax <- 50
ymin <- -6
ymax <- 4
plot(x = sig_gene_count_means, y = limma_results[limma_results$genes %in% sig_genes,"AveExpr"], xlim = c(0, 50), cex = 0.2, pch = 18, xlab = "Avg. Count", ylab = "Avg. Expression", yaxt = "n")
axis(2, at = ymin:ymax)

plot(x = sig_gene_count_means, y = limma_results[limma_results$genes %in% sig_genes,"AveExpr"], cex = 0.2, pch = 18, xlab = "Avg. Count", ylab = "Avg. Expression", yaxt = "n")



# plot average counts in CFS

sig_gene_count_means <- colMeans(dct_raw_cts[dct_raw_cts$group=="CFS",sig_genes])
all_gene_count_means <- colMeans(dct_raw_cts[dct_raw_cts$group=="CFS",3:length(dct_raw_cts)])
breaks <- nclass.FD(sig_gene_count_means)
sig_mean_counts_hist <- hist(sig_gene_count_means, breaks = length(sig_gene_count_means), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")

# plot average counts in CTL

sig_gene_count_means <- colMeans(dct_raw_cts[dct_raw_cts$group=="CTL",sig_genes])
all_gene_count_means <- colMeans(dct_raw_cts[dct_raw_cts$group=="CTL",3:length(dct_raw_cts)])
breaks <- nclass.FD(sig_gene_count_means)
sig_mean_counts_hist <- hist(sig_gene_count_means, breaks = length(sig_gene_count_means), xlab = "Average Counts", main = "Ave. Counts, Significant Genes")


par(mfrow = c(1,1))
# define indices of significant genes (rows) from whole blood limma results
sig_gene_idx <- which(limma_results$adj.P.Val< 0.20)
# define genes above significance threshold and above average expression threshold
sig_ave_gene_idx <- which(limma_results$adj.P.Val < 0.20 & limma_results$AveExpr > -1.5)
# define significant genes from vector of indices
sel_genes <- limma_results$genes[sig_ave_gene_idx]
# plot
plot(x = limma_results$AveExpr, y = limma_results$logFC, col = ifelse(limma_results$genes %in% sel_genes, "red", "black"), pch = ifelse(limma_results$genes %in% sel_genes, 19, 18), cex = ifelse(limma_results$genes %in% sel_genes, 0.6, 0.2), ylab = bquote("log" [2] ~ "fold change"), xlab = bquote("log" [2] ~ "expression"))
segments(x0 = -1.5, y0 = -1000, x1 = -1.5, y1 = 1000)


# end of NK Cell Analysis


##################

### Histogram of Average Expression


breaks <- nclass.FD(NK_limma_results_incl_outliers$AveExpr)
hist1 <- hist(NK_limma_results_incl_outliers$AveExpr, breaks = breaks, freq = F, xlab = "Average Expression", main = "Histogram of Average Expression (Limma)")
hist2 <- hist(NK_limma_results_incl_outliers$AveExpr, breaks = breaks, freq = T, xlab = "Average Expression", main = "Histogram of Average Expression (Limma)")


hist1 <- hist(NK_limma_results_incl_outliers$AveExpr, breaks = breaks, freq = F, ylim = c(0,max(hist1$density)), xlab = "Average Expression", main = "Histogram of Average Expression in NK Cells (Limma)")
lines(density(NK_limma_results_incl_outliers$AveExpr))
hist2 <- hist(NK_limma_results_incl_outliers$AveExpr, breaks = breaks, freq = T, ylim = c(0,max(hist2$counts)), xlab = "Average Expression", main = "Histogram of Average Expression in NK Cells (Limma)")

max(hist1$density)
max(hist2$counts)


### Function
histograms_KH <- function (df) {
  
  breaks <- nclass.FD(df$AveExpr)
  hist1 <- hist(df$AveExpr, breaks = breaks, freq = F, xlab = "Average Expression", main = "Histogram of Average Expression (Limma)")
  hist2 <- hist(df$AveExpr, breaks = breaks, freq = T, xlab = "Average Expression", main = "Histogram of Average Expression (Limma)")
  

  hist1 <- hist(df$AveExpr, breaks = breaks, freq = F, ylim = c(0,max(hist1$density)), xlab = "Average Expression", main = "Histogram of Average Expression (Limma)")
  lines(density(df$AveExpr))
  hist2 <- hist(df$AveExpr, breaks = breaks, freq = T, ylim = c(0,max(hist2$counts)), xlab = "Average Expression", main = "Histogram of Average Expression (Limma)")
  
}

par(mfrow=c(1,2))
histograms_KH(WB_limma_results)
histograms_KH(NK_limma_results)


par(mfrow=c(length(lst_of_limma),2))
lapply(lst_of_limma, histograms_KH)


