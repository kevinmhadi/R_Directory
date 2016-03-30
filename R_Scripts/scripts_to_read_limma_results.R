NK_limma_results_excl_outliers <- read.table(file = "Documents/NKcell_all_limmaresults.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

NK_limma_results_incl_outliers <- read.table(file = "Documents/NKcell_all_limmaresults_incl_outliers.tsv", sep = "\t", stringsAsFactors = FALSE , header = TRUE)

# colnames(NK_limma_results)
par(mfrow = c(1,2))

hist(NK_limma_results_excl_outliers$AveExpr)
hist(NK_limma_results_incl_outliers$AveExpr)