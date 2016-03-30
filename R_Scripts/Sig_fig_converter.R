print_vector <- c(1.27, 0.26)


for (i in print_vector) {
  cat(log2(i))
  cat("\n")
}

ctl_baseline <- c(9.02, 2.00, 2.33, 1.23, 5.42, 3.25, 1.13, 4.54)
cfs_baseline <- c(7.79, 1.86, 2.43, 1.41, 7.01, 7.27, 9.98, 4.34)

expr_matrix <- rbind(ctl_baseline, cfs_baseline)
dimnames(expr_matrix) <- list(c("CTL","CFS"),c("ASIC3","P2X4","P2X5","TRPV1","ADRA2A","ADRB1","ADRB2",
                                               "TLR4"))

for (i in 1:ncol(expr_matrix)) {
  foldchange <- expr_matrix[2,i]/expr_matrix[1,i]
  log2_FC <- signif(log2(foldchange), digits = 2)
  cat(paste0(colnames(expr_matrix)[i],"\t",log2_FC,"\n"))
}
### not in log2
sig_vector <- "0.431494871
0.476482682
0.396642494
0.446471752
0.136174868
0.386990347
0.821632974
0.674630106
0.827320713
0.179965725"

library("stringr")
splitstr <- strsplit(str_replace_all(sig_vector, "[\r\n]", ","), split = ",")
sig_fig_vector <- formatC(signif(as.numeric(unlist(splitstr)),digits = 2),digits = 2, format="fg",flag="#")

cat(paste0(sig_fig_vector, "\n"))


### in log2

sig_vector <- "4.1
5.6
3.6
6
3.9
4.8
0.29
4.1
3.8
7.2"

library("stringr")
splitstr <- strsplit(str_replace_all(sig_vector, "[\r\n]", ","), split = ",")
numeric_vec <- as.numeric(unlist(splitstr))
numeric_vec <- log2(numeric_vec)
sig_fig_vector <- formatC(signif(numeric_vec,digits = 2),digits=2,format="fg",flag="#")

cat(paste0(sig_fig_vector, "\n"))


