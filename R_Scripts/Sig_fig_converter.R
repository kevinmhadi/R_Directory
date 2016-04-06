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

###

sig_vector <- "4.77760892	2.008514413	6.770168028	8.768801448	6.351368712	2.892322756	5.849206181	3.232798611	3.68911472	4.835226632	4.409582733	8.637660758	5.28191896	7.310961106	4.537822808	6.83679154	5.441025736	1.000332923	7.853684796	3.719033183	5.864733191	6.139193109	2.325170546	4.282689213	3.386255659	4.926862738	1.930606893	7.658609009	4.122419307	7.915243836	7.255831732	5.974453212	4.666419514	5.6162333	4.070360058	4.794238917	7.633291566	4.157650919	6.23379946	5.806956892	3.698140812	1.884611628	0.321929485	2.417343469	5.276580266	7.575154627	8.560294151	6.032836414	0.007027486	6.751378558	4.169838854	7.350148482	5.7127924	4.129047858	6.343473311	6.965036984	6.306667026	5.585262628	3.3588926	5.733782713	4.067915783	8.819143143	7.196162146	7.947947329	5.532178585	6.983260621	0.455390733"
library("stringr")

splitstr <- strsplit(str_replace_all(sig_vector, "[\t]", ","), split = ",")
numeric_vec <- as.numeric(unlist(splitstr))
mean(numeric_vec)
