table_of_comparison_of_gene_expr_in_literature <- read.table(file = "Documents/CFS_WBanalysis_Literature_comparison_files/expression_of_genes_WB_literature_comparison.tsv", sep = "\t",stringsAsFactors = FALSE,header = TRUE)

a <- table_of_comparison_of_gene_expr_in_literature

p_val_current <- a$P.Value.of.Log2.Fold.Change.from.Current.Study


library(graphics)

###
idx_for_lab <- which(a[,"Gene"] %in% c("POLR2G","DEFA4"))


df_to_add <- data.frame(Citation_shorthand = rep(NA,nrow(a)))
citations_to_add <- c("Kaushik 2005", "Kerr 2008", "Gow 2009")
df_to_add[idx_for_lab,"Citation_shorthand"] <- citations_to_add

###

### original scatterplots (x and y scales set) see below for z-scored data

current <- a$Log2.Fold.Change.in.Current.Study..Baseline.
citation <- a$Log2.Fold.Change.in.Citation..Baseline.
Gene <- a$Gene

plot(current,citation, pch = 19, xlim = c(-2.5,3.5), ylim = c(-2.5,3.5), xlab = "Fold Changes of Genes in Current Study", ylab = "Fold Changes of Genes in Citations", main = "Differential Expression in Current Study vs. Citations")
## abline(lm(citation~current), col = "red")

text(current[idx_for_lab], citation[idx_for_lab],labels = Gene[idx_for_lab], cex = 0.6, adj = c(-0.1, -1.2))
text(current[idx_for_lab], citation[idx_for_lab], labels = bquote(.(p_val_current[idx_for_lab])), cex = 0.6, adj = c(-0.5, 0))
text(current, citation, labels = df_to_add$Citation_shorthand, cex = 0.6, adj = c(-0.1,1.2))
segments(-1000, -1000, 1000, 1000)

# r <- cor.test(current,citation, method = "spearman", exact = F)$estimate
# r2 <- r^2
# r2 <- formatC(signif(r2,digits = 3),digits=3,format="fg",flag="#")
# 
# p_val <- cor.test(current,citation, method = "spearman", exact = F)$p.value
# p_val <- formatC(signif(p_val,digits = 3),digits=3,format="fg",flag="#")
# 
# mtext(bquote(r^2 == .(r2)),side = 3,adj = 0, cex = 0.9)
# mtext(bquote(p == .(p_val)),side = 3, adj = 1, cex = 0.9)





###z-Scored scatterplots
z_score_current <- scale(a$Log2.Fold.Change.in.Current.Study..Baseline.,center = T, scale = T)
z_score_citation <- scale (a$Log2.Fold.Change.in.Citation..Baseline.,center = T,scale = T)

plot(z_score_current,z_score_citation, pch = 19, xlim = c(-3.5,3.5), ylim = c(-3.5,3.5), xlab = "Z-Scored Fold Changes of Genes in Current Study", ylab = "Z-Scored Fold Changes of Genes in Citations", main = "Differential Expression in Current Study vs. Citations")
abline(lm(z_score_citation~z_score_current), col = "red")

text(z_score_current[idx_for_lab], z_score_citation[idx_for_lab],labels = Gene[idx_for_lab], cex = 0.7, pos = 3, offset = 1.2)
text(z_score_current[idx_for_lab], z_score_citation[idx_for_lab],labels = bquote(.(p_val_current[idx_for_lab])), cex = 0.7, pos = 3)




r <- cor.test(z_score_current,z_score_citation, method = "spearman", exact = F)$estimate
r2 <- r^2
r2 <- formatC(signif(r2,digits = 3),digits=3,format="fg",flag="#")

p_val <- cor.test(z_score_current,z_score_citation, method = "spearman", exact = F)$p.value
p_val <- formatC(signif(p_val,digits = 3),digits=3,format="fg",flag="#")

mtext(bquote(r^2 == .(r2)),side = 3,adj = 0, cex = 0.9)
mtext(bquote(p == .(p_val)),side = 3, adj = 1, cex = 0.9)


cor.test(a$Log2.Fold.Change.in.Current.Study..Baseline.,a$Log2.Fold.Change.in.Citation..Baseline., method = "spearman", exact = F)


cor.test(a$Log2.Fold.Change.in.Current.Study..Baseline.,a$Log2.Fold.Change.in.Citation..Baseline.)