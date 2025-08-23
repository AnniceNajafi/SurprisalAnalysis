#Let us test the normalization methods here
#pseudocount vs log1p

#read.csv("~/Downloads/pone.0195142.s013.csv")->testing.case.1


#' Plot Bland Altman and density plots of the difference between the two zero handling
#' methods
#'
#' @param aa the gene expression matrix
#' @param pc the pseudocount
#'
#' @return two plots, the Bland Altman plot and the density plot
plot.norm.diffs <- function(aa, pc){


  mat_logpc <- log(X + pc)
  mat_log1p <- log1p(X)

  avg_vals  <- (mat_logpc + mat_log1p) / 2
  diff_vals <- mat_logpc - mat_log1p

  df_ba <- data.frame(
    Avg  = as.vector(avg_vals),
    Diff = as.vector(diff_vals)
  )


  p_ba <- ggplot(df_ba, aes(x = Avg, y = Diff)) +
    geom_point(size = 0.8, color = "darkred") +
    geom_hline(yintercept = 0, linetype = 2, color = "navy") +
    labs(
      title = "Bland–Altman Plot",
      subtitle = "log(x+1e-6) vs log1p(x)",
      x = "Average of two transforms",
      y = "Difference (log(x+1e-6) - log1p(x))"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text = element_text(size = 15),
      text = element_text(size = 18)
    )


  df_diff <- data.frame(Diff = as.vector(diff_vals))

  p_density <- ggplot(df_diff, aes(x = Diff)) +
    geom_density(aes(y = after_stat(scaled)), fill = "steelblue", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, color = "navy") +
    labs(
      title = "Distribution of Differences",
      x = "Difference (log(x+1e-6) - log1p(x))",
      y = "Density"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text = element_text(size = 15),
      text = element_text(size = 18)
    )

  return(p_ba + p_density)
}

#pseudocount
pc <- 1e-6

plot.norm.diffs(testing.case.1, pc)
plot.norm.diffs(testing.case.2, pc)


plt.A/plt.B




read.csv("~/Downloads/GSE60450_Lactation-GenewiseCounts.csv")->testing.case.1

rownames(testing.case.1)<- testing.case.1[,1]

testing.case.1[,2:ncol(testing.case.1)]->testing.case.1



df <- read.table("~/Downloads/GSE60450_Lactation-GenewiseCounts (1).txt",
                 header = TRUE, sep = "\t", check.names = FALSE)


stopifnot(all(c("EntrezGeneID","Length") %in% names(df)))


gene_ids <- df$EntrezGeneID
len_bp   <- df$Length
counts   <- as.matrix(df[ , -(1:2)])


keep <- is.finite(len_bp) & (len_bp > 0)
counts <- counts[keep, , drop = FALSE]
len_bp <- len_bp[keep]
gene_ids <- gene_ids[keep]
rownames(counts) <- gene_ids


len_kb <- len_bp / 1000

#FPKM
#FPKM_ij = counts_ij / (len_kb_i * libsize_j_in_millions)
libsize        <- colSums(counts, na.rm = TRUE)
libsize_mill   <- libsize / 1e6

fpkm <- sweep(counts, 1, len_kb, "/")
fpkm <- sweep(fpkm, 2, libsize_mill, "/")

#TPM
# TPM_ij = (counts_ij / len_kb_i) / sum_i(counts_ij / len_kb_i) * 1e6
rpk  <- sweep(counts, 1, len_kb, "/")
scale <- colSums(rpk, na.rm = TRUE)
tpm <- sweep(rpk, 2, scale, "/") * 1e6

testing.case.2 <- fpkm
testing.case.3 <- tpm

plot.norm.diffs(testing.case.1, pc)

plot.norm.diffs(testing.case.2, pc)

plot.norm.diffs(testing.case.3, pc)



baseline <- data.frame(raw = surprisal_analysis(testing.case.1)[[1]][,1],
           fpkm = surprisal_analysis(testing.case.2)[[1]][,1],
           tpm = surprisal_analysis(testing.case.3)[[1]][,1])

baseline.pseudo <-
  ggplot(baseline)+
  geom_line(aes(x=1:nrow(baseline), y=raw, colour="Raw"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=raw, colour="Raw"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=fpkm, color = "FPKM Normalized"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=fpkm, colour = "FPKM Normalized"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), size=1,)+
  geom_point(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), shape=8, stroke=1)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )+
  scale_color_manual(values=c("Raw"="#9ECAD6", "FPKM Normalized" = "#4A9782","TPM Normalized"="#004030" ))+
  labs(x="Sample", title=TeX("$\\lambda_0$"), y="Value",colour="Normalization")




lambda_1 <- data.frame(raw = surprisal_analysis(testing.case.1)[[1]][,2],
                       fpkm = surprisal_analysis(testing.case.2)[[1]][,2],
                       tpm = surprisal_analysis(testing.case.3)[[1]][,2])

lambda.1.pseudo <-
  ggplot(lambda_1)+
  geom_line(aes(x=1:nrow(baseline), y=raw, colour="Raw"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=raw, colour="Raw"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=fpkm, color = "FPKM Normalized"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=fpkm, colour = "FPKM Normalized"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), size=1,)+
  geom_point(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), shape=8, stroke=1)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )+
  scale_color_manual(values=c("Raw"="#9ECAD6", "FPKM Normalized" = "#4A9782","TPM Normalized"="#004030" ))+
  labs(x="Sample", title=TeX("$\\lambda_1$"), y="Value",colour="Normalization")




lambda_2 <- data.frame(raw = surprisal_analysis(testing.case.1)[[1]][,3],
                       fpkm = surprisal_analysis(testing.case.2)[[1]][,3],
                       tpm = surprisal_analysis(testing.case.3)[[1]][,3])

lambda_2.pseudo <-
  ggplot(lambda_2)+
  geom_line(aes(x=1:nrow(baseline), y=raw, colour="Raw"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=raw, colour="Raw"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=fpkm, color = "FPKM Normalized"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=fpkm, colour = "FPKM Normalized"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), size=1,)+
  geom_point(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), shape=8, stroke=1)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )+
  scale_color_manual(values=c("Raw"="#9ECAD6", "FPKM Normalized" = "#4A9782","TPM Normalized"="#004030" ))+
  labs(x="Sample", title=TeX("$\\lambda_2$"), y="Value",colour="Normalization")





baseline <- data.frame(raw = surprisal_analysis(testing.case.1, "log1p")[[1]][,1],
                       fpkm = surprisal_analysis(testing.case.2, "log1p")[[1]][,1],
                       tpm = surprisal_analysis(testing.case.3, "log1p")[[1]][,1])

baseline.log1p <-
  ggplot(baseline)+
  geom_line(aes(x=1:nrow(baseline), y=raw, colour="Raw"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=raw, colour="Raw"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=fpkm, color = "FPKM Normalized"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=fpkm, colour = "FPKM Normalized"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), size=1,)+
  geom_point(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), shape=8, stroke=1)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )+
  scale_color_manual(values=c("Raw"="#9ECAD6", "FPKM Normalized" = "#4A9782","TPM Normalized"="#004030" ))+
  labs(x="Sample", title=TeX("$\\lambda_0$"), y="Value",colour="Normalization")




lambda_1 <- data.frame(raw = surprisal_analysis(testing.case.1, "log1p")[[1]][,2],
                       fpkm = surprisal_analysis(testing.case.2, "log1p")[[1]][,2],
                       tpm = surprisal_analysis(testing.case.3, "log1p")[[1]][,2])

lambda.1.log1p <-
  ggplot(lambda_1)+
  geom_line(aes(x=1:nrow(baseline), y=raw, colour="Raw"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=raw, colour="Raw"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=fpkm, color = "FPKM Normalized"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=fpkm, colour = "FPKM Normalized"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), size=1,)+
  geom_point(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), shape=8, stroke=1)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )+
  scale_color_manual(values=c("Raw"="#9ECAD6", "FPKM Normalized" = "#4A9782","TPM Normalized"="#004030" ))+
  labs(x="Sample", title=TeX("$\\lambda_1$"), y="Value",colour="Normalization")





lambda_2 <- data.frame(raw = surprisal_analysis(testing.case.1, "log1p")[[1]][,3],
                       fpkm = surprisal_analysis(testing.case.2, "log1p")[[1]][,3],
                       tpm = surprisal_analysis(testing.case.3, "log1p")[[1]][,3])

lambda_2.log1p<-
  ggplot(lambda_2)+
  geom_line(aes(x=1:nrow(baseline), y=raw, colour="Raw"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=raw, colour="Raw"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=fpkm, color = "FPKM Normalized"), size=1)+
  geom_point(aes(x=1:nrow(baseline), y=fpkm, colour = "FPKM Normalized"), shape=8, stroke=1)+
  geom_line(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), size=1,)+
  geom_point(aes(x=1:nrow(baseline), y=tpm, colour="TPM Normalized"), shape=8, stroke=1)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )+
  scale_color_manual(values=c("Raw"="#9ECAD6", "FPKM Normalized" = "#4A9782","TPM Normalized"="#004030" ))+
  labs(x="Sample", title=TeX("$\\lambda_2$"), y="Value",colour="Normalization")




(baseline.pseudo+baseline.log1p)/(lambda.1.pseudo+lambda.1.log1p)/(lambda_2.pseudo+lambda_2.log1p)






running.test.case.1<-data.frame(Sample = as.numeric(rownames(testing.case.1)), testing.case.1)

running.test.case.2<-data.frame(Sample = as.numeric(rownames(testing.case.2)), testing.case.2)

running.test.case.3<-data.frame(Sample = as.numeric(rownames(testing.case.3)), testing.case.3)

testing.case.2.log1p <- GO_analysis_surprisal_analysis(surprisal_analysis(running.test.case.2, "log1p")[[2]], 0.85, 2, key_type = "ENTREZID", flip = FALSE, species.db.str =  "org.Mm.eg.db", top_GO_terms=15)

testing.case.2.pseudo <- GO_analysis_surprisal_analysis(surprisal_analysis(running.test.case.2, "pseudocount")[[2]], 0.85, 2, key_type = "ENTREZID", flip = FALSE, species.db.str =  "org.Mm.eg.db", top_GO_terms=15)





ggplot(testing.case.2.pseudo, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(
    title = "Top 15 Enriched GO Terms",
    x = "GO Term",
    y = "Gene Count",
    fill = "Adj. p-value"
  ) +
  theme_minimal(base_size = 14)





#Jaccard similarity


extreme_pct_lambda <- function(sa_mat, lambda_col = "lambda_1", p = 0.05) {
  stopifnot(lambda_col %in% colnames(sa_mat))
  v <- sa_mat[, lambda_col]
  n <- length(v)
  k <- max(1, ceiling(p * n))


  top_ids <- names(sort(v, decreasing = TRUE))[seq_len(k)]
  bot_ids <- names(sort(v, decreasing = FALSE))[seq_len(k)]

  return(sort(unique(c(top_ids, bot_ids))))
}


set1_log1p <- extreme_pct_lambda(sa1_log1p, "lambda_1", 0.05)
set1_pc    <- extreme_pct_lambda(sa1_pc,    "lambda_1", 0.05)

set2_log1p <- extreme_pct_lambda(sa2_log1p, "lambda_1", 0.05)
set2_pc    <- extreme_pct_lambda(sa2_pc,    "lambda_1", 0.05)

set3_log1p <- extreme_pct_lambda(sa3_log1p, "lambda_1", 0.05)
set3_pc    <- extreme_pct_lambda(sa3_pc,    "lambda_1", 0.05)

sets <- list(
  `Raw-log1p`        = set1_log1p,
  `Raw-pseudocount`  = set1_pc,
  `FPKM-log1p`       = set2_log1p,
  `FPKM-pseudocount` = set2_pc,
  `TPM-log1p`        = set3_log1p,
  `TPM-pseudocount`  = set3_pc
)


pairwise_matrix <- function(sets, fun) {
  n <- length(sets)
  M <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
  for (i in seq_len(n)) for (j in seq_len(n)) M[i, j] <- fun(sets[[i]], sets[[j]])
  M
}

overlap_fun <- function(a, b) length(intersect(a, b))
jaccard_fun <- function(a, b) {
  inter <- length(intersect(a, b)); uni <- length(union(a, b))
  if (uni == 0) return(NA_real_) else inter / uni
}

overlap_mat <- pairwise_matrix(sets, overlap_fun)
diag(overlap_mat) <- vapply(sets, length, integer(1))
jaccard_mat <- pairwise_matrix(sets, jaccard_fun)


pheatmap::pheatmap(
  jaccard_mat,
  main = TeX("Jaccard similarity (Top & Bottom 5% of $\\lambda_0$)"),
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_row = 10, fontsize_col = 10, border_color = NA,
  color = colorRampPalette(c("#1B3C53", "#456882", "#91C8E4"))(100)
)


set1_log1p <- extreme_pct_lambda(sa1_log1p, "lambda_2", 0.05)
set1_pc    <- extreme_pct_lambda(sa1_pc,    "lambda_2", 0.05)

set2_log1p <- extreme_pct_lambda(sa2_log1p, "lambda_2", 0.05)
set2_pc    <- extreme_pct_lambda(sa2_pc,    "lambda_2", 0.05)

set3_log1p <- extreme_pct_lambda(sa3_log1p, "lambda_2", 0.05)
set3_pc    <- extreme_pct_lambda(sa3_pc,    "lambda_2", 0.05)

sets <- list(
  `Raw-log1p`        = set1_log1p,
  `Raw-pseudocount`  = set1_pc,
  `FPKM-log1p`       = set2_log1p,
  `FPKM-pseudocount` = set2_pc,
  `TPM-log1p`        = set3_log1p,
  `TPM-pseudocount`  = set3_pc
)


overlap_mat <- pairwise_matrix(sets, overlap_fun)
diag(overlap_mat) <- vapply(sets, length, integer(1))
jaccard_mat <- pairwise_matrix(sets, jaccard_fun)


pheatmap::pheatmap(
  jaccard_mat,
  main = TeX("Jaccard similarity (Top & Bottom 5% of $\\lambda_1$)"),
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_row = 10, fontsize_col = 10, border_color = NA,
  color = colorRampPalette(c("#1B3C53", "#456882", "#91C8E4"))(100)
)



set1_log1p <- extreme_pct_lambda(sa1_log1p, "lambda_3", 0.05)
set1_pc    <- extreme_pct_lambda(sa1_pc,    "lambda_3", 0.05)

set2_log1p <- extreme_pct_lambda(sa2_log1p, "lambda_3", 0.05)
set2_pc    <- extreme_pct_lambda(sa2_pc,    "lambda_3", 0.05)

set3_log1p <- extreme_pct_lambda(sa3_log1p, "lambda_3", 0.05)
set3_pc    <- extreme_pct_lambda(sa3_pc,    "lambda_3", 0.05)

sets <- list(
  `Raw-log1p`        = set1_log1p,
  `Raw-pseudocount`  = set1_pc,
  `FPKM-log1p`       = set2_log1p,
  `FPKM-pseudocount` = set2_pc,
  `TPM-log1p`        = set3_log1p,
  `TPM-pseudocount`  = set3_pc
)


overlap_mat <- pairwise_matrix(sets, overlap_fun)
diag(overlap_mat) <- vapply(sets, length, integer(1))
jaccard_mat <- pairwise_matrix(sets, jaccard_fun)


pheatmap::pheatmap(
  jaccard_mat,
  main = TeX("Jaccard similarity (Top & Bottom 5% of $\\lambda_2$)"),
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_row = 10, fontsize_col = 10, border_color = NA,
  color = colorRampPalette(c("#1B3C53", "#456882", "#91C8E4"))(100)
)



































