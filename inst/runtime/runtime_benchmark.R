# Load relevant libraries
library(SurprisalAnalysis)
library(ggplot2)
library(peakRAM)


#setting seed to ensure reproducibility
set.seed(123)

#params
n_genes       <- 1000
sample_counts <- c(5, 10, 20, 50, 100)
n_reps        <- 100


#' Function records run time and RAM peak usage
#'
#' @param n_genes number of genes
#' @param n_samples
#'
#' @return run time and RAM peak usage
#' @export
time_ram_single_run <- function(n_genes, n_samples) {
  elapsed <- NA_real_
  pr <- peakRAM({

    expr_matrix <- matrix(
      rlnorm(n = n_genes * n_samples, meanlog = 2, sdlog = 1),
      nrow = n_genes,
      ncol = n_samples
    )
    expr_df <- as.data.frame(expr_matrix)
    rownames(expr_df) <- paste0("Gene", seq_len(n_genes))
    colnames(expr_df) <- paste0("Sample", seq_len(n_samples))

    t0 <- Sys.time()
    res <- surprisal_analysis(expr_df)
    t1 <- Sys.time()
    elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  })
  list(
    time_sec = elapsed,
    peak_ram_mib = pr$Peak_RAM_Used_MiB[1]
  )
}


#store info

results_df <- data.frame(
  n_samples               = integer(length(sample_counts)),
  mean_time               = numeric(length(sample_counts)),
  sd_time                 = numeric(length(sample_counts)),
  ci_lower                = numeric(length(sample_counts)),
  ci_upper                = numeric(length(sample_counts)),
  mean_peak_ram_mib       = numeric(length(sample_counts)),
  sd_peak_ram_mib         = numeric(length(sample_counts)),
  ci_lower_peak_ram_mib   = numeric(length(sample_counts)),
  ci_upper_peak_ram_mib   = numeric(length(sample_counts)),
  stringsAsFactors = FALSE
)


ci_halfwidth_fun <- function(x, n) {
  sdev <- sd(x)
  se   <- sdev / sqrt(n)
  qt(0.975, df = n - 1) * se
}


for (i in seq_along(sample_counts)) {
  s <- sample_counts[i]
  rep_times <- numeric(n_reps)
  rep_ram   <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    gc()
    tr <- time_ram_single_run(n_genes = n_genes, n_samples = s)
    rep_times[r] <- tr$time_sec
    rep_ram[r]   <- tr$peak_ram_mib
  }

  mu_time <- mean(rep_times)
  sd_time <- sd(rep_times)
  hw_time <- ci_halfwidth_fun(rep_times, n_reps)

  mu_ram  <- mean(rep_ram)
  sd_ram  <- sd(rep_ram)
  hw_ram  <- ci_halfwidth_fun(rep_ram, n_reps)

  results_df[i, ] <- list(
    s,
    mu_time, sd_time, mu_time - hw_time, mu_time + hw_time,
    mu_ram,  sd_ram,  mu_ram  - hw_ram,  mu_ram  + hw_ram
  )
  print(i)
  message(sprintf(
    "n_samples = %3d : time = %.4f s (sd=%.4f, 95%% CI=[%.4f, %.4f]);  peak RAM = %.2f MiB (sd=%.2f, 95%% CI=[%.2f, %.2f])",
    s, mu_time, sd_time, mu_time - hw_time, mu_time + hw_time,
    mu_ram, sd_ram, mu_ram - hw_ram, mu_ram + hw_ram
  ))
}

print(results_df)

model_cubic <- lm(mean_time ~ poly(n_samples, 3), data = results_df)
summary(model_cubic)

results_df$pred_cubic <- predict(model_cubic, newdata = results_df)

ggplot(results_df, aes(x = n_samples, y = mean_time)) +
  geom_ribbon(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time),
              fill = "steelblue", alpha = 0.2) +
  geom_point(color = "steelblue4", size = 2, shape=8, stroke=1.5) +
  geom_line(color = "steelblue4") +
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 3),
              se = FALSE,
              color = "firebrick",
              linetype = "dashed",
              size = 1) +
  labs(
    x = "Number of samples",
    y = "Average elapsed time (seconds)",
    title = sprintf("SurprisalAnalysis Runtime vs. Number of Samples\n(gene count = %d; mean ± 95%% CI over %d runs)", n_genes, n_reps), tag="A"
  ) +
  theme_minimal(base_size = 12)->plt.A


is_bytes <- max(results_df$mean_peak_ram_mib, na.rm = TRUE) > 1e5
scale_factor <- if (is_bytes) 1/(1024^2) else 1

use_ci <- all(c("ci_lower_peak_ram_mib", "ci_upper_peak_ram_mib") %in% names(results_df))

df_plot <- within(results_df, {
  mean_ram_plot  <- mean_peak_ram_mib * scale_factor
  lower_ram_plot <- (if (use_ci) ci_lower_peak_ram_mib else mean_peak_ram_mib - sd_peak_ram_mib) * scale_factor
  upper_ram_plot <- (if (use_ci) ci_upper_peak_ram_mib else mean_peak_ram_mib + sd_peak_ram_mib) * scale_factor
})


ggplot(df_plot, aes(x = n_samples, y = mean_ram_plot)) +
  geom_ribbon(aes(ymin = lower_ram_plot, ymax = upper_ram_plot),
              fill = "steelblue", alpha = 0.2) +
  geom_point(color = "steelblue4", size = 2, shape=8, stroke = 1.5) +
  geom_line(color = "steelblue4") +
  scale_x_continuous(breaks = sample_counts) +
  labs(
    title = sprintf("SurprisalAnalysis Peak RAM vs. Number of Samples\n(gene count = %d; mean ± 95%% CI over %d runs)", n_genes, n_reps),
    x = "Number of Samples",
    y = "Peak RAM (MiB)", tag="B"
  ) +
  theme_minimal(base_size = 12)->plt.B



plt.A+plt.B


write.csv(results_df, "~/Downloads/runtime_altering_samples_try_2.csv")

####Altering the number of genes


#params
n_samples_fixed <- 50
gene_counts     <- c(100, 500, 1000, 2000, 5000)
n_reps          <- 100


results_genes_df <- data.frame(
  n_genes                 = integer(length(gene_counts)),
  mean_time               = numeric(length(gene_counts)),
  sd_time                 = numeric(length(gene_counts)),
  ci_lower                = numeric(length(gene_counts)),
  ci_upper                = numeric(length(gene_counts)),
  mean_peak_ram_mib       = numeric(length(gene_counts)),
  sd_peak_ram_mib         = numeric(length(gene_counts)),
  ci_lower_peak_ram_mib   = numeric(length(gene_counts)),
  ci_upper_peak_ram_mib   = numeric(length(gene_counts)),
  stringsAsFactors = FALSE
)

ci_halfwidth_fun <- function(x, n) {
  sdev <- sd(x); se <- sdev / sqrt(n); qt(0.975, df = n - 1) * se
}


for (i in seq_along(gene_counts)) {
  g <- gene_counts[i]
  rep_times <- numeric(n_reps)
  rep_ram   <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    gc()
    tr <- time_ram_single_run(n_genes = g, n_samples = n_samples_fixed)
    rep_times[r] <- tr$time_sec
    rep_ram[r]   <- tr$peak_ram_mib
  }

  mu_time <- mean(rep_times); sd_time <- sd(rep_times); hw_time <- ci_halfwidth_fun(rep_times, n_reps)
  mu_ram  <- mean(rep_ram);   sd_ram  <- sd(rep_ram);  hw_ram  <- ci_halfwidth_fun(rep_ram,   n_reps)

  results_genes_df[i, ] <- list(
    g,
    mu_time, sd_time, mu_time - hw_time, mu_time + hw_time,
    mu_ram,  sd_ram,  mu_ram  - hw_ram,  mu_ram  + hw_ram
  )

  message(sprintf(
    "n_genes = %5d : time = %.4f s (sd=%.4f, 95%% CI=[%.4f, %.4f]);  peak RAM = %.2f MiB (sd=%.2f, 95%% CI=[%.2f, %.2f])",
    g, mu_time, sd_time, mu_time - hw_time, mu_time + hw_time,
    mu_ram, sd_ram, mu_ram - hw_ram, mu_ram + hw_ram
  ))
}

print(results_genes_df)

ggplot(results_genes_df, aes(x = n_genes, y = mean_time)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "steelblue", alpha = 0.2) +
  geom_point(color = "steelblue4", size = 2, shape = 8, stroke = 1.5) +
  geom_line(color = "steelblue4") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", size = 1) +
  scale_x_continuous(breaks = gene_counts) +
  labs(
    x = "Number of genes",
    y = "Average elapsed time (seconds)",
    title = sprintf(
      "SurprisalAnalysis Runtime vs. Number of Genes\n(sample count = %d; mean ± 95%% CI over %d runs)",
      n_samples_fixed, n_reps
    ),
    tag = "C"
  ) +
  theme_minimal(base_size = 12) -> plt.C


is_bytes_g <- max(results_genes_df$mean_peak_ram_mib, na.rm = TRUE) > 1e5
scale_factor_g <- if (is_bytes_g) 1/(1024^2) else 1
use_ci_g <- all(c("ci_lower_peak_ram_mib", "ci_upper_peak_ram_mib") %in% names(results_genes_df))

df_plot_genes <- within(results_genes_df, {
  mean_ram_plot_g  <- mean_peak_ram_mib * scale_factor_g
  lower_ram_plot_g <- (if (use_ci_g) ci_lower_peak_ram_mib else mean_peak_ram_mib - sd_peak_ram_mib) * scale_factor_g
  upper_ram_plot_g <- (if (use_ci_g) ci_upper_peak_ram_mib else mean_peak_ram_mib + sd_peak_ram_mib) * scale_factor_g
})


ggplot(df_plot_genes, aes(x = n_genes, y = mean_ram_plot_g)) +
  geom_ribbon(aes(ymin = lower_ram_plot_g, ymax = upper_ram_plot_g),
              fill = "steelblue", alpha = 0.2) +
  geom_point(color = "steelblue4", size = 2, shape = 8, stroke = 1.5) +
  geom_line(color = "steelblue4") +
  scale_x_continuous(breaks = gene_counts) +
  labs(
    title = sprintf("SurprisalAnalysis Peak RAM vs. Number of Genes\n(sample count = %d; mean ± 95%% CI over %d runs)", n_samples_fixed, n_reps),
    x = "Number of genes",
    y = "Peak RAM (MiB)",
    tag = "D"
  ) +
  theme_minimal(base_size = 12) -> plt.D


plt.C + plt.D


write.csv(results_genes_df, "~/Downloads/runtime_altering_genes_try_2.csv")



(plt.A+plt.B)/(plt.C+plt.D+ylim(c(0, 12)))







