#Code for benchmarking package against WGCNA
#@author Annice


#Load relevant libraries

library(WGCNA)
library(ggplot2)
library(pheatmap)
library(SurprisalAnalysis)
library(data.table)
library(bapred)
library(tidyverse)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#Step I: Apply Surprisal Analysis 

#Read data for surprisal analysis
expr_df <- read.csv("C:\\Users\\annic\\Downloads\\helper_T_cell_0_test.csv")


sa <- surprisal_analysis(expr_df)
holds <- sa[[1]]        
alph  <- sa[[2]]      






#Step II: Apply WGCNA
#substeps:
#'
#'Consider the entire data, 
#'Extract differentially expressed genes through GeoR
#'Perform batch adjustment
#'Separate out the particular sample and run WGCNA on it
#'

#Batch adjustment is not necessary for Surprisal analysis since we log normalize the data
#AND if batch effects happen we would have a constraint related to the batch
#For WGCNA, however, it is necessary

#Load data first

##Set paths
gse_id       <- "GSE200250"
out_dir      <- "C:/Users/annic/Downloads/test_data"
annot_path   <- "C:/Users/annic/Downloads/MouseWG-6_V2_0_R3_11278593_A.csv"
deg_path     <- "C:/Users/annic/Downloads/GSE200250.top.table.csv"



#Download GSE files
gse_obj <- getGEO(gse_id, destdir = out_dir)
gse_obj <- if (is.list(gse_obj)) gse_obj[[1]] else gse_obj
test.data <- pData(gse_obj)  


gsm_ids   <- test.data$geo_accession
gsm_tbls  <- lapply(gsm_ids, function(id) GEOquery::Table(getGEO(id)))
names(gsm_tbls) <- gsm_ids


invisible(mapply(function(tb, id) {
  write.csv(tb, file.path(out_dir, sprintf("%s.csv", id)), row.names = FALSE)
}, gsm_tbls, gsm_ids))



annotation <- read.csv(annot_path, stringsAsFactors = FALSE)
probe_ids_last <- gsm_tbls[[length(gsm_tbls)]]$ID_REF

gene_by_probe <- annotation$ILMN_Gene[ match(probe_ids_last, annotation$Probe_Id) ]
tx_by_probe   <- annotation$Transcript[ match(probe_ids_last, annotation$Probe_Id) ]
gene_by_probe[is.na(gene_by_probe)] <- tx_by_probe[is.na(gene_by_probe)]


keep_idx <- !duplicated(gene_by_probe)
my.rows  <- gene_by_probe[keep_idx]


sample_titles <- test.data$title
df_mat <- matrix(NA_real_, nrow = length(my.rows), ncol = length(sample_titles),
                 dimnames = list(my.rows, sample_titles))


for (i in seq_along(gsm_ids)) {
  
  tb <- gsm_tbls[[i]]

  g <- annotation$ILMN_Gene[ match(tb$ID_REF, annotation$Probe_Id) ]
  t <- annotation$Transcript[ match(tb$ID_REF, annotation$Probe_Id) ]
  g[is.na(g)] <- t[is.na(g)]
  
  dt  <- data.table(gene = g, value = tb$VALUE)
  agg <- dt[, .(meanmean = mean(value, na.rm = TRUE)), by = gene]
  
  df_mat[ , i] <- agg$meanmean[ match(my.rows, agg$gene) ]
}

df <- as.data.frame(df_mat)


time.points <- vapply(sample_titles, function(s) strsplit(s, "_")[[1]][2], character(1))


ref.deg <- read.csv(deg_path)
ref.deg <- dplyr::filter(ref.deg, P.Value < 0.1)
df <- df[ rownames(df) %in% toupper(ref.deg$Gene.symbol), , drop = FALSE]


t.df <- t(data.frame(time = time.points))
colnames(t.df) <- colnames(df)
df <- rbind(df, t.df)


th.labels <- colnames(df)

#batch adjustment

th.combatted <- combatba( t(data.matrix(df[-nrow(df), , drop = FALSE])),
                          as.factor(time.points) )
adjusted.th.data <- th.combatted$xadj


#Typical WGCNA steps

soft_power_grid <- c(1:10, seq(12, 50, by = 2))
allowWGCNAThreads()

#IMPORTANT - Separate out only the sample we are testing - Th0 replicate 1.

filt.test.dat <- adjusted.th.data[1:9,]


R.power.table <- pickSoftThreshold(filt.test.dat,
                                   powerVector = soft_power_grid,
                                   networkType = "signed")[[2]]

colors <- colorRampPalette(c("#0B2447", "#00809D"))(nrow(R.power.table))
p.n.1 <- ggplot(R.power.table, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(color = colors, stroke = 1.5, shape=8) +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "#0B2447") +
  labs(x = "Power", y = "Scale Free Topology Model Fit, signed R^2", tag = "A") +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(size = 15),
        text = element_text(size = 18))

p.n.2 <- ggplot(R.power.table, aes(Power, mean.k., label = Power)) +
  geom_point(color = colors, stroke = 1.5, shape = 8) +
  geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean Connectivity", tag = "B") +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(size = 15),
        text = element_text(size = 18))

p.n.1 + p.n.2






soft.thresh <- 12
bwmod <- blockwiseModules(filt.test.dat, maxBlockSize = 4780, 
                          TOMType = "signed", power=soft.thresh, 
                          mergeCutHeight = 0.75, numericLabels = FALSE,
                          RandomSeed = 1234, verbose = 3)


mod.eigengenes <- bwmod$MEs
table(bwmod$colors)

plotDendroAndColors(bwmod$dendrograms[[1]], cbind(bwmod$unmergedColors[bwmod$blockGenes[[1]]],
                                                  bwmod$colors[bwmod$blockGenes[[1]]]), c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang=0.03, guideHang=0.05)





module.gene.mapping <- as.data.frame(bwmod$colors)

module.gene.mapping %>% filter(`bwmod$colors`=="blue") %>% rownames()->blue.genes
module.gene.mapping %>% filter(`bwmod$colors`=="turquoise") %>% rownames()->turquoise.genes
module.gene.mapping %>% filter(`bwmod$colors`=="grey") %>% rownames()->grey.genes
module.gene.mapping %>% filter(`bwmod$colors`=="yellow") %>% rownames()->yellow.genes
module.gene.mapping %>% filter(`bwmod$colors`=="brown") %>% rownames()->brown.genes
module.gene.mapping %>% filter(`bwmod$colors`=="green") %>% rownames()->green.genes




entrez_ids <- mapIds(org.Mm.eg.db, keys=str_to_title(blue.genes), column="ENTREZID", keytype="SYMBOL")

GO_results <- enrichGO(gene = entrez_ids[!is.na(entrez_ids)], OrgDb = "org.Mm.eg.db", keyType="ENTREZID", ont = "BP")
fit <- plot(barplot(GO_results, showCategory = 15))

head(GO_results@result, 15)->Go.top.15

ggplot(Go.top.15, aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+scale_fill_gradient(low = "#790915", high = "#062c5c")+theme_minimal()+
  
  theme(
    # Remove panel border
    panel.border=element_blank(),  
    #plot.border = element_blank(),
    # Remove panel grid lines
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black"),
    #axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.text = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    #axis.text = element_text(size = 15),
    
    text = element_text(size=18)
  ) +coord_flip()





percentile_GO <- 0.95 #change based on your preference
lambda_no <- 3 #change based on your preference
GO_analysis_surprisal_analysis(alph, percentile_GO, lambda_no, key_type = "SYMBOL", flip = FALSE, species.db.str =  "org.Mm.eg.db", top_GO_terms=15)

ggplot(Go.top.15, aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+scale_fill_gradient(low = "#790915", high = "#062c5c")+theme_minimal()+
  
  theme(
    # Remove panel border
    panel.border=element_blank(),  
    #plot.border = element_blank(),
    # Remove panel grid lines
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black"),
    #axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.text = element_blank(),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=20),
    #axis.text = element_text(size = 15),
    
    text = element_text(size=18)
  ) +coord_flip()







library(pheatmap)
library(dplyr)


module_colors <- bwmod$colors  
module_sets   <- split(names(module_colors), module_colors)


get_sa_set <- function(G_mat, lambda_idx = 2, top_frac = 0.05, combine_top_bottom = TRUE) {
  g <- G_mat[, lambda_idx]
  if (combine_top_bottom) {
    k <- max(1, ceiling(top_frac * length(g)))
    names(sort(abs(g), decreasing = TRUE))[seq_len(k)]
  } else {
    k <- max(1, ceiling(top_frac * length(g)))
    list(top = names(sort(g, decreasing = TRUE))[seq_len(k)],
         bot = names(sort(g, decreasing = FALSE))[seq_len(k)])
  }
}


sa_sets <- list(
  `SA top30% (lambda_1)` = get_sa_set(alph, lambda_idx = 2, top_frac = 0.1),
  `SA top30% (lambda_2)` = get_sa_set(alph, lambda_idx = 3, top_frac = 0.1),
  `SA top30% (lambda_3)` = get_sa_set(alph, lambda_idx = 4, top_frac = 0.1)
)


universe <- Reduce(intersect, list(rownames(alph), names(module_colors)))
sa_sets_u     <- lapply(sa_sets,     function(s) intersect(s, universe))
module_sets_u <- lapply(module_sets, function(s) intersect(s, universe))


jaccard <- function(a, b) {
  u <- union(a, b); if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}




Jac <- sapply(module_sets_u, function(mod) vapply(sa_sets_u, jaccard, b = mod, numeric(1)))



pheatmap(Jac,
         main = "Jaccard similarity: Surprisal sets vs WGCNA modules",
         color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
         display_numbers = TRUE, number_format = "%.2f", border_color = NA)









