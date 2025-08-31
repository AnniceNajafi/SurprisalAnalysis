#' This function performs surprisal analysis on transcriptomics data
#'
#' @param input.data transcriptomics data stores as dataframe
#' @param zero.handling zero handling method. Can be either 'pseudocount' or
#' 'log1p'. By default it is set to 'pseudocount'
#' @return a list containing two matrix array objects, first one holding the
#' lambda values representing the constraints or Lagrange multipliers and the
#' second one holding the corresponding weights of transcripts stored (G matrix)
#' @importFrom matlib inv
#' @examples
#' expr.df <- data.frame(gene_id = paste0("Gene", 1:6),
#' S1 = c(0, 12, 3, 0, 50, 7),
#' S2 = c(5, 0, 2, 9, 0, 4),
#' S3 = c(8, 15, 0, 1, 25, 0),
#' S4 = c(0, 7, 6, 0, 40, 3),
#' check.names = FALSE)
#' surprisal_analysis(expr.df, zero.handling = "pseudocount")
#'
#' @export surprisal_analysis
surprisal_analysis <- function(input.data, zero.handling = 'pseudocount'){


  #input.data[input.data==0]<-0.000001 #replace zero values to avoid undefined issues when taking the log later
  if(zero.handling == 'pseudocount'){

    my.mat.filtered <- input.data
    my.mat.filtered[,2:ncol(my.mat.filtered)] <- input.data[,2:ncol(input.data)] + 1e-6

    log(my.mat.filtered[,2:ncol(my.mat.filtered)])->my.mat.sa #take the log


  }else if(zero.handling == 'log1p'){

    input.data->my.mat.filtered

    my.mat.sa <- log1p(my.mat.filtered[, 2:ncol(my.mat.filtered)])
  }else{

    message("Error. Please choose a valid zero handling method.")

  }


  hold.mat <- matrix(NA, nrow = ncol(my.mat.sa), ncol = ncol(my.mat.sa)) #create matrix to hold lambda results

  for(i in 1:ncol(my.mat.sa)){

    for(j in 1:ncol(my.mat.sa)){

      hold <-as.vector(my.mat.sa[,i])*as.vector(my.mat.sa[,j])

      hold[is.na(hold)]<-0

      sum(hold)->hold.mat[i, j]

    }
  }

  C <-hold.mat

  # Compute eigenvalue decomposition
  eigen_decomp <- eigen(C)

  P <- eigen_decomp$vectors

  # Eigenvalues vector D
  D <- diag(sqrt(eigen_decomp$values))


  Y <- as.matrix(my.mat.sa)

  #hold all lambda pattern results
  alph_lst<-list()

  for(alph in 1:ncol(Y)){

    (colSums(inv(D))[alph]*Y%*%P[,alph])->G
    #make.names(toupper(annot$Gene.symbol), unique=TRUE) -> rownames(G)

    assign(paste0("alpha_", alph), G)

    alph_lst[[alph]]<-G

  }

  alph_all <- do.call(cbind, alph_lst)

  hold_lst<-list()

  for(u in 1:ncol(Y)){

    hold <- P[,u]*D[u,u]

    hold_lst[[u]]<-hold

  }

  holds <- do.call(cbind, hold_lst)

  #change column names

  paste("lambda", seq(1, ncol(holds), 1), sep = "_")->colnames(holds)

  colnames(holds)->rownames(holds)

  paste("lambda", seq(1, ncol(alph_all), 1), sep = "_")->colnames(alph_all)

  rownames(alph_all)<-input.data[,1]

  return(list(holds, alph_all))
}



#' Perform Gene ontology analysis on a pattern of interest
#'
#' @param transcript_weights a dataframe containing the weight of transcripts in
#' each pattern
#' @param percentile_GO the percentile of transcript to be used for GO analysis,
#' for example 95 will run GO on transcripts in the 95th percentile and above
#' @param lambda_no the lambda pattern the user is interested in analyzing
#' @param key_type type of transcripts which can be either SYMBOL, ENTREZID,
#' ENSEMBL, or PROBEID
#' @param flip a boolean variable which can either be true or false, if it is
#' set to true, the lambda values will be multiplied by -1
#' @param species.db.str the type of species used for GO analysis, by default set to
#' Homo sapiens, can be either 'org.Hs.eg.db' or 'org.Mm.eg.db'
#' @param ont the ontology term for GO enrichment analysis. Can be either "BP", "MF" or "CC".
#' They stand for "Biological Process", "Molecular Function" or "Cellular Component". Set to "BP" by default
#' @param pAdjustMethod multiple testing correction method. Could be one of "BH", "bonferroni",
#' "holm", "hochberg", "hommel", "BY", or "none". The default setting is "BH"
#' @param top_GO_terms number of GO terms returns, by default set to 15
#'
#' @return dataframe, the important GO terms related to a lambda gene pattern
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi mapIds
#' @importFrom clusterProfiler enrichGO
#' @importFrom stats quantile
#' @importFrom utils head
#' @examples
#'
#' csv.path <- system.file(
#'   "extdata", "helper_T_cell_0_test.csv",
#'   package = "SurprisalAnalysis"
#' )
#'
#' expr.df <- utils::read.csv(csv.path, check.names = FALSE)
#'
#' sa.res <- surprisal_analysis(expr.df, zero.handling = "log1p")
#' alph.all <- sa.res[[2]]
#'
#' go_top <- GO_analysis_surprisal_analysis(
#'     transcript_weights = alph.all,
#'     percentile_GO      = 95,
#'     lambda_no          = "lambda_1",
#'     key_type           = "SYMBOL",
#'     flip               = FALSE,
#'     species.db.str     = "org.Hs.eg.db",
#'     ont                = "BP",
#'     pAdjustMethod      = "BH",
#'     top_GO_terms       = 15
#'     )
#'
#' @export GO_analysis_surprisal_analysis
GO_analysis_surprisal_analysis <- function(transcript_weights, percentile_GO, lambda_no, key_type = "SYMBOL", flip = FALSE, species.db.str =  "org.Hs.eg.db",
                                           ont = "BP", pAdjustMethod = "BH", top_GO_terms=15){

  transcript_weights->alph_all

  if(flip==TRUE){ #flip lambda values if interested

    -1*alph_all[,lambda_no]->alph_all[,lambda_no]

  }

  percentile_int <- quantile(alph_all[,lambda_no], percentile_GO*0.01)

  #Extract values in X column that are higher than the 95th percentile of alpha_1
  values_above_percentile_int <- toupper(rownames(alph_all))[alph_all[,lambda_no] > percentile_int]


  if(species.db.str == "org.Hs.eg.db"){

  species.db <- org.Hs.eg.db

  }else if(species.db.str == 'org.Mm.eg.db'){

  species.db <- org.Mm.eg.db

  }else{
    stop('Please enter either "org.Hs.eg.db" or "org.Mm.eg.db".')
  }

  entrez_ids <- tryCatch(mapIds(species.db, keys=values_above_percentile_int,column="ENTREZID",keytype=key_type,multiVals="first"),error=function(e)NULL)

  GO_results <-enrichGO(gene=entrez_ids, OrgDb=species.db,keyType="ENTREZID",ont=ont, pAdjustMethod = pAdjustMethod)

  head(GO_results@result, top_GO_terms)->Go.top

  return(Go.top)
}

