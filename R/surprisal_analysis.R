#' This function performs surprisal analysis on transcriptomics data
#'
#' @param input.data transcriptomics data stores as dataframe
#'
#' @return a list containing the lambda values and corresponding weights of
#' transcripts stored
#' @export
surprisal_analysis <- function(input.data){

  input.data[input.data==0]<-0.000001 #replace zero values to avoid undefined issues when taking the log later

  input.data->my.mat.filtered

  log(my.mat.filtered)->my.mat.sa #take the log

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


  Y <- as.matrix(log(input.data))

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

  paste("lambda", seq(1, ncol(holds), 1), sep = "_")->holds

  paste("lambda", seq(1, ncol(alph_all), 1), sep = "_")->alph_all

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
#' @param species.db the type of species used for GO analysis, by default set to
#' Homo sapiens, can be either org.Hs.eg.db or org.Mm.eg.db
#' @param top_GO_terms number of GO terms returns, by default set to 15
#'
#' @return the important GO terms related to a lambda gene pattern
#' @importFrom AnnotationDbi
#' @export
GO_analysis_surprisal_analysis <- function(transcript_weights, percentile_GO, lambda_no, key_type = "SYMBOL", flip = FALSE, species.db =  org.Hs.eg.db, top_GO_terms=15){

  transcript_weights->alph_all

  if(flip==TRUE){ #flip lambda values if interested

    -1*alph_all[,lambda_no]->alph_all[,lambda_no]

  }

  percentile_int <- quantile(alph_all[,lambda_no], percentile_GO*0.01)

  #Extract values in X column that are higher than the 95th percentile of alpha_1
  values_above_percentile_int <- toupper(rownames(alph_all))[alph_all[,lambda_no] > percentile_int]



  if(species.db == org.Hs.eg.db){

  species.db.str <- "org.Hs.eg.db"

  }else{

  species.db.str <- "org.Mm.eg.db"

  }


  GO_results <- enrichGO(gene = entrez_ids[!is.na(entrez_ids)], OrgDb = species.db.str, keyType=key_type, ont = "BP")

  head(GO_results@result, top_GO_terms)->Go.top

  return(Go.top)
}

