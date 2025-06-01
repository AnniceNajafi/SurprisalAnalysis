<h1>SurprisalAnalysis R package guidelines</h1>

<h3>Installation</h3>

To install the R package:

```
install.packages('devtools')
devtools::install_github('AnniceNajafi/SurprisalAnalysis')
```

<h3>Usage</h3>


To use the R package you should follow the steps below:

I. Store gene expression data in a csv file with the first row holding the sample names and the first column holding the gene names.</li>


II. Read the csv file and run the following code:
  <br>

  
```
  input.data <- read.csv('expression_data.csv')
  results <- surprisal_analysis(input.data)
  
```

III. To run GO analysis on the patterns simply use the code below:

```

results[[2]]-> transcript_weights
percentile_GO <- 0.95 #change based on your preference
lambda_no <- 1 #change based on your preference
GO_analysis_surprisal_analysis(transcript_weights, percentile_GO, lambda_no, key_type = "SYMBOL", flip = FALSE, species.db.str =  "org.Hs.eg.db", top_GO_terms=15)

```

<h3>Use GUI from R package</h3>

Simply run the following code:

```
runSurprisalApp()
```


<h3>Web-based application</h3>

A web-based application based on the above has been deployed on <a href = "https://najafiannice.shinyapps.io/surprisal_analysis_app/">this link</a>.








