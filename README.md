# Surprisal Analysis

This is an R package for applying surprisal analysis to microarray and RNAsequencing data.

 Please follow the instructions below to run the R package:

<h4><strong>I.</strong> Install the R package</h4>

```
devtools::install_github("AnniceNajafi/SurprisalAnalysis")
```

<h4><strong>II.</strong> Load example dataset</h4>

Example data set can be downloaded from <a href="https://drive.google.com/file/d/1exoPw_Cnn_vNJACea68oSMTJ4Fg7DNN3/view?usp=drive_link">this link</a>.
```
data.th <- read.csv("~/Downloads/helper_T_cell_0_test.csv", header=TRUE)

make.names(data.th$Time, unique=TRUE) ->rownames(data.th)
data.th$Time <- NULL
as.numeric(sub("^X", "", colnames(data.th)))->colnames(data.th)

```
<h4><strong>III.</strong> Apply surprisal analysis to data</h4>

```
surprisal_analysis(data.th)->res
res[[1]] ->lambda_values
res[[2]] ->transcript_weights
```

<h4><strong>IV.</strong> Plot lambdas over time</h4>

```
plot_lambda(lambda_values, 2, colnames(data.th))
```
<h4><strong>V.</strong> Apply GO analysis on pattern</h4>

```
GO_analysis_surprisal_analysis(transcript_weights, 95, 2)
```



