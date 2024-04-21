# SurprisalAnalysis

This is an R package for applying surprisal analysis to microarray and RNAsequencing data.

 Please follow the instructions below to run the R package:

<h4>Install the R package</h4>

```
devtools::install_github("AnniceNajafi/SurprisalAnalysis")
```

<h4>Load example dataset</h4>

```
data.th <- read.csv("~/Downloads/helper_T_cell_0_test.csv", header=TRUE)

make.names(data.th$Time, unique=TRUE) ->rownames(data.th)
data.th$Time <- NULL
as.numeric(sub("^X", "", colnames(data.th)))->colnames(data.th)

```
<h4>Apply surprisal analysis to data</h4>

```
surprisal_analysis(data.th)->res
res[[1]] ->lambdas_values
res[[2]] ->transcript_weights
```

<h4>Plot lambdas over time</h4>

```
plot_lambda(lambda_values, 2)
```

