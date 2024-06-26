% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_diff_expr.R
\name{run_diff_expr_analysis}
\alias{run_diff_expr_analysis}
\title{Differential expression analysis of single study}
\usage{
run_diff_expr_analysis(
  se,
  assay,
  assay.voom.weights = NULL,
  var,
  covar = NULL,
  var.id = NULL,
  var.ref.level = NULL,
  res.file
)
}
\arguments{
\item{se}{[SummarizedExperiment object] study data.}

\item{assay}{[character(1)] name of assay with expression data for the 
differential expression analysis.}

\item{assay.voom.weights}{[character(1)] name of assay with weights from 
\code{\link[limma]{voom}} transformed (might be applied to RNAseq data).}

\item{var}{[character(1)] name of variable of interest (e.g. group 
definition); needs to be available in colData() of se.}

\item{covar}{[character(n)] additional covariates to be used in the linear
(mixed) model; need to be available in colData() of se.}

\item{var.id}{[character(1)] name of variable with subject identifiers; 
needs to be available in colData() of se. If given, a linear mixed model 
with subject identifers as random effect is fitted.}

\item{var.ref.level}{[character(1)] name of reference category for variable
of interest.}

\item{res.file}{[character(1)] name of file for saving results.}
}
\description{
Perform differential expression analysis of a single study. If a variable 
with subject identifiers is provided a linear mixed model with subject 
identifier as random effect is used based on the function
\code{\link[variancePartition]{dream}}. Otherwise a linear model as 
implemented in \code{\link[limma]{lmFit}} is fitted.
}
\details{
This function writes a file with the following results of the differential
expression analysis for each gene (see \code{\link[limma]{topTable}} for 
more information):
\itemize{
\item gene = gene name (corresponding to rownames() of se)
\item logFC = log2 fold change corresponding to beta estimate of var
\item AveExpr = mean expression value
\item t = moderated t statistic corresponding to var
\item P.Value = P-value corresponding to var
\item adj.P.Value = P-value corresponding to var after adjustment for 
multiple testing using the Benjamini-Hochberg approach
\item SE = standard error corresponding to var
\item mean.<group1> mean expression value in group 1 (if applicable)
\item mean.<group2> mean expression value in group 2 (if applicable)
}
+ additional information available in rowData() of se
}
