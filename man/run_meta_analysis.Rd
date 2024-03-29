% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_meta_analysis.R
\name{run_meta_analysis}
\alias{run_meta_analysis}
\title{Meta analysis of multiple studies}
\usage{
run_meta_analysis(
  res.studies,
  min.no.studies = 3,
  fixed = TRUE,
  random = TRUE,
  res.file,
  anno = NULL
)
}
\arguments{
\item{res.studies}{[list] results of study specific differential expression 
analyses (e.g. loaded by \code{\link{load_study_results}}).}

\item{min.no.studies}{[numeric(1)] minimal number of studies for which gene
level results need to be available in order to be included in the meta
analysis.}

\item{fixed}{[logical(1)] should fixed effect meta analysis be performed? 
(default: TRUE).}

\item{random}{[logical(1)] should random effect meta analysis be performed? 
(default: TRUE).}

\item{res.file}{[character(1)] name of file for saving results.}

\item{anno}{[data.frame] data.frame including annotation information 
with gene identifiers as row names (default: NULL).}
}
\description{
Perform fixed and random effect meta analysis for each gene across multiple
studies using the \code{\link[metafor]{rma.uni}} function.
}
\details{
This function writes a file with the following results of the meta analysis 
for each gene (see \code{\link[metafor]{rma.uni}} for 
more information):
\itemize{
\item gene = gene name
\item number.studies = number of studies with differential expression 
analysis results
\item estimate.fixed = estimated coefficient of the fixed effect meta 
analysis
\item ci.l.fixed = lower bound of confidence interval for the coefficient of
the fixed effect meta analysis
\item ci.u.fixed = upper bound of confidence interval for the coefficient of
the fixed effect meta analysis
\item z.fixed = test statistic of the fixed effect meta analysis
\item pvalue.fixed = P-value of the fixed effect meta analysis
\item estimate.random = estimated coefficient of the random effect meta 
analysis
\item ci.l.random = lower bound of confidence interval for the coefficient of
the random effect meta analysis
\item ci.u.random = upper bound of confidence interval for the coefficient of
the random effect meta analysis
\item z.random = test statistic of the random effect meta analysis
\item pvalue.random = P-value of the random effect meta analysis
\item tau2 = estimated amount of (residual) heterogeneity
\item I2 = value of \eqn{I^2}
\item Q = test statistic for the test of (residual) heterogeneity
\item pvalue.Q = P-value for the test of (residual) heterogeneity
\item freq.estimate.pos = realtive frequency of studies with positive beta
coefficients
\item pvalue.fixed.adj = P-value of the fixed effect meta analysis after 
adjustment for multiple testing using the Benjamini-Hochberg approach
\item pvalue.random.adj = P-value of the random effect meta analysis after 
adjustment for multiple testing using the Benjamini-Hochberg approach
}
}
