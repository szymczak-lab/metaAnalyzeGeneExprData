% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_qc.R
\name{qc_estimate}
\alias{qc_estimate}
\title{QC on estimates}
\usage{
qc_estimate(res.studies, all = FALSE)
}
\arguments{
\item{res.studies}{[list] study specific
results (as generated by \code{\link{load_study_results}}.}

\item{all}{[logical(1)] should all genes be plotted? If FALSE, only 
genes with nominal P-value < 0.05 are used. (default: FALSE).}
}
\value{
data.frame with gene level statistics for each study:
\itemize{
\item id = unique study identifier
\item median.cor.estimate = median number of correlation of estimates with
other studies
\item median.cor.estimate.outlier = outlier based on
median.cor.estimate (1: yes, 0: no, NA: not calculated due to small
number of studies)
}
}
\description{
Calculate correlation of estimates in each pair of studies
based on overlapping genes.
}
