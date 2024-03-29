
#' Meta analysis of multiple studies
#'
#' Perform fixed and random effect meta analysis for each gene across multiple
#' studies using the \code{\link[metafor]{rma.uni}} function.
#' 
#' This function writes a file with the following results of the meta analysis 
#' for each gene (see \code{\link[metafor]{rma.uni}} for 
#' more information):
#' \itemize{
#' \item gene = gene name
#' \item number.studies = number of studies with differential expression 
#' analysis results
#' \item estimate.fixed = estimated coefficient of the fixed effect meta 
#' analysis
#' \item ci.l.fixed = lower bound of confidence interval for the coefficient of
#' the fixed effect meta analysis
#' \item ci.u.fixed = upper bound of confidence interval for the coefficient of
#' the fixed effect meta analysis
#' \item z.fixed = test statistic of the fixed effect meta analysis
#' \item pvalue.fixed = P-value of the fixed effect meta analysis
#' \item estimate.random = estimated coefficient of the random effect meta 
#' analysis
#' \item ci.l.random = lower bound of confidence interval for the coefficient of
#' the random effect meta analysis
#' \item ci.u.random = upper bound of confidence interval for the coefficient of
#' the random effect meta analysis
#' \item z.random = test statistic of the random effect meta analysis
#' \item pvalue.random = P-value of the random effect meta analysis
#' \item tau2 = estimated amount of (residual) heterogeneity
#' \item I2 = value of \eqn{I^2}
#' \item Q = test statistic for the test of (residual) heterogeneity
#' \item pvalue.Q = P-value for the test of (residual) heterogeneity
#' \item freq.estimate.pos = realtive frequency of studies with positive beta
#' coefficients
#' \item pvalue.fixed.adj = P-value of the fixed effect meta analysis after 
#' adjustment for multiple testing using the Benjamini-Hochberg approach
#' \item pvalue.random.adj = P-value of the random effect meta analysis after 
#' adjustment for multiple testing using the Benjamini-Hochberg approach
#' }
#' 
#' @param res.studies [list] results of study specific differential expression 
#' analyses (e.g. loaded by \code{\link{load_study_results}}).
#' @param min.no.studies [numeric(1)] minimal number of studies for which gene
#' level results need to be available in order to be included in the meta
#' analysis.
#' @param fixed [logical(1)] should fixed effect meta analysis be performed? 
#' (default: TRUE).
#' @param random [logical(1)] should random effect meta analysis be performed? 
#' (default: TRUE).
#' @param res.file [character(1)] name of file for saving results.
#' @param anno [data.frame] data.frame including annotation information 
#' with gene identifiers as row names (default: NULL).
#' 
#' @export
run_meta_analysis <- function(res.studies,
                              min.no.studies = 3,
                              fixed = TRUE,
                              random = TRUE,
                              res.file,
                              anno = NULL) {

  if (!fixed & !random) {
    stop("at least one of fixed or random needs to be set to TRUE")
  }
  
  ## extract genes to be used
  info.gene.names = unlist(lapply(res.studies, function(x) {x$gene}))
  tab.name = table(info.gene.names)
  genes = sort(names(tab.name)[tab.name >= min.no.studies])

  ## gene level meta analysis
  res.all = NULL
  for (g in genes) {
    dat = stats::na.omit(data.frame(
      estimate = sapply(res.studies, function(x) {x[g, "estimate"]}),
      se = sapply(res.studies, function(x) {x[g, "se"]})))

    if (fixed) {
      res.fixed = metafor::rma.uni(
        yi = dat$estimate,
        sei = dat$se,
        method = "FE")
    } else {
      res.fixed = NULL
    }
    
    if (random) {
      res.random = metafor::rma.uni(
        yi = dat$estimate,
        sei = dat$se,
        method = "DL")
    } else {
      res.random = NULL
    }
    
    res.g = c(
      number.studies = res.fixed$k,
      estimate.fixed = res.fixed$beta,
      ci.l.fixed = res.fixed$ci.lb,
      ci.u.fixed = res.fixed$ci.ub,
      z.fixed = res.fixed$zval,
      pvalue.fixed = res.fixed$pval,
      estimate.random = res.random$beta,
      ci.l.random = res.random$ci.lb,
      ci.u.random = res.random$ci.ub,
      z.random = res.random$zval,
      pvalue.random = res.random$pval,
      tau2 = res.random$tau2,
      I2 = res.random$I2,
      Q = res.random$QE,
      pvalue.Q = res.random$QEp,
      freq.estimate.pos = sum(dat$estimate > 0) / res.fixed$k)
    res.all = rbind(res.all, res.g)
  }
  res.all = data.frame(gene = genes,
                       res.all,
                       stringsAsFactors = FALSE)

  ## adjust P-values
  if (fixed) {
    res.all$pvalue.fixed.adj = stats::p.adjust(
      res.all$pvalue.fixed,
      method = "BH")
  }
  if (random) {
    res.all$pvalue.random.adj = stats::p.adjust(
      res.all$pvalue.random,
      method = "BH")
  }
  
  ## add annotation information
  if (!is.null(anno)) {
    if (length(intersect(rownames(anno), res.all$gene)) == 0) {
      stop(paste("gene identifiers in annotation information and study results",
                 "do not match"))
    }
    if (length(setdiff(res.all$gene, rownames(anno))) > 0) {
      warning(paste("some gene identifiers not available in annotation",
                    "information"))
    }
    
    res.all = merge(
      res.all, 
      anno,
      all.x = TRUE, 
      by.x = "gene",
      by.y = "row.names")
  }
  
  rio::export(res.all,
              file = res.file)
}
