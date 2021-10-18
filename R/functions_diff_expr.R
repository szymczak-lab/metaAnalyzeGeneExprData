
#' Differential expression analysis of single study
#'
#' Perform differential expression analysis of a single study. If a variable 
#' with subject identifiers is provided a linear mixed model with subject 
#' identifier as random effect is used based on the function
#' \code{\link[variancePartition]{dream}}. Otherwise a linear model as 
#' implemented in \code{\link[limma]{lmFit}} is fitted.
#' 
#' This function writes a file with the following results of the differential
#' expression analysis for each gene (see \code{\link[limma]{topTable}} for 
#' more information):
#' \itemize{
#' \item gene = gene name (corresponding to rownames() of se)
#' \item logFC = log2 fold change corresponding to beta estimate of var
#' \item AveExpr = mean expression value
#' \item t = moderated t statistic corresponding to var
#' \item P.Value = P-value corresponding to var
#' \item adj.P.Value = P-value corresponding to var after adjustment for 
#' multiple testing using the Benjamini-Hochberg approach
#' \item SE = standard error corresponding to var
#' \item mean.<group1> mean expression value in group 1 (if applicable)
#' \item mean.<group2> mean expression value in group 2 (if applicable)
#' }
#' + additional information available in rowData() of se
#' 
#' @param se [SummarizedExperiment object] study data.
#' @param assay [character(1)] name of assay with expression data for the 
#' differential expression analysis.
#' @param assay.voom.weights [character(1)] name of assay with weights from 
#' \code{\link[limma]{voom}} transformed (might be applied to RNAseq data).
#' @param var [character(1)] name of variable of interest (e.g. group 
#' definition); needs to be available in colData() of se.
#' @param covar [character(n)] additional covariates to be used in the linear
#' (mixed) model; need to be available in colData() of se.
#' @param var.id [character(1)] name of variable with subject identifiers; 
#' needs to be available in colData() of se. If given, a linear mixed model 
#' with subject identifers as random effect is fitted.
#' @param res.file [character(1)] name of file for saving results.
#' 
#' @export
run_diff_expr_analysis <- function(
  se,
  assay,
  assay.voom.weights = NULL,
  var,
  covar = NULL,
  var.id = NULL,
  res.file) {
  
  ## define formulas for linear (mixed) model
  form = paste0("~", var)
  if (!is.null(covar)) {
    form = paste0(form, "+",
                  paste(covar, collapse = "+"))
  }
  if (!is.null(var.id)) {
    form.random = paste0("(1|", var.id, ")")
    form = paste(form, form.random, sep = " + ")
  }
  
  ## extract expression data
  if (!(assay %in% names(SummarizedExperiment::assays(se)))) {
    stop(paste(assay, "not available in assays!"))
  }
  expr = as.matrix(SummarizedExperiment::assays(se)[[assay]])
  if (!is.null(assay.voom.weights)) {
    if (!(assay.voom.weights %in% names(SummarizedExperiment::assays(se)))) {
      stop(paste(assay, "not available in assays!"))
    }
    weights = SummarizedExperiment::assays(se)[[assay.voom.weights]]
    expr = methods::new("EList",
                        list(E = expr,
                             weights = weights))
    
  }
  
  ## extract phenotype data
  pheno = as.data.frame(
    SummarizedExperiment::colData(se)[, c(var, covar, var.id), drop = FALSE])
  
  ## differential expression analysis
  # linear model
  if (is.null(var.id)) {
    design = stats::model.matrix(
      stats::as.formula(form),
      data = pheno)
    fit = limma::lmFit(expr[, rownames(design)], design)
    fit = limma::eBayes(fit)
  } else {
    fit = variancePartition::dream(
      exprObj = expr,
      formula = stats::as.formula(form),
      data = pheno)
  }
  
  ## extract results and adjust for multiple testing
  res = limma::topTable(
    fit,
    coef = 2,
    number = nrow(expr),
    adjust.method = "BH",
    sort.by = "none")
  res$B = NULL
  
  ## calculate SE
  res$SE = res$logFC / res$t
  
  ## add mean expression per group
  if (length(unique(pheno[, var])) == 2) {
    res = data.frame(
      res,
      estimate_means(
        gr = pheno[, var],
        expr = expr))
  }
  
  ## store results in data.frame
  res = data.frame(
    gene = rownames(res),
    res,
    stringsAsFactors = FALSE)
  
  ## add information about genes
  res = merge(
    res,
    as.data.frame(SummarizedExperiment::rowData(se)),
    by.x = "gene",
    by.y = "row.names")

  ## save in file
  rio::export(res,
              file = res.file)
  return(res)
  
}

#' internal function
#' @keywords internal
estimate_means <- function(gr, expr) {
  if (is.factor(gr) || is.character(gr) || length(unique(gr)) == 2) {
    if (is.matrix(expr)) {
      x = expr
    } else {
      x = expr$E
    }
    expr.l = split(data.frame(t(x)), gr)
    info.mean = sapply(expr.l, colMeans)
    colnames(info.mean) = paste0("mean.", colnames(info.mean))
    sum = apply(info.mean, 2, function(x) {sum(!(is.nan(x)))})
    info.mean = info.mean[, sum > 0]
  } else {
    info.mean = NULL
  }
  return(info.mean)
}
