
## Information about normalization of counts, svaseq and limma from:
## https://github.com/ben-laufer/RNA-seq/blob/main/04-limma-voom.R

#' Differential expression analysis of single study
#'
#' Perform differential expression analysis of a single study. If a variable 
#' with subject identifiers is provided a linear mixed model with subject 
#' identifier as random effect is used based on the function
#' \code{\link[variancePartition]{dream}}. Otherwise a linear model as 
#' implemented in \code{\link[limma]{lmFit}} is fitted.
#' 
#' This function writes a file with the following results of the 
#' differential expression analysis for each gene (see 
#' \code{\link[limma]{topTable}} for more information):
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
#' differential expression analysis. Needs to be "counts" for RNASeq data.
#' @param var [character(1)] name of variable of interest (e.g. group 
#' definition); needs to be available in colData() of se.
#' @param covar [character(n)] additional covariates to be used in the linear
#' (mixed) model; need to be available in colData() of se.
#' @param var.id [character(1)] name of variable with subject identifiers; 
#' needs to be available in colData() of se. If given, a linear mixed model 
#' with subject identifers as random effect is fitted.
#' @param var.ref.level [character(1)] name of reference category for variable
#' of interest.
#' @param sv [logical(1)] should surrogate variables be estimated and used as 
#' (additional) covariates? (default: FALSE).
#' @param sv.file [character(1)] name of file for saving surrogate variables
#' (default: NULL so that surrogate variables are not saved). Note that no file
#' will be created if no surrogate variables could be identified.
#' @param standardize [logical(1)] should expression values be standardized (centered
#' and scaled) before analysis? (default: TRUE).
#' @param res.file [character(1)] name of file for saving gene level results.
#' @param BPPARAM [bpparamClass] parameters for parallel evaluation (see
#' \code{\link[BiocParallel]{bpparam}} for more information).
#' 
#' @export
run_diff_expr_analysis <- function(
  se,
  assay,
  var,
  covar = NULL,
  var.id = NULL,
  var.ref.level = NULL,
  sv = FALSE,
  sv.file = NULL,
  standardize = TRUE,
  res.file,
  BPPARAM = BiocParallel::bpparam()) {
  
  if (!is.null(covar)) {
    ## check covariates in sva and voom 
    stop("covariate adjustment not yet implemented for RNASeq!")
  }
  
  ## extract phenotype data
  pheno = as.data.frame(
    SummarizedExperiment::colData(se)[, c(var, covar, var.id), drop = FALSE])
  
  ## remove individuals with missing values
  pheno = stats::na.omit(pheno)
  if (nrow(pheno) < 3) {
    stop("less than 3 observations for analysis!")
  }
  
  ## set reference level
  if (!is.null(var.ref.level)) {
    if (length(unique(pheno[, var])) > 5) {
      warning(paste("variable", var, "has more than 5 different values\n"))
    }
    pheno[, var] = stats::relevel(factor(pheno[, var]),
                                  ref = var.ref.level)
  }
  
  ## extract expression data
  if (!(assay %in% names(SummarizedExperiment::assays(se)))) {
    stop(paste(assay, "not available in assays!"))
  }
  expr = as.matrix(
    SummarizedExperiment::assays(se)[[assay]])[, rownames(pheno)]
  
  ## RNASeq data:
  ## convert to DGEList object and calculate normalization factors
  if (assay == "counts") {
    expr = edgeR::DGEList(counts = expr)
    expr = edgeR:: calcNormFactors(
      object = expr)
  }

  ## surrogate variables
  if (sv) {
    info.sv = estimate_surrogate_var(
      expr = expr,
      pheno = pheno,
      var = var)
    n.sv = info.sv$n.sv
    if (n.sv > 0) {
       pheno = info.sv$pheno
       covar = c(covar, 
                 paste0("sv_", 1:n.sv))
    }
  } else {
    n.sv = 0
  }

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
  
  ## voom transformation for RNASeq data
  if (assay == "counts") {
    mod = stats::model.matrix(
      stats::as.formula(form),
      data = pheno)
    expr = limma::voom(
      counts = expr,
      design = mod);
  }
  
  ## standardize per gene
  if (standardize) {
    if (assay == "counts") { ## RNASeq
      expr$E = t(apply(expr$E, 1, function(x) {
        (x - mean(x, nar.rm = TRUE)) / stats::sd(x, na.rm = TRUE)}))
    } else {
      expr = t(apply(expr, 1, function(x) {
        (x - mean(x, nar.rm = TRUE)) / stats::sd(x, na.rm = TRUE)}))
    }
  }

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
      data = pheno, 
      BPPARAM = BPPARAM)
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
    n.sv = n.sv,
    stringsAsFactors = FALSE)
  
  ## add information about genes
  res = merge(
    res,
    as.data.frame(SummarizedExperiment::rowData(se)),
    by.x = "gene",
    by.y = "row.names")

  ## save gene level results in file
  rio::export(res,
              file = res.file)
  
  ## save surrogate variables in file
  if (n.sv > 0 & !is.null(sv.file)) {
    rio::export(pheno[, grep("sv_", colnames(pheno)), drop = FALSE],
                file = sv.file,
                row.names = TRUE)
  }

}

#' internal function
#' @keywords internal
estimate_surrogate_var <- function(
    expr,
    pheno,
    var) {
  
  ## model with variable of interest
  mod = stats::model.matrix(
    stats::as.formula(paste0("~",  var)),
    data = pheno)
  
  ## null model
  mod0 = stats::model.matrix(~1,
                      data = pheno) 
  
  ## estimate normalized counts for RNASeq data
  if (methods::is(expr, "DGEList")) {
    expr = edgeR::cpm(expr, log = FALSE)
  }

  ## estimate surrogate variables
  n.sv = sva::num.sv(
    dat = expr,
    mod = mod,
    method = "leek")
  
  if (n.sv > 0) {
    
    if (n.sv > 5) {
      warning(
        paste("number of estimated surrogate variables",
              n.sv,
              "is restricted to 5!"))
      n.sv = 5
    }
    if (methods::is(expr, "DGEList")) { # RNASeq
      svobj = sva::svaseq(
        dat = expr,
        mod = mod,
        mod0 = mod0,
        n.sv = n.sv)
    } else {
      svobj = sva::sva(
        dat = expr,
        mod = mod,
        mod0 = mod0,
        n.sv = n.sv)
    }
    sv = svobj$sv
    colnames(sv) = paste0("sv_", 1:ncol(sv))
    
    if (length(intersect(colnames(pheno), colnames(sv))) > 0) {
      stop("surrogate variables already available in colData")
    }
    
    pheno = data.frame(pheno, sv)
  } else {
    warning("number of estimated surrogate variables is zero")
  }
  return(list(n.sv = n.sv,
              pheno = pheno))
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
