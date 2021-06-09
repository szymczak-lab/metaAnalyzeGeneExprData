
#' @export
run_diff_expr_analysis_dream <- function(
  se,
  assay,
  assay.voom.weights = NULL,
  var,
  covar = NULL,
  var.id,
  res.file) {

  ## define formulas for linear mixed model
  form = paste0("~", var)
  if (!is.null(covar)) {
    form = paste0(form, "+",
                  paste(covar, collapse = "+"))
  }
  form.random = paste0("(1|", var.id, ")")
  form.final = paste(form, form.random, sep = " + ")

  ## extract expression data
  if (!(assay %in% names(SummarizedExperiment::assays(se)))) {
    stop(paste(assay, "not available in assays!"))
  }
  expr = SummarizedExperiment::assays(se)[[assay]]
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
    SummarizedExperiment::colData(se)[, c(var, covar, var.id)])

  ## differential expression analysis
  fitmm = variancePartition::dream(
    exprObj = expr,
    formula = stats::as.formula(form.final),
    data = pheno)

  ## extract results and adjust for multiple testing
  res = limma::topTable(
    fitmm,
    coef = 2,
    number = nrow(expr),
    adjust.method = "BH",
    sort.by = "none")

  ## calculate SE
  res$SE = res$logFC / res$t

  ## add mean expression per group
  info.mean = estimate_means(gr = pheno[, var],
                             expr = expr)

  ## store results in data.frame
  res = data.frame(
    gene = rownames(se),
    res,
    info.mean,
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE)

  ## save in file
  rio::export(res,
              file = res.file)
  return(res)

}

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
