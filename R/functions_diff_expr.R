
#' @export
run_diff_expr_analysis_dream <- function(
  se,
  assay,
  assay.voom.weights = NULL,
  var,
  covar = NULL,
  var.id,
  res.file) {

  ## define formulas
  form = paste0("~", var)
  if (!is.null(covar)) {
    form = paste0(form, "+",
                  paste(covar, collapse = "+"))
  }
  form.random = paste0("(1|", var.id, ")")
  form.final = paste(form, form.random, sep = " + ")

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

  pheno = as.data.frame(
    SummarizedExperiment::colData(se)[, c(var, covar, var.id)])

  fitmm = variancePartition::dream(
    exprObj = expr,
    formula = stats::as.formula(form.final),
    data = pheno)
  res = limma::topTable(
    fitmm,
    coef = 2,
    number = nrow(expr),
    adjust.method = "BH",
    sort.by = "none")
  res = data.frame(
    SummarizedExperiment::rowData(se[rownames(res), ]),
    res,
    stringsAsFactors = FALSE)

  ## add mean expression per group
  info.mean = estimate_means(gr = pheno[, var],
                             expr = expr)
  if (!is.null(info.mean)) {
    res = data.frame(res,
                     info.mean)
  }

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
