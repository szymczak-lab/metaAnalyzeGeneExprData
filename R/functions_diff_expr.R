
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
#' @param n.sv [numeric(1)] number of surrogate variables to be estimated
#' @param sv.file [character(1)] name of file for saving surrogate variables
#' (default: NULL so that surrogate variables are not saved).
#' @param sv.var.check [character(n)] name of variables in colData() of se that
#' should be checked for association with surrogate variables (default: NULL).
#' @param title [character(1)] title used in plots for surrogate variables, 
#' i.e. only used if sv = TRUE (default: NULL).
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
  n.sv = 5,
  sv.file = NULL,
  sv.var.check = NULL,
  title = NULL,
  standardize = TRUE,
  res.file,
  BPPARAM = BiocParallel::bpparam()) {
  
  if (!is.null(covar) & (!is.null(var.id))) {
    stop("covariate adjustment not yet implemented for mixed model!")
  }
  
  ## extract phenotype data
  if (!is.null(sv.var.check)) {
    sv.var.check = intersect(sv.var.check,
                             colnames(SummarizedExperiment::colData(se)))
  }
  pheno = as.data.frame(
    SummarizedExperiment::colData(se)[, unique(c(
      var, covar, sv.var.check, var.id)), 
      drop = FALSE])
  ind.rm.l = apply(
    pheno[, c(var, covar), drop = FALSE], 2, function(x) {
    rownames(pheno)[which(is.na(x))]})
  ind.rm = unique(unlist(ind.rm.l))
  if (length(ind.rm) > 0) {
    warning(paste(
      "removed", length(ind.rm), "individuals because of missing 
      values in variable of interest or covariate(s)"))
    pheno = pheno[setdiff(rownames(pheno), ind.rm), ]
    se = se[, rownames(pheno)]
  }

  ## remove variables with only NAs
  all.na = apply(pheno, 2, function(x) {all(is.na(x))})
  pheno = pheno[, !all.na, drop = FALSE]
  if (!is.null(sv.var.check)) {
    sv.var.check = setdiff(colnames(pheno), var)
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
    expr = edgeR::calcNormFactors(
      object = expr)
  }

  ## surrogate variables
  if (sv) {
    pheno = estimate_surrogate_var(
      expr = expr,
      pheno = pheno,
      var = var,
      covar = covar,
      n.sv = n.sv,
      sv.var.check = sv.var.check,
      title = title)
    covar = c(covar, 
              grep("sv_", colnames(pheno), value = TRUE))
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
    ## Note: cannot use formula with random effect in stats::model.matrix()
    if (!is.null(var.id)) {
      form.use = stats::as.formula(paste0("~", var))
      if (!is.null(covar)) {
        form.use = paste0(form.use, "+",
                          paste(covar, collapse = "+"))
      }
    } else {
      form.use = form
    }
    mod = stats::model.matrix(
      stats::as.formula(form.use),
      data = pheno)
    expr = limma::voom(
      counts = expr,
      design = mod)
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
  if (sv & !is.null(sv.file)) {
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
    var,
    covar = NULL,
    n.sv = 5,
    sv.var.check,
    title) {

  ## always estimate 5 surrogate variables
  #n.sv = 5
    
  ## model with variable of interest
  mod = stats::model.matrix(
    stats::as.formula(paste0(
      "~",  
      paste(
        c(var, covar), 
        collapse = "+"))),
    data = pheno)
  
  ## null model
  mod0 = stats::model.matrix(~1,
                             data = pheno) 
  
  ## estimate normalized counts for RNASeq data
  if (methods::is(expr, "DGEList")) {
    type = "rnaseq"
    expr = edgeR::cpm(expr, log = FALSE)
    
    ## remove genes with total count = 0
    sum = apply(expr, 1, sum)
    expr = expr[sum > 0, ]
  } else {
    type = "array"
  }

  ## estimate surrogate variables
  #n.sv = sva::num.sv(
  #  dat = expr,
  #  mod = mod,
  #  method = "leek")
  
  if (type == "rnaseq") { # RNASeq
    svaobj = sva::svaseq(
      dat = expr,
      mod = mod,
      mod0 = mod0,
      n.sv = n.sv)
  } else {
    svaobj = sva::sva(
      dat = expr,
      mod = mod,
      mod0 = mod0,
      n.sv = n.sv)
  }
  
  plot_variance_explained_surrogate_var(
    svaobj = svaobj,
    expr = expr,
    title = title)
  
  sv = svaobj$sv
  colnames(sv) = paste0("sv_", 1:ncol(sv))

  if (length(intersect(colnames(pheno), colnames(sv))) > 0) {
    stop("surrogate variables already available in colData")
  }
  pheno = data.frame(pheno, sv)

  if (!is.null(sv.var.check) & length(sv.var.check) > 0) {
    plot_association_surrogate_var(
      pheno = pheno,
      sv.var.check = sv.var.check,
      title = title)
  }
  
  return(pheno)
}



#' internal function
#' 
#' based on code from
#' https://support.bioconductor.org/p/88553/
#' 
#' @keywords internal
plot_variance_explained_surrogate_var <- function(
    svaobj,
    expr, 
    title = "") {
  
  pprob = svaobj$pprob.gam*(1 - svaobj$pprob.b)
  dats = expr * pprob
  dats = dats - base::rowMeans(dats)
  uu = base::eigen(t(dats) %*% dats)
  uu_val = uu$values / sum(uu$values)
  
  x_val = 1:svaobj$n.sv
  expl_var_plot = data.frame(
    no = 1:svaobj$n.sv, 
    var.expl = uu_val[1:svaobj$n.sv])
  
  g = ggpubr::ggdotchart(
    expl_var_plot,
    x = "no",
    y = "var.expl",
    xlab = "surrogate variable",
    ylab = "Proportion of variance explained",
    dot.size = 2,
    title = title,
    sorting = "descending",
    ylim = c(0, max(expl_var_plot$var.expl))) +
    ggpubr::rotate_x_text(0)
  print(g)
}


#' internal function
#' 
#' @keywords internal
plot_association_surrogate_var <- function(
  pheno,
  sv.var.check,
  title) {
  
  if (is.null(sv.var.check)) {
    stop("sv.var.check needs to be defined")
  }
  r.sq = NULL
  var.used = NULL
  for (v in sv.var.check) {
    
    if (v %in% colnames(pheno) & 
        !all(is.na(pheno[, v])) &
        length(unique(pheno[, v])) > 1) {
#      print(paste("analyzing variable", v))
      
      group = pheno[, v]
      
      n.sv = length(grep("^sv_", colnames(pheno)))
      r.sq.v = NULL
      for (j in 1:n.sv) {
        x = pheno[, paste0("sv_", j)]
      
        ## using code from check_batch_effects() in QCnormSE
        fit = stats::lm(x ~ group)
        s = summary(fit)
        r.sq.v[j] = s$adj.r.squared
      }
      r.sq = rbind(r.sq, r.sq.v)
      var.used = c(var.used, v)
      
    }
  }
  dimnames(r.sq) = list(var.used, 
                        paste0("sv_", 1:ncol(r.sq)))

  g = plot_heatmap(
    matrix = abs(stats::na.omit(r.sq)), 
    title = title)
  print(g)

}

#' internal function
#' @keywords internal
plot_heatmap <- function(
    matrix,
    title = NULL) {
  
    col = circlize::colorRamp2(
      c(0, 1),
      c("white", "blue"))
    name = "|r^2|"
    
  hm = ComplexHeatmap::Heatmap(
    matrix,
    rect_gp = grid::gpar(col = "black"),
    name = name,
    col = col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(border = "black"),
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 45,
    column_title = title)
  
  return(hm)
  
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
