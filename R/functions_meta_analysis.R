
## res.studies: list of results of differential expression analysis
## min.no.studies: only genes available in min.no.studies will be used
## res.file: text file for saving results of fixed and random effect meta analysis
#' @export
run_meta_analysis <- function(res.studies,
                              min.no.studies = 3,
                              res.file) {

  ## extract genes to be used
  info.gene.names = unlist(lapply(res.studies, function(x) {x$gene}))
  tab.name = table(info.gene.names)
  genes = sort(names(tab.name)[tab.name >= min.no.studies])

  ## gene level meta analysis
  res.all = NULL
  for (g in genes) {
    dat = stats::na.omit(data.frame(
      estimate = sapply(res.studies, function(x) {x[g, "estimate"]}),
      se = sapply(res.studies, function(x) {x[g, "se"]}),
      pvalue = sapply(res.studies, function(x) {x[g, "pvalue"]})))

    res.fixed = metafor::rma.uni(
      yi = dat$estimate,
      sei = dat$se,
      method = "FE")
    res.random = metafor::rma.uni(
      yi = dat$estimate,
      sei = dat$se,
      method = "DL")
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
  res.all$pvalue.fixed.adj = stats::p.adjust(
    res.all$pvalue.fixed,
    method = "BH")
  res.all$pvalue.random.adj = stats::p.adjust(
    res.all$pvalue.random,
    method = "BH")

  rio::export(res.all,
              file = res.file)
}
