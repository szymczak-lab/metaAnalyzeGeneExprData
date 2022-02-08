
## todo:
## mapping to correct HGNC symbol?

#' Load study results
#'
#' Load results of study specific differential expression analysis.
#'
#' The data.frame info.studies has to contain the following required and 
#' optional columns:
#' \itemize{
#' \item id = unique study identifier (required)
#' \item number.subjects = number of study subjects
#' \item number.samples = number of study samples (might be identical to
#' number.subjects)
#' \item technology = technology (e.g. "RNA-seq" or "array")
#' \item platform = specific array or sequencing platform
#' \item file = path to text file with study results (required)
#' \item column.gene = name of column with gene identifier (required)
#' \item column.estimate = name of column with estimate (e.g. beta coefficient) 
#' (required)
#' \item column.se = name of column with standard error (required)
#' \item column.pvalue = name of column with unadjusted P value
#' \item column.mean.expr = name of column with mean expression value
#' \item column.mean.gr.1 = name of column with mean expression value in group 1
#' \item column.mean.gr.2 = name of column with mean expression value in group 2
#' }
#'
#' @param info.studies [data.frame] information about each
#' study (see Details).
#' @param verbose [logical(1)] should additional information be printed?
#' (default: TRUE).
#'
#' @return list of data.frames with study results
#'
#' @export
load_study_results <- function(
  info.studies,
  verbose = TRUE) {
  
  ## check column names
  cols.req = c("id",
               "file",
               paste("column", 
                     c("gene", "estimate", "se"),
                     sep = "."))
  cols.opt = c(paste("column",
                     c("pvalue", "mean.expr",
                       "mean.gr.1", "mean.gr.2"),
                     sep = "."))
  cols.not.found = setdiff(cols.req, colnames(info.studies))
  if (length(cols.not.found) > 0) {
    stop(paste0("Some columns not found in info.studies:\n",
                   paste(cols.not.found, collapse = ", "), "\n"))
  }
  cols.use = grep("^column", c(cols.req, cols.opt), value = TRUE)
  
  k = nrow(info.studies)
  res.studies = vector("list", length = k)
  
  if (any(duplicated(info.studies$id))) {
    stop("id in info.studies must be unique!")
  }
  names(res.studies) = info.studies$id
  
  for (i in seq_len(k)) {
    info.i = info.studies[i, ]
    if (verbose)
      print(paste("loading study", info.i$id, "..."))
    if (!file.exists(info.i$file)) {
      stop(paste("file", info.i$file, "does not exist!"))
    }
    dat = rio::import(file = info.i$file)
    
    dat = check_gene_ids(dat = dat)
    
    res = NULL
    for (c in cols.use) {
      if (!(c %in% colnames(info.i))) {
        warning(paste("column", c, "not available in info.studies!\n"))
        info.col = rep(NA, nrow(dat))
      } else {
        c.use = info.i[, c]
        if (!(c.use %in% colnames(dat))) {
          warning(paste("column", c.use, "not available in dat\n!"))
          info.col = rep(NA, nrow(dat))
        } else {
          info.col = dat[, c.use]
        }
      }
      if (is.null(res)) {
        res = data.frame(
          info.col,
          stringsAsFactors = FALSE)
      } else {
        res = data.frame(
          res,
          info.col,
          stringsAsFactors = FALSE)
      }
    }
    colnames(res) = gsub("column.", "", cols.use)
    rownames(res) = res$gene
    
    res$pvalue.adj = stats::p.adjust(res$pvalue, method = "BH")
    
    res.studies[[info.i$id]] = res
  }
  
  return(res.studies)
}

#' internal function
#' @keywords internal
check_gene_ids <- function(dat, verbose = TRUE) {
  genes = dat$gene
  index.bad = grep(" |;|,|\\(|\\|", genes)
  if (length(index.bad) > 0) {
    if (verbose)
    print(paste("removing", length(index.bad), "genes"))
    dat = dat[-index.bad, ]
  }
  if (any(duplicated(dat$gene))) {
    stop("gene ids must be unique!")
  }
  return(dat)
}
