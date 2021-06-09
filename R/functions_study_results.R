
## todo:
## column mean.expr facultative
## mapping to correct HGNC symbol?

#' Load study results
#'
#' Load results of study specific differential expression analysis.
#'
#' The data.frame info.studies needs to contain the following columns:
#' \itemize{
#' \item id = unique study identifier
#' \item number.subjects = number of study subjects
#' \item number.samples = number of study samples (might be identical to
#' number.subjects)
#' \item technology = technology (e.g. "RNA-seq" or "array")
#' \item platform = specific array or sequencing platform
#' \item file = path to text file with study results
#' \item column.gene = name of column with gene identifier
#' \item column.estimate = name of column with estimate (e.g. beta coefficient)
#' \item column.se = name of column with standard error
#' \item column.pvalue = name of column with unadjusted P value
#' \item column.mean.expr = name of column with mean expression value
#' }
#'
#' @param info.studies data.frame with information about each
#' study (see Details).
#' @param verbose Logical. Should additional information be printed? (default:
#' TRUE).
#'
#' @return list of data.frames with study results
#'
#' @export
load_study_results <- function(info.studies,
                               verbose = TRUE) {

    ## check column names
    cols = c("id",
             "number.subjects",
             "number.samples",
             "technology",
             "platform",
             "file",
             paste("column",
                   c("gene", "estimate", "se", "pvalue", "mean.expr"),
                   sep = "."))
    for (c in cols) {
        if (!(c %in% colnames(info.studies))) {
            stop(paste(c, "not found in info.studies"))
        }
    }

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

        res = data.frame(
            gene = extract_column(dat = dat,
                                  column = info.i$column.gene),
            estimate = extract_column(dat = dat,
                                      column = info.i$column.estimate),
            se = extract_column(dat = dat,
                                column = info.i$column.se),
            pvalue = extract_column(dat = dat,
                                    column = info.i$column.pvalue),
            mean.expr = extract_column(dat = dat,
                                       column = info.i$column.mean.expr),
            stringsAsFactors = FALSE
        )
        res$pvalue.adj = stats::p.adjust(res$pvalue, method = "BH")
        rownames(res) = res$gene
        res.studies[[info.i$id]] = res
    }

    return(res.studies)
}


#' internal function
#' @keywords internal
extract_column <- function(dat, column) {
    if (!(column %in% colnames(dat))) {
        stop(paste("column", column, "not available!"))
    }
    return(dat[, column])
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