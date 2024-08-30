
# metaAnalyzeGeneExprData 0.1.4

## bug fix

run_diff_expr_analysis()
use formula without random effect for function voom()

# metaAnalyzeGeneExprData 0.1.3

## new features

run_diff_expr_analysis()
added fixing of covariates for estimation of surrogate variable 

# metaAnalyzeGeneExprData 0.1.2

## new features

run_diff_expr_analysis()
adjusting for covariates (not yet implemented for surrogate variables or mixed model)

# metaAnalyzeGeneExprData 0.1.1

## bug fixes

estimate_surrogate_var() (internal function)  
- use svaseq for count data and not erroneously sva  
- remove genes with total count = 0  
- remove variables that should be checked for association with surrogate 
variables if they contain only missing values

# metaAnalyzeGeneExprData 0.1.0

## modifications

run_diff_expr_analysis()
- change output to list

## new features

run_diff_expr_analysis()
- estimating and adjusting for surrogate variables (including plots)

# metaAnalyzeGeneExprData 0.0.0.9001

## bug fix

qc_genes() and qc_diff_expr_genes()
- no error message if some pairs of studies have no overlap or no 
differentially expressed genes

## new features

new function plot_nr_de_genes()

# metaAnalyzeGeneExprData 0.0.0.9000

first version
