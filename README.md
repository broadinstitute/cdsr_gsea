cdsrgsea
================

cdsrgsea contains lightweight wrappers around useful gene set enrichment
functions

## Install

``` r
library(devtools)
devtools::install_github("broadinstitute/cdsr_gsea")
```

The package can then be loaded by calling

``` r
library(cdsrgsea)
```

## load\_gene\_sets

Loads a list of gene sets in term2gene format from taiga. These gene
sets are drawn from the
[Enrichr](https://amp.pharm.mssm.edu/Enrichr/#stats) and
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp)
    libraries.

``` r
gene_sets <- cdsrgsea::load_gene_sets()
```

``` r
gene_sets %>% names() %>% head()
```

    ## [1] "BioCarta"       "BioPlanet"      "BioPlex"        "Cancer_Modules"
    ## [5] "Canonical"      "ChEA"

``` r
gene_sets$Hallmark %>% head()
```

    ## # A tibble: 6 x 2
    ##   term                             gene   
    ##   <chr>                            <chr>  
    ## 1 HALLMARK_TNFA_SIGNALING_VIA_NFKB JUNB   
    ## 2 HALLMARK_TNFA_SIGNALING_VIA_NFKB CXCL2  
    ## 3 HALLMARK_TNFA_SIGNALING_VIA_NFKB ATF3   
    ## 4 HALLMARK_TNFA_SIGNALING_VIA_NFKB NFKBIA 
    ## 5 HALLMARK_TNFA_SIGNALING_VIA_NFKB TNFAIP3
    ## 6 HALLMARK_TNFA_SIGNALING_VIA_NFKB PTGS2

## run\_hyper

Runs overrepresentation analysis based on the hypergeometric
distribution. The function excepts either a small list of significant
genes or a table with gene level stats as input. If a table is provided
the `gene_var`, `rank_var`, and `n_genes` parameters are used to define
a small list of significant genes.

1.  Small list of significant genes

<!-- end list -->

``` r
genes <- diff_expr %>% arrange(logFC) %>% tail(100) %>% .[["Gene"]]
genes %>% head()
```

    ## [1] "MYEOV2" "RRM2B"  "ANKRA2" "RPL35A" "RGS20"  "ORAI3"

``` r
hyper_res <- cdsrgsea::run_hyper(genes,gene_sets$Hallmark)
hyper_res
```

    ## # A tibble: 50 x 8
    ##    term         p_value p_adjust odds_ratio direction  size overlap_size overlap
    ##    <chr>          <dbl>    <dbl>      <dbl> <lgl>     <int>        <int> <list> 
    ##  1 HALLMARK_P… 5.52e-21 2.76e-19      25.8  NA          200           24 <chr […
    ##  2 HALLMARK_A… 3.33e- 5 8.32e- 4       6.70 NA          161            9 <chr […
    ##  3 HALLMARK_T… 1.81e- 4 3.02e- 3       5.28 NA          200            9 <chr […
    ##  4 HALLMARK_O… 9.52e- 4 1.19e- 2       4.55 NA          200            8 <chr […
    ##  5 HALLMARK_M… 1.71e- 2 1.71e- 1       3.20 NA          200            6 <chr […
    ##  6 HALLMARK_P… 9.67e- 2 6.05e- 1       2.90 NA          105            3 <chr […
    ##  7 HALLMARK_H… 5.66e- 2 4.04e- 1       2.59 NA          200            5 <chr […
    ##  8 HALLMARK_X… 5.66e- 2 4.04e- 1       2.59 NA          200            5 <chr […
    ##  9 HALLMARK_C… 1.75e- 1 7.96e- 1       2.17 NA          138            3 <chr […
    ## 10 HALLMARK_U… 1.91e- 1 7.96e- 1       2.08 NA          144            3 <chr […
    ## # … with 40 more rows

1.  A table with gene level
stats

<!-- end list -->

``` r
ora_res <- cdsrgsea::run_hyper(diff_expr,gene_sets$Hallmark,gene_var = "Gene",rank_var = "logFC",dir = "pos")
ora_res
```

    ## # A tibble: 50 x 8
    ##    term         p_value p_adjust odds_ratio direction  size overlap_size overlap
    ##    <chr>          <dbl>    <dbl>      <dbl> <chr>     <int>        <int> <list> 
    ##  1 HALLMARK_P… 5.52e-21 2.76e-19      25.8  pos         200           24 <chr […
    ##  2 HALLMARK_A… 3.33e- 5 8.32e- 4       6.70 pos         161            9 <chr […
    ##  3 HALLMARK_T… 1.81e- 4 3.02e- 3       5.28 pos         200            9 <chr […
    ##  4 HALLMARK_O… 9.52e- 4 1.19e- 2       4.55 pos         200            8 <chr […
    ##  5 HALLMARK_M… 1.71e- 2 1.71e- 1       3.20 pos         200            6 <chr […
    ##  6 HALLMARK_P… 9.67e- 2 6.05e- 1       2.90 pos         105            3 <chr […
    ##  7 HALLMARK_H… 5.66e- 2 4.04e- 1       2.59 pos         200            5 <chr […
    ##  8 HALLMARK_X… 5.66e- 2 4.04e- 1       2.59 pos         200            5 <chr […
    ##  9 HALLMARK_C… 1.75e- 1 7.96e- 1       2.17 pos         138            3 <chr […
    ## 10 HALLMARK_U… 1.91e- 1 7.96e- 1       2.08 pos         144            3 <chr […
    ## # … with 40 more rows

## run\_gsea

Runs gene set enrichment analysis using the `fgseaMultilevel` method
from the `fgsea` package. The function excepts a table with gene level
stats. The `dir` parameter indicates whether to return positive gene
sets, negative gene sets, or
both.

``` r
gsea_res <- cdsrgsea::run_gsea(diff_expr,gene_sets$Hallmark,gene_var = "Gene",rank_var = "logFC",dir = "neg")
gsea_res
```

    ## # A tibble: 39 x 8
    ##    term               p_value p_adjust     ES   NES direction  size leading_edge
    ##    <chr>                <dbl>    <dbl>  <dbl> <dbl> <chr>     <int> <list>      
    ##  1 HALLMARK_E2F_TAR… 1.41e-23 3.53e-22 -0.646 -2.15 neg         179 <chr [118]> 
    ##  2 HALLMARK_G2M_CHE… 2.00e-19 3.33e-18 -0.634 -2.10 neg         163 <chr [88]>  
    ##  3 HALLMARK_MYC_TAR… 1.28e-14 1.60e-13 -0.574 -1.92 neg         196 <chr [93]>  
    ##  4 HALLMARK_MYC_TAR… 4.67e- 4 2.59e- 3 -0.548 -1.65 neg          51 <chr [25]>  
    ##  5 HALLMARK_PANCREA… 3.18e- 2 1.11e- 1 -0.692 -1.49 neg           8 <chr [4]>   
    ##  6 HALLMARK_MITOTIC… 6.71e- 4 3.35e- 3 -0.442 -1.45 neg         144 <chr [83]>  
    ##  7 HALLMARK_INTERFE… 3.33e- 2 1.11e- 1 -0.466 -1.40 neg          52 <chr [21]>  
    ##  8 HALLMARK_SPERMAT… 7.74e- 2 1.94e- 1 -0.447 -1.33 neg          45 <chr [19]>  
    ##  9 HALLMARK_DNA_REP… 2.61e- 2 1.00e- 1 -0.405 -1.32 neg         124 <chr [41]>  
    ## 10 HALLMARK_ESTROGE… 5.03e- 2 1.57e- 1 -0.399 -1.28 neg         106 <chr [33]>  
    ## # … with 29 more rows
