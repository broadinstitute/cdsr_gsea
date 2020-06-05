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

## run\_ora

Runs overrepresentation analysis with Fisher’s exact test. The function
excepts either a small list of significant genes or a table with gene
level stats as input. If a table is provided the `gene_var`, `rank_var`,
and `n_genes` parameters are used to define a small list of significant
genes.

1.  Small list of significant genes

<!-- end list -->

``` r
genes <- diff_expr %>% arrange(logFC) %>% tail(100) %>% .[["Gene"]]
genes %>% head()
```

    ## [1] "MYEOV2" "RRM2B"  "ANKRA2" "RPL35A" "RGS20"  "ORAI3"

``` r
ora_res <- cdsrgsea::run_ora(genes,gene_sets$Hallmark)
ora_res
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
ora_res <- cdsrgsea::run_ora(diff_expr,gene_sets$Hallmark,gene_var = "Gene",rank_var = "logFC",dir = "pos")
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
gsea_res <- cdsrgsea::run_gsea(diff_expr,gene_sets$Hallmark,gene_var = "Gene",rank_var = "logFC",dir = "pos")
gsea_res
```

    ## # A tibble: 11 x 8
    ##    term                p_value p_adjust    ES   NES direction  size leading_edge
    ##    <chr>                 <dbl>    <dbl> <dbl> <dbl> <chr>     <int> <list>      
    ##  1 HALLMARK_P53_PATH… 6.65e-34 3.32e-32 0.672 3.66  pos         142 <chr [70]>  
    ##  2 HALLMARK_TNFA_SIG… 4.29e-11 4.29e-10 0.467 2.35  pos         124 <chr [53]>  
    ##  3 HALLMARK_HYPOXIA   3.02e-10 2.51e- 9 0.450 2.29  pos         116 <chr [28]>  
    ##  4 HALLMARK_APOPTOSIS 4.29e- 7 3.06e- 6 0.416 2.03  pos         109 <chr [30]>  
    ##  5 HALLMARK_INFLAMMA… 3.81e- 3 1.59e- 2 0.337 1.67  pos          74 <chr [14]>  
    ##  6 HALLMARK_COAGULAT… 3.45e- 3 1.57e- 2 0.383 1.66  pos          50 <chr [14]>  
    ##  7 HALLMARK_EPITHELI… 1.68e- 3 8.40e- 3 0.304 1.54  pos         108 <chr [27]>  
    ##  8 HALLMARK_KRAS_SIG… 4.88e- 2 1.51e- 1 0.390 1.41  pos          18 <chr [3]>   
    ##  9 HALLMARK_MYOGENES… 1.33e- 1 3.17e- 1 0.269 1.37  pos          81 <chr [15]>  
    ## 10 HALLMARK_KRAS_SIG… 9.67e- 2 2.42e- 1 0.258 1.29  pos          76 <chr [26]>  
    ## 11 HALLMARK_COMPLEME… 1.00e+ 0 1.00e+ 0 0.176 0.883 pos          88 <chr [21]>
