require(fgsea)
require(plyr)
require(tidyverse)

#' Load gene sets from taiga
#'
#' @return A list of gene set in term2gene format
#' @export load_gene_sets
#'
#'
load_gene_sets <- function() {
  read_rds(download.raw.from.taiga(data.name='msigdb-gene-set-collections-8453',data.file='gsc_data_term2gene'))
}


#' Performs gene set enrichment analysis with fgsea on a table of gene stats
#'
#' @param gene_stats a table of gene level stats.
#' @param term2gene a table of gene_sets in term2gene format.
#' @param gene_var name of the column containing gene names.
#' @param rank_var name of the column to rank genes by.
#' @param dir direction(s) to consider. Options are ("both","pos","neg")
#' @param method GSEA methods. Options are ("fgsea").
#' @param min_size minimal size of a gene set to test.
#' @param max_size maximal size of a gene set to test. Defaults to 500.
#'
#'
#' @return A table with the following columns:
#' \describe{
#' term: gene set name. \cr
#' p_value: enrichment p-value. \cr
#' p_adjust: BH-adjusted p-value. \cr
#' ES: enrichment score. \cr
#' NES: enrichment score normalized to mean enrichment of random samples of the same size. \cr
#' direction: enrichment direction. \cr
#' size: size of the gene set. \cr
#' leading_edge: vector with indexes of leading edge genes that drive the enrichment. \cr
#' }
#'
#' @export run_gsea
run_gsea <- function(gene_stats, term2gene, gene_var = "Gene", rank_var = "logFC",dir = "both",
                     method = "fgsea",min_size = 1,max_size = 500) {
  #Check inputs
  if(!gene_var %in% names(gene_stats)){stop("gene_var must be a column in gene_stats")}
  if(!rank_var %in% names(gene_stats)){stop("rank_var must be a column in gene_stats")}
  if(!dir %in% c("pos","neg","both")){stop("type must be in ('pos','neg','both')")}
  if(!method %in% c("fgsea")){stop("dir must be in ('fgsea')")}
  if(length(intersect(gene_stats[[gene_var]],term2gene$gene)) < 5) {
    stop(str_c("expecting genes in format ",term2gene$gene[1]," ",term2gene$gene[2],"..."))
  }

  #run fgsea
  if(method == "fgsea") {
    gene_set_list <- plyr::dlply(term2gene,"term",function(df) {return(df[["gene"]])})
    stats_vec <- gene_stats[[rank_var]]
    names(stats_vec) <- gene_stats[[gene_var]]
    stats_vec <- rev(sort(stats_vec))
    fgsea_res <- suppressWarnings(fgsea::fgseaMultilevel(pathways = gene_set_list,stats = stats_vec,
                 minSize = min_size,maxSize = max_size))
    res <- fgsea_res %>% as_tibble() %>%
      dplyr::transmute(term=pathway,p_value=pval,p_adjust=padj,ES=ES,NES=NES,direction = ifelse(NES > 0,"pos","neg"),
      size=size,leading_edge = leadingEdge) %>% arrange(-abs(NES))
    if(dir != "both") {
      res <- res %>% filter(direction == dir)
    }
  }

  return(res)
}

fisher_test <- function(genes,term2gene,p_adjust_method,min_size,max_size,dir) {

  res <- term2gene %>% dplyr::group_by(term) %>%
    dplyr::summarise(size = n(),overlap_size = length(intersect(gene,genes)),
                     overlap = list(intersect(gene,genes))) %>%
    dplyr::filter((size > min_size) & (size < max_size))

  #hypergeometric test
  #k <- length(genes)
  m <- res$size
  x <- res$overlap_size
  n <- rep(length(unique(term2gene$gene)),length(m)) - m
  k <- intersect(genes,term2gene$gene) %>% length()
  n <- rep(length(unique(term2gene$gene)),length(m)) - m
  args_df <- tibble(q=x-1,m=m,n=n,k=k)
  p_value <- apply(args_df, 1, function(n)
    phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
  )

  #odds ratio
  prob_in_set <- x/m
  prob_not_in_set <- (k-x)/n
  odds_ratio <- (prob_in_set / (1 - prob_in_set)) / (prob_not_in_set / (1 - prob_not_in_set))

  res <- res %>% mutate(p_value = p_value,p_adjust = p.adjust(p_value, method=p_adjust_method),odds_ratio = odds_ratio,direction = dir) %>%
    select(term,p_value,p_adjust,odds_ratio,direction,size,overlap_size,overlap) %>%
    arrange(-odds_ratio)

  return(res)
}


#' Performs over representation analysis with fisher's exact test on a gene list or table of gene stats
#'
#' @param gene_stats a small list of significant genes or table of gene level stats.
#' @param term2gene a table of gene_sets in term2gene format.
#' @param gene_var if gene_stats is table, name of the column containing gene names.
#' @param rank_var if gene_stats is table, name of the column to rank genes by.
#' @param n_genes if gene_stats is table, number of genes to select after ranking by rank_var
#' @param dir if gene_stats, direction(s) to consider. Options are ("both","pos","neg")
#' @param p_adjust_method Method for p-value adjustment. Defaults to BH
#' @param min_size minimal size of a gene set to test.
#' @param max_size maximal size of a gene set to test.
#'
#'
#' @return A table with the following columns:
#' \describe{
#' term: gene set name. \cr
#' p_value: enrichment p-value. \cr
#' p_adjust: adjusted p-value. \cr
#' odds_ratio: the odds ratio the genes being in the set vs not in the set. \cr
#' direction: enrichment direction. \cr
#' size: size of the gene set. \cr
#' overlap_size: size of the overlap between the gene set and significant genes list. \cr
#' overlap: genes in both the gene set and significant genes list. \cr
#' }
#'
#' @export run_ora
run_ora <- function(gene_stats = NULL,term2gene,gene_var = "Gene",rank_var = "logFC",
                    n_genes = 100,dir = "both",p_adjust_method = "BH",min_size = 1,max_size = Inf) {

  #Small list of genes
  if(is.null(ncol(gene_stats))) {
    if(length(intersect(gene_stats,term2gene$gene)) < 5) {
      stop(str_c("expecting genes in format ",term2gene$gene[1]," ",term2gene$gene[2],"..."))
    }
    res <- fisher_test(gene_stats,term2gene,p_adjust_method,min_size,max_size,dir = NA)
    return(res)
  }

  #Check inputs
  if(!gene_var %in% names(gene_stats)){stop("gene_var must be a column in gene_stats")}
  if(!rank_var %in% names(gene_stats)){stop("rank_var must be a column in gene_stats")}
  if(!is.na(dir) & !dir %in% c("pos","neg","both")){stop("dir must be in ('pos','neg','both')")}

  #Table of gene stats
  stats_vec <- gene_stats[[rank_var]]
  names(stats_vec) <- gene_stats[[gene_var]]
  if(length(intersect(gene_stats[[gene_var]],term2gene$gene)) < 5) {
    stop(str_c("expecting genes in format ",term2gene$gene[1]," ",term2gene$gene[2],"..."))
  }
  genes_pos <- sort(stats_vec) %>% names() %>% tail(n_genes)
  genes_neg <- sort(stats_vec) %>% names() %>% head(n_genes)
  if(!is.null(gene_stats) & dir == "pos"){
    res <- fisher_test(genes_pos,term2gene,p_adjust_method,min_size,max_size,dir = dir)
    return(res)
  }
  if(!is.null(gene_stats) & dir == "neg"){
    res <- fisher_test(genes_neg,term2gene,p_adjust_method,min_size,max_size,dir = dir)
    return(res)
  }
  if(!is.null(gene_stats) & dir == "both"){
    res_pos <- fisher_test(genes_pos,term2gene,p_adjust_method,min_size,max_size,dir = "pos")
    res_neg <- fisher_test(genes_neg,term2gene,p_adjust_method,min_size,max_size,dir = "neg")
    return(rbind(res_pos,res_neg) %>% arrange(-odds_ratio))
  }
}