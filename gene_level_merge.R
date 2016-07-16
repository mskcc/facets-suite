#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# Merging gene-level FACETS calls from purity, hisens runs
# TODO: Extend to handle 3+ Cvals
##########################################################################################
##########################################################################################

merge_calls <- function(hisens, purity) {
  
  hisens[, id := basename(gsub("_hisens.cncf.txt", "", Tumor_Sample_Barcode))]
  purity[, id := basename(gsub("_purity.cncf.txt", "", Tumor_Sample_Barcode))]
  
  hisens_fails <- setdiff(unique(purity$id), unique(hisens$id))
  purity <- purity[id %in% hisens_fails]
  merged <- rbind(hisens, purity)
  merged[, id := NULL]
  
  return(merged)
  
}

if( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse')
  lapply(pkgs, require, character.only = T)
  
  parser=ArgumentParser()
  parser$add_argument('-s', '--hisens', type='character', help='FACETS hisens gene-level calls')
  parser$add_argument('-p', '--purity', type='character', help='FACETS purity gene-level calls')
  parser$add_argument('-o', '--outfile', type='character', help='Output filename.')
  args=parser$parse_args()
  
  hisens <- fread(args$hisens)
  purity <- fread(args$purity)
  outfile <- args$outfile
  
  merged_calls <- merge_calls(hisens, purity)
  write.table(merged_calls, outfile, 
              quote = F, 
              col.names = T, 
              row.names = F,
              sep = "\t")
   
}
