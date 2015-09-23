#!/usr/bin/env Rscript

#### usage ./get_gene_level_calls.R output_file.txt *_cncf.txt
  
library(data.table)

write.text <- function (...) {
  write.table(..., quote = F, col.names = T, row.names = F, 
              sep = "\t")
}

getSDIR <- function(){
  args=commandArgs(trailing=F)
  TAG="--file="
  path_idx=grep(TAG,args)
  SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
  if(length(SDIR)==0) {
    return(getwd())
  } else {
    return(SDIR)
  }
}

### get IMPACT341 loci and gene names
msk_impact_341 <- scan(file.path(getSDIR(), "data", "cv3_hg19_gene_impact341.list"), what="")
IMPACT341_targets <- fread(paste0('grep -v "^@" ', getSDIR(), '/data/cv3_hg19_picard_targets.interval_list'))
setnames(IMPACT341_targets, c("chr", "start", "end", "strand", "name"))
setkey(IMPACT341_targets, chr, start, end)


get_gene_level_calls <- function(cncf_files,
                                 amp_threshold = 5, ### total copy number greater than this value for an amplification
                                 mean_chrom_threshold = 0, ### total copy number also greater than this value multiplied by the chromosome mean for an amplification
                                 fun.rename = function(filename){filename}
){
  cncf_txt_list <- lapply(cncf_files, fread)
  names(cncf_txt_list) <- cncf_files
  concat_cncf_txt <- rbindlist(cncf_txt_list, idcol = "filename")
  concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
  setkey(concat_cncf_txt, chrom, loc.start, loc.end)
  concat_cncf_txt[,Tumor_Sample_Barcode := fun.rename(filename)]
  concat_cncf_txt[, filename := NULL]
  concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
  concat_cncf_txt[chrom == "23", chrom := "X"]
  setkey(concat_cncf_txt, chrom, loc.start, loc.end)
  
  ### estimate WGD
  concat_cncf_txt[, frac_elev_major_cn := sum(
    as.numeric((tcn - lcn) >= 2) * 
      as.numeric(loc.end-loc.start), na.rm = T) / 
      sum(as.numeric(loc.end-loc.start)
      ), 
    by=Tumor_Sample_Barcode]  
  
  ### Extract integer copy number for each probe from concat_cncf_txt
  fo_impact <- foverlaps(IMPACT341_targets, concat_cncf_txt, nomatch=NA)
  fo_impact[,Hugo_Symbol:=gsub("_.*$", "", name)]

  ### Summarize copy number for each gene
  gene_level <- fo_impact[!Hugo_Symbol %in% c("Tiling", "FP"), 
                          list(chr = unique(chr), 
                               seg.start=unique(loc.start), 
                               seg.end=unique(loc.end), 
                               start=unique(start), 
                               end=unique(end), 
                               frac_elev_major_cn=unique(frac_elev_major_cn),
                               Nprobes = .N), 
                          keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn.em, lcn=lcn.em)]
  
  ### Call TCGA-style copy number (-2, -1, 0, 1, 2)
  gene_level[, FACETS_CNA := ifelse(tcn == 0, -2, #HOM DEL
                                    ifelse(tcn == 1 | (tcn <= 2 & !is.na(lcn) & lcn == 0), -1, #HET LOSS
                                           ifelse(tcn >= 6, 2, #AMP
                                                  ifelse(tcn >= 4, 1, 0)#GAIN
                                           )
                                    )
  )
  ]

  ### focality requirement
  ### get (weighted) mean total copy number for the chromosome
  mean_tcn_per_chr <- concat_cncf_txt[, list(mean_tcn = 
                                       sum(tcn * as.numeric(loc.end - loc.start)) / 
                                       sum(as.numeric(loc.end - loc.start))), 
                              keyby=list(Tumor_Sample_Barcode, chrom)]
  mean_tcn_per_chr_v <- with(mean_tcn_per_chr, 
                             structure(mean_tcn, 
                                       .Names = paste(Tumor_Sample_Barcode, 
                                                      chrom)))
  gene_level[, mean_tcn_per_chr := mean_tcn_per_chr_v[paste(Tumor_Sample_Barcode, chr)]]

  
  ### More advanced copy number calling, including CNLOH for example
  gene_level[, FACETS_CALL := factor(
    levels = c("HOMDEL", "HETLOSS", "CNLOH", "DIPLOID", "GAIN", "AMP", "OTHER"),
    ifelse(tcn == 0, "HOMDEL", #HOM DEL
           ifelse(tcn == 1 & !is.na(lcn) & lcn == 0, "HETLOSS", #HET LOSS
                  ifelse(tcn == 2 & !is.na(lcn) & lcn == 0, "CNLOH", #CN LOH
                         ifelse(tcn == 2 & !is.na(lcn) & lcn == 1, "DIPLOID", #DIPLOID
                                ifelse(tcn > amp_threshold & 
                                         tcn > mean_tcn_per_chr * mean_chrom_threshold, ### factor elevation above chrom average
                                       "AMP", #AMP
                                       ifelse(tcn > amp_threshold - 2, "GAIN", "OTHER"
                                       )
                                )
                         )
                  )
           )
    )
  )
  ]
  
  ### remove duplicate entries for partial deletions
  ### the lower value is chosen on the basis that a 
  ### partial deletion is a deletion but a 
  ### partial amplification is not an amplification
  
  gene_level <- gene_level[order(FACETS_CNA, -Nprobes)]
  gene_level <- gene_level[!duplicated(gene_level, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]

  ### Sort & set column classes
  gene_level <- gene_level[order(Tumor_Sample_Barcode, Hugo_Symbol)]
  gene_level[, FACETS_CNA := factor(FACETS_CNA, levels = c(-2:2))]
  gene_level[, chr := factor(chr, levels = c(1:22, "X", "Y"))]
  gene_level
}


if(!interactive()){
  args <- commandArgs(TRUE)
  #### usage ./get_gene_level_calls.R output_file.txt *_cncf.txt
  output_file <- args[1]; args <- args[-1]
  filenames <- args
  gene_level_calls <- get_gene_level_calls(filenames)
  write.text(gene_level_calls, output_file)
}
