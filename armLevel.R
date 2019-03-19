#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript

#### usage ./get_gene_level_calls.R output_file.txt *_cncf.txt

library(argparse)
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
arm_definitions <- suppressWarnings(fread(file.path(getSDIR(), 'data/CytobandTableRed.Ensemblgrch37.txt')))
setnames(arm_definitions, c("chr", "start", "end", "arm"))
arm_definitions <- arm_definitions[chr != "Y"]
arm_definitions[, arm := factor(arm, levels = unique(arm))]
arm_definitions[, chr := factor(chr, levels = unique(chr))]
setkey(arm_definitions, chr, start, end)


#############################################
### definition of copy number calls in WGD
FACETS_CALL_table <- fread(file.path(getSDIR(), 'data/FACETS_CALL_table.tsv'))
setkey(FACETS_CALL_table, WGD, mcn, lcn)
### lowest value of tcn for AMP
AMP_thresh_tcn <- 6
#############################################



annotate_integer_copy_number <- function(gene_level){
  ### Portal/TCGA copy number calls as well as
  ### more advanced copy number calling, including CNLOH etc.

  ### assumes that WGD variable has already been set
  input_key <- key(gene_level)
  input_cols <- names(gene_level)
  gene_level[, mcn := tcn - lcn]
  setkey(gene_level, WGD, mcn, lcn)

  if("FACETS_CNA" %in% names(gene_level)) gene_level[, FACETS_CNA := NULL]
  if("FACETS_CALL" %in% names(gene_level)) gene_level[, FACETS_CALL := NULL]

  gene_level <- merge(gene_level, FACETS_CALL_table, sort = F, all.x = T)

  gene_level[tcn >= AMP_thresh_tcn, FACETS_CNA := 2]
  gene_level[tcn >= AMP_thresh_tcn, FACETS_CALL := "AMP"]

  output_cols <- c(input_cols, c("FACETS_CNA", "FACETS_CALL"))
  setcolorder(gene_level, output_cols)
  setkeyv(gene_level, input_key)
  gene_level
}

get_gene_level_calls <- function(cncf_files,
                                 method,
                                 WGD_threshold = 0.5, ### least value of frac_elev_major_cn for WGD
                                 amp_threshold = 5, ### total copy number greater than this value for an amplification
                                 mean_chrom_threshold = 0, ### total copy number also greater than this value multiplied by the chromosome mean for an amplification
                                 fun.rename = function(filename){filename}
){

  ### check version of data.table
  if(packageVersion("data.table") < "1.9.6"){
    stop("please update data.table to v1.9.6")
  }

  ### concatenate input files
  cncf_txt_list <- lapply(cncf_files, fread)
  names(cncf_txt_list) <- cncf_files
  concat_cncf_txt <- rbindlist(cncf_txt_list, idcol = "filename")

  ### format concat_cncf_txt segment table
  concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
  concat_cncf_txt[,Tumor_Sample_Barcode := fun.rename(filename)]
  concat_cncf_txt[, filename := NULL]
  concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
  concat_cncf_txt[chrom == "23", chrom := "X"]
  setkey(concat_cncf_txt, chrom, loc.start, loc.end)

  if (!("tcn" %in% names(concat_cncf_txt))) {
    concat_cncf_txt[, c("tcn", "lcn") := list(tcn.em, lcn.em)]
  }

  ### estimate fraction of the genome with more than one copy from a parent
  ### a large fraction implies whole genome duplication
  concat_cncf_txt[, frac_elev_major_cn := sum(
    as.numeric((tcn - lcn) >= 2) *
      as.numeric(loc.end-loc.start), na.rm = T) /
      sum(as.numeric(loc.end-loc.start)
      ),
    by=Tumor_Sample_Barcode]
  
  ### Extract integer copy number for each probe from concat_cncf_txt
  fo_impact <- foverlaps(arm_definitions, concat_cncf_txt, nomatch=NA)
  ### Truncate segments that span two arms
  fo_impact[, loc.start := ifelse(loc.start < start, start, loc.start)]
  fo_impact[, loc.end := ifelse(loc.end > end, end, loc.end)]
  fo_impact[, arm := factor(arm, levels = levels(arm_definitions$arm))]
  #fo_impact[,Hugo_Symbol:=gsub("_.*$", "", name)]

  ### Summarize copy number for each gene
  if(method == 'cncf'){

      gene_level <- fo_impact[,
                       list(frac_elev_major_cn=unique(frac_elev_major_cn),
                            Nsegments = .N,
                            length_CN = sum(as.numeric(loc.end - loc.start))),
                          by=list(Tumor_Sample_Barcode, arm, tcn=tcn, lcn=lcn)]

  }
  if(method == 'em'){

      gene_level <- fo_impact[,
                          list(frac_elev_major_cn=unique(frac_elev_major_cn),
                               Nsegments = .N,
                               length_CN = sum(as.numeric(loc.end - loc.start))),
                          by=list(Tumor_Sample_Barcode, arm, tcn=tcn.em, lcn=lcn.em)]
  }
  
  setkey(gene_level, Tumor_Sample_Barcode, arm)
  ### for each CN status, calculate fraction of arm covered
  setkey(arm_definitions, arm)
  gene_level[, arm_frac:= round(length_CN / arm_definitions[as.character(gene_level$arm), as.numeric(end - start)], digits = 4)]
  setkey(arm_definitions, chr, start, end)
  ### ignore CN status valuesif present in less than 10% of arm
  ## gene_level <- gene_level[arm_frac >= 0.1 ]
  gene_level <- gene_level[order(Tumor_Sample_Barcode, arm, -arm_frac)]
  ### estimate WGD from frac_elev_major_cn
  gene_level[, WGD := factor(ifelse(frac_elev_major_cn > WGD_threshold, "WGD", "no WGD"))]

  ### fix bug where lcn == NA even when tcn is 1
  gene_level[tcn == 1 & is.na(lcn), lcn := 0]


  ### annotate integer copy number
  gene_level <- annotate_integer_copy_number(gene_level)

  ### remove duplicate entries for partial deletions
  ### the lower value is chosen on the basis that a
  ### partial deletion is a deletion but a
  ### partial amplification is not an amplification
#   gene_level <- gene_level[order(FACETS_CNA, -Nprobes)]
#   gene_level <- gene_level[!duplicated(gene_level, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]

  gene_level <- gene_level[order(Tumor_Sample_Barcode, arm, -arm_frac)]
  gene_level[,primary := as.integer(!duplicated(gene_level, by=c("Tumor_Sample_Barcode", "arm")))]

  setkey(gene_level, Tumor_Sample_Barcode, arm)
  gene_level[, tcn_summary := as.character(tcn)]
  gene_level[duplicated(gene_level) | duplicated(gene_level, fromLast = T),
             tcn_summary := paste(collapse = " ",
                                  paste(sep = ":",
                                        paste0(
                                          round(arm_frac*100, 0), "%"),
                                        tcn)),
             by = list(Tumor_Sample_Barcode, arm)]
  gene_level
}

#####################################################################################
#####################################################################################


if(!interactive()){

  parser = ArgumentParser()
  parser$add_argument('-f', '--filenames', type='character', nargs='+', help='list of filenames to be processed.')
  parser$add_argument('-o', '--outfile', type='character', help='Output filename.')
  parser$add_argument('-m', '--method', type='character', default='cncf', help='Method used to calculate integer copy number. Allowed values cncf or em')

  args=parser$parse_args()

  filenames = args$filenames
  outfile = args$outfile
  method = args$method
  
  #### usage ./get_gene_level_calls.R output_file.txt *_cncf.txt
  gene_level_calls = get_gene_level_calls(filenames, method)
  write.text(gene_level_calls, outfile)
}



