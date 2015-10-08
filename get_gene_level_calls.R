#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript

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
msk_impact_341 <- scan('/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv3/genelist', what="", quiet = TRUE)
IMPACT341_targets <- suppressWarnings(fread(paste0('grep -v "^@" /ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv3/picard_targets.interval_list')))
setnames(IMPACT341_targets, c("chr", "start", "end", "strand", "name"))
setkey(IMPACT341_targets, chr, start, end)


#############################################
### definition of copy number calls in WGD
WGD_CALL_table <- fread("
WGD	mcn	lcn	FACETS_CNA	FACETS_CALL
no WGD	0	0	-2	HOMDEL
no WGD	1	0	-1	HETLOSS
no WGD	2	0	-1	CNLOH
no WGD	3	0	1	CNLOH & GAIN
no WGD	4	0	1	CNLOH & GAIN
no WGD	5	0	2	AMP (LOH)
no WGD	6	0	2	AMP (LOH)
no WGD	1	1	0	DIPLOID
no WGD	2	1	1	GAIN
no WGD	3	1	1	GAIN
no WGD	4	1	2	AMP
no WGD	5	1	2	AMP
no WGD	6	1	2	AMP
no WGD	2	2	1	TETRAPLOID
no WGD	3	2	2	AMP
no WGD	4	2	2	AMP
no WGD	5	2	2	AMP
no WGD	6	2	2	AMP
no WGD	3	3	2	AMP (BALANCED)
no WGD	4	3	2	AMP
no WGD	5	3	2	AMP
no WGD	6	3	2	AMP
WGD	0	0	-2	HOMDEL
WGD	1	0	-1	LOSS BEFORE & AFTER
WGD	2	0	-1	LOSS BEFORE
WGD	3	0	-1	CNLOH BEFORE & LOSS
WGD	4	0	-1	CNLOH BEFORE
WGD	5	0	1	CNLOH BEFORE & GAIN
WGD	6	0	2	AMP (LOH)
WGD	1	1	-1	DOUBLE LOSS AFTER
WGD	2	1	-1	LOSS AFTER
WGD	3	1	-1	CNLOH AFTER
WGD	4	1	1	LOSS & GAIN
WGD	5	1	2	AMP
WGD	6	1	2	AMP
WGD	2	2	0	TETRAPLOID
WGD	3	2	1	GAIN
WGD	4	2	2	AMP
WGD	5	2	2	AMP
WGD	6	2	2	AMP
WGD	3	3	2	AMP (BALANCED)
WGD	4	3	2	AMP
WGD	5	3	2	AMP
WGD	6	3	2	AMP
")
setkey(WGD_CALL_table, WGD, mcn, lcn)
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

  gene_level <- merge(gene_level, WGD_CALL_table, sort = F, all.x = T)

  gene_level[tcn >= AMP_thresh_tcn, FACETS_CNA := 2]
  gene_level[tcn >= AMP_thresh_tcn, FACETS_CALL := "AMP"]

  output_cols <- c(input_cols, c("FACETS_CNA", "FACETS_CALL"))
  setcolorder(gene_level, output_cols)
  setkeyv(gene_level, input_key)
  gene_level
}

get_gene_level_calls <- function(cncf_files,
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

  ### fix bug where lcn == NA even when tcn is 1
  concat_cncf_txt[tcn == 1 & is.na(lcn), lcn := 0]

  ### format concat_cncf_txt segment table
  concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
  concat_cncf_txt[,Tumor_Sample_Barcode := fun.rename(filename)]
  concat_cncf_txt[, filename := NULL]
  concat_cncf_txt$chrom <- as.character(concat_cncf_txt$chrom)
  concat_cncf_txt[chrom == "23", chrom := "X"]
  setkey(concat_cncf_txt, chrom, loc.start, loc.end)

  ### estimate fraction of the genome with more than one copy from a parent
  ### a large fraction implies whole genome duplication
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
#                                start=unique(start),   ### with these uncommented, the per-gene summarization is broken (??)
#                                end=unique(end),
                               frac_elev_major_cn=unique(frac_elev_major_cn),
                               Nprobes = .N),
                          keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn.em, lcn=lcn.em)]

  ### apply WGD threshold
  gene_level[, WGD := factor(ifelse(frac_elev_major_cn > WGD_threshold, "WGD", "no WGD"))]
  setkey(gene_level, Tumor_Sample_Barcode, Hugo_Symbol)

  ### focality requirement
  ### get (weighted) mean total copy number for the chromosome
#   mean_tcn_per_chr <- concat_cncf_txt[, list(mean_tcn =
#                                                sum(tcn * as.numeric(loc.end - loc.start)) /
#                                                sum(as.numeric(loc.end - loc.start))),
#                                       keyby=list(Tumor_Sample_Barcode, chrom)]
#   mean_tcn_per_chr_v <- with(mean_tcn_per_chr,
#                              structure(mean_tcn,
#                                        .Names = paste(Tumor_Sample_Barcode,
#                                                       chrom)))
#   gene_level[, mean_tcn_per_chr := mean_tcn_per_chr_v[paste(Tumor_Sample_Barcode, chr)]]

  ### annotate integer copy number
  gene_level <- annotate_integer_copy_number(gene_level)

  ### remove duplicate entries for partial deletions
  ### the lower value is chosen on the basis that a
  ### partial deletion is a deletion but a
  ### partial amplification is not an amplification
  gene_level <- gene_level[order(FACETS_CNA, -Nprobes)]
  gene_level <- gene_level[!duplicated(gene_level, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]

  ### Sort & set column classes
  gene_level[, FACETS_CNA := factor(FACETS_CNA, levels = c(-2:2))]
  gene_level[, chr := factor(chr, levels = c(1:22, "X", "Y"))]
  setkey(gene_level, Tumor_Sample_Barcode, Hugo_Symbol)
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
