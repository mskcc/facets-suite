#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript

# Generate gene-level TCN/LCN and CNA call, given cncf.txt files from a doFacets run

library(argparse)
library(data.table)

write.text <- function (...) {
    write.table(..., quote = F, col.names = T, row.names = F, sep = "\t")
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

# Get IMPACT341 loci and gene names
#msk_impact_341 <- scan('/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv3/genelist', what="", quiet = TRUE)
IMPACT341_targets <- suppressWarnings(fread('grep -v ^@ /ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv3/picard_targets.interval_list'))
setnames(IMPACT341_targets, c("chr", "start", "end", "strand", "name"))
setkey(IMPACT341_targets, chr, start, end)

# Get IMPACT410 loci and gene names
#msk_impact_410 <- scan('/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv5/genelist', what="", quiet = TRUE)
IMPACT410_targets <- suppressWarnings(fread('grep -v ^@ /ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv5/picard_targets.interval_list'))
setnames(IMPACT410_targets, c("chr", "start", "end", "strand", "name"))
setkey(IMPACT410_targets, chr, start, end)

#############################################
### definition of copy number calls in WGD
FACETS_CALL_table <- fread(paste0(getSDIR(), "/FACETS_CALL_table.tsv"))
setkey(FACETS_CALL_table, WGD, mcn, lcn)
#############################################

annotate_integer_copy_number <- function(gene_level, amp_threshold){
  ### Portal/TCGA copy number calls as well as
  ### more advanced copy number calling, including CNLOH etc.

  ### lowest value of tcn for AMP
  AMP_thresh_tcn <- amp_threshold + 1

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
                                 method=c('em','cncf'),
                                 gene_targets,
                                 WGD_threshold = 0.5, ### least value of frac_elev_major_cn for WGD
                                 amp_threshold = 5, ### total copy number greater than this value for an amplification
                                 mean_chrom_threshold = 0, ### total copy number also greater than this value multiplied by the chromosome mean for an amplification
                                 fun.rename = function(filename){filename}){

  ### check version of data.table
  if(packageVersion("data.table") < "1.9.6"){stop("please update data.table to v1.9.6")}

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
  fo_impact <- foverlaps(gene_targets, concat_cncf_txt, nomatch=NA)
  fo_impact <- fo_impact[!is.na(ID)]
  fo_impact[,Hugo_Symbol:=gsub("_.*$", "", name)]

  ### Summarize copy number for each gene
  if(method == 'cncf'){

      gene_level <- fo_impact[!Hugo_Symbol %in% c("Tiling", "FP", "intron"),
                        list(chr = unique(chr),
                             seg.start=unique(loc.start),
                             seg.end=unique(loc.end),
#                                start=unique(start),   ### with these uncommented, the per-gene summarization is broken (??)
#                                end=unique(end),
                             frac_elev_major_cn=unique(frac_elev_major_cn),
                             Nprobes = .N),
                        keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn, lcn=lcn)]
  }

  if(method == 'em'){

      gene_level <- fo_impact[!Hugo_Symbol %in% c("Tiling", "FP", "intron"),
                        list(chr = unique(chr),
                             seg.start=unique(loc.start),
                             seg.end=unique(loc.end),
#                                start=unique(start),   ### with these uncommented, the per-gene summarization is broken (??)
#                                end=unique(end),
                             frac_elev_major_cn=unique(frac_elev_major_cn),
                             Nprobes = .N),
                        keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn=tcn.em, lcn=lcn.em)]
  }


  ### fix bug where lcn == NA even when tcn is 1
  gene_level[tcn == 1 & is.na(lcn), lcn := 0]

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
  gene_level <- annotate_integer_copy_number(gene_level, amp_threshold)

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

#####################################################################################
#####################################################################################


if(!interactive()){

  parser = ArgumentParser()
  parser$add_argument('-f', '--filenames', type='character', nargs='+', help='list of filenames to be processed.')
  parser$add_argument('-o', '--outfile', type='character', help='Output filename.')
  parser$add_argument('-m', '--method', type='character', default='cncf', help='Method used to calculate integer copy number. Allowed values cncf or em [default cncf]')
  parser$add_argument('-t', '--targetFile', type='character', default='IMPACT341', help="IMPACT341, IMPACT410, or a Picard interval list file of gene target coordinates [default IMPACT341]")
  args=parser$parse_args()

  filenames = args$filenames
  outfile = args$outfile
  method = args$method

  if(args$targetFile=="IMPACT341") {

    geneTargets=IMPACT341_targets

  } else if(args$targetFile=="IMPACT410") {

    geneTargets=IMPACT410_targets

  } else {

    # Note the target file needs to not only be in the PICARD interval list format
    # But the names must match the regex: /GENESYMBOL_.*/ (e.g. TP53_target_02)
    geneTargets <- suppressWarnings(fread(paste0('grep -v "^@" ',args$targetFile)))
    setnames(geneTargets, c("chr", "start", "end", "strand", "name"))
    setkey(geneTargets, chr, start, end)

  }

  gene_level_calls = get_gene_level_calls(filenames, method, geneTargets)
  write.text(gene_level_calls, outfile)

}
