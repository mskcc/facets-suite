#!/usr/bin/env Rscript

################################################################################################################################
################################################################################################################################

estimated_af_life_history = function(purity, ns, nw, m, M, copies=1, limit=TRUE){

  if(any(is.na(c(purity, ns, nw, m, M)))){return(NA)}

  if(!copies %in% c(1, "M")) stop('copies must be 1 or "M"')
  if(copies == "M"){r = M } ## number of copies present
  else {r = 1}

  allele_fraction = ns / (ns + nw)

  if(limit){frac = min(1, allele_fraction * (purity * (M+m) + 2*(1-purity)) / purity / r)}

  else{frac = allele_fraction * (purity * (M+m) + 2*(1-purity)) / purity / r}
  as.numeric(frac)
}


################################################################################################################################
################################################################################################################################

integer_cn_table = function(out, fit, em=FALSE){

  df = out$IGV
  n.xchr <- nrow(df[df$chrom == 23,])
  if(n.xchr > 0) {
    df[df$chrom == 23,]$chrom = "X"
  }
  df$chrom = factor(df$chrom)
  if(em==TRUE){
    dt = data.table(df,
                    cf=fit$cncf$cf.em,
                    tcn=fit$cncf$tcn.em,
                    mcn=fit$cncf$tcn.em - fit$cncf$lcn.em,
                    lcn=fit$cncf$lcn.em,
                    ploidy=fit$ploidy,
                    purity=fit$purity,
                    dipLogR=out$dipLogR)
  }
  if(em==FALSE){
    dt = data.table(df,
                    cf=fit$cncf$cf,
                    tcn=fit$cncf$tcn,
                    mcn=fit$cncf$tcn - fit$cncf$lcn,
                    lcn=fit$cncf$lcn,
                    ploidy=fit$ploidy,
                    purity=fit$purity,
                    dipLogR=out$dipLogR)
  }
  setkey(dt, chrom, loc.start, loc.end)
  dt
}


################################################################################################################################
################################################################################################################################

annotate_maf_with_facets_cf_tcn_lcn = function(maf, out, fit, iTumor_Sample_Barcode=NULL){

  maf = as.data.table(maf)
  maf_cols = colnames(maf)
  maf$Chromosome = factor(maf$Chromosome)
  setkey(maf,Chromosome,Start_Position,End_Position)
  dt = integer_cn_table(out, fit)

  ### check for duplicate columns
  if(any(duplicated(names(maf)))){
    warning("duplicate columns removed from maf file")
    maf[, which(duplicated(names(maf))) := NULL, with = F]
  }

  if(is.null(iTumor_Sample_Barcode)){maf_ann = foverlaps(maf, dt, mult="first",nomatch=NA)}
  else{maf_ann = foverlaps(maf[Tumor_Sample_Barcode == iTumor_Sample_Barcode], dt, mult="first",nomatch=)}

  maf_ann[,c(maf_cols, 'dipLogR', 'seg.mean', 'cf', 'tcn', 'lcn', 'purity', 'ploidy'), with=F]
}


################################################################################################################################
################################################################################################################################

ccf.likelihood = function(purity, absCN, alt_allele, coverage, copies){

  #From McGranahan_and_Swanton_2015

  CCFs = seq(0.001,1,0.001)
  vac.ccf  = function(CCF, purity, absCN){purity * CCF * copies / (2*(1 - purity) + purity * absCN)}
  probs = sapply(CCFs, function(c){dbinom(alt_allele, coverage, vac.ccf(c, purity, absCN))})
  probs = probs/sum(probs)

  ccf.max = which.max(probs)
  ccf.gt.half.max = which(probs > max(probs)/2)
  ccf.lower = max(ccf.gt.half.max[1] - 1, 1) ### closest ccf value before half-max range (within 0-1 range)
  ccf.upper = min(ccf.gt.half.max[length(ccf.gt.half.max)] + 1, length(CCFs)) ### closest ccf value after half-max range (within 0-1 range)
  if(is.na(purity)){ccf.upper=NA}
  ccf.max = ccf.max/length(CCFs)
  ccf.lower = ccf.lower/length(CCFs)
  ccf.upper = ccf.upper/length(CCFs)
  prob.95 = sum(probs[950:1000])
  prob.90 = sum(probs[900:1000])
  #if(is.na(purity)){ccf.upper=NA}
  list(ccf.max,ccf.lower,ccf.upper,prob.95,prob.90)
}


################################################################################################################################
################################################################################################################################

main = function(maf,facets_files){

  maf = as.data.table(maf)
  maf_Tumor_Sample_Barcodes = unique(maf$Tumor_Sample_Barcode)

  not.in.maf = setdiff(names(facets_files),maf_Tumor_Sample_Barcodes)
  no.facets = setdiff(maf_Tumor_Sample_Barcodes, names(facets_files))

  no.facets.data = maf[maf$Tumor_Sample_Barcode %in% no.facets,]
  maf_Tumor_Sample_Barcodes = maf_Tumor_Sample_Barcodes[!maf_Tumor_Sample_Barcodes %in% no.facets]

  write(paste('Missing facets data:', no.facets), stderr())
  write(paste('Not in MAF:', not.in.maf), stderr())

  idi = intersect(names(facets_files), maf_Tumor_Sample_Barcodes)
  maf = maf[maf$Tumor_Sample_Barcode %in% idi]
  maf_list = lapply(idi, function(x){load(facets_files[x]);
                                     maf = annotate_maf_with_facets_cf_tcn_lcn(maf, out, fit, x)})

  if(length(no.facets)){maf_list = c(maf_list, list(no.facets.data))}
  maf = rbindlist(maf_list,fill=T)

  maf[,c("ccf_Mcopies", "ccf_Mcopies_lower", "ccf_Mcopies_upper", "ccf_Mcopies_prob95", "ccf_Mcopies_prob90"):=ccf.likelihood(purity,
                                                                                                                    tcn,
                                                                                                                    t_alt_count,
                                                                                                                    (t_alt_count + t_ref_count),
                                                                                                                    copies=(tcn-lcn)), by= 1:nrow(maf)]

  maf[,c("ccf_1copy", "ccf_1copy_lower", "ccf_1copy_upper", "ccf_1copy_prob95", "ccf_1copy_prob90"):=ccf.likelihood(purity,
                                                                                                                    tcn,
                                                                                                                    t_alt_count,
                                                                                                                    (t_alt_count + t_ref_count),
                                                                                                                    copies=1), by= 1:nrow(maf)]
   maf
}


################################################################################################################################
################################################################################################################################

suppressPackageStartupMessages(library(data.table))
library(argparse)


if(!interactive()){

  parser=ArgumentParser()
  parser$add_argument('-m','--maf', type='character', help='file name of maf file to be annotated.')
  parser$add_argument('-f','--facets_files', type='character',
                      help='Mapping of "Tumor_Sample_Barcode" from maf and "Rdata_filename" from FACETS (tab-delimited with header)')
  parser$add_argument('-o','--out_maf', type='character', help='file name of CN annotated maf.')
  parser$add_argument('-c','--save_comments', action = 'store_true', default = FALSE, help = 'Retains comments in file header')
  args=parser$parse_args()

  maf_file = args$maf
  facets_samples_file = args$facets_files
  output_maf_file = args$out_maf
  save_comments = args$save_comments

  maf = fread(paste0('grep -v "^#" ', maf_file))
  facets_samples = fread(facets_samples_file)

  facets_files = with(facets_samples, structure(Rdata_filename, .Names = Tumor_Sample_Barcode))

  maf = main(maf, facets_files)

  if (save_comments) {
		header = readLines(maf_file, n = 10)
		header = header[unlist(lapply(header, function(x) substr(x,1,1)=='#'))]
		out_file = file(output_maf_file, open = 'wt')
		for (i in 1:length(header)) cat(header[i], '\n', file = out_file, append = T)
		write.table(maf, file = out_file, quote = F, col.names = T, row.names = F, sep = "\t")
		close(out_file)
  } else {
		  write.table(maf, file = output_maf_file,
              quote = F, col.names = T, row.names = F, sep = "\t")
  }
}

