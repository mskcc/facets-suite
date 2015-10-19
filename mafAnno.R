#!/usr/bin/env Rscript



################################################################################################################################
################################################################################################################################

estimated_af_life_history = function(purity, ns, nw, m, M, copies=1, limit=TRUE){
  #Alex's Penson's code and notes
  # taken from The life history of 21 breast cancers page S1.
  # but assuming the mutation is present at M copies rather than 1

  #m = minor allele#; M = major allele#; ns = number of success (t_alt_count); nw = number of failures (t_ref_count)

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
  #Alex Penson's function and notes with a couple of minor edits from me.
  ### replace facets chr 23 with "X"
  df = out$IGV
  df[df$chrom == 23,]$chrom = "X"
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
  #print(dt)
  dt
}


################################################################################################################################
################################################################################################################################

annotate_maf_with_facets_cf_tcn_lcn = function(maf, out, fit, iTumor_Sample_Barcode=NULL){

  #' Alex Penson's function and notes
  #' believe it or not, the most elegant way i could find
  #' to assign NA to mutations that fall outside of segmented
  #' regions is to fill in the gaps in the GRanges object
  #' this requires pulling the chromosome lengths and also

  maf = as.data.table(maf)
  maf_cols = colnames(maf)
  maf$Chromosome = factor(maf$Chromosome)
  setkey(maf,Chromosome,Start_Position,End_Position)

  ## #seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg19) = "NCBI" ### without "chr", "MT" for mitochondrial
  ## sl = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  ## sl = sl[!grepl("_", names(sl))]
  #names(sl) <- gsub("chr", "", names(sl))

  dt = integer_cn_table(out, fit)

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

  #Adapted from Alex's code
  maf = as.data.table(maf)

  maf_Tumor_Sample_Barcodes = unique(maf$Tumor_Sample_Barcode)
  cat(sort(maf_Tumor_Sample_Barcodes))
  cat('\n')
  cat(sort(names(facets_files)))
  cat('\n')
  cat(paste(c("Tumor_Sample_Barcodes missing in maf:", setdiff(maf_Tumor_Sample_Barcodes, names(facets_files))), collapse=" "))
  cat('\n')
  # facets_files = facets_files[names(facets_files) %in% maf_Tumor_Sample_Barcodes]
  no.facets = setdiff(maf_Tumor_Sample_Barcodes, names(facets_files))
  no.facets.data = maf[maf$Tumor_Sample_Barcode %in% no.facets,]
  idi = intersect(names(facets_files), maf_Tumor_Sample_Barcodes)
  maf = maf[maf$Tumor_Sample_Barcode %in% idi]
  maf_list = lapply(idi, function(x){load(facets_files[x]);
                                     maf = annotate_maf_with_facets_cf_tcn_lcn(maf, out, fit, x)})
  if(length(no.facets)){
      maf_list <- c(maf_list, no.facets.data)
  }
  maf = rbindlist(maf_list,fill=T)

  ## maf[,ccf_1copy_:=as.numeric(estimated_af_life_history(purity,
  ##                                                      t_alt_count,
  ##                                                      t_ref_count,
  ##                                                      lcn,
  ##                                                      tcn-lcn,
  ##                                                      copies=1,
  ##                                                      limit=TRUE)),by = 1: nrow(maf)]

  ## maf[,ccf_Mcopies_:=as.numeric(estimated_af_life_history(purity,
  ##                                                        t_alt_count,
  ##                                                        t_ref_count,
  ##                                                        lcn,
  ##                                                        tcn-lcn,
  ##                                                        copies='M',
  ##                                                        limit=TRUE)), by = 1:nrow(maf)]


  #mine
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
  parser$add_argument('-f','--facets_files', type='character', help='file name contains names of facets files.')
  parser$add_argument('-o','--out_maf', type='character', help='file name of CN annotated maf.')
  args=parser$parse_args()

  maf_file = args$maf
  facets_samples_file = args$facets_files
  output_maf_file = args$out_maf

  maf = fread(paste0('grep -v "^#" ', maf_file))
  facets_samples = fread(facets_samples_file)

  facets_files = with(facets_samples, structure(Rdata_filename, .Names = Tumor_Sample_Barcode))

  maf = main(maf, facets_files)

  write.table(maf, file = output_maf_file,
              quote = F, col.names = T, row.names = F, sep = "\t")
}


#facets_files = Sys.glob('TCGA*/*100*.Rdata')
#names(facets_files) = matrix(unlist(strsplit(facets_files,'/')),byrow=T,nc=2)[,1]
#maf = fread('AKT1_UCEC.maf')
#maf$Tumor_Sample_Barcode = maf$pid
#maf_ = main(maf,facets_files)

#nCCF = 1000
#CCFs = seq(0.001,1,0.001)
#alt_allele = 52   #4     #52
#ref_allele = 13   #29    #13
#coverage = 65     #33    #65
#absCN = 3         #2     #3
#m=0
#M=3
#r=3
#purity= 0.693104435 #0.410028


