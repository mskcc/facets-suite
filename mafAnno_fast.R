#!/usr/bin/env Rscript

fread_rbind <- function(filenames, fn = fread){
  rbindlist(idcol = "filename",
            sapply(USE.NAMES = T,
                   simplify = F,
                   filenames,
                   fn))
}

out_file_table <- function(filenames){
  out <- fread_rbind(
    filenames,
    fn = function(filename)
      fread(
        paste0("(sed 's/^#//' | grep -v '^$' | grep -v '^ $') < ", filename),
        sep = "=", header = F, fill = T
      )
  )
  setnames(out,c("Tumor_Sample_Barcode", "variable", "value"))
  out <- out[!variable %in% c("INPUT PARAMETERS GIVEN",
                              "LOADED MODULE INFO",
                              "FACETS OUTPUT")]
  out
}

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
  dt = integer_cn_table(out, fit, em = F)

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

annotate_maf_with_facets_cf_tcn_lcn_cncf = function(maf, cncf){

  maf = as.data.table(maf)
  maf <- maf[Tumor_Sample_Barcode %in% cncf_files$Tumor_Sample_Barcode]
  maf_cols = colnames(maf)
  maf[, Chromosome := factor(as.character(Chromosome))]
  setkey(maf,
         Tumor_Sample_Barcode,
         Chromosome,
         Start_Position,
         End_Position)
  cncf[, chrom := factor(as.character(chrom))]
  cncf[chrom == 23, chrom := "X"]

  setkey(cncf,
         Tumor_Sample_Barcode,
         chrom,
         loc.start,
         loc.end)

  #dt = integer_cn_table(out, fit, em = F)

  ### check for duplicate columns
  if(any(duplicated(names(maf)))){
    warning("duplicate columns removed from maf file")
    maf[, which(duplicated(names(maf))) := NULL, with = F]
  }

  # if(is.null(iTumor_Sample_Barcode)){
  #   maf_ann = foverlaps(maf, dt, mult="first",nomatch=NA)
  # }else{
  maf_ann = foverlaps(maf, cncf,
                      by.x = c("Tumor_Sample_Barcode",
                               "Chromosome",
                               "Start_Position",
                               "End_Position"),
                      mult="first",
                      nomatch=NA)
  # }

  setcolorder(maf_ann, c(maf_cols, setdiff(names(maf_ann), maf_cols)))
  copy(maf_ann)
  #maf_ann[,c(maf_cols, 'dipLogR', 'seg.mean', 'cf', 'tcn', 'lcn', 'purity', 'ploidy'), with=F]
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


ccf.likelihood.purity.fit = function(purity, absCN, alt_allele, coverage, copies){

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


main = function(maf, cncf_files, purity_fit_filename = ""){

  ### read cncf
  cncf <- rbindlist(idcol = "Tumor_Sample_Barcode", fill = T,
                    sapply(USE.NAMES = TRUE, simplify = FALSE,
                           cncf_files[file.exists(cncf_filename), structure(cncf_filename, .Names = Tumor_Sample_Barcode)],
                           fread))
  cncf_files[, out_filename := gsub(".cncf.txt$", ".out", cncf_filename)]
  out <- out_file_table(
    cncf_files[file.exists(out_filename),
               structure(out_filename,
                         .Names = Tumor_Sample_Barcode)]
  )
  out_table <- dcast.data.table(
    out[variable %in% c("Purity", "Ploidy", "dipLogR", "dipt", "loglik")],
    Tumor_Sample_Barcode ~ variable,
    value.var = "value"
  )
  if("Purity" %in% names(out_table)) setnames(out_table, "Purity", "purity")
  if("Ploidy" %in% names(out_table)) setnames(out_table, "Ploidy", "ploidy")
  out_table[,purity := as.numeric(purity)]
  out_table[,ploidy := as.numeric(ploidy)]
  out_table[,dipLogR := as.numeric(dipLogR)]
  out_table[,loglik := as.numeric(loglik)]

  cncf_out <- merge(cncf, out_table, by = "Tumor_Sample_Barcode", all = T, allow.cartesian=TRUE)


  maf = as.data.table(maf)
  # maf_Tumor_Sample_Barcodes = unique(maf$Tumor_Sample_Barcode)
  #
  # not.in.maf = setdiff(names(cncf_files),maf_Tumor_Sample_Barcodes)
  # no.facets = setdiff(maf_Tumor_Sample_Barcodes, names(cncf_files))
  #
  # no.facets.data = maf[maf$Tumor_Sample_Barcode %in% no.facets,]
  # maf_Tumor_Sample_Barcodes = maf_Tumor_Sample_Barcodes[!maf_Tumor_Sample_Barcodes %in% no.facets]
  #
  # #write(paste('Missing facets data:', no.facets), stderr())
  # write(paste('Not in MAF:', not.in.maf), stderr())
  #
  # idi = intersect(names(cncf_files), maf_Tumor_Sample_Barcodes)
  # maf = maf[maf$Tumor_Sample_Barcode %in% idi]
  maf = annotate_maf_with_facets_cf_tcn_lcn_cncf(maf, cncf_out)

  # if(length(no.facets)){maf_list = c(maf_list, list(no.facets.data))}
  # maf = rbindlist(maf_list,fill=T)

  maf[,t_alt_count := as.numeric(t_alt_count)]
  maf[,t_ref_count := as.numeric(t_ref_count)]
  maf[, paste0("ccf_Mcopies",
              c("", "_lower", "_upper", "_prob95", "_prob90")) :=
                ccf.likelihood(purity,
                               tcn,
                               t_alt_count,
                               (t_alt_count + t_ref_count),
                               copies=(tcn-lcn)), by= 1:nrow(maf)]

  maf[, paste0("ccf_1copy",
              c("", "_lower", "_upper", "_prob95", "_prob90")) :=
                ccf.likelihood(purity,
                               tcn,
                               t_alt_count,
                               (t_alt_count + t_ref_count),
                               copies=1), by= 1:nrow(maf)]

  if (purity_fit_filename != ""){
    purity_fit <- readRDS(purity_fit_filename)
    sink(file = "/dev/null")
    maf[!is.infinite(purity) & !is.na(purity), c("mu", "sigma", "nu", "tau") := predictAll(purity_fit, newdata = data.frame(purity = purity))]
    sink()
    maf[!is.infinite(purity) & !is.na(purity),
        c("purity_max", "purity_lower", "purity_upper") :=
        {
          suppressWarnings(rm(purity_pdf))
          purity_pdf <- Vectorize(
            function(x) {
              dBCTtr(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
            }
          )
          optimum <-
            optimize(
              purity_pdf,
              c(0, 1),
              maximum = T,
              tol = tolerance)
          purity_max <- optimum$maximum
          purity_objective <- optimum$objective
          purity_lower <- uniroot(function(z){purity_pdf(z) - optimum$objective / 2}, c(0, optimum$maximum), tol = tolerance)$root
          purity_upper <- tryCatch(
            ### if the purity is very high, an upper bound on the purity equal to one leads to an error
            uniroot(function(z){purity_pdf(z) - optimum$objective / 2}, c(optimum$maximum, 1), tol = tolerance)$root,
            error = function(e) {
              return(1)
            })
          list(purity_max = purity_max, purity_lower = purity_lower, purity_upper = purity_upper)
          #          list(purity_max = purity_max, purity_objective = purity_objective)
        },
        .(mu, sigma, nu, tau)]

    maf[, paste0("ccf_Mcopies_loP",
                 c("", "_lower", "_upper", "_prob95", "_prob90")) :=
          ccf.likelihood(purity,
                         tcn,
                         t_alt_count,
                         (t_alt_count + t_ref_count),
                         copies=(tcn-lcn)), by= 1:nrow(maf)]
    maf[, paste0("ccf_Mcopies_hiP",
                 c("", "_lower", "_upper", "_prob95", "_prob90")) :=
          ccf.likelihood(purity,
                         tcn,
                         t_alt_count,
                         (t_alt_count + t_ref_count),
                         copies=(tcn-lcn)), by= 1:nrow(maf)]

    maf[, paste0("ccf_1copy_loP",
                 c("", "_lower", "_upper", "_prob95", "_prob90")) :=
          ccf.likelihood(purity,
                         tcn,
                         t_alt_count,
                         (t_alt_count + t_ref_count),
                         copies=1), by= 1:nrow(maf)]
    maf[, paste0("ccf_1copy_hiP",
                 c("", "_lower", "_upper", "_prob95", "_prob90")) :=
          ccf.likelihood(purity,
                         tcn,
                         t_alt_count,
                         (t_alt_count + t_ref_count),
                         copies=1), by= 1:nrow(maf)]


  }
  maf
}


################################################################################################################################
################################################################################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))


if(!interactive()){

  parser=ArgumentParser()
  parser$add_argument('-m','--maf', type='character', help='file name of maf file to be annotated.')
  parser$add_argument('-f','--cncf_files', type='character',
                      help='Mapping of "Tumor_Sample_Barcode" from maf and "cncf_filename" from FACETS (tab-delimited with header)')
  parser$add_argument('-o','--out_maf', type='character', help='file name of CN annotated maf.')
  parser$add_argument('-p','--purity_fit_filename', type='character', default = "", help='file name of CN annotated maf.')
  parser$add_argument('-c','--save_comments', action = 'store_true', default = FALSE, help = 'Retains comments in file header')
  args=parser$parse_args()

  maf_file = args$maf
  cncf_file_list = args$cncf_files
  output_maf_file = args$out_maf
  purity_fit_filename = args$purity_fit_filename
  save_comments = args$save_comments

  maf = fread(paste0('grep -v "^#" ', maf_file))
  cncf_files = fread(cncf_file_list)

  #cncf_files = with(cncf_files_table, structure(cncf_filename, .Names = Tumor_Sample_Barcode))

  maf = main(maf, cncf_files, purity_fit_filename = purity_fit_filename)

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

