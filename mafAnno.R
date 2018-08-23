#!/usr/bin/env Rscript

# Functions ------------------------------------------------------------------------------------------------------------
integer_cn_table = function(diplogr, fit) {

    df = fit$cncf
    if (!'loc.start' %in% names(df) | !'loc.end' %in% names(df)) { # this should capture all differently formatted outputs
        if (all(c('start', 'end') %in% names(fit$cncf))) {
                df$loc.start = df$start
                df$loc.end = df$end
        } else {
                df$loc.start = fit$start
                df$loc.end = fit$end
        }
    }
    df$chrom[which(df$chrom == 23)] = 'X'
    df$chrom = factor(df$chrom)
    dt = data.table(df,
    #                 cf = fit$cncf$cf,
    #                 tcn = fit$cncf$tcn,
                    mcn = fit$cncf$tcn - fit$cncf$lcn,
    #                 lcn = fit$cncf$lcn,
    #                 cf.em = fit$cncf$cf.em,
    #                 tcn.em = fit$cncf$tcn.em,
                    mcn.em = fit$cncf$tcn.em - fit$cncf$lcn.em,
    #                 lcn.em = fit$cncf$lcn.em,
                    ploidy = fit$ploidy,
                    purity = fit$purity,
                    dipLogR = diplogr)

    setkey(dt, chrom, loc.start, loc.end)
    dt
}

annotate_maf_with_facets_cf_tcn_lcn = function(maf, diplogr, fit, iTumor_Sample_Barcode = NULL) {

    maf_cols = colnames(maf)
    maf$Chromosome = factor(maf$Chromosome)
    setkey(maf, Chromosome, Start_Position, End_Position)
    dt = integer_cn_table(diplogr, fit)

    # check for duplicate columns
    if (any(duplicated(names(maf)))) {
        warning("Duplicate columns removed from MAF file")
        maf[, which(duplicated(names(maf))) := NULL, with = F]
    }

    if (is.null(iTumor_Sample_Barcode)) {
        maf_ann = foverlaps(maf, dt, mult = "first", nomatch = NA)
    } else {
        maf_ann = foverlaps(maf[Tumor_Sample_Barcode == iTumor_Sample_Barcode], dt, mult = "first", nomatch = NA)
    }

    maf_ann = maf_ann[, c(maf_cols, 'dipLogR', 'cf', 'tcn', 'lcn', 'cf.em', 'tcn.em', 'lcn.em', 'purity', 'ploidy'), with=F]
    maf_ann$purity = fit$purity
    maf_ann
}

ccf.likelihood = function(purity, absCN, alt_allele, coverage, copies) { #From McGranahan_and_Swanton_2015

    CCFs = seq(0.001,1,0.001)
    vac.ccf  = function(CCF, purity, absCN){purity * CCF * copies / (2*(1 - purity) + purity * absCN)}
    probs = sapply(CCFs, function(c){dbinom(alt_allele, coverage, vac.ccf(c, purity, absCN))})
    probs = probs/sum(probs)

    ccf.max = which.max(probs)
    if (identical(ccf.max, integer(0))) ccf.max = NA
    ccf.gt.half.max = which(probs > max(probs)/2)
    ccf.lower = max(ccf.gt.half.max[1] - 1, 1) ### closest ccf value before half-max range (within 0-1 range)
    ccf.upper = min(ccf.gt.half.max[length(ccf.gt.half.max)] + 1, length(CCFs)) ### closest ccf value after half-max range (within 0-1 range)
    if(is.na(purity)){ ccf.upper = NA }
    ccf.max = ccf.max/length(CCFs)
    ccf.lower = ccf.lower/length(CCFs)
    ccf.upper = ccf.upper/length(CCFs)
    prob.95 = sum(probs[950:1000])
    prob.90 = sum(probs[900:1000])
    list(ccf.max,ccf.lower,ccf.upper,prob.95,prob.90)
}

main = function(maf, facets_files, file_type = 'Rdata'){

    maf = as.data.table(maf)
    maf_Tumor_Sample_Barcodes = unique(maf$Tumor_Sample_Barcode)

    not.in.maf = setdiff(names(facets_files), maf_Tumor_Sample_Barcodes)
    no.facets = setdiff(maf_Tumor_Sample_Barcodes, names(facets_files))

    no.facets.data = maf[maf$Tumor_Sample_Barcode %in% no.facets,]
    maf_Tumor_Sample_Barcodes = maf_Tumor_Sample_Barcodes[!maf_Tumor_Sample_Barcodes %in% no.facets]

    if (length(no.facets) > 0) write(paste('Missing facets data:', no.facets), stderr())
    if (length(not.in.maf) > 0) write(paste('Not in MAF:', not.in.maf), stderr())

    idi = intersect(names(facets_files), maf_Tumor_Sample_Barcodes)
    maf = maf[maf$Tumor_Sample_Barcode %in% idi]

    if (file_type == 'Rdata') {
        maf_list = mclapply(idi, function(x) {
            load(facets_files[x])
            annotate_maf_with_facets_cf_tcn_lcn(maf, out$dipLogR, fit, x)
        }, mc.cores = detectCores())
    } else if (file_type == 'cncf') {
        maf_list = mclapply(idi, function(x) { 
            fit = list()
            fit$cncf = fread(facets_files[x])
            out = system(paste('grep -P "Purity|dipLogR|Ploidy"',
                         gsub('.cncf.txt', '.out', facets_files[x])), intern = T)
            names(out) = str_extract(out, '(Purity|dipLogR|Ploidy)')
            out = lapply(out, function(x) str_extract(x, '(?<=Purity = |dipLogR = |Ploidy = )[0-9\\-\\.]+'))
            out$dipLogR = as.numeric(out$dipLogR)
            fit$purity = as.numeric(out$Purity)
            fit$ploidy = as.numeric(out$Ploidy)
            annotate_maf_with_facets_cf_tcn_lcn(maf, out$dipLogR, fit, x)
        }, mc.cores = detectCores())
    }

    if (length(no.facets)) { maf_list = c(maf_list, list(no.facets.data)) }
    maf = rbindlist(maf_list,fill = T)

    maf[, t_alt_count := as.numeric(t_alt_count)]
    maf[, t_ref_count := as.numeric(t_ref_count)]

    maf[, c("ccf_Mcopies", "ccf_Mcopies_lower", "ccf_Mcopies_upper", "ccf_Mcopies_prob95", "ccf_Mcopies_prob90") :=
            ccf.likelihood(purity, # CCF estimates if mutation on mcn, cncf algorithm
                                         tcn,
                                         t_alt_count,
                                         (t_alt_count + t_ref_count),
                                         copies=(tcn-lcn)),
        by = 1:nrow(maf)]

    maf[, c("ccf_Mcopies", "ccf_Mcopies_lower", "ccf_Mcopies_upper", "ccf_Mcopies_prob95", "ccf_Mcopies_prob90") :=
            ccf.likelihood(purity, # CCF estimates if mutation in 1 copy, cncf algorithm
                                         tcn,
                                         t_alt_count,
                                         (t_alt_count + t_ref_count),
                                         copies = 1),
        by= 1:nrow(maf)]

    maf[, c("ccf_Mcopies_em", "ccf_Mcopies_lower_em", "ccf_Mcopies_upper_em", "ccf_Mcopies_prob95_em", "ccf_Mcopies_prob90_em") :=
            ccf.likelihood(purity, # ditto for em algorithm, mcn copies
                                       tcn.em,
                                       t_alt_count,
                                       (t_alt_count + t_ref_count),
                                       copies=(tcn.em-lcn.em)),
        by= 1:nrow(maf)]

    maf[, c("ccf_Mcopies_em", "ccf_Mcopies_lower_em", "ccf_Mcopies_upper_em", "ccf_Mcopies_prob95_em", "ccf_Mcopies_prob90_em") :=
            ccf.likelihood(purity, # ditto for em algorithm, 1 copy
                                       tcn.em,
                                       t_alt_count,
                                       (t_alt_count + t_ref_count),
                                       copies=1),
        by= 1:nrow(maf)]

   maf
}

# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(data.table)
    library(argparse)
    library(stringr)
    library(parallel)
    })

if (!interactive()) {

    parser=ArgumentParser()
    parser$add_argument('-m', '--maf', type='character', help='file name of maf file to be annotated')
    parser$add_argument('-f', '--facets_files', type = 'character',
                        help = 'Mapping of "Tumor_Sample_Barcode" from maf to "Rdata_filename" or "cncf_filename" from FACETS (tab-delimited with header)')
    parser$add_argument('-o', '--out_maf', type = 'character', help = 'file name of CN annotated maf')
    parser$add_argument('-c', '--save_comments', action = 'store_true', default = TRUE, help = 'Retains comments in file header')
    args=parser$parse_args()

    maf_file = args$maf
    facets_samples_file = args$facets_files
    output_maf_file = args$out_maf
    save_comments = args$save_comments

    maf = fread(maf_file, skip = 'Hugo_Symbol')
    facets_samples = fread(facets_samples_file)

    if ('Rdata_filename' %in% names(facets_samples)) { # Rdata input
        facets_files = with(facets_samples, structure(Rdata_filename, .Names = Tumor_Sample_Barcode))
        maf = main(maf, facets_files, 'Rdata')
    } else if ('cncf_filename' %in% names(facets_samples)) { # cncf input, requires that *.out file named similarily
        facets_files = with(facets_samples, structure(cncf_filename, .Names = Tumor_Sample_Barcode))
        maf = main(maf, facets_files, 'cncf')
    }

    if (save_comments) {
        header = readLines(maf_file, n = 10)
        header = header[unlist(lapply(header, function(x) substr(x,1,1)=='#'))]
        out_file = file(output_maf_file, open = 'wt')
        if (length(header) > 0) for (i in 1:length(header)) cat(header[i], '\n', file = out_file, append = T)
        write.table(maf, file = out_file, quote = F, col.names = T, row.names = F, sep = "\t")
        close(out_file)
    } else {
        write.table(maf, file = output_maf_file,
                    quote = F, col.names = T, row.names = F, sep = "\t")
    }
}
