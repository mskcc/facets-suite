#' Estimate CCFs of somatic mutations
#'
#' Based on FACETS data, infer cancer-cell fraction (CCF) for somatic mutations in a sample.
#'
#' @param maf Input MAF file.
#' @param segs FACETS segmentation output.
#' @param purity Sample purity estimate.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#' 
#' @importFrom data.table setDT foverlaps :=
#' @importFrom stats dbinom
#'
#' @return MAF file annotated with clonality estimates for each mutation, where the following column prefixes are used:
#' \itemize{
#'   \item{\code{ccf_Mcopies*}:} {Inferred CCF if mutation is on the major allele.}
#'   \item{\code{ccf_1copy*}:} {Inferred CCF if mutation exists in one copy.}
#'   \item{\code{ccf_expected_copies*}:} {Inferred CCF if mutation exists in number of copies expected from observed VAF and local ploidy.}
#' }
#' 
#' @export
ccf_annotate_maf = function(maf,
                            segs,
                            purity,
                            algorithm = c('em', 'cncf')) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    if (!all(c('t_alt_count', 't_ref_count') %in% names(maf))) {
        stop('Columns t_alt_count and t_depth required.')
    }
    
    # Find segments overlapping mutation loci
    data.table::setDT(maf, key = c('Chromosome', 'Start_Position', 'End_Position'))
    original_cols = names(maf)
    maf <- maf[, Chromosome:=as.character(Chromosome)]
    
    segs$chrom[segs$chrom == 23] = 'X'
    if (algorithm == 'em') {
        segs$tcn = segs$tcn.em
        segs$lcn = segs$lcn.em
    }
    segs = segs[, c('chrom', 'start', 'end', 'tcn', 'lcn', 'cf')]
    data.table::setDT(segs, key = c('chrom', 'start', 'end'))
    
    maf = data.table::foverlaps(maf, segs,
                                by.x = c('Chromosome', 'Start_Position', 'End_Position'),
                                by.y = c('chrom', 'start', 'end'),
                                type = 'within', mult = 'first', nomatch = NA)
    maf = maf[, c(original_cols, 'tcn', 'lcn'), with = FALSE]
    
    # Calculate CCFs
    maf[, `:=` (
        purity = purity,
        t_alt_count = as.numeric(t_alt_count),
        t_ref_count = as.numeric(t_ref_count)
    )]
    maf[, t_depth := t_alt_count + t_ref_count]
    maf[, t_var_freq := t_alt_count / t_depth]
    maf[, expected_alt_copies := expected_mutant_copies(t_var_freq, tcn, purity), by = seq_len(nrow(maf))]
    maf[, c('ccf_Mcopies', 'ccf_Mcopies_lower', 'ccf_Mcopies_upper', 'ccf_Mcopies_prob95', 'ccf_Mcopies_prob90') :=
            estimate_ccf(purity, tcn, tcn - lcn, t_alt_count, t_depth), by = seq_len(nrow(maf))]
    maf[, c('ccf_1copy', 'ccf_1copy_lower', 'ccf_1copy_upper', 'ccf_1copy_prob95', 'ccf_1copy_prob90') :=
            estimate_ccf(purity, tcn, tcn - lcn, t_alt_count, t_depth), by = seq_len(nrow(maf))]
    maf[, c('ccf_expected_copies', 'ccf_expected_copies_lower', 'ccf_expected_copies_upper',
            'ccf_expected_copies_prob95', 'ccf_expected_copies_prob90') :=
            estimate_ccf(purity, tcn, tcn - lcn, t_alt_count, t_depth), by = seq_len(nrow(maf))]
    as.data.frame(maf)
}

#' Estimate CCFs of somatic mutations from cncf.txt file. 
#'
#' Based on FACETS data, infer cancer-cell fraction (CCF) for somatic mutations in a sample.
#'
#' @param maf Input MAF file.
#' @param cncf_txt_file .cncf.txt file created with legacy output of facets-suite.
#' @param purity Sample purity estimate.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#' 
#' @importFrom data.table setDT foverlaps :=
#' @importFrom stats dbinom
#'
#' @return MAF file annotated with clonality estimates for each mutation, where the following column prefixes are used:
#' \itemize{
#'   \item{\code{ccf_Mcopies*}:} {Inferred CCF if mutation is on the major allele.}
#'   \item{\code{ccf_1copy*}:} {Inferred CCF if mutation exists in one copy.}
#'   \item{\code{ccf_expected_copies*}:} {Inferred CCF if mutation exists in number of copies expected from observed VAF and local ploidy.}
#' }
#' 
#' @export
ccf_annotate_maf_legacy <- function(maf, 
                                    cncf_txt_file,
                                    purity,
                                    algorithm = c('em', 'cncf')) {
    segs <-
        fread(cncf_txt_file) %>%
        select(chrom, seg, num.mark, nhet, cnlr.median,
               mafR, segclust, cnlr.median.clust, mafR.clust, 
               start = loc.start, end = loc.end, cf.em, tcn.em, lcn.em, cf, tcn, lcn)
    
    ccf_annotate_maf(maf, segs, purity, algorithm)    
}

# Estimate most likely CCF given observed VAF, purity and local ploidy
# Based on PMID 25877892
estimate_ccf = function(purity,
                        total_copies,
                        mutant_copies,
                        t_alt_count,
                        t_depth) {
    
    
    ccfs = seq(0.001, 1, 0.001)
    expected_vaf  = function(ccf, purity, total_copies) {
        purity * ccf * mutant_copies / (2 * (1 - purity) + purity * total_copies)
    }
    
    probs = sapply(ccfs, function(c) {
        stats::dbinom(t_alt_count, t_depth, expected_vaf(c, purity, total_copies))
    })
    probs = probs / sum(probs)
    
    ccf_max = which.max(probs)
    if (identical(ccf_max, integer(0))) ccf_max = NA
    ccf_half_max = which(probs > max(probs) / 2)
    ccf_lower = max(ccf_half_max[1] - 1, 1) # closest ccf value before half-max range (within 0-1 range)
    ccf_upper = min(ccf_half_max[length(ccf_half_max)] + 1, length(ccfs)) # closest ccf value after half-max range (within 0-1 range)
    if (is.na(purity)) ccf.upper = NA 
    ccf_max = ccf_max / length(ccfs)
    ccf_lower = ccf_lower / length(ccfs)
    ccf_upper = ccf_upper / length(ccfs)
    prob95 = sum(probs[950:1000])
    prob90 = sum(probs[900:1000])
    
    list(ccf_max, ccf_lower, ccf_upper, prob95, prob90)
}

# Estimate mutant copy number, given observed VAF, purity, and local ploidy
# Based on PMID 28270531
expected_mutant_copies = function(t_var_freq,
                                  total_copies,
                                  purity) {
    
    if (is.na(total_copies)) {
        NA_real_
    } else {
        if (total_copies == 0) total_copies = 1
        mu = t_var_freq * (1 / purity) * (purity * total_copies + (1 - purity) * 2)
        alt_copies = ifelse(mu < 1, 1, abs(mu)) # mu < 1 ~ 1, mu >= 1 ~ abs(mu)
        round(alt_copies)
    }
}
