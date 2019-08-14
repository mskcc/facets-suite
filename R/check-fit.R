#' Sample QC
#' 
#' Generate QC metrics for sample. Can use mutation calls in MAF file.
#'
#' @param facets_output Output from \code{run_facets}.
#' @param maf A MAF file with somatic mutations from sample (optional).
#' @param algorithm Choose assessing the fit from the \code{em} or \code{cncf} algorithm.
#'
#' @importFrom dplyr distinct
#' @import data.table

check_fit = function(facets_output,
                     maf = NULL,
                     genome = c('hg19', 'hg18', 'hg38'),
                     algorithm = c('em', 'cncf')) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = F)
    genome_choice = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = F))
    
    # Set variables
    segs = facets_output$segs
    snps = facets_output$snps
    diplogr = facets_output$diplogr
    alballogr = as.numeric(facets_output$alballogr[, 1])
    mafr_thresh = facets_output$mafr_thresh
    fcna_output = calculate_fraction_cna(segs, facets_output$ploidy, genome, algorithm)
    wgd = fcna_output$genome_doubled
    
    segs = parse_segs(segs, algorithm)
    segs = as.data.table(segs)
    setkey(segs, chrom, start, end)
    snps = as.data.table(snps)
    setkey(snps, chrom, maploc)

    # Check diplogr
    diplogr_flag = abs(facets_output$diplogr) > 1
    
    # Check alternative diplogr
    alt_diplogr_flag = length(setdiff(alballogr, diplogr) > 0)
    
    # Check for mutations in homdels
    if (!is.null(maf)) {
        
    }
    
    # Balanced diploid regions in WGD case
    if (wgd) {
        balanced_segs = segs[cnlr.median.clust == diplogr & mafR.clust < mafr_thresh]
        balanced_segs = paste0(unique(balanced_segs$chrom), collapse = ',')
    } else {
        balanced_segs = NA_character_
    } 
    
    # Number of amplifications and deletions
    n_amps = nrow(segs[tcn.em >= 10])
    n_dels = nrow(segs[tcn.em == 0])
    
    # Number of unique copy-number states
    n_cn_states = nrow(distinct(segs, tcn, lcn))
    n_segs = nrow(segs)
    
    # Number of NA lcns
    n_lcn_na = nrow(segs[is.na(lcn)])
    
    # SNP couns
    n_snps = nrow(snps)
    n_het_snps = nrow(snps[het == 1])
    frac_het_snps = n_het_snps/n_snps
    
    # Output
    list(
        diplogr_flag = diplogr_flag,
        balanced_segs = balanced_segs,
        amplifications = n_amps,
        deletions = n_dels,
        distinct_states = n_cn_states,
        segments = n_segs,
        lcn_is_na = n_lcn_na,
        snps = n_snps,
        het_snps = n_het_snps,
        fraction_het_snps = signif(frac_het_snps, 2)
    )
}
    
