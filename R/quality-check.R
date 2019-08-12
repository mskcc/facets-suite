#' Sample QC
#' 
#' Generate QC metrics for sample. Can use mutation calls in MAF file.
#'
#' @param facets_output 
#' @param maf 
#' @param algorithm 
#'
#' @importFrom dplyr distinct

quality_check = function(facets_output,
                         maf = NULL,
                         genome = c('hg19', 'hg18', 'hg38'),
                         algorithm = c('em', 'cncf')) {
    
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = F)
    genome = match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = F)
    
    fcna_output = calculate_fraction_cna(facets_output$segs, facets_output$ploidy, genome, algorithm)
    
    genome = get(genome)
    
    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm)
    
    # Check diplogr
    diplogr_flag = abs(facets_output$diplogr) > 1
    
    
    # Check for mutations in homdels
    if (!is.null(maf)) {
        
    } else {
        
    }
    
    
    # Balanced diploid regions in WGD case
    if (wgd) {
        balanced_segs = facets_output$segs[which(facets_output$segs$cnlr.median.clust == facets_output$diplogr &
                                                facets_output$segs$mafR.clust < facets_output$mafr_thresh), ]
        balanced_segs = paste0(unique(balanced_segs$chrom), collapse = ',')
    } else {
        balanced_segs = NA_character_
    } 
    
    # Number of amplifications and deletions
    n_amps = nrow(segs[, tcn >= 10])
    n_dels = nrow(segs[, tcn == 10])
    
    # Number of unique copy-number states
    n_cn_states = nrow(distinct(segs, tcn, lcn))
    n_segs = nrow(segs)
    
    # Number of NA lcns
    n_lcn_na = nrow(segs[is.na(lcn), ])
    
    # SNP couns
    n_snps = nrow(facets_output$snps)
    n_het_snps = nrow(facets_output$snps[het == 1, ])
    frac_het_snps = n_het_snps/n_snps
    
    # Output
    list(
        diplogr_flag = diplogr_flag,
        balanced_segs = balanced_segs,
        amplifications = n_amp,
        deletions = n_dels,
        distinct_states = n_cn_states,
        segments = n_segs,
        lcn_is_na = n_lcn_na,
        snps = n_snps,
        het_snps = n_het_snps,
        fraction_het_snps = frac_het_snps
    )
}
    
