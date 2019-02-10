#' Arm-level changes 
#' 
#' Get the altered chromosome arms in sample. Does not include the acrocentric p arms of chromosomes 12, 14, 15, 31, and 22.
#'  
#' @param segs FACETS segmentation output.
#' @param ploidy Sample ploidy.
#' @param genome Genome build.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#'
#' @return List with one or more values from function.
#' 
#' @importFrom dplyr left_join filter summarize select %>%
#' @importFrom purrr map_dfr map_lgl discard
#' @importFrom tidyr gather
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/29622463}

#' @export 
arm_level_changes = function(segs,
                             ploidy,
                             genome = c('hg19', 'hg18', 'hg38'),
                             algorithm = c('em', 'cncf')) {
    
    genome = match.arg(genome)
    algorithm = match.arg(algorithm)
    
    # Get WGD status
    fcna_output = calculate_fcna(segs, ploidy, genome, algorithm)
    
    # Centromere locations
    genome = get(genome)
    
    # Create chrom_info for sample
    segs = parse_segs(segs, algorithm)
    sample_chrom_info = get_sample_genome(segs, genome)

    # Altered arms
    segs = left_join(segs, select(sample_chrom_info, chr, centromere), by = c('chrom' = 'chr'))
    altered_arms = sapply(unique(segs$chrom), function(x) {
        seg_p = segs[which(segs$chrom == x & segs$start < segs$centromere), ]
        seg_p$end = ifelse(seg_p$end > seg_p$centromere, seg_p$centromere, seg_p$end)
        seg_p$length = seg_p$end - seg_p$start
        seg_q = segs[which(segs$chrom == x & segs$end > segs$centromere), ]
        seg_q$start = ifelse(seg_q$start < seg_q$centromere, seg_q$centromere, seg_q$start)
        seg_q$length = seg_q$end - seg_q$start
        
        if (!fcna_output$wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn == 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn == 1)])
        } else if (fcna_output$wgd) {
            seg_p_unaltered = sum(seg_p$length[which(seg_p$tcn == sample_ploidy & seg_p$lcn >= 1)])
            seg_q_unaltered = sum(seg_q$length[which(seg_q$tcn == sample_ploidy & seg_q$lcn >= 1)])
        }
        paste0(paste0(x, c('p', 'q'))[c((sample_chrom_info$plength[sample_chrom_info$chr == x] - 
                                             seg_p_unaltered) / sample_chrom_info$plength[sample_chrom_info$chr == x] > .8,
                                        (sample_chrom_info$qlength[sample_chrom_info$chr == x] -
                                             seg_q_unaltered) / sample_chrom_info$qlength[sample_chrom_info$chr == x] > .8) %>%
                                          map_lgl(~ifelse(is.na(.), FALSE, .))], collapse = ',')
    }
    ) %>% strsplit(., ',') %>%
        unlist %>%
        discard(. %in% c('', '13p', '14p', '15p', '21p', '22p'))
    
    # Weighted fraction copy-number altered
    frac_altered_w = select(sample_chrom_info, chr, p = plength, q = qlength) %>%
        gather(arm, length, -chr) %>%
        filter(!paste0(chr, arm) %in% c('13p', '14p', '15p', '21p', '22p')) %>%
        summarize(sum(length[paste0(chr, arm) %in% altered_arms]) / sum(length)) %>%
        as.numeric()
    
    list(
        genome_doubled = fcna_output$wgd,
        fcna = fcna_output$fcna,
        weighted_fcna = fcna_output$weighted_fcna,
        aneuploidy_score = length(altered_arms),
        altered_arms = paste0(altered_arms, collapse = ',')
    )
}
