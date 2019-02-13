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
#' @importFrom dplyr left_join filter summarize select %>% mutate_at case_when group_by
#' @importFrom purrr map_dfr map_lgl map_chr discard
#' @importFrom tidyr gather separate_rows
#' @importFrom plyr mapvalues
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/29622463}

#' @export 
arm_level_changes = function(segs,
                             ploidy,
                             genome = c('hg19', 'hg18', 'hg38'),
                             algorithm = c('em', 'cncf')) {
    
    genome = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Get WGD status
    fcna_output = calculate_fraction_cna(segs, ploidy, genome, algorithm)
    wgd = fcna_output$genome_doubled
    sample_ploidy = ifelse(wgd, round(ploidy), 2)
    
    # Create chrom_info for sample
    sample_chrom_info = get_sample_genome(segs, genome)
    segs = parse_segs(segs, algorithm) %>% 
        left_join(., select(sample_chrom_info, chr, centromere), by = c('chrom' = 'chr'))
    
    # Find altered arms
    # Split centromere-spanning segments
    segs = rowwise(segs) %>% 
        mutate(
            arm = case_when(
                start < centromere & end <= centromere ~ 'p',
                start >= centromere ~ 'q',
                TRUE ~ 'span'),
            start = ifelse(arm == 'span', paste(c(start, centromere), collapse = ','), as.character(start)),
            end = ifelse(arm == 'span', paste(c(centromere, end), collapse = ','), as.character(end))
        ) %>% 
        separate_rows(start, end, sep = ',') %>% 
        mutate(start = as.numeric(start),
               end = as.numeric(end),
               arm = case_when(
                   start < centromere & end <= centromere ~ paste0(chrom, 'p'),
                   start >= centromere ~ paste0(chrom, 'q')),
               length = end - start)
    
    # Find distinct copy-number states 
    copy_number_states = mutate(copy_number_states, map_string = paste(wgd, mcn, lcn, sep = ':'))
    altered_arms = group_by(segs, arm, tcn, lcn) %>% 
        summarize(cn_length = sum(length)) %>% 
        group_by(arm) %>% 
        mutate(arm_length = sum(cn_length),
               majority = cn_length >= 0.8 * arm_length,
               cn_state = mapvalues(paste(wgd, tcn-lcn, lcn, sep = ':'),
                                    copy_number_states$map_string, copy_number_states$call,
                                    warn_missing = FALSE)) %>% 
        filter(majority == TRUE,
               cn_state != 'DIPLOID',
               !arm %in% c('13p', '14p', '15p', '21p', '22p')) # acrocentric chromsomes

    # Weighted fraction copy-number altered
    frac_altered_w = select(sample_chrom_info, chr, p = plength, q = qlength) %>%
        gather(arm, length, -chr) %>%
        filter(!paste0(chr, arm) %in% c('13p', '14p', '15p', '21p', '22p')) %>%
        summarize(sum(length[paste0(chr, arm) %in% altered_arms$arm]) / sum(length)) %>%
        as.numeric()
    
    list(
        genome_doubled = fcna_output$genome_doubled,
        fcna = fcna_output$fraction_cna,
        weighted_fcna = frac_altered_w,
        aneuploidy_score = length(altered_arms),
        altered_arms = paste(paste0(altered_arms$arm, ':', altered_arms$cn_state), collapse = ',')
    )
}
