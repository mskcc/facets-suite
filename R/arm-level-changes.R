#' Arm-level changes 
#' 
#' Get the altered chromosome arms in sample. Does not include the acrocentric p arms of chromosomes 12, 14, 15, 31, and 22.
#'  
#' @param segs FACETS segmentation output.
#' @param ploidy Sample ploidy.
#' @param genome Genome build.
#' @param algorithm Choice between FACETS \code{em} and \code{cncf} algorithm.
#'
#' @return List of items, containing:
#' @return \code{data.frame} for all genes mapping onto a segment in the output segmentation, with the columns:
#' \itemize{
#'     \item{\code{genome_doubled}:} {Boolean indicating whether sample genome is doubled.}
#'     \item{\code{fraction_cna}:} {Fraction of genome altered.}
#'     \item{\code{weighted_fraction_cna}:} {A weighted version of \code{fraction_cna} where only altered chromosomes are counted and weighted according to their length relative to total genome.}
#'     \item{\code{aneuploidy_scores}:} {Count of the number of altered arms, see source URL.}
#'     \item{\code{full_output}:} {Full per-arm copy-number status.}
#' }
#'
#' @importFrom dplyr left_join filter summarize select %>% mutate_at case_when group_by rowwise arrange
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
    
    genome_choice = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)
    
    # Get WGD status
    fcna_output = calculate_fraction_cna(segs, ploidy, genome, algorithm)
    wgd = fcna_output$genome_doubled
    sample_ploidy = ifelse(wgd, round(ploidy), 2)
    
    # Create chrom_info for sample
    sample_chrom_info = get_sample_genome(segs, genome_choice)
    segs = parse_segs(segs, algorithm) %>% 
        left_join(., select(sample_chrom_info, chr, centromere), by = c('chrom' = 'chr'))
    
    # Find altered arms
    # Split centromere-spanning segments
    # Remove segments where lcn is NA
    segs = filter(segs, !is.na(lcn)) %>% 
        rowwise() %>% 
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
    # Requires that >=80% exist at given copy-number state
    acro_arms = c('13p', '14p', '15p', '21p', '22p') # acrocentric chromsomes
    chrom_arms = setdiff(paste0(rep(unique(test_facets_output$segs$chrom), each = 2), c('p', 'q')), acro_arms)
    
    segs = group_by(segs, arm, tcn, lcn) %>% 
        summarize(cn_length = sum(length)) %>% 
        group_by(arm) %>% 
        mutate(arm_length = sum(cn_length),
               majority = cn_length >= 0.8 * arm_length,
               frac_of_arm = signif(cn_length/arm_length, 2),
               cn_state = mapvalues(paste(wgd, tcn-lcn, lcn, sep = ':'),
                                    copy_number_states$map_string, copy_number_states$call,
                                    warn_missing = FALSE)) %>% 
        ungroup() %>% 
        filter(majority == TRUE, arm %in% chrom_arms) %>% 
        select(-majority) %>% 
        mutate(arm = factor(arm, chrom_arms, ordered = T)) %>% 
        arrange(arm)
    
    altered_arms = filter(segs, cn_state != 'DIPLOID')
    
    # Weighted fraction copy-number altered
    frac_altered_w = select(sample_chrom_info, chr, p = plength, q = qlength) %>%
        gather(arm, length, -chr) %>%
        filter(paste0(chr, arm) %in% chrom_arms) %>%
        summarize(sum(length[paste0(chr, arm) %in% altered_arms$arm]) / sum(length)) %>%
        as.numeric()
    
    list(
        genome_doubled = fcna_output$genome_doubled,
        fraction_cna = fcna_output$fraction_cna,
        weighted_fraction_cna = frac_altered_w,
        aneuploidy_score = length(altered_arms),
        full_output = segs
    )
}

#' @export 
arm_level_changes_dmp = function(segs,
                                 ploidy,
                                 genome = c('hg19', 'hg18', 'hg38'),
                                 algorithm = c('em', 'cncf'), 
                                 diplogr) {
  
    genome_choice = get(match.arg(genome, c('hg19', 'hg18', 'hg38'), several.ok = FALSE))
    algorithm = match.arg(algorithm, c('em', 'cncf'), several.ok = FALSE)

    # Get WGD status
    fcna_output = calculate_fraction_cna(segs, ploidy, genome, algorithm)
    wgd = fcna_output$genome_doubled
    sample_ploidy = ifelse(wgd, round(ploidy), 2)
    
    
    sample_chrom_info = get_sample_genome(segs, genome_choice)

    segs = parse_segs(segs, algorithm) %>%
        left_join(., select(sample_chrom_info, chr, centromere), by = c('chrom' = 'chr'))
    
  

     segs = segs %>% 
         rowwise() %>% 
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
               length = end - start,
               cnlr.adj = cnlr.median - diplogr,
               phet = nhet / num.mark)
     
     arm_lengths = segs %>% 
         group_by(arm) %>% 
         summarize(arm_length = sum(length))
     
     segs = left_join(segs, arm_lengths)
     
    # Find distinct copy-number states 
    # Requires that >=50% exist at given copy-number state
    acro_arms = c('13p', '14p', '15p', '21p', '22p') # acrocentric chromsomes
    chrom_arms = setdiff(paste0(rep(unique(test_facets_output$segs$chrom), each = 2), c('p', 'q')), acro_arms)
    

    segs = segs %>% 
        mutate(tcn = case_when(algorithm=="em" ~ tcn.em, TRUE ~ tcn),
               lcn = case_when(algorithm=="em" ~ lcn.em, TRUE ~ lcn),
               cf = case_when(algorithm=="em" ~ cf.em, TRUE ~ cf))
    
    #annotate sample sex based on proportion of heterozygous SNPs on chrX
   chrx =  segs %>% filter(chrom==23) %>% group_by(chrom) %>% summarize(prop_het = max(phet))
   sex = case_when(chrx$prop_het > 0.01 ~ "Female", TRUE ~"Male")
   
   #correct NAs for high confidence CNLOH
   #theoretical values cnloh
   phis = seq(0, 0.9, by = 0.01)
   cnlr = function(phi, m = 0, p = 1) {
         log2((2 * (1 - phi) + (m + p) * phi) / 2)
     }

     valor = function(phi, m = 0, p = 1) {
           abs(log((m * phi + 1 - phi) / (p * phi + 1 - phi)))
       }
   cnloh_line = data.frame(
        phi = phis,
         cnlr = sapply(
               phis,
               function(phi){cnlr(phi, m = 2, p = 0)}
           ),
         valor = sapply(
               phis,
               function(phi){valor(phi, m = 2, p = 0)}
           ))
   
   hetloss_line = data.frame(
       phi = phis,
       cnlr = sapply(
           phis,
           function(phi){cnlr(phi, m = 0, p = 1)}
       ),
       valor = sapply(
           phis,
           function(phi){valor(phi, m = 0, p = 1)}
       ))
   
   #calculate estimated cellular fraction with obvious CNLOH misses
   segs = segs %>%
       mutate(
        #correct CNLOH to het loss for obvious errors
        tcn = case_when(
            tcn==2 & cnlr.adj < -0.5 & mafR>0.5 & lcn==0 & chrom !=23 ~ 1, 
            TRUE ~ as.numeric(tcn)),
        
        cf = case_when(
            is.na(lcn) & tcn ==2 & mafR>1 & nhet >= 5  & abs(cnlr.adj)<0.2 ~ cnloh_line[findInterval(abs(mafR), cnloh_line$valor),]$phi, 
            lcn==1 & tcn==2 & mafR>1 & abs(cnlr.adj)<0.2  & nhet>=5 ~ cnloh_line[findInterval(abs(mafR), cnloh_line$valor),]$phi,
            tcn==2 & cnlr.adj < -0.5 & mafR>0.5 & lcn==0 & chrom !=23 ~ hetloss_line[findInterval(abs(mafR), hetloss_line$valor),]$phi,
            TRUE ~ as.numeric(cf)),
        
       lcn = case_when(
            is.na(lcn) & tcn ==2 & mafR>1 & nhet>=5 & abs(cnlr.adj)<0.2 ~ 0, 
            lcn==1 & tcn ==2 & mafR>1 & abs(cnlr.adj) <0.2 & nhet>=5  ~ 0,
       TRUE ~ as.numeric(lcn)),
       arm_fraction = length / arm_length
   )

   
    #annotate cn_state
    segs = segs %>% 
        mutate(
            cn_state = mapvalues(paste(wgd, tcn-lcn, lcn, sep = ':'),
                       copy_number_states$map_string, 
                       copy_number_states$call,
                       warn_missing = FALSE),
               
            cn_state_num = mapvalues(paste(wgd, tcn-lcn, lcn, sep = ':'),
                           copy_number_states$map_string, 
                           copy_number_states$numeric_call,
                           warn_missing = FALSE)
               )
    
    #correct X for patient sex
    segs = segs %>%
        mutate(cn_state = case_when(chrom==23 & sex=="Male" & tcn==1 ~ "DIPLOID", TRUE ~cn_state),
               cn_state = case_when(chrom==23 & sex=="Male" & tcn > 1 ~ "GAIN", TRUE ~cn_state),
               cn_state_num =  case_when(chrom==23 & sex=="Male" & tcn==1 ~ 0, TRUE ~ as.numeric(cn_state_num)),
               cn_state_num =  case_when(chrom==23 & sex=="Male" & tcn > 1 ~ 1, TRUE ~ as.numeric(cn_state_num)),
               lcn = case_when(chrom==23 & sex=="Male" ~ 0, TRUE ~ as.numeric(lcn))
        )
    
    segs.full = segs %>%
        mutate(chrom = gsub('23', 'X', chrom), 
               arm = gsub('23', 'X', arm))
    
     # val_arms = c("3q_Loss", "5q_Loss", "7q_Loss",
     #              "8p_Gain", "8q_Gain", "11q_Gain", "11q_Loss",
     #              "12p_Loss", "12p_Gain", "12q_Gain","13q_Loss",
     #              "17p_Loss", "19p_Gain", "19q_Gain", "20q_Loss")
    
    segs.full = segs.full %>%
        mutate(Class = case_when(cn_state=="DIPLOID" ~ 'Diploid', #exclude 1,0 males from being called loss
                                 chrom==23 & sex=="Male" & cn_state =='GAIN' ~ 'Gain',
                                 tcn==2 & lcn==0  ~ 'Loss', #CNLOH
                                 # **gain/loss definitions do not account for WGD, e.g. 3,0 | 4,0 states**
                                 tcn>2 ~ 'Gain', 
                                 lcn<2 ~ 'Loss',
                                 ),
               arm_change = paste(arm, Class, sep = "_"))
    
    maxarm = segs.full %>% 
        #filter(arm_change %in% val_arms) %>% 
        group_by(arm) %>%
        summarize(max_arm_len= max(length)) %>%
        mutate(key = paste(arm, max_arm_len, sep = "_")) 
    
    chrom_arms = gsub('23', 'X', chrom_arms)
    
    segs.filt = segs.full %>% 
        mutate(key = paste(arm, length, sep = "_")) %>% 
        filter(key %in% maxarm$key, arm %in% chrom_arms) %>% 
      dplyr::select(arm,tcn, lcn, cf, arm_fraction, cn_state, Class)
    
    list(
        full_output = segs.full,
        filtered_output = segs.filt
    )
    
}

