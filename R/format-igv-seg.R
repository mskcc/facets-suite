#' Format `.seg` file
#'
#' Format Facets output for viewing in IGV.
#'
#' @param facets_output Facets output, from `run_facets`.
#' @param sample_id Sample name.
#' @param normalize Adjust copy-number log-ratio by dipLogR.
#' 
#' @importFrom dplyr group_by left_join summarize mutate select ungroup
#'
#' @return Segmentation output formatted for IGV.

#' @export
format_igv_seg = function(facets_output, sample_id, normalize = TRUE) {
    
    if (!all(c('snps', 'segs', 'diplogr') %in% names(facets_output))) {
        stop(paste0('Input is missing segs, snps or diplogr ojbect.'), call. = FALSE)
    }
    
    seg = group_by(facets_output$snps, chrom, seg) %>% 
        summarize(loc.start = min(maploc),
                  loc.end = max(maploc)) %>% 
        ungroup() %>% 
        left_join(., select(facets_output$segs, chrom, seg, num.mark, seg.mean = cnlr.median),
                  by = c('chrom', 'seg')) %>% 
        mutate(ID = sample_id) %>% 
        select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
    
    if (normalize) { seg = mutate(seg, seg.mean = seg.mean - facets_output$diplogr) }
    data.frame(seg)
}
