#' Format `.seg` file
#'
#' Format Facets output for viewing in IGV.
#'
#' @param facets_data Facets output, from `run_facets`.
#' @param sample_id Sample name.
#' @param normalize Adjust copy-number log-ratio by dipLogR.
#' 
#' @importFrom dplyr group_by left_join summarize mutate select
#'
#' @return Segmentation output formatted for IGV.

#' @export
format_igv_seg = function(facets_data, sample_id, normalize = TRUE) {
    
    seg = group_by(facets_data$snps, chrom, seg) %>% 
        summarize(loc.start = min(maploc),
                  loc.end = max(maploc)) %>% 
        left_join(., select(facets_data$segs, chrom, seg, num.mark, seg.mean = cnlr.median),
                  by = c('chrom', 'seg')) %>% 
        mutate(ID = sample_id) %>% 
        select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
    
    if (normalize) {
        mutate(seg, seg.mean = seg.mean - facets_data$diplogr)
    } else {
        seg
    }
}
