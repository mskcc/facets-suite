#' Format .seg file
#'
#' Format Facets output for viewing in IGV.
#'
#' @param out Facets out object.
#' @param sample_id Sample name.
#' 
#' @importFrom dplyr group_by left_join summarize mutate select
#'
#' @return Segmentation output formatted for IGV.

#' @export
format_seg = function(facets_output, sample_id, normalize = TRUE) {
    
    seg = group_by(facets_output$snps, chrom, seg) %>% 
        summarize(loc.start = min(maploc),
                  loc.end = max(maploc)) %>% 
        left_join(., select(facets_output$segs, chrom, seg, num.mark, seg.mean = cnlr.median),
                  by = c('chrom', 'seg')) %>% 
        mutate(ID = sample_id) %>% 
        select(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
    
    if (normalize) {
        mutate(seg, seg.mean = seg.mean - facets_output$diplogr)
    } else {
        seg
    }
}
