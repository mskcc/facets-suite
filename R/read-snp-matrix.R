#' Read counts file
#'
#' Reads a gzipped output file from \code{snp-pileup}.
#'
#' @param input_file Path to input file.
#' @param err.thresh Threshold for errors at locus.
#' @param del.thresh Threshold for deletions at locus.
#' 
#' @source \code{snp-pileup is part of} \url{www.github.com/mskcc/facets}
#'
#' @return Count matrix.
#' 
#' @importFrom data.table fread :=

#' @export
read_snp_matrix = function(input_file,
                           err.thresh = 10,
                           del.thresh = 10) {
    
    read_counts = data.table::fread(cmd = paste('gunzip -c', input_file), key = c('Chromosome', 'Position'))
    read_counts = read_counts[File1E <= err.thresh & File2E <= err.thresh &
                              File1D <= del.thresh & File2D <= del.thresh &
                              !Chromosome %in% c('MT', 'chrM', 'Y', 'chrY')]
    
    read_counts[, `:=`(
        NOR.DP = File1R + File1A,
        TUM.DP = File2R + File2A,
        NOR.RD = File1R,
        TUM.RD = File2R,
        Chromosome = gsub('chr', '', Chromosome)
        )][, ('Chromosome') := factor(get('Chromosome'), levels = c(1:22, 'X'))]

    read_counts[order(Chromosome, Position)][, list(Chromosome, Position, NOR.DP, TUM.DP, NOR.RD, TUM.RD)]
}
