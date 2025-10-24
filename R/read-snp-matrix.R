
read_snp_matrix = function(input_file,
                           err.thresh = 10,
                           del.thresh = 10) {
    #' Reads a gzipped output file from \code{snp-pileup}.
    #'
    #' @param input_file Path to input file.
    #' @param err.thresh Threshold for errors at locus.
    #' @param del.thresh Threshold for deletions at locus.
    #' 
    #' @source \code{snp-pileup} is part of \href{www.github.com/mskcc/facets}{FACETS}.
    #'
    #' @return A SNP count matrix, with the following columns:
    #' \itemize{
    #'       \item{\code{Chromosome} and \code{Position}:} {Chromosomal position of SNP.}
    #'       \item{\code{NOR.DP}:} {Total read depth in normal sample.}
    #'       \item{\code{TUM.DP}:} {Total read depth in tumor sample.}
    #'       \item{\code{NOR.RD}:} {Reference allele read depth in normal sample.}
    #'       \item{\code{TUM.RD}:} {Reference allele read depth in tumor sample.}
    #' }
    #' 
    #' @import data.table
    #' @export
    
    read_counts = data.table::fread(cmd = paste('gunzip -c', input_file), key = c('Chromosome', 'Position'))
    
    if (nrow(read_counts) == 0) { # necessary since fread command doesn't throw errors for certain cases
        stop(paste(input_file, 'does not exist or cannot be read properly.'), call. = F)
    }
    
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

read_snp_matrix_facets2n = function(filename, err.thresh=Inf, del.thresh=Inf, perl.pileup=FALSE,
                                    MandUnormal=FALSE, spanT=0.2, spanA=0.2, spanX=0.2, gbuild="hg19",
                                    ReferencePileupFile=NULL, ReferenceLoessFile=NULL, MinOverlap=0.9, useMatchedX=FALSE, refX=FALSE, unmatched=FALSE, donorCounts=FALSE){
    #' #' Read in the snp-pileup generated SNP read count matrix file for facets2n processing
    #' @importFrom utils read.csv
    #' @param filename counts file from snp-pileup
    #' @param skip (character) Skip n number of lines in the input file.
    #' @param err.thresh (numeric) Error threshold to be used to filter snp-pileup data frame.
    #' @param del.thresh (numeric) Deletion threshold to be used to filter snp-pileup data frame.
    #' @param perl.pileup (logical) Is the pileup data generated using perl pileup tool?
    #' @param MandUnormal (logical) Is CNLR analysis to be peformed using unmatched reference normals?
    #' @param unmatched (logical) is the tumor being analyzed unmatched
    #' @param spanT (numeric) Default span value to be used for loess normalization in tumor sample.
    #' @param spanA (numeric) Default span value to be used for loess normalization across autosomal chromosomes in the normal sample.
    #' @param spanX (numeric) Default span value to be used for loess normalization in Chr X in the normal sample.
    #' @param gbuild (character) Genome build (Default: hg19).
    #' @param ReferencePileupFile (character) Filepath to an optional snp-pileup generated pileup data of one or more reference normals.
    #' @param ReferenceLoessFile (character) Filepath to an optional loess data, generated using the facets2n package, of one or more reference normals. The number of normals in this data should match that in the ReferencePileupFile, and should be in the same order.
    #' @param MinOverlap (numeric) Mininum overlap fraction of loci between a tumor pileup and reference pileup data.
    #' @param useMatchedX (logical) Is the matched normal to be used for ChrX normalization?
    #' @param refX (logical) Use sex matched reference normal for chrX normalization
    #' @param donorCounts snp read count matrix for donor sample(s). Required columns: Chromosome Position Ref Alt and for each donor sample,i: RefDonoriR RefDonoriA RefDonoriE RefDonoriD RefDonoriDP
    #' @return A dataframe of pileup depth values for Tumor and Matched Normal if MandUnormal is FALSE. Else, a list of data frame with pileup depth values of Tumor, matched Normal, and a best unmatched normal, and the associated span values.
    #' @source \code{snp-pileup} is part of \href{www.github.com/mskcc/facets}{FACETS}.
    #' @import facets2n
    #' @importFrom pctGCdata getGCpct
    #' @export
    if (MandUnormal){
        facets2n::readSnpMatrix(filename, MandUnormal = MandUnormal, ReferencePileupFile = ReferencePileupFile, ReferenceLoessFile = ReferenceLoessFile, useMatchedX = useMatchedX, refX = refX, unmatched=unmatched, donorCounts = donorCounts)
    }else{
        facets2n::readSnpMatrix(filename, MandUnormal = MandUnormal, unmatched=unmatched, donorCounts = donorCounts)
    }
}
