#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript
### Reads a counts.dat.gz FACETS input file and 
### returns counts.norm_normal_depth.dat.gz.
###
### The tumor depth and allele counts are normalized
### to remove the dependence of the log-ratio on
### normal depth, using a lowess fit.
###
### This is designed to remove the "waterfall" effect
### produced when degraded tumor DNA with small
### fragment size means that SNPs adjacent to target exons
### are consistently covered less well in the tumor
### than in the normal.
###
### Usage: 
### ./norm_normal_depth.R countsMerged.dat.gz | gzip -c > countsMerged.norm_normal_depth.dat.gz
### 
### Alex Penson <pensona@mskcc.org>
###


library(data.table)
library(argparse)


normalize_facets_depth <- function(args){
  countsMerged <- fread(paste0('gunzip --stdout ', args$file))
  countsMerged$Chrom <- factor(countsMerged$Chrom, levels=c(1:22, "X", "Y"))
  dt <- with(countsMerged[TUM.DP>0][order(NOR.DP)], 
             data.table(Chrom, 
                        Pos, 
                        TUM.DP,
                        NOR.DP,
                        M=log2(NOR.DP), 
                        A=log2(TUM.DP/NOR.DP)
             )
  )

  
  lx <- with(dt, lowess(M, A)) ### variables in "MA plot"
  dt[,A.lowess := lx$y]
  A.mean <- mean(dt$A, na.rm=T)
  dt[,A.norm := A - A.lowess + A.mean]
  dt[,norm.factor := NOR.DP*2^A.norm / TUM.DP]
  
  countsMerged_merge <- merge(countsMerged,
                              dt[,c("Chrom", "Pos", "norm.factor"), with=F],
                              by=c("Chrom", "Pos"),
                              all.x=T)
  for (i in 5:13){ ### loop over TUM.* columns
    set(countsMerged_merge, NULL, i, countsMerged_merge[[i]] * countsMerged_merge$norm.factor) ### mulitply by norm.factor
    countsMerged_merge[[i]][is.na(countsMerged_merge[[i]])] <- 0 ### retain TUM.DP == 0
  }
  countsMerged_merge[,norm.factor := NULL]

  #### DIRECT GZIP OUTPUT DOESN'T WORK ??
  #### GZIP SEEMS TO BE CORRUPT
  ## output_filename <- gsub(".dat.gz$", ".norm_normal_depth.dat.gz", filename)
  ## gz = gzfile(output_filename,"w");
  ## write.maf(countsMerged_merge, gz)

##  output_filename <- gsub(".dat.gz$", ".norm_normal_depth.dat", filename)
  GZOUT = gzfile(args$outfile)
  write.table(countsMerged_merge, file=GZOUT, 
              quote = F, 
              col.names = T, 
              row.names = F, 
              sep = "\t")
###  system(paste0("gzip -f ", output_filename))
}

if(!interactive()){
  parser = ArgumentParser()
  parser$add_argument('-f', '--file', type='character', help='Name of counts file to be normalized.')
  parser$add_argument('-o', '--outfile', type='character', help='Normalized counts file')
  args=parser$parse_args()
  normalize_facets_depth(args)
}
