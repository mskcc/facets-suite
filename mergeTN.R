#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript
library(data.table)

args <- commandArgs(TRUE)
TUMOR <- args[1]; args <- args[-1]
NORMAL <- args[1]; args <- args[-1]
MINCOV_NORMAL=25

read_counts <- function(file){
  dt <- fread(paste0("gunzip --stdout ", file), showProgress=FALSE)
  setkeyv(dt, c("Chrom", "Pos", "Ref", "Alt"))
  dt[,Refidx := NULL]
  dt[,TOTAL_depth := NULL]
  dt[,MAPQ_depth := NULL]
  dt[,ID := NULL]
  dt[,INS := NULL]
  dt[,DEL := NULL]
  dt[Ref != "N" & nchar(Ref) == 1 & nchar(Alt) == 1]
}

write("Reading normal ...", stderr())
NORMAL_dt <- read_counts(NORMAL)
NORMAL_dt <- NORMAL_dt[BASEQ_depth >= MINCOV_NORMAL]
write("done ...", stderr())

write("Reading tumor ...", stderr())
TUMOR_dt <- read_counts(TUMOR)
write("done ...", stderr())

mergeTN <- merge(TUMOR_dt, NORMAL_dt, suffixes = c(".TUM", ".NOR"))

setnames(mergeTN,
         c("Chrom", "Pos", "Ref", "Alt", "TUM.DP", "TUM.Ap", "TUM.Cp", 
           "TUM.Gp", "TUM.Tp", "TUM.An", "TUM.Cn", "TUM.Gn", "TUM.Tn", "NOR.DP", 
           "NOR.Ap", "NOR.Cp", "NOR.Gp", "NOR.Tp", "NOR.An", "NOR.Cn", "NOR.Gn", 
           "NOR.Tn"))

mergeTN[, Chrom := factor(Chrom, levels=c(c(1:22, "X", "Y", "MT"), paste0("chr", c(1:22, "X", "Y", "M"))))]
mergeTN <- mergeTN[order(Chrom, Pos)]

write.table(mergeTN, file=stdout(), 
            quote = F, 
            col.names = T, 
            row.names = F, 
            sep = "\t")
