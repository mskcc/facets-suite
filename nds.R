###
# nds funcs
#
cc <- function(...) {paste(...,sep='_')}
len <- function(x) {length(x)}
write.xls <- function(dd,filename,row.names=T,col.names=NA) {
  if (!is.data.frame(dd)) {
    dd <- data.frame(dd,check.names=F)
  }
  if(!row.names) {
    col.names=T
  }
  write.table(dd,file=filename,sep="\t",quote=FALSE,
              col.names=col.names,row.names=row.names)
}


###
