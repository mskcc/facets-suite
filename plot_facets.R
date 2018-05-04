#!/opt/common/CentOS_6-dev/R/R-3.2.2/bin/Rscript

library(ggplot2)
library(Cairo)
library(argparse)
library(stringr)
getSDIR = function(){
    args=commandArgs(trailing=F)
    TAG="--file="
    path_idx=grep(TAG,args)
    SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
    if(length(SDIR)==0) {
        return(getwd())
    } else {
        return(SDIR)
    }
}

source(file.path(getSDIR(),"fPlots_ggplot2.R"))

args = commandArgs(TRUE)
if (length(args) < 1) stop('Usage: plot-facets.R sample.Rdata optional-name')

rdata = args[1]
if (length(args)==2) { name = args[2] } else { name = NULL }
load(rdata)
outfile = gsub('.Rdata', '.out', rdata)
outfile = readLines(outfile)

DIRECTORY=dirname(rdata)
TAG=gsub('.Rdata', '', basename(rdata))
if (!is.null(name)) {
		DIRECTORY = getwd()
		suffix = str_extract(TAG, 'purity|hisens')
		if (is.na(suffix)) { TAG = name } else { TAG = str_c(name, '_', suffix) }
}

filename = paste0(DIRECTORY, "/", TAG,".CNCF")

fit$dipLogR = tryCatch({round(fit$dipLogR,2)},error=function(cond){return(fit$dipLogR)})
fit$ploidy = tryCatch({round(fit$ploidy,2)},error=function(cond){return(fit$ploidy)})
fit$purity = tryCatch({round(fit$purity,2)},error=function(cond){return(fit$purity)})
TAG = strsplit(TAG,'__')[[1]]

if (grepl('purity', rdata)) {
    CVAL = stringr::str_extract(grep('# purity_cval',  outfile, value = T), '[0-9]+')
} else { CVAL = stringr::str_extract(grep('# cval',  outfile, value = T), '[0-9]+') }


main = paste(TAG, ' | cval=', CVAL, ' | purity=', fit$purity, ' | ploidy= ', fit$ploidy, ' | dipLogR=', fit$dipLogR, sep='')

plot.facets.all.output(out, fit, type='png', main=main, plotname=filename, em.plot = T)
