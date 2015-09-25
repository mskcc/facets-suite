#!/opt/common/CentOS_6-dev/R/R-3.1.3/bin/Rscript

### modified version of doFacets.R
###
### runs FACETS twice with two sets of input parameters
### dipLogR from the first iteration is used for the second
### output files for both iterations are retained


getSDIR <- function(){
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

source(file.path(getSDIR(),"funcs.R"))
source(file.path(getSDIR(),"nds.R"))

#buildData=installed.packages()["facets",]
#cat("#Module Info\n")
#for(fi in c("Package","LibPath","Version","Built")){
#    cat("#",paste(fi,":",sep=""),buildData[fi],"\n")
#}
#version=buildData["Version"]
#cat("\n")

library(argparse)
parser=ArgumentParser()

parser$add_argument("-c","--cval",type="integer",default=50,help="critical value for segmentation")
parser$add_argument("-s","--snp_nbhd",type="integer",default=250,help="window size")
parser$add_argument("-n","--ndepth",type="integer",default=35,help="threshold for depth in the normal sample")
parser$add_argument("-m","--min_nhet",type="integer",default=25,
    help="minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segments")

parser$add_argument("-pc","--purity_cval",type="integer",default=-99,help="critical value for segmentation")
parser$add_argument("-ps","--purity_snp_nbhd",type="integer",default=250,help="window size")
parser$add_argument("-pn","--purity_ndepth",type="integer",default=35,help="threshold for depth in the normal sample")
parser$add_argument("-pm","--purity_min_nhet",type="integer",default=25,
    help="minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segments")

parser$add_argument("-d","--dipLogR",type="double",default=-99,help="diploid log ratio")
parser$add_argument("--genome",type="character",default="hg19",help="Genome of counts file")
parser$add_argument("counts_file",nargs=1,help="Paired Counts File")
parser$add_argument("TAG",nargs=1,help="output prefix")
parser$add_argument("directory",nargs=1,help="output prefix")
parser$add_argument("-r", "--R_lib", type="character", default="latest", help="Which version of FACETs to load into R")
parser$add_argument("-C", "--single_chrom", type="character", default='F',help="Perform analysis on single chromosome")
parser$add_argument("-g", "--ggplot2", type="character", default='T', help="Plots using  ggplot2")
args=parser$parse_args()

CVAL=args$cval
SNP_NBHD=args$snp_nbhd
NDEPTH=args$ndepth
MIN_NHET=args$min_nhet

PURITY_CVAL=args$purity_cval
PURITY_SNP_NBHD=args$purity_snp_nbhd
PURITY_NDEPTH=args$purity_ndepth
PURITY_MIN_NHET=args$purity_min_nhet

COUNTS_FILE=args$counts_file
TAG=args$TAG
DIRECTORY=args$directory
DIPLOGR=args$dipLogR

RLIB_PATH=args$R_lib
if(exists('RLIB_PATH')){
    library(facets, lib.loc=RLIB_PATH)
} else{
    library(facets)
}
RLIB_VERSION = packageVersion('facets')
write(paste('Path to facets library loaded:', RLIB_PATH), stdout())
write(paste('Version of facets library installed:', RLIB_VERSION), stdout())

## BASE=basename(FILE)
## BASE=gsub("countsMerged____","",gsub(".dat.*","",BASE))

## sampleNames=gsub(".*recal_","",strsplit(BASE,"____")[[1]])
## tumorName=sampleNames[1]
## normalName=sampleNames[2]
## projectName=gsub("_indel.*","",strsplit(BASE,"____")[[1]])[1]

SINGLE_CHROM = args$single_chrom
GENOME=args$genome
GGPLOT=args$ggplot2

facets_iteration <- function(COUNTS_FILE = COUNTS_FILE,
                             TAG = TAG,
                             DIRECTORY = DIRECTORY,
                             CVAL = CVAL,
                             DIPLOGR = DIPLOGR,
                             NDEPTH = NDEPTH,
                             SNP_NBHD = SNP_NBHD,
                             MIN_NHET = MIN_NHET,
                             GENOME = GENOME,
                             SINGLE_CHROM = SINGLE_CHROM){

    if(PURITY_CVAL==-99){
        PURITY_CVAL=NULL
    }
    if(DIPLOGR==-99){
        DIPLOGR=NULL
    }

    cat("COUNTS_FILE =", COUNTS_FILE, "\n")
    cat("TAG =", TAG, "\n")
    cat("DIRECTORY =", DIRECTORY, "\n")
    cat("CVAL =", CVAL, "\n")
    cat("DIPLOGR =", DIPLOGR, "\n")
    cat("NDEPTH =", NDEPTH, "\n")
    cat("SNP_NBHD =", SNP_NBHD, "\n")
    cat("MIN_NHET =", MIN_NHET, "\n")
    cat("GENOME =", GENOME, "\n")
    
    
    buildData=installed.packages()["facets",]
    cat("#Module Info\n")
    for(fi in c("Package","LibPath","Version","Built")){
        cat("#",paste(fi,":",sep=""),buildData[fi],"\n")
    }
    version=buildData["Version"]
    cat("\n")

    switch(args$genome,
    hg19={
        data(hg19gcpct)
        if(SINGLE_CHROM == 'F'){chromLevels=c(1:22, "X")}
        if(SINGLE_CHROM != 'F'){chromLevels=SINGLE_CHROM}
    },
    mm9={
        data(mm9gcpct)
        if(SINGLE_CHROM == 'F'){chromLevels=c(1:19)}
        if(SINGLE_CHROM != 'F'){chromLevels=SINGLE_CHROM}
    },
    {
        stop(paste("Invalid Genome",args$genome))
    }
)


    pre.CVAL=CVAL
    dat=preProcSample(COUNTS_FILE,snp.nbhd=SNP_NBHD,cval=pre.CVAL,chromlevels=chromLevels,ndepth=NDEPTH)
    out=procSample(dat,cval=CVAL,min.nhet=MIN_NHET,dipLogR=DIPLOGR)
    CairoPNG(file=paste0(DIRECTORY, "/", TAG,".BiSeg.png"),height=1000,width=800)
    plotSample(out,chromlevels=chromLevels)
    text(-.08,-.08,paste(TAG,"cval =",CVAL),xpd=T,pos=4)
    dev.off()
                                        #fit=emcncf(out$jointseg,out$out,dipLogR=out$dipLogR) OLD
    fit=emcncf(out) # NEW
    out$IGV=formatSegmentOutput(out,TAG)
    save(out,fit,file=paste0(DIRECTORY, "/", TAG,".Rdata"),compress=T)
    write.table(out$IGV,file=paste0(DIRECTORY, "/", TAG,'.seg'),row.names=F,quote=F,sep="\t") #NEW
    
    ff=paste0(DIRECTORY, "/", TAG,".out")
    cat("# TAG ="        ,TAG                   ,"\n" ,file=ff)
    cat("# Facets version ="    ,as.character(RLIB_VERSION) ,"\n" ,file=ff ,append=T)
    cat("# Input ="      ,basename(COUNTS_FILE) ,"\n" ,file=ff ,append=T)
    cat("# snp.nbhd ="   ,SNP_NBHD              ,"\n" ,file=ff ,append=T)
    cat("# cval ="       ,CVAL                  ,"\n" ,file=ff ,append=T)
    cat("# min.nhet ="   ,MIN_NHET              ,"\n" ,file=ff ,append=T)
    cat("# genome ="     ,args$genome           ,"\n" ,file=ff ,append=T)
    ## cat("# Project =" ,projectName           ,"\n" ,file=ff ,append=T)
    ## cat("# Tumor ="   ,tumorName             ,"\n" ,file=ff ,append=T)
    ## cat("# Normal ="  ,normalName            ,"\n" ,file=ff ,append=T)
    cat("# Purity ="     ,fit$purity            ,"\n" ,file=ff ,append=T)
    cat("# Ploidy ="     ,fit$ploidy            ,"\n" ,file=ff ,append=T)
    cat("# dipLogR ="    ,fit$dipLogR           ,"\n" ,file=ff ,append=T)
    cat("# dipt ="       ,fit$dipt              ,"\n" ,file=ff ,append=T)
    cat("# loglik ="     ,fit$loglik            ,"\n" ,file=ff ,append=T)
    cat("# output flags\n"                            ,file=ff ,append=T)
    cat(paste0("# ", out$flags, "\n"), sep=""         ,file=ff ,append=T)    
    
    write.xls(cbind(out$IGV[,1:4],fit$cncf[,2:ncol(fit$cncf)]),
              paste0(DIRECTORY, "/", TAG,".cncf.txt"),row.names=F)

    if(SINGLE_CHROM == 'F'){

        filename = paste0(DIRECTORY, "/", TAG,".CNCF.png")
        h = 1100
        w = 850
        if(purity_cval == -99){main = paste(TAG, ' | cval=', CVAL, ' | purity=', round(fit$purity,2), ' | ploidy= ', round(fit$ploidy,2), ' | dipLogR=', round(fit$dipLogR,2), sep='')}
        if(purity_cval != -99){main = paste(TAG, ' | purity_cval=', PURITY_CVAL, ' | cval=', CVAL, ' | purity=', round(fit$purity,2), ' | ploidy= ', round(fit$ploidy,2), ' | dipLogR=', round(fit$dipLogR,2), sep='')}

        if(ggplot2 == 'F'){
            library(Cairo)
            source(file.path(getSDIR(),"fPlots.R"))
            base graphics version
            CairoPNG(file=filename, height=h, width=w)
            plotSampleCNCF.custom(out$jointseg, out$out,fit, main=main)
            dev.off()
        }

        if(ggplot2 == 'T'){
            source(file.path(getSDIR(),"fPlots_ggplot2.R"))
            plot.facets.all.output(out, fit, type='png', main=main, plotname=filename)
        }
        
    }

    return(fit$dipLogR)
}

if(!is.null(PURITY_CVAL)){
### if "PURITY_CVAL" is specified, run FACETS twice. Take dipLogR from the first run...
    estimated_dipLogR <- facets_iteration(COUNTS_FILE, paste0(TAG, "_purity"), DIRECTORY, PURITY_CVAL,           DIPLOGR, PURITY_NDEPTH, PURITY_SNP_NBHD, PURITY_MIN_NHET, GENOME, SINGLE_CHROM)
    estimated_dipLogR <- ifelse(is.null(estimated_dipLogR), -99, estimated_dipLogR) ### if the fit fails (dipLogR is null), then take a second attempt
### ... and use it for a second run of FACETS
                         facets_iteration(COUNTS_FILE, paste0(TAG, "_hisens"), DIRECTORY,        CVAL, estimated_dipLogR,        NDEPTH,        SNP_NBHD,        MIN_NHET, GENOME, SINGLE_CHROM)
} else {
### if "PURITY_CVAL" is not specified, run FACETS only once, with dipLogR taken from the arguments
                         facets_iteration(COUNTS_FILE, TAG, DIRECTORY,        CVAL,            DIPLOGR,       NDEPTH,        SNP_NBHD,        MIN_NHET, GENOME, SINGLE_CHROM)
}
    

