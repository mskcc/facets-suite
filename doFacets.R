#!/opt/common/CentOS_6-dev/R/R-3.2.2/bin/Rscript

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


##########################################################################################
##########################################################################################

print_run_details <- function(out, fit, COUNTS_FILE, TAG, DIRECTORY, CVAL,
                              DIPLOGR, NDEPTH, SNP_NBHD,MIN_NHET, GENOME,
                              GGPLOT, SINGLE_CHROM, SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL, unmatched){

    if(is.null(GIVE_PCVAL)){GIVE_PCVAL = 'NA'}
    ff=paste0(DIRECTORY, "/", TAG,".out")
    cat("# INPUT PARAMETERS GIVEN\n", file=ff)
    cat("# TAG ="        , TAG                   , "\n" , file=ff, append=T)
    cat("# Facets version ="    , as.character(RLIB_VERSION), "\n", file=ff, append=T)
    cat("# Input ="      , basename(COUNTS_FILE) , "\n", file=ff, append=T)
    cat("# Seed ="       , SEED                  , "\n", file=ff, append=T)
    cat("# snp.nbhd ="   , SNP_NBHD              , "\n", file=ff, append=T)
    cat("# ndepth ="     , NDEPTH                , "\n", file=ff, append=T)
    cat("# purity_cval = "  , GIVE_PCVAL             , "\n", file=ff, append=T)
    cat("# unmatched = "  , unmatched            , "\n", file=ff, append=T)
    cat("# cval ="       , CVAL                  , "\n", file=ff, append=T)
    cat("# min.nhet ="   , MIN_NHET              , "\n", file=ff, append=T)
    cat("# genome ="     , GENOME                , "\n", file=ff, append=T)

    cat("\n# LOADED MODULE INFO\n", file=ff, append=T)
    pv = packageVersion('facets')
    cat("# Facets version =", as.character(pv), "\n" , file=ff, append=T)
    cat("\n", file=ff, append=T)
    #buildData=installed.packages()["facets",]
    #version=buildData["Version"]
    #for(fi in c("Package","LibPath","Version","Built")){
    #   cat("#",paste(fi,":",sep=""), buildData[fi], "\n", file=ff, append=T)
    #}

    cat("\n# FACETS OUTPUT\n", file=ff, append=T)
    cat("# Purity ="     , fit$purity            ,"\n", file=ff, append=T)
    cat("# Ploidy ="     , fit$ploidy            ,"\n", file=ff, append=T)
    cat("# dipLogR ="    , fit$dipLogR           ,"\n", file=ff, append=T)
    cat("# dipt ="       , fit$dipt              ,"\n", file=ff, append=T)
    cat("# loglik ="     , fit$loglik            ,"\n", file=ff, append=T)
    cat("# output flags\n"                           , file=ff, append=T)
    cat(paste0("# ", out$flags, "\n"), sep=""        , file=ff,append=T)

    ## cat("# Project =" ,projectName           ,"\n", file=ff ,append=T)
    ## cat("# Tumor ="   ,tumorName             ,"\n", file=ff ,append=T)
    ## cat("# Normal ="  ,normalName            ,"\n", file=ff ,append=T)
}


##########################################################################################
##########################################################################################

select_genome <- function(GENOME, SINGLE_CHROM){

  if(GENOME == 'hg19'){
    data(hg19gcpct)
    if(SINGLE_CHROM == 'F'){chromLevels=c(1:22, "X")}
    if(SINGLE_CHROM != 'F'){chromLevels=SINGLE_CHROM}
  }

  if(GENOME == 'hg18'){ #Is hg18gcpct available?
    data(hg18gcpct)
    if(SINGLE_CHROM == 'F'){chromLevels=c(1:22, "X")}
    if(SINGLE_CHROM != 'F'){chromLevels=SINGLE_CHROM}
  }

  if(GENOME == 'mm9'){
    data(mm9gcpct)
    if(SINGLE_CHROM == 'F'){chromLevels=c(1:19)}
    if(SINGLE_CHROM != 'F'){chromLevels=SINGLE_CHROM}
  }

  if(!GENOME %in% c('hg19', 'hg18', 'mm9')){
    stop(paste("Invalid Genome",GENOME))
  }
  chromLevels
}


##########################################################################################
##########################################################################################

seed_setting <- function(SEED){

  if(is.null(SEED)){

    machine_name_and_date_string <- paste(Sys.info()['nodename'], format(Sys.time(), "%a %b %d %Y %H:%M:%OS3"))
    multiply_mod_integer_max <- function(a, b){(a*b) %% .Machine$integer.max}

    SEED <- Reduce(multiply_mod_integer_max, as.numeric(charToRaw(machine_name_and_date_string)))
  }
  set.seed(seed = SEED)
}


##########################################################################################
##########################################################################################

seg_figure <- function(out, DIRECTORY, TAG, chromLevels, CVAL){

  CairoPNG(file=paste0(DIRECTORY, "/", TAG,".BiSeg.png"),height=1000,width=800)
  plotSample(out, chromlevels=chromLevels)
  text(-.08, -.08, paste(TAG,"cval =",CVAL), xpd=T, pos=4)
  dev.off()
}


##########################################################################################
##########################################################################################

write_output <- function(out, fit, DIRECTORY, TAG){

  out$IGV=formatSegmentOutput(out, TAG)
  save(out, fit, file=paste0(DIRECTORY, "/", TAG,".Rdata"), compress=T)
  write.table(out$IGV,file=paste0(DIRECTORY, "/", TAG,'.seg'), row.names=F, quote=F, sep="\t") #NEW

  out_cols = c('ID','chrom','loc.start','loc.end')
  cncf_cols = c("seg","num.mark","nhet","cnlr.median","mafR","segclust","cnlr.median.clust","mafR.clust","cf","tcn",
    "lcn","cf.em","tcn.em","lcn.em")
  write.xls(cbind(out$IGV[,out_cols],
    fit$cncf[,cncf_cols]), paste0(DIRECTORY, "/", TAG,".cncf.txt"), row.names=F)
}


##########################################################################################
##########################################################################################

results_figure <- function(out, fit, DIRECTORY, TAG, CVAL, GGPLOT, SINGLE_CHROM, GIVE_PCVAL, EM_PLOT=FALSE){

    if(SINGLE_CHROM == 'F'){

        filename = paste0(DIRECTORY, "/", TAG,".CNCF")
        h = 1100
        w = 850
        fit$dipLogR = tryCatch({round(fit$dipLogR,2)},error=function(cond){return(fit$dipLogR)})
        fit$ploidy = tryCatch({round(fit$ploidy,2)},error=function(cond){return(fit$ploidy)})
        fit$purity = tryCatch({round(fit$purity,2)},error=function(cond){return(fit$purity)})
        TAG = tryCatch({unlist(strsplit(TAG,'__'))[1]},error=function(cond){return(TAG)})

        if(is.null(GIVE_PCVAL)){
          main = paste(TAG, ' | cval=', CVAL, ' | purity=', fit$purity, ' | ploidy= ', fit$ploidy, ' | dipLogR=', fit$dipLogR, sep='')
        }
        if(!is.null(GIVE_PCVAL)){
          main = paste(TAG, ' | purity_cval=', PURITY_CVAL, ' | cval=', CVAL, ' | purity=', fit$purity, ' | ploidy= ', fit$ploidy, ' | dipLogR=', fit$dipLogR, sep='')
        }

        if(GGPLOT == 'F'){
            source(file.path(getSDIR(),"fPlots.R"))
            ## base graphics version
            CairoPNG(file=filename, height=h, width=w)
            plotSampleCNCF.custom(out$jointseg, out$out,fit, main=main)
            dev.off()
        }

        if(GGPLOT == 'T'){
            source(file.path(getSDIR(),"fPlots_ggplot2.R"))
            plot.facets.all.output(out, fit, type='png', main=main, plotname=filename, em.plot = EM_PLOT)
        }
    }
}


##########################################################################################
##########################################################################################

extract_six_column_counts_matrix <- function(COUNTS_FILE){
  mat <- fread(paste0("gunzip --stdout ", COUNTS_FILE))
  mat[, TUM.RD := get(paste0("TUM.", Ref, "p")) + get(paste0("TUM.", Ref, "n")),
      1:nrow(mat)]
  mat[, NOR.RD := get(paste0("NOR.", Ref, "p")) + get(paste0("NOR.", Ref, "n")),
      1:nrow(mat)]

  mat[, c("Ref", "Alt",
          "TUM.Ap", "TUM.Cp", "TUM.Gp", "TUM.Tp",
          "TUM.An", "TUM.Cn", "TUM.Gn", "TUM.Tn",
          "NOR.Ap", "NOR.Cp", "NOR.Gp", "NOR.Tp",
          "NOR.An", "NOR.Cn", "NOR.Gn", "NOR.Tn") := NULL]
  setnames(mat, 1, "Chromosome") ## must be called Chromosome
  setnames(mat, 2, "Position") ## must be called Position
  setcolorder(mat,
              c("Chromosome", "Position",
                "NOR.DP", "NOR.RD",
                "TUM.DP", "TUM.RD"))

  mat <- as.data.frame(mat)
  mat
}

##########################################################################################
##########################################################################################

facets_iteration <- function(COUNTS_FILE, TAG, DIRECTORY, CVAL, DIPLOGR, NDEPTH,
                             SNP_NBHD, MIN_NHET, GENOME, GGPLOT, SINGLE_CHROM,
                             SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL, unmatched){

    if (RLIB_VERSION >= '0.5.2') {


      rcmat = readSnpMatrix(COUNTS_FILE, err.thresh = 10, del.thresh = 10)

      dat = preProcSample(rcmat, ndepth = NDEPTH, het.thresh = 0.25, snp.nbhd = SNP_NBHD, cval = CVAL,
        gbuild = GENOME, hetscale = TRUE, unmatched = unmatched, ndepthmax = 1000)

      out = procSample(dat, cval = CVAL, min.nhet = MIN_NHET, dipLogR = DIPLOGR)
      fit = emcncf(out)

      fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)

      write_output(out, fit, DIRECTORY, TAG)
      print_run_details(out, fit, COUNTS_FILE, TAG, DIRECTORY, CVAL, DIPLOGR, NDEPTH, SNP_NBHD,
                        MIN_NHET, GENOME, GGPLOT, SINGLE_CHROM, SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL, unmatched)

      results_figure(out, fit, DIRECTORY, TAG, CVAL, GGPLOT, SINGLE_CHROM, GIVE_PCVAL, EM_PLOT=TRUE)

      return(fit$dipLogR)

    } else {

      chromLevels = select_genome(GENOME, SINGLE_CHROM)

      dat=preProcSample(COUNTS_FILE,snp.nbhd=SNP_NBHD,cval=CVAL,chromlevels=chromLevels,ndepth=NDEPTH)
      out=procSample(dat,cval=CVAL,min.nhet=MIN_NHET,dipLogR=DIPLOGR)

      #seg_figure(out, DIRECTORY, TAG, chromLevels, CVAL) no need for this

      fit=emcncf(out) #fit=emcncf(out$jointseg,out$out,dipLogR=out$dipLogR) OLD

      write_output(out, fit, DIRECTORY, TAG)
      print_run_details(out, fit, COUNTS_FILE, TAG, DIRECTORY, CVAL, DIPLOGR, NDEPTH, SNP_NBHD,
                        MIN_NHET, GENOME, GGPLOT, SINGLE_CHROM, SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL, unmatched)

      results_figure(out, fit, DIRECTORY, TAG, CVAL, GGPLOT, SINGLE_CHROM, GIVE_PCVAL)

      return(fit$dipLogR)
    }

}


##########################################################################################
##########################################################################################

source(file.path(getSDIR(),"funcs.R"))
source(file.path(getSDIR(),"nds.R"))

library(ggplot2)
library(Cairo)
library(argparse)


parser=ArgumentParser()
parser$add_argument("-c", "--cval",type="integer",default=50,help="critical value for segmentation")
parser$add_argument("-s", "--snp_nbhd",type="integer",default=250,help="window size")
parser$add_argument("-n", "--ndepth",type="integer",default=35,help="threshold for depth in the normal sample")
parser$add_argument("-m", "--min_nhet",type="integer",default=15,
    help="minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segments")

parser$add_argument("-pc", "--purity_cval",type="integer", help="critical value for segmentation")
parser$add_argument("-ps", "--purity_snp_nbhd",type="integer",default=250,help="window size")
parser$add_argument("-pn", "--purity_ndepth",type="integer",default=35,help="threshold for depth in the normal sample")
parser$add_argument("-pm", "--purity_min_nhet",type="integer",default=15,
    help="minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segments")

parser$add_argument("-d", "--dipLogR",type="double",help="diploid log ratio")
parser$add_argument("-g", "--genome",type="character",default="hg19",help="Genome of counts file")
parser$add_argument("-f", "--counts_file", type="character", required=T, help="paired Counts File")
parser$add_argument("-t", "--TAG", type="character", required=T, help="output prefix")
parser$add_argument("-D", "--directory", type="character", required=T, help="output prefix")
parser$add_argument("-r", "--R_lib", type="character", default='latest', help="Which version of FACETs to load into R")
parser$add_argument("-C", "--single_chrom", type="character", default='F',help="Perform analysis on single chromosome")
parser$add_argument("-G", "--ggplot2", type="character", default='T', help="Plots using  ggplot2")
parser$add_argument("--seed", type="integer", help="Set the seed for reproducibility")
parser$add_argument("-u", "--unmatched", type="character", default='F', help="run using a pooled normal")
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
if(RLIB_PATH != "latest"){
    library(facets, lib.loc=RLIB_PATH)
    #library(facets)
} else{
    library(facets)
}
RLIB_VERSION = packageVersion('facets')
print(RLIB_VERSION)

SINGLE_CHROM = args$single_chrom
GENOME=args$genome
GGPLOT=args$ggplot2
SEED=args$seed
unmatched=args$unmatched
if(unmatched == 'F'){
unmatched = FALSE
} else {
unmatched = TRUE
}
seed_setting(SEED)

if(!is.null(PURITY_CVAL)){

    ### if "PURITY_CVAL" is specified, run FACETS twice. Take dipLogR from the first run...
    estimated_dipLogR = facets_iteration(COUNTS_FILE, paste0(TAG, "_purity"), DIRECTORY, PURITY_CVAL, DIPLOGR, NDEPTH, SNP_NBHD, PURITY_MIN_NHET, GENOME, GGPLOT, SINGLE_CHROM, SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL=NULL, unmatched = unmatched)

    ### ... and use it for a second run of FACETS
    facets_iteration(COUNTS_FILE, paste0(TAG, "_hisens"), DIRECTORY, CVAL, estimated_dipLogR, NDEPTH, SNP_NBHD, MIN_NHET, GENOME, GGPLOT, SINGLE_CHROM, SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL=PURITY_CVAL, unmatched = unmatched)
}
if(is.null(PURITY_CVAL)){

    ### if "PURITY_CVAL" is not specified, run FACETS only once, with dipLogR taken from the arguments
    facets_iteration(COUNTS_FILE, TAG, DIRECTORY, CVAL, DIPLOGR, NDEPTH, SNP_NBHD, MIN_NHET, GENOME, GGPLOT, SINGLE_CHROM, SEED, RLIB_PATH, RLIB_VERSION, GIVE_PCVAL=NULL, unmatched = unmatched)
}



## BASE=basename(FILE)
## BASE=gsub("countsMerged____","",gsub(".dat.*","",BASE))

## sampleNames=gsub(".*recal_","",strsplit(BASE,"____")[[1]])
## tumorName=sampleNames[1]
## normalName=sampleNames[2]
## projectName=gsub("_indel.*","",strsplit(BASE,"____")[[1]])[1]
