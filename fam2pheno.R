# To use this file properly:
# -have fam2pheno.R in the same directory as folders containing the PLINK binary .fam files
#    -fam (.fam) file format: https://www.cog-genomics.org/plink2/formats#fam
#    -pheno (.csv) file format: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#pheno
# -need to have phenotypes in association with our plate ID (e.g. Wisc12H01), and labeled "ID_2"
#    -these are sometimes located in existing SNPtest sample files or original Excel/csv files
# -set working directory to a level AT OR ABOVE where fam2pheno.R is located
fam2pheno <- function(famdir="fams4pheno_addition",
                      famkeepcols=c(1:2),
                      phenodir="phenos2add",
                      phenokeepcols=c(2:ncol(phenodf)),
                      studyname="WSC",
                      phenoname="PLMD",
                      mergeon="ID_2",
                      GenoArray="other") {
    #get to the directory in which fam2pheno resides
    #should be at least a directory above directories containing fam and pheno files
    dirpattern <- "fam2pheno.R"
    myfunction.path <- list.files(getwd(),recursive=TRUE,full.names=TRUE,pattern=dirpattern)
    setwd(paste(gsub(dirpattern,'\\1',myfunction.path),famdir,sep="/"))
    #create sub-directory for placement of new pheno files
    if(!dir.exists(paste(studyname,phenoname,"phenofiles",sep="_"))) {
        dir.create(paste(studyname,phenoname,"phenofiles",sep="_"))
    }
    phenoout.path <- paste(getwd(),paste(studyname,phenoname,"phenofiles",sep="_"),sep="/")
    #file list of fam files into which you will add phenotypes
    famlist <- list.files(pattern='.fam$')
    setwd(paste(gsub(dirpattern,'\\1',myfunction.path),paste(phenodir,phenoname,sep="/"),sep=""))
    #file list of pheno files from which you will add phenotype
    #-need to be in csv format
    #-need to have plateID (e.g. Wisc12H01) in second column, which has heading "ID_2"
    phenolist <- list.files(pattern='.csv$')
    for(j in famlist) {
        #read in fam file
        famfile <- read.table(paste(paste(gsub(dirpattern,'\\1',myfunction.path),famdir,sep=""),j,sep="/"), header=F)
        #identify which chip pseudocohort the fam file belongs to
        #ideally the fam files have the format pseudocohortname_whatever_else.fam
        famcohort <- as.character(gsub('([0-9]|[a-zA-Z]+)_.*','\\1',j))
        #create a code for fam chip pseudocohorts (for comparison later): 1 - 500K, 2 - ADD, 3 - Affy6
        if(grepl("500[Kk]|[Ff][Ii][Vv][Ee][Kk]",famcohort)){
            famcode<-1
            famcohort<-"500K"
        } else if(grepl("[Aa][Dd][Dd]",famcohort)) {
            famcode<-2
            famcohort<-"ADD"
        } else if(grepl("[Aa][Ff][Ff][Yy]",famcohort)) {
            famcode<-3
            famcohort<-"Affy6"
        } else {
            famcode<-4
            famcohort<-GenoArray
        }
        #read in pheno file
        phenodf <- read.csv(phenolist,header=T)
        #subset fam file
        famsub <- famfile[,famkeepcols]
        names(famsub) <- c("ID_1",mergeon)
        #subset phenodf
        phenosub <- phenodf[,phenokeepcols]
        #merge INTO famsub based on "mergeon" column, don't sort, don't add suffixes based on origin df
        pham <- merge.data.frame(famsub,phenosub,by=mergeon,sort=F,suffixes=NULL)
        phamreor <- pham[,c(2,1,3:ncol(pham))]
        names(phamreor)[1] <- "FID"
        names(phamreor)[2] <- "IID"
        write.table(phamreor, file=paste(phenoout.path,paste(paste(famcohort,"pheno",sep="_"),"txt",sep="."),sep="/"),quote=F,row.names=F)
    }
    setwd(gsub(dirpattern,'\\1',myfunction.path))
}
