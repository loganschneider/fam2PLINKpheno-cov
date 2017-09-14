# To use this file properly:
# -have fam2cov.R in the same directory as folders containing the PLINK binary .fam files
#    -fam (.fam) file format: https://www.cog-genomics.org/plink2/formats#fam
#    -cov (.txt) file format: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#cov
# -need to have cov types in association with our plate ID (e.g. Wisc12H01), and labeled "ID_2"
#    -these are sometimes located in existing SNPtest sample files or original Excel/csv files
# -set working directory to a level AT OR ABOVE where fam2cov.R is located

###TODO: need to generate cov columns for plate #, age, optimal number of MDS dimensions (pulled from .mds file)
fam2cov <- function(famdir="fams4cov_addition",
                    famkeepcols=c(1:2),
                    covdir="covs2add",
                    covkeepcols=c(2:ncol(covdf)),
                    studyname="WSC",
                    phenoname="PLMD",
                    mergeon="ID_2",
                    GenoArray="other") {
    #get to the directory in which fam2cov resides
    #should be at least a directory above directories containing fam and cov files
    dirpattern <- "fam2cov.R"
    myfunction.path <- list.files(getwd(),recursive=TRUE,full.names=TRUE,pattern=dirpattern)
    setwd(paste(gsub(dirpattern,'\\1',myfunction.path),famdir,sep="/"))
    #create sub-directory for placement of new cov files
    if(!dir.exists(paste(studyname,phenoname,"covfiles",sep="_"))) {
        dir.create(paste(studyname,phenoname,"covfiles",sep="_"))
    }
    covout.path <- paste(getwd(),paste(studyname,phenoname,"covfiles",sep="_"),sep="/")
    #file list of fam files into which you will add covtypes
    famlist <- list.files(pattern='.fam$')
    setwd(paste(gsub(dirpattern,'\\1',myfunction.path),paste(covdir,phenoname,sep="/"),sep=""))
    #file list of cov files from which you will add covtype
    #-need to be in csv format
    #-need to have plateID (e.g. Wisc12H01) in second column, which has heading "ID_2"
    covlist <- list.files(pattern='.csv$')
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
        #read in cov file
        covdf <- read.csv(covlist,header=T)
        #subset fam file
        famsub <- famfile[,famkeepcols]
        names(famsub) <- c("ID_1",mergeon)
        #request columns to keep
        display <- data.frame(seq(1:ncol(covdf)),colnames(covdf))
        colnames(display) <- c("Number","Column Name")
        print(display[3:nrow(display),])
        covkeepcols <- readline("List columns to keep, separated by commas, WITHOUT spaces: ")
        print(covkeepcols)
        covkeepcols <- paste("2",covkeepcols,sep=",")
        covcols2keep <- as.numeric(unlist(strsplit(covkeepcols, ",")))
        #subset covdf
        covsub <- covdf[,c(covcols2keep)]
        #merge INTO famsub based on "mergeon" column, don't sort, don't add suffixes based on origin df
        cam <- merge.data.frame(famsub,covsub,by=mergeon,sort=F,suffixes=NULL)
        camreor <- cam[,c(2,1,3:ncol(cam))]
        camadd <- camreor
        camadd$Array_psuedocohort <- famcohort
        for(k in 1:nrow(camadd)){
            extract <- gregexpr('[0-9]+',camadd[k,mergeon])
            camadd[k,"plate"] <- as.integer(regmatches(camadd[k,mergeon],extract)[[1]][1])
            #camadd[k,"plate"] <- gsub('.*[a-zA-Z]([0-9]+)[a-zA-Z].*','\\1',camadd[k,mergeon]) ###Doesn't work for things like NA06984
        }
        #add in MDS dimensions that were chosen (assuming there are any chosen)
        MDS <- read.table(paste(famcohort,"optimalMDScovs.txt",sep="_"),header=T)
        if (ncol(MDS) > 2){
            names(MDS)[1] <- c("ID_1")
            names(MDS)[1] <- c(mergeon)
            camMDS <- merge.data.frame(camadd,MDS[,2:ncol(MDS)],by=mergeon,sort=F,suffixes=NULL)
            camadd <- cam[,c(2,1,3:ncol(camMDS))]
        }
        write.table(camadd,
                    file=paste(covout.path,paste(paste(famcohort,studyname,"cov",sep="_"),"txt",sep="."),sep="/"),
                    quote=F,
                    row.names=F)
    }
    setwd(gsub(dirpattern,'\\1',myfunction.path))
}