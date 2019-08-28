library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)
library(Biobase)

getgCSI <-
  function (verbose=FALSE,
            nthread=1){
            
options(stringsAsFactors = FALSE)
z <- list()

matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}

load("/pfs/downloadrnagCSI/gCSI_2017_molecprofile.RData")
load("/pfs/gcsi2017ProfilesAssemble/profiles.RData")
load("/pfs/getgCSI2017/gcsidrugpost.RData")
load("/pfs/gcsi2017raw/sensitivity.RData")

    
recomputed_2017 <- res
sensitivityProfiles_2017 <- data.frame("AAC" = NA, "IC50"=NA, "ic50_published"=NA, "meanviability_published"=NA, "HS"=NA, "E_inf"=NA, "EC50"=NA)
sensitivityProfiles_2017[nrow(sensitivityProfiles_2017)+6454,] <- NA
sensitivityProfiles_2017[,"AAC"] <- as.numeric(recomputed_2017[,"AAC"])
sensitivityProfiles_2017[,"IC50"] <- as.numeric((recomputed_2017[,"IC50"]))
sensitivityProfiles_2017[,"HS"] <- as.numeric((recomputed_2017[,"HS"]))
sensitivityProfiles_2017[,"E_inf"] <- as.numeric((recomputed_2017[,"E_inf"]))
sensitivityProfiles_2017[,"EC50"] <- as.numeric((recomputed_2017[,"EC50"]))
sensitivityProfiles_2017[,"ic50_published"] <- gcsi.ic50[rownames(raw.sensitivity),"IC50"]
sensitivityProfiles_2017[,"meanviability_published"] <- gcsi.mv[rownames(raw.sensitivity),"Mean Viability"]
rownames(sensitivityProfiles_2017) <- rownames(recomputed_2017)
sensitivityProfiles_2017 <- sensitivityProfiles_2017[rownames(raw.sensitivity),]

cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))


curationCell <- cell_all[which(!is.na(cell_all[ , "gCSI.cellid"])),]
curationTissue <- cell_all[which(!is.na(cell_all[ , "gCSI.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "gCSI.cellid")]
curationTissue <- curationTissue[ , c("unique.tissueid", "gCSI.tissueid")]

rownames(curationTissue) <- curationCell[ , "unique.cellid"]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "gCSI.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "gCSI.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]


sensitivityInfo_2017 <- sensitivityInfo_2017[rownames(raw.sensitivity),]

emptyEset <- ExpressionSet()
annotation(emptyEset) <- "issue with cell annotation, will fix later"

gCSI_2017 <- PharmacoSet(molecularProfiles=list("rnaseq" = emptyEset),
                       name="gCSI",
                       cell=curationCell,
                       drug=curationDrug,
                       sensitivityInfo=sensitivityInfo_2017,
                       sensitivityRaw=raw.sensitivity,
                       sensitivityProfiles=sensitivityProfiles_2017,
                       curationCell=curationCell,
                       curationDrug=curationDrug,
                       curationTissue=curationTissue,
                       datasetType="sensitivity")

save(gCSI_2017, file="/pfs/out/gCSI_2017.RData")
    
return (gCSI_2017)

}

getgCSI(verbose=FALSE, nthread=1)
