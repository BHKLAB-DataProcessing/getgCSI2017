library(PharmacoGx)

library(PharmacoGx)
# library(devtools)

# install.packages("BiocManager")
# library(BiocManager)
# message("trying install BiocManager")
# message(.libPaths())
# install(c("multtest", "NMF", "rmarkdown", "RColorBrewer"))

# library(devtools)
# message("trying install")
# install_url("http://research-pub.gene.com/gCSI-cellline-data/compareDrugScreens_current.tar.gz")


# library(compareDrugScreens)

library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)
library(Biobase)
library(reshape2)
install.packages("abind")
library(abind)
verbose=FALSE
nthread=1
            
options(stringsAsFactors = FALSE)
z <- list()

dir.prefix <- "/pfs"

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
        if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}

# load(file.path(dir.prefix, "downloadrnagCSI/gCSI_2017_molecprofile.RData"))
load(file.path(dir.prefix, "gcsi2017ProfilesAssemble/profiles.RData"))
# load(file.path(dir.prefix, "getgCSI2017/gcsidrugpost.RData"))
load(file.path(dir.prefix, "gcsi2017RawSensitivity/raw.sensitivity.RData"))
load(file.path(dir.prefix, "gcsi2017RawSensitivity/gCSI_molData.RData"))



    
recomputed_2017 <- res
sensitivityProfiles_2017 <- data.frame("aac_recomputed" = NA, "ic50_recomputed"=NA, "ic50_published"=NA, "meanviability_published"=NA, "HS"=NA, "E_inf"=NA, "EC50"=NA)
sensitivityProfiles_2017[nrow(published.profiles),] <- NA
rownames(sensitivityProfiles_2017) <- rownames(published.profiles)

sensitivityProfiles_2017[rownames(recomputed_2017),"aac_recomputed"] <- as.numeric(recomputed_2017[,"AAC"])
sensitivityProfiles_2017[rownames(recomputed_2017),"ic50_recomputed"] <- as.numeric((recomputed_2017[,"IC50"]))
sensitivityProfiles_2017[rownames(recomputed_2017),"HS"] <- as.numeric((recomputed_2017[,"HS"]))
sensitivityProfiles_2017[rownames(recomputed_2017),"E_inf"] <- as.numeric((recomputed_2017[,"E_inf"]))
sensitivityProfiles_2017[rownames(recomputed_2017),"EC50"] <- as.numeric((recomputed_2017[,"EC50"]))
sensitivityProfiles_2017[,"ic50_published"] <- published.profiles$ic50_published
sensitivityProfiles_2017[,"meanviability_published"] <- published.profiles$mean.viability_published

### Need to add empty rows to the raw array because not all of the experiments had raw values.

new.rows <- setdiff(rownames(sensitivity.info), rownames(raw.sensitivity))
new.rows <- array(NA_real_, dim=c(length(new.rows), ncol(raw.sensitivity),2), dimnames = list(new.rows, colnames(raw.sensitivity), dimnames(raw.sensitivity)[[3]]))

library(abind)

raw.sensitivity <- abind(raw.sensitivity, new.rows, along=1)

raw.sensitivity <- raw.sensitivity[rownames(sensitivity.info),,]

cell_all <- read.csv(file.path(dir.prefix, "downAnnotations/cell_annotation_all.csv"), na.strings=c("", " ", "NA"))
drug_all <- read.csv(file.path(dir.prefix, "downAnnotations/drugs_with_ids.csv"), na.strings=c("", " ", "NA"))


curationCell <- cell_all[which(!is.na(cell_all[ , "gCSI.cellid"])),]
curationTissue <- cell_all[which(!is.na(cell_all[ , "gCSI.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "gCSI.cellid")]
curationTissue <- curationTissue[ , c("unique.tissueid", "gCSI.tissueid")]

rownames(curationTissue) <- curationCell[ , "unique.cellid"]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "gCSI.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "gCSI.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]

## Only doing this for data added to the pset.

reps <- matchToIDTable(mut$cellid, curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
mut$cellid <- reps
mut$tissueid <- curationTissue[mut$cellid, "unique.tissueid"]

reps <- matchToIDTable(rnaseq$cellid, curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
rnaseq$cellid <- reps
rnaseq$tissueid <- curationTissue[rnaseq$cellid, "unique.tissueid"]


reps <- matchToIDTable(cnv$cellid, curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
cnv$cellid <- reps
cnv$tissueid <- curationTissue[cnv$cellid, "unique.tissueid"]




# curationCell[803,] <- c("U-266", "U-266")
# rownames(curationCell)[803] <- "U-266"

mapInfo <- data.frame(gCSI.cellid = unique(sensitivity.info$cellid), 
                      unique.cellid = matchToIDTable(unique(sensitivity.info$cellid), curationCell, "gCSI.cellid", "unique.cellid"))
stopifnot(!anyNA(mapInfo[,2]))
sensitivity.info$cellid <- mapInfo[match(sensitivity.info$cellid, mapInfo[,1]),2]


mapInfo <- data.frame(gCSI.drugid = unique(sensitivity.info$drugid), 
                      unique.drugid = matchToIDTable(unique(sensitivity.info$drugid), curationDrug, "gCSI.drugid", "unique.drugid"))
stopifnot(!anyNA(mapInfo[,2]))
sensitivity.info$drugid <- mapInfo[match(sensitivity.info$drugid, mapInfo[,1]),2]



cellInfo <- as.data.frame(cellInfo)

cellinall <- unionList(rnaseq$cellid, sensitivity.info$cellid, cnv$cellid, mut$cellid)

newCells <- setdiff(cellinall, rownames(cellInfo))

newrows <- cellInfo[newCells,]
rownames(newrows) <- newCells
cellInfo <- rbind(cellInfo, newrows)

cellInfo$tissueid <- curationTissue[rownames(cellInfo), "unique.tissueid"]



z <- list()
z <- c(z,c(
  "rnaseq"=rnaseq,
  "cnv"=cnv,
  "mutation"=mut)
)

gCSI_2017 <- PharmacoSet(molecularProfiles=z,
                       name="gCSI",
                       cell=curationCell,
                       drug=curationDrug,
                       sensitivityInfo=sensitivity.info,
                       sensitivityRaw=raw.sensitivity,
                       sensitivityProfiles=sensitivityProfiles_2017,
                       curationCell=curationCell,
                       curationDrug=curationDrug,
                       curationTissue=curationTissue,
                       datasetType="sensitivity")

save(gCSI_2017, file=file.path(dir.prefix, "out/gCSI_2017.RData"))
    