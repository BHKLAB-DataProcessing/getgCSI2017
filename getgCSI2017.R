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
#library(CoreGx)
library(Biobase)
library(reshape2)
library(abind)
library(data.table)
library(parallel)
library(CoreGx)
library(SummarizedExperiment)

verbose=FALSE
nthread=1
            
options(stringsAsFactors = FALSE)

myDirPrefix <- "/pfs/"
args = commandArgs(trailingOnly=TRUE)
rnaseq_select <- args
print(rnaseq_select)
rnaseq_results <- list()
ORCESTRA_ID = tail(rnaseq_select, n=1)

cnv_select <-  grep('cnv', rnaseq_select)
mutation_select <-  grep('mutation', rnaseq_select)
	  
tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
tools <- rnaseq_select[tools]
tools <- gsub("-", "_", tools)
transcriptome <- grep(pattern = 'Gencode|Ensembl', x = rnaseq_select)
transcriptome <- rnaseq_select[transcriptome]
tool_path = expand.grid(a = tools,b = transcriptome)
tool_path = paste0(tool_path$a, "_",tool_path$b)
	  
print(tool_path)

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
load(file.path(dir.prefix, "gcsi2017ProfileAssemble/profiles.RData"))
# load(file.path(dir.prefix, "getgCSI2017/gcsidrugpost.RData"))
load(file.path(dir.prefix, "gcsi2017RawSensitivity/raw.sensitivity.RData"))
load(file.path(dir.prefix, "gcsi2017RawSensitivity/gCSI_molData.RData"))
#load(file.path(dir.prefix, "downloadrnagCSI/gCSI_2017_molecprofile.RData"))



    
recomputed_2017 <- res
sensitivityProfiles_2017 <- data.frame("aac_recomputed" = NA, "ic50_recomputed"=NA, "ic50_published"=NA, "meanviability_published"=NA, "HS"=NA, "E_inf"=NA, "EC50"=NA)
sensitivityProfiles_2017[nrow(published.profiles),] <- NA
rownames(sensitivityProfiles_2017) <- rownames(published.profiles)

sensitivityProfiles_2017[rownames(recomputed_2017),"aac_recomputed"] <- as.numeric(recomputed_2017[,"aac_recomputed"])
sensitivityProfiles_2017[rownames(recomputed_2017),"ic50_recomputed"] <- as.numeric((recomputed_2017[,"ic50_recomputed"]))
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
drug_all <- drug_all[ , c("unique.drugid", "gCSI.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

curationCell <- cell_all[apply(!is.na(cell_all[,c("gCSI.cellid", "GNE.cellid")]),1,any),]
curationTissue <- cell_all[apply(!is.na(cell_all[,c("gCSI.cellid", "GNE.cellid")]),1,any),]
curationCell <- curationCell[ , c("unique.cellid", "gCSI.cellid", "GNE.cellid")]
curationTissue <- curationTissue[ , c("unique.tissueid", "gCSI.tissueid", "GNE.tissueid")]

rownames(curationTissue) <- curationCell[ , "unique.cellid"]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "gCSI.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "gCSI.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]

drugInfo <- data.frame("drugid"=curationDrug$unique.drugid)
rownames(drugInfo) <- drugInfo$drugid
drug_all <- drug_all[rownames(drugInfo),]
drugInfo[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]

## Only doing this for data added to the pset.

reps <- matchToIDTable(mut$cellid, curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
mut$cellid <- reps
mut$tissueid <- curationTissue[mut$cellid, "unique.tissueid"]


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


summarizeRnaSeq <- function (dir, 
                             features_annotation,
                             samples_annotation,
			      method) {
  library(Biobase)
  library(readr)
  library(tximport)
  
  load(features_annotation)
    
  tx2gene <- as.data.frame(cbind("transcript"=tx2gene$transcripts, "gene"=tx2gene$genes))
  
  files <- list.files(dir, recursive = TRUE, full.names = T)
  if(method=="kallisto"){
  resFiles <- grep("abundance.h5", files)
  }else{
  resFiles <- grep("quant.sf", files)
  }
  resFiles <- files[resFiles]
  length(resFiles)
  names(resFiles) <- basename(dirname(resFiles))
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  } else{
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)	  
  }
	  
  head(txi$counts[,1:5])
  dim(txi$counts)
	  
  xx <- txi$abundance
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- features_gene[featureNames(gene.exp),]
  pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
  annotation(gene.exp) <- "rnaseq"
  
  xx <- txi$counts
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- features_gene[featureNames(gene.count),]
  pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
  annotation(gene.count) <- "rnaseq"
  
  txii <- tximport(resFiles, type=method, txOut=T)
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  #remove non-coding transcripts in ensembl 	  
  rownames(txii$abundance) <-  gsub("\\..*","",rownames(txii$abundance))
  txii$abundance[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  txii$abundance <- txii$abundance [-which(rownames(txii$abundance) %in% missing_transcript),]
  }
  	  
  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  if(features_annotation == "/pfs/downAnnotations/Gencode.v33.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  featureNames(transcript.exp) <- gsub("\\|.*","",featureNames(transcript.exp))
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }else{
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
  annotation(transcript.exp) <- "isoform"
  
	  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  #remove non-coding transcripts in ensembl
  rownames(txii$counts) <-  gsub("\\..*","",rownames(txii$counts))
  txii$counts <- txii$counts [-which(rownames(txii$counts) %in% missing_transcript),]	  
  }	  
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  if(features_annotation == "/pfs/downAnnotations/Gencode.v33.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  featureNames(transcript.count) <- gsub("\\|.*","",featureNames(transcript.count))
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }else{
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
  annotation(transcript.count) <- "isoform"
  
  return(list("rnaseq"=gene.exp, 
              "rnaseq.counts"=gene.count, 
              "isoforms"=transcript.exp, 
              "isoforms.counts"=transcript.count))
}


rnaseq.sampleinfo <- read.csv(file="/pfs/downAnnotations/gCSI_rnaseq_meta.csv", stringsAsFactors=FALSE, row.names=1)
rnaseq.sampleinfo[ , "cellid"] <-  matchToIDTable(ids=rnaseq.sampleinfo[ , "Cell_line"], tbl=curationCell, column = "GNE.cellid", returnColumn = "unique.cellid")


 for (r in 1:length(tool_path)){
  print(tool_path[r])
  if (length(grep(pattern = 'Kallisto', x = tool_path[r])) > 0){
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    tdir = paste0("gcsi_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")  
    rnatool="kallisto"	  
  } else {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    tdir = paste0("gcsi_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
    rnatool="salmon"	  
  }
  
  
  if (length(grep(pattern = 'lift37', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"
  } else if (length(grep(pattern = 'v33', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33.annotation.RData"
  } else {
    annot = "/pfs/downAnnotations/Ensembl.v99.annotation.RData"
  }
    print(annot)
  
  print(tdir)
  rnaseq <- summarizeRnaSeq(dir=file.path(paste0(myDirPrefix, tdir, tool_path[r])),
                            features_annotation=annot,
                            samples_annotation=rnaseq.sampleinfo,
			    method = rnatool)
	 
  	 
  reps <- matchToIDTable(rnaseq$rnaseq$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$rnaseq$cellid <- reps
	 
  reps <- matchToIDTable(rnaseq$rnaseq.counts$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$rnaseq.counts$cellid <- reps

  reps <- matchToIDTable(rnaseq$isoforms$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$isoforms$cellid <- reps
  
  reps <- matchToIDTable(rnaseq$isoforms.counts$Cell_line, curationCell, "GNE.cellid", "unique.cellid")
  stopifnot(!anyNA(reps))
  rnaseq$isoforms.counts$cellid <- reps
  
  rnaseq_results <- c(rnaseq_results,c(
    rnaseq <- setNames(rnaseq,  paste0(tool,".", names(rnaseq)))
  )
  )
}

rnaseq_cellid_all <- pData(rnaseq_results[[1]])[,"cellid"]
reps <- matchToIDTable(rownames(cellInfo), curationCell, "gCSI.cellid", "unique.cellid")
stopifnot(!anyNA(reps))
rownames(cellInfo) <- reps


cellInfo <- as.data.frame(cellInfo)

cellinall <- CoreGx::.unionList(rnaseq_cellid_all, sensitivity.info$cellid, cnv$cellid, mut$cellid)

newCells <- setdiff(cellinall, rownames(cellInfo))

newrows <- cellInfo[newCells,]
rownames(newrows) <- newCells
cellInfo <- rbind(cellInfo, newrows)

cellInfo$tissueid <- curationTissue[rownames(cellInfo), "unique.tissueid"]
cellInfo$cellid <- rownames(cellInfo)

curationTissue <- curationTissue[rownames(cellInfo),]
curationCell <- curationCell[rownames(cellInfo),]

if (length(cnv_select) > 0){
  cnv_cells_id <- cnv$cellid
} else {
  cnv_cells_id <- c()
  cnv <- ExpressionSet()
  pData(cnv)$cellid <- character()
  pData(cnv)$batchid <- character()
  fData(cnv)$BEST <- vector()
  fData(cnv)$Symbol <- character()
  annotation(cnv) <- "CNV data was not selected for on ORCESTRA"
}
		 
if (length(mutation_select) > 0){
  mutation_cells_id <- mut$cellid
} else {
  mutation_cells_id <- c()
  mut <-  ExpressionSet()
  pData(mut)$cellid <- character()
  pData(mut)$batchid <- character()
  fData(mut)$BEST <- vector()
  fData(mut)$Symbol <- character()
  annotation(mut) <- "Mutation data was not selected for on ORCESTRA"
}


z <- list()

z <- c(z,c(
  rnaseq_results,
  "cnv"=cnv,
  "mutation" = mut)
)

sensitivity.info <- as.data.frame(sensitivity.info)

drugInfo <- drugInfo[unique(sensitivity.info$drugid),]
curationDrug <- curationDrug[rownames(drugInfo),]

		 		 
#add cellosaurus disease type to cell-info

disease <- cell_all$Cellosaurus.Disease.Type[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$Cellosaurus.Disease.Type <- disease
		 
#add cellosaurus assession to cell-info
assession <- cell_all$Cellosaurus.Accession.id[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$Cellosaurus.Accession.id <- assession
		 
#add pharmacodb id to cell-info
pdb <- cell_all$PharmacoDB.id[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$PharmacoDB.id <- pdb

#add study tissue id to cell_info
study_tissue <- cell_all$unique.tissueid.fromstudies[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$unique.tissueid.fromstudies <- study_tissue
		 
#add study cell-line type to cell_info
cell_type <- cell_all$CellLine.Type[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$CellLine.Type <- cell_type
		 
#add metastatic info to cell_info		 
metastatic <- cell_all$Metastatic[match(cellInfo$cellid, cell_all$unique.cellid)]
cellInfo$Metastatic <- metastatic
	

.converteSetToSE <- function(eSets) {
  
  SEfinal <- lapply(eSets,
         function(eSet){
             # Change rownames from probes to EnsemblGeneId for rna data type
             if (grepl("^rna$", Biobase::annotation(eSet))) {
               rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
             }
             
             # Build summarized experiment from eSet
             SE <- SummarizedExperiment::SummarizedExperiment(
               ## TODO:: Do we want to pass an environment for better memory efficiency?
               assays=S4Vectors::SimpleList(as.list(Biobase::assayData(eSet))
               ),
               # Switch rearrange columns so that IDs are first, probes second
               rowData=S4Vectors::DataFrame(Biobase::fData(eSet),
                                            rownames=rownames(Biobase::fData(eSet)) 
               ),
               colData=S4Vectors::DataFrame(Biobase::pData(eSet),
                                            rownames=rownames(Biobase::pData(eSet))
               ),
               metadata=list("experimentData" = eSet@experimentData, 
                             "annotation" = Biobase::annotation(eSet), 
                             "protocolData" = Biobase::protocolData(eSet)
               )
             )
             ## TODO:: Determine if this can be done in the SE constructor?
             # Extract names from expression set
             SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
             mDataType <- Biobase::annotation(eSet)
             eSets[[mDataType]] <- SE
         })
  #setNames(pSet@molecularProfiles, names(eSets))
  return(SEfinal)
}
		 
z <- .converteSetToSE(z)

cells_keep <- unique(c(rnaseq_cellid_all, sensitivity.info$cellid, cnv_cells_id, mutation_cells_id))
		 
cellInfo <- cellInfo[cells_keep,]
curationCell <- curationCell[cells_keep,]
curationTissue <- curationTissue[cells_keep,]

gCSI_2017 <- PharmacoGx::PharmacoSet(molecularProfiles=z,
                       name="gCSI",
                       cell=cellInfo,
                       drug=drugInfo,
                       sensitivityInfo=sensitivity.info,
                       sensitivityRaw=raw.sensitivity,
                       sensitivityProfiles=sensitivityProfiles_2017,
                       curationCell=curationCell,
                       curationDrug=curationDrug,
                       curationTissue=curationTissue,
                       datasetType="sensitivity")
                              

gCSI_2017@annotation$version <- 2
saveRDS(gCSI_2017, file="/pfs/out/gCSI_2017.rds")
		 
dataset <- "gCSI"	
		 
#output ORCESTRA_ID and Pachyderm commit id
write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F) 

    
