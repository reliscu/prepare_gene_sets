library(GSEABase)
library(parallel)

source("/home/rebecca/code/map_identifiers/map_identifiers_function.R")

broad_info <- function(idx, sets) {
  working_set <- sets[[idx]]
  legend <- c(
    SetID=GSEABase::setIdentifier(working_set),
    SetName=GSEABase::setName(working_set),
    Category=GSEABase::collectionType(working_set)@category,
    SetSize=length(GSEABase::geneIds(working_set)),
    Species=GSEABase::organism(working_set),
    PubMed=paste(GSEABase::pubMedIds(working_set), collapse=" | "),
    Description=GSEABase::description(working_set)
  )
  return(legend)
}

prep_broad_sets <- function(xml_file, version, out_dir){

  bsets <- GSEABase::getBroadSets(xml_file)
  broad_legend <- lapply(1:length(bsets), broad_info, bsets)
  broad_legend <- data.frame(do.call(rbind, broad_legend))
  broad_sets <- lapply(bsets, GSEABase::geneIds)
  names(broad_sets) <- broad_legend$SetName
  broad_sets <- broad_sets[broad_legend$SetSize>0]
  broad_legend <- broad_legend[broad_legend$SetSize>0,]
  save(broad_sets, broad_legend, file=paste0(file.path(out_dir, "broad_sets_v"), version, ".RData"))
  
}

prep_oldham_sets <- function(set_dir, oldham_legend, out_dir) {
  
  setwd(set_dir)
  set_files <- list.files()[list.files()%in%paste0(oldham_legend$SetID, ".csv")]
  set_id <- gsub(".csv", "", set_files)
  oldham_legend <- oldham_legend[match(set_id, oldham_legend$SetID),]
  if(!identical(set_id, oldham_legend$SetID)) {
    stop("Some sets in the legend file were not found in the gene set directory, or vice versa.")
  }
  oldham_sets <- lapply(1:length(set_files), function(idx) read.csv(set_files[idx], header=F)[,1])
  names(oldham_sets) <- oldham_legend$SetName
  save(oldham_sets, oldham_legend, file=file.path(out_dir, "oldham_sets.RData"))
  
}

map_sets <- function(sets, legend, mapping_tables_dir, out_dir, n_threads=15){
  
  ## For multiple identfiers, select first value:
  sets <- lapply(sets, function(x) sapply(strsplit(x, "///", fixed=T), "[", 1))
  sets <- lapply(sets, function(x) sapply(strsplit(x, ";", fixed=T), "[", 1))
  
  ## Clean up:
  sets <- lapply(sets, function(x) x[x!="---"])
  sets <- lapply(sets, trimws)
  
  ## Map to most recent gene symbol:
  sets_mapped <- mclapply(sets, function(x) {
    mapAlias2Symbol(features=data.frame(x), unique_id_col=1, mapping_tables_dir, keep_all_mappings=F, fill_NAs=T)[,2]
  }, mc.cores=n_threads)
  
  return(sets_mapped)
  
} ## prep_oldam_sets=function(){
