# 
# annotated_samples <- readRDS("annotation_samples.rds")
# gse_gpls <- strsplit(names(annotated_samples), ".", fixed = T)
# orgAbbr <- c("Rattus norvegicus" = "rnorvegicus_gene_ensembl",
#              "Mus musculus" = "mmusculus_gene_ensembl")
# ensembl_release = "92"
# hkgDataDir <- "hkg-data"
# 
# #### Load Datasets ####
library(stringi)
library(data.table)
library(GEOquery)
# 
# all_gsets <- list()
# for (gse_gpl in gse_gpls) {
#   fns <- list.files(path = hkgDataDir, 
#              pattern = paste0("^", gse_gpl[1],".*_series_matrix.txt.gz$"),
#              full.names = T)
#   gsets <- list()
#   for (fn in fns) {
#     gpl <- stri_extract(fn, regex="GPL\\d+")
#     gset_raw <- getGEO(filename = fn, 
#                        GSEMatrix=TRUE, AnnotGPL=T, getGPL = T, destdir = hkgDataDir)
#     if (is.na(gpl))
#       gpl <- gset_raw@annotation
#     gsets[[gpl]] <- gset_raw
#   }
#   all_gsets[[gse_gpl[1]]] <- gsets
# }
# 
# #### Load Annotations ####
library(AnnotationHub)
library(ensembldb)
library(biomaRt)
library(limma)
library(splitstackshape)
library(genefilter)
# 
# ah <- query(AnnotationHub(), "EnsDb")
# orgAnno <- lapply(
#   names(orgAbbr),
#   function(o) {
#     queryResult <- query(x=ah, pattern=o)
#     queryResult[[grep(paste("Ensembl",ensembl_release,"EnsDb"), mcols(queryResult)$title, fixed=T)]]
#   }
# )
# names(orgAnno) <- tolower(names(orgAbbr))
# 
# orgMart <- lapply(
#   orgAbbr,
#   function(orgAbbr){
#     useEnsembl("ensembl", dataset = orgAbbr, version = ensembl_release)
#   }
# )
# names(orgMart) <- tolower(names(orgAbbr))


#### Analysis ####
source("housekeepr-profiling/housekeepr-functions-comparison.R")
load("housekeepr-profiling/comparison-preparation.RData")
ah <- AnnotationHub() #reload ah because cache paths may be corrupted
data.table::setDTthreads(1)

orgAbbr <- c("Rattus norvegicus" = "rnorvegicus_gene_ensembl",
             "Mus musculus" = "mmusculus_gene_ensembl")

orgAnno <- lapply(
  names(orgAbbr),
  function(o) {
    queryResult <- query(x=ah, pattern=o)
    queryResult[[grep(paste("Ensembl",ensembl_release,"EnsDb"), mcols(queryResult)$title, fixed=T)]]
  }
)
names(orgAnno) <- tolower(names(orgAbbr))

housekeepr_profiles <- list()

num_hkg <- 50
bootstrap_sample_size <- 12
bootstrap_replications <- 10000

# top tables of all data sets
tT_list <- list()

# keep track of organisms of each data set
dataset_organisms <- list()

# mapping of gene id -> paralog gene id
all_paralogs <- data.table(gene_id=character(0), paralog_gene_id=character(0))

# mapping of gene id -> "paralog group"
paralogs_groups <- data.table(gene_id=character(0), paralog_group=numeric(0), key="gene_id")

template_organism <- NULL

gset_dims <- list()

for (i in 1:length(annotated_samples)) {
  gse_gpl <- names(annotated_samples)[[i]]
 
  gse_gpls <- strsplit(gse_gpl, ".", fixed = T)[[1]]
  gse <- gse_gpls[1]
  gpl <- gse_gpls[2]
  
  print(sprintf("%s", gse_gpl))
  
  gset_raw <- all_gsets[[gse]]
  gsetEx <- getExpressions(gset_raw, gpl=gpl)[[1]]
  
  gset <- gsetEx$gset
  currentOrganism <- tolower(identifyOrganism(gset))
  
  dataset_organisms[gse] <- currentOrganism
  
  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  if (is.null(template_organism) | i == 1) {
    template_organism <- tolower(dataset_organisms[[1]])
    targetOrg <- template_organism
    target_orgDb <- orgAnno[[targetOrg]]
    target_org_mart <- orgMart[[targetOrg]]
  }
  
  if (currentOrganism != targetOrg) {
    source_orgDb <- orgAnno[[currentOrganism]]
    source_org_mart <- orgMart[[currentOrganism]]
  } else {
    source_orgDb <- target_orgDb
    source_org_mart <- target_org_mart
  }
  
  start.time <- Sys.time()
  tT <- calculateTopTable(annotated_samples, gse_gpl, gset)
  tT <- ensureGeneIdColumn(tT, currentOrganism, ah, ensembl_release, orgDb=source_orgDb)
  
  tT <- removeWithMissingEntrezId(tT, "gene_id")
  tT <- subset(tT, select=c("ID","adj.P.Val","P.Value", "logAveExpr","AveExpr", "variances","t","B","logFC", "gene_id"))
  
  tT <- unique(tT)
  res <- extractParalogs(tT, source_org_mart, all_paralogs, paralogs_groups, "gene_id")
  
  tT.paralogs <- res$tT.paralogs
  all_paralogs <- res$all_paralogs
  paralogs_groups <- res$paralogs_groups
  
  tT_target <- mapToTargetOrganismEntrez(tT.paralogs, currentOrganism, targetOrg, target_orgDb, source_org_mart, target_org_mart)
  
  if (currentOrganism != targetOrg) {
    # clean up
    tT_target <- removeWithMissingEntrezId(tT_target, "target_entrez")
    
    res <- extractParalogs(tT_target, target_org_mart, all_paralogs, paralogs_groups, "target_entrez")
    tT_target <- res$tT.paralogs
    all_paralogs <- res$all_paralogs
    paralogs_groups <- res$paralogs_groups
  }
  
  tT_aggr <- aggregateByParalogGroup(tT_target)
  tT_aggr <- calculateGeneScoresAndRanks(tT_aggr)
  
  housekeepr_profiles[[gse_gpl]] <- Sys.time() - start.time
  print(housekeepr_profiles[[gse_gpl]])
  gset_dims[[gse_gpl]] <- c(nrow(tT), ncol(exprs(gset)))
  # save all datasets to same list
  tT_list[[gse_gpl]] <- tT_aggr
}
saveRDS(gset_dims, "housekeepr-profiling/gset_dims.rds")
print("### Calculating the ranking of genes in each individual dataset")

res <- constructRankingMatrices(tT_list)

ranking_matrix <- res$ranking_matrix
ranking_matrix_long <- res$ranking_matrix_long

print("### Running bootstrapping with the following parameters:")
print(paste("- Number of HKG genes:", num_hkg))
print(paste("- Bootstrapping sample size:", bootstrap_sample_size))
print(paste("- Bootstrapping replications:", bootstrap_replications))

# genes_of_interest <- c(Gapdh = 'ENSRNOG00000018630', Sdha = 'ENSRNOG00000013331', B2m = 'ENSRNOG00000017123')
# par_groups_of_interest <- sapply(genes_of_interest, function(g) {
#   unique(sapply(tT_list, function(gs)gs[target_entrez == g, paralog_group]))
# })
# par_groups_of_interest
# Gapdh  Sdha   B2m 
# 12365   429 10440 

start.time <- Sys.time()
bootstrapRes <- calculateBootstrapMatrix(ranking_matrix, bootstrap_replications, bootstrap_sample_size, num_hkg)#, print_paralog_groups = par_groups_of_interest)
# paralog_group na_ratio num_na_factor rank_mean   rank_var total_rank
# 1:         12365        0             0   69.3161   26.13037         16
# 2:           429        0             0  842.4812  688.40753        233
# 3:         10440        0             0 2605.0518 1974.60652        774
bootstrap_matrix <- bootstrapRes$bootstrap_matrix
bootstrap_matrix[,merge_idx:=1:nrow(bootstrap_matrix)]
bootstrap_matrix_long <- bootstrapRes$bootstrap_matrix_long

annotateTargetParalogGroupsWithSymbols(ranking_matrix_long, paralogs_groups, "gene_id", target_orgDb,
                                       paralogue_groups_to_be_annotated = bootstrap_matrix[1:num_hkg,paralog_group])
annotateTargetParalogGroupsWithSymbols(bootstrap_matrix, paralogs_groups, "gene_id", target_orgDb,
                                       paralogue_groups_to_be_annotated = bootstrap_matrix[1:num_hkg,paralog_group])
annotateTargetParalogGroupsWithSymbols(bootstrap_matrix_long, paralogs_groups, "gene_id", target_orgDb,
                                       paralogue_groups_to_be_annotated = bootstrap_matrix[1:num_hkg,paralog_group])

candidate_hkg_ranking <- bootstrap_matrix[, .(paralog_group, symbols, na_ratio,
                                              rank_mean, rank_var,
                                              total_rank)]
housekeepr_profiles[["Bootstrapping"]] <- Sys.time() - start.time
print(housekeepr_profiles[["Bootstrapping"]])

#check that number of hkg given as input is not above the total number of genes
num_hkg <- ifelse(num_hkg>nrow(candidate_hkg_ranking),
                  nrow(candidate_hkg_ranking), num_hkg)

# data preparations for plotting
top_hkg <- candidate_hkg_ranking$paralog_group[1:num_hkg]

res <- list(bootstrap_ranking=bootstrap_matrix[1:num_hkg], 
            bootstrap_ranking_long=bootstrap_matrix_long[paralog_group%in%top_hkg], 
            ranking_long=ranking_matrix_long[paralog_group%in%top_hkg])

saveRDS(housekeepr_profiles, file = "housekeepr-profiling/´housekeepr_profiles.rds")