
# #### Load Datasets ####
library(stringi)
library(data.table)
library(GEOquery)
# 
# 
# #### Load Annotations ####
library(AnnotationHub)
library(ensembldb)
library(biomaRt)
library(limma)
library(splitstackshape)
library(genefilter)
# 


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

num_hkg <- 100
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

saveRDS(tT_list, 'housekeepr-profiling/tT_list.rds')

# OLD----

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
# housekeepr_profiles[["Bootstrapping"]] <- Sys.time() - start.time
# print(housekeepr_profiles[["Bootstrapping"]])

#check that number of hkg given as input is not above the total number of genes
num_hkg <- ifelse(num_hkg>nrow(candidate_hkg_ranking),
                  nrow(candidate_hkg_ranking), num_hkg)

# data preparations for plotting
top_hkg <- candidate_hkg_ranking$paralog_group[1:num_hkg]

res_old <- list(bootstrap_ranking=bootstrap_matrix[1:num_hkg],
            bootstrap_ranking_long=bootstrap_matrix_long[paralog_group%in%top_hkg],
            ranking_long=ranking_matrix_long[paralog_group%in%top_hkg])


# CV_FC----

# now for new score
invisible(lapply(tT_list, function(x) setnames(x, old = c('rank_CV_FC', 'gene_rank'), c('gene_rank', 'rank_old'))))

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
# housekeepr_profiles[["Bootstrapping"]] <- Sys.time() - start.time
# print(housekeepr_profiles[["Bootstrapping"]])

#check that number of hkg given as input is not above the total number of genes
num_hkg <- ifelse(num_hkg>nrow(candidate_hkg_ranking),
                  nrow(candidate_hkg_ranking), num_hkg)

# data preparations for plotting
top_hkg <- candidate_hkg_ranking$paralog_group[1:num_hkg]

res_CV_FC <- list(bootstrap_ranking=bootstrap_matrix[1:num_hkg],
                bootstrap_ranking_long=bootstrap_matrix_long[paralog_group%in%top_hkg],
                ranking_long=ranking_matrix_long[paralog_group%in%top_hkg])

# rank_product----

invisible(lapply(tT_list, function(x) setnames(x, old = c('gene_rank', 'rank_product'), c('rank_CV_FC', 'gene_rank'))))

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
# housekeepr_profiles[["Bootstrapping"]] <- Sys.time() - start.time
# print(housekeepr_profiles[["Bootstrapping"]])

#check that number of hkg given as input is not above the total number of genes
num_hkg <- ifelse(num_hkg>nrow(candidate_hkg_ranking),
                  nrow(candidate_hkg_ranking), num_hkg)

# data preparations for plotting
top_hkg <- candidate_hkg_ranking$paralog_group[1:num_hkg]

res_rank_product <- list(bootstrap_ranking=bootstrap_matrix[1:num_hkg],
                bootstrap_ranking_long=bootstrap_matrix_long[paralog_group%in%top_hkg],
                ranking_long=ranking_matrix_long[paralog_group%in%top_hkg])

repeated_boot <- lapply(1:10, function(i){
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
  
  return(copy(bootstrap_matrix[1:num_hkg][, run := i]))
})
repeated_boot_dt <- rbindlist(repeated_boot)
library(ggplot2)
ggplot(repeated_boot_dt, aes(x=reorder(symbols, total_rank), y=total_rank)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_x_discrete(breaks=repeated_boot_dt$symbols[Reduce(`|`, sapply(experimentally_validated, grepl, x = repeated_boot_dt$symbols, simplify = F))])

saveRDS(res_old, 'housekeepr-profiling/res_old_score.rds')
saveRDS(res_CV_FC, 'housekeepr-profiling/res_CV_FC.rds')
saveRDS(res_rank_product, 'housekeepr-profiling/res_rank_product.rds')

# Check Results ----

library(data.table)
tT_list <- readRDS('housekeepr-profiling/tT_list.rds')

tT_list$GSE100235.GPL17117[AveExpr < 0 & variances > 0][order(gene_rank)][1, .(target_entrez, variances, variances_unscaled, AveExpr, AveExpr_unscaled, score_CV_FC, rank_CV_FC, rank_product, score, gene_rank)]

res_old <- readRDS('housekeepr-profiling/res_old_score.rds')
res_CV_FC <- readRDS('housekeepr-profiling/res_CV_FC.rds')
res_rank_product <- readRDS('housekeepr-profiling/res_rank_product.rds')

results_bootstrap <- rbindlist(list(old = res_old$bootstrap_ranking, CV_FC = res_CV_FC$bootstrap_ranking, rank_product = res_rank_product$bootstrap_ranking), idcol = 'rank')

results_bootstrap[paralog_group %in% results_bootstrap[,.(keep = .N == 1), by=paralog_group][keep == TRUE, paralog_group], .(rank, total_rank, symbols)]

experimentally_validated <- c('Ubb', 'Fth1', 'Fau', 'Ppial4d', 'Cst3', 'Rps23', 'Actb', 'Tuba1b', 'Rps3', 'Rplp2')
experimentally_validated_symbols <- unique(results_bootstrap$symbols[Reduce(`|`, sapply(experimentally_validated, grepl, x = results_bootstrap$symbols, simplify = F))])

library(ggplot2)
library(ggrepel)
# results_bootstrap[, .(rank_diff = abs(ifelse(length(diff(total_rank))==0, NA, diff(total_rank)))), by=.(symbols)]
ggplot(results_bootstrap[, .(rank_diff = abs(diff(total_rank))), by=.(symbols)][Reduce(`|`, sapply(experimentally_validated, grepl, x = symbols, simplify = F))], aes(y=reorder(symbols, rank_diff), x= rank_diff)) +
  geom_point()

res_merged <- merge(res_old$bootstrap_ranking, res_new$bootstrap_ranking, by=c('paralog_group', 'symbols', 'gene_ids'), suffixes=c('.old', '.new'))
ggplot(res_merged, aes(x=total_rank.new, y=total_rank.old, label=ifelse(Reduce(`|`, sapply(experimentally_validated, grepl, x = symbols, simplify = F)), symbols, ""))) +
  geom_point() + geom_text_repel(min.segment.length = 0, size=7) + theme(text=element_text(size=20)) + labs(title="Genes shared in the Top100")

first50_no_overlap <- results_bootstrap[paralog_group %in% results_bootstrap[,.(keep = .N == 1), by=paralog_group][keep == TRUE, paralog_group], 
                                        .(rank, total_rank, symbols)][total_rank<=50]

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#CC79A7")
ggplot(results_bootstrap,
       aes(y=total_rank, x=reorder(symbols, total_rank), color=rank, shape=rank, label=ifelse(Reduce(`|`, sapply(experimentally_validated, grepl, x = results_bootstrap$symbols, simplify = F)), symbols, ""))) + 
  geom_point(size=5, alpha=.8) + 
  # geom_text_repel(min.segment.length = unit(0, 'lines'), size=5) +
  labs(title="Ranking Strategies", x='Gene', y='Total Rank') +
  theme_bw() + 
  theme(text=element_text(size=14), axis.text.x = element_text(angle = 40, hjust=1)) + #, axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  scale_x_discrete(breaks=experimentally_validated_symbols)

results_wo_Mt <- results_bootstrap[!startsWith(symbols, 'Mt-')]
results_wo_Mt[, total_rank_wo_MT := rank(total_rank), by =rank]
ggplot(results_wo_Mt[symbols %in% experimentally_validated_symbols], aes(x=factor(symbols, levels=results_wo_Mt[symbols %in% experimentally_validated_symbols & rank == 'old'][order(total_rank_wo_MT), symbols]), y =total_rank_wo_MT, color = rank, shape=rank)) + 
  geom_point(size=5, alpha=.8) + 
  # geom_text_repel(min.segment.length = unit(0, 'lines'), size=5) +
  labs(title="Ranks for Validated Genes (Computed after Removing Mitochondrial Gene)", x='Gene', y='Total Rank') +
  theme_bw() + 
  theme(text=element_text(size=14), axis.text.x = element_text(angle = 40, hjust=1)) + #, axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_colour_manual(values=cbPalette)

dcast(results_wo_Mt[, .(rank, paralog_group, gene_ids, symbols, total_rank_wo_MT)], ... ~ rank, value.var = 'total_rank_wo_MT')[symbols %in% experimentally_validated_symbols]

