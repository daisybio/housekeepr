#### HELPER FUNCTIONS #####
#------------------------------------------------------------------------------#
# Standardization from 0 to 1 #
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Calculating the variances of all samples
RowVar <- function (x,na.rm = TRUE) 
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}
#------------------------------------------------------------------------------#

downloadFromGEO <- function(gse, 
                            GSEMatrix=TRUE, 
                            AnnotGPL=T, 
                            destdir="hkg-data",
                            attempts=10, 
                            sleepBetween=2) {
  gset_raw <- NULL
  attempt <- 1
  timed_out <- TRUE
  while(timed_out && attempt <= attempts ) {
    if (attempt > 1)
      print(paste("Retrying download:", paste(attempt,attempts,sep="/")))
    Sys.sleep(sleepBetween)
    attempt <- attempt + 1
    tryCatch({
      withTimeout({
        gset_raw <- getGEO(GEO=gse, GSEMatrix=GSEMatrix, AnnotGPL=AnnotGPL, destdir = destdir)
      }, timeout=30*60, onTimeout="error")
      timed_out <- FALSE
    }, TimeoutException=function(e) {
      print(sprintf("Caught error downloading %s: %s", gse, e$message))
      timed_out <- TRUE
    })
  }
  list(gset_raw=gset_raw, timed_out=timed_out)
}

getExpressions <- function(gset_raws, gpls) {
  res <- lapply(gpls, function(gpl) {
    gset_raw <- gset_raws[[gpl]]
    if(class(gset_raw) == "GDS") {
      gset <- GDS2eSet(gset_raw)
      ex <- exprs(gset)
    } else if (class(gset_raw) == "ExpressionSet") {
      gset <- gset_raw
      ex <- exprs(gset)
    } else {
      ex <- NULL
      print(paste('Unknown data type'))
    }
    
    if (!is.null(ex) && nrow(ex) == 0) {
      print(paste('Could not load expression values'))
      ex <- NULL
    }
    list(gset=gset, gpl=gpl, ex=ex)
  })
  names(res) <- gpls
  res
}


chooseGPLandGetExpression <- function(gset_raw) {
  if (is.vector(gset_raw)) {
    if (length(gset_raw) == 0)
      return(NULL)
    if (class(gset_raw[[1]]) != "ExpressionSet")
      return(NULL)
    
    # check if more than 1 GPL is available, and take the first one
    gpllist <- as.character(lapply( gset_raw , function(x) `@`(x, annotation)))
    if (length(gset_raw) > 1) idx <- grep(gpllist[1], attr(gset_raw, "names")) else idx <- 1
    gset <- gset_raw[[idx]]
    
    #get expression matrix
    ex <- exprs(gset)
  }
  else if(class(gset_raw) == "ExpressionSet") {
    gset <- gset_raw
    #get expression matrix
    ex <- exprs(gset)
  }
  else if(class(gset_raw) == "GDS") {
    gset <- GDS2eSet(gset_raw,do.log2=T)
    ex <- exprs(gset)
  } else {
    print(paste('Skipping: Unknown data type'))
    return(NULL)
  }
  
  if (nrow(ex) == 0) {
    print(paste('Skipping: Could not load expression values'))
    return(NULL)
  }
  
  return(list(gset=gset, ex=ex))
}

identifyOrganism <- function(gset_raw) {
  # find organism of data set
  currentOrganism <- NULL
  if (is.list(gset_raw) && length(gset_raw) > 0)
    gset <- gset_raw[[1]]
  else if (class(gset_raw) == "ExpressionSet")
    gset <- gset_raw
  
  if (is.null(gset)) {
    print("Failed to find organism of data set. Skipping")
    return(NULL)
  }
  if (.hasSlot(gset, "header"))
    currentOrganism <- gset@header$sample_organism
  else if (.hasSlot(gset, "phenoData") && .hasSlot(gset@phenoData, "data"))
    currentOrganism <- as.character(gset@phenoData@data$organism_ch1[1])
  
  if (is.null(currentOrganism)) {
    print("Failed to find organism of data set. Skipping")
  }
  return(currentOrganism)
}

calculateTopTable <- function(GSE_samples_annotations, GSE, gset) {
  print("calculateTopTable")
  # group names for all samples
  gsms <- GSE_samples_annotations[[GSE]][sampleNames(gset)]
  sml <- unlist(gsms)
  
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(detail='Formatting and normalizing data set')
  # Make sure to unlog the expression before calculating the variance
  qx <- as.numeric(quantile(exprs(gset), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { exprs(gset)[which(exprs(gset) <= 0)] <- NaN
  exprs(gset) <- log2(exprs(gset)) }
  variances <- RowVar(exprs(gset))
  # if (max(exprs(gset), na.rm = T) > 10000){
  #   variances <- RowVar(exprs(gset))
  # } else {
  #   exprs(gset) <- 2^exprs(gset)
  #   variances <- RowVar(exprs(gset))
  # }
  # # Reapply log transform to the expression for the diff analysis afterwards
  # exprs(gset) <- log2(exprs(gset))
  
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(detail='Identifying differential genes and calculating fold changes')
  # prepare data for differential analysis
  # if (isolate(reactive_vals$use_normfinder)) {
  #   norm_frame <- rbind(as.data.frame(exprs(gset)), as.numeric(sml))
  #   rownames(norm_frame)[nrow(norm_frame)] <- "group "
  #   norm_res <- Normfinder(norm_frame) 
  #   print("done Normfinder")
  #   browser()
  # } else {
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(diffExpr=G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- as.data.table(topTable(fit2, adjust="fdr", sort.by="A", number=Inf), keep.rownames=TRUE)
    tT[, variances:=variances]
    # The AveExpr is the log mean expression, so renaming the col, and calculating the
    # mean expression without the log
    setnames(tT, old = "AveExpr", new = "logAveExpr")
    # tT[, AveExpr:= 2^logAveExpr]
  # }
  return(tT)
}

findEnsemblGeneIdColumn <- function(tT, currentOrganism, orgDb, targetGeneIdColName) {
  print("findEnsemblGeneIdColumn")
  #TODO: check if more than one column makes sense   
  #TODO: look for SPOT_ID for ensembl transcript id
  #for (candidateColNamePattern in c("mrna_assignment", "ensembl")) {
    
  col.idx <- grep("mrna_assignment", tolower(colnames(tT)))
  if (length(col.idx) == 0)
    return(list(tT=tT, gene.id.col=targetGeneIdColName))
#      next()
  else
    col.idx <- colnames(tT)[col.idx]
    #changing + to * in regex to ensure also human ens ids
  ensemblIds <- as.data.table(stri_extract_all(tT[, get(col.idx)], 
                                               regex="ENS[A-Z]*G\\d{11}", 
                                               simplify=T))
  ensemblIds[, ID:=tT[, merge_idx]]
  ensemblIds <- melt(ensemblIds, direction="long", variable.factor=F, id.vars="ID")
  setnames(ensemblIds, "value", "ENSEMBL_ID")
  ensemblIds[, variable:=NULL]
  ensemblIds <- ensemblIds[!is.na(ENSEMBL_ID) & ENSEMBL_ID != ""]

  # verify that these ensembl gene ids are still valid and not "retired" (according to ensembl)
  tryCatch({
    mapped <- ensembldb::select(orgDb, keys=ensemblIds[!is.na(ENSEMBL_ID),ENSEMBL_ID], keytype="GENEID", columns="SYMBOL") %>%
      dplyr::distinct(GENEID) %>%
      dplyr::pull(GENEID)
    ensemblIds <- ensemblIds[ENSEMBL_ID%in%mapped]
  }, error=function(e) {
    print(e)
    ensemblIds <- data.table(ID=character(0), ENSEMBL_ID=character(0))
  })
  
  
  tT <- merge(tT, ensemblIds, by.x="merge_idx", by.y="ID", suffixes=c(".x", ""), all.x=TRUE, allow.cartesian=TRUE)
  
#  if (targetGeneIdColName %in% names(tT)) {
#      tT[is.na(get(targetGeneIdColName)) & !is.na(ENSEMBL_ID) & ENSEMBL_ID != "", (targetGeneIdColName):=ENSEMBL_ID]
#      tT[, ENSEMBL_ID:=NULL]
#    } else {
  setnames(tT, "ENSEMBL_ID", targetGeneIdColName)
#    }
  
  print(paste(length(unique(ensemblIds$ID)),"probe ids ->", length(unique(ensemblIds$ENSEMBL_ID)) , "ensemble gene ids" ,
              paste0("(",countNumberProbesAnnotated(tT, "ID", targetGeneIdColName),
                     "/", length(unique(tT$ID)), " probe ids annotated)")
  ))
#  }
  list(tT=tT, gene.id.col=targetGeneIdColName)
}

countNumberProbesAnnotated <- function(tT, idColName, targetGeneIdColName) {
  nrow(tT[, length(na.omit(get(targetGeneIdColName))), by=idColName][V1 > 0])
}

mapEnsemblTranscriptToGeneId <- function(tT, currentOrganism, orgDb, targetGeneIdColName, mergeIdxWithMissingGeneIds=NULL) {
  print("mapEnsemblTranscriptToGeneId")
#  for (candidateColNamePattern in c("accession_string", "mrna_assignment")) {

  if(!is.null(mergeIdxWithMissingGeneIds) & length(mergeIdxWithMissingGeneIds) == 0)
    return(list(tT=tT, gene.id.col=targetGeneIdColName))
  
  if (is.null(mergeIdxWithMissingGeneIds)) {
    mergeIdxWithMissingGeneIds <- tT$merge_idx
    mergeIdxRows <- mergeIdxWithMissingGeneIds
  } else {
    mergeIdxRows <- match(mergeIdxWithMissingGeneIds, tT$merge_idx)
  }
  
  # find ensembl transcript column
  ensembltrans.id.col <- grep("accession_string", tolower(colnames(tT)))
  if (length(ensembltrans.id.col) == 0){
    ensembltrans.id.col <- grep("mrna_assignment", tolower(colnames(tT)))
    if (length(ensembltrans.id.col) == 0)
      return(list(tT=tT, gene.id.col=targetGeneIdColName))
    else 
      ensembltrans.id.colName <- colnames(tT)[ensembltrans.id.col]
  } else
    ensembltrans.id.colName <- colnames(tT)[ensembltrans.id.col]
  #changing + to * in regex to ensure also human ens ids
  ensemblIds <- as.data.table(stri_extract_all(tT[mergeIdxRows, get(ensembltrans.id.colName)], 
                                            regex="ENS[A-Z]*T\\d{11}", 
                                            simplify=T))
  ensemblIds[, ID:=tT[mergeIdxRows, ID]]
  ensemblIds <- melt(ensemblIds, direction="long", variable.factor=F, id.vars="ID")
  setnames(ensemblIds, "value", "ENSEMBL_ID")
  ensemblIds[, variable:=NULL]
  ensemblIds <- ensemblIds[!is.na(ENSEMBL_ID) & ENSEMBL_ID != ""]
  
  
  
  tT <- merge(tT, ensemblIds, by.x="ID", by.y="ID", suffixes=c(".x", ""), all.x=TRUE)
  
  # map ensembl transcript to ensembl gene id
  tryCatch({
    mapped <- ensembldb::select(orgDb, keys=tT[!is.na(ENSEMBL_ID),ENSEMBL_ID], keytype="TXID", columns="GENEID")
    tT <- merge(tT, mapped, by.x="ENSEMBL_ID", by.y="TXID", allow.cartesian=TRUE, all.x=TRUE)
    
    if (targetGeneIdColName %in% names(tT)) {
      tT[is.na(get(targetGeneIdColName)) & !is.na(GENEID) & GENEID != "", (targetGeneIdColName):=GENEID]
      tT[, GENEID:=NULL]
    } else {
      setnames(tT, "GENEID", targetGeneIdColName)
      gene.id.col <- targetGeneIdColName
    }
    
    print(paste(length(unique(ensemblIds$ID)), "probe ids ->", length(unique(ensemblIds$ENSEMBL_ID)),"ensembl transcript ids;", 
                length(unique(mapped$TXID)),"ensemble transcript ids ->", length(unique(mapped$GENEID)),"ensemble gene ids",
                paste0("(",countNumberProbesAnnotated(tT, "ID", targetGeneIdColName),
                       "/", length(unique(tT$ID)), " probe ids annotated)")))
    
  }, error=function(e){print(e)})
#  }
  return(list(tT=tT, gene.id.col=targetGeneIdColName))
}

mapEntrezIdToEnsemblGeneId <- function(tT, currentOrganism, orgDb, targetGeneIdColName, mergeIdxWithMissingGeneIds=NULL) {
  print("mapEntrezIdToEnsemblGeneId")
  
  for (candidateColNamePattern in c("^gene.*id$", "^gene$", "^locus.*link.*id$", "^gene_assignment$", "^entrez.*gene.*id$")) {# ".*entrez.*gene.*id.*")) {
    
    if(!is.null(mergeIdxWithMissingGeneIds) & length(mergeIdxWithMissingGeneIds) == 0)
      return(list(tT=tT, gene.id.col=targetGeneIdColName))
    
    if (is.null(mergeIdxWithMissingGeneIds)) {
      mergeIdxWithMissingGeneIds <- tT$merge_idx
      mergeIdxRows <- mergeIdxWithMissingGeneIds
    } else {
      mergeIdxRows <- match(mergeIdxWithMissingGeneIds, tT$merge_idx)
    }
    
    col.idx <- grep(candidateColNamePattern,tolower(colnames(tT)))
    
    if (length(col.idx) == 0)
      next()
    else
      col.name <- colnames(tT)[col.idx]
    
    # TODO: for gene_assignment this parsing did not always work
    # old sep="\\s*//+\\s*", 
    # new sep = "/"
    entrezIds <- cSplit(
      indt=na.omit(tT[mergeIdxRows, mget(c("ID",col.name))], cols="ID"),
      splitCols=col.name,
      sep="\\s*//+\\s*",
      type.convert=F,
      fixed=F,
      direction="long"
    )[grep("\\d+", get(col.name))]
    entrezIds <- unique(na.omit(entrezIds))
    
    colnames(entrezIds) <- c("ID", "ENTREZ_ID")
    
    # map entrez IDs to ensembl gene ids
    if (!nrow(entrezIds) == 0){
    tryCatch({
      mapped <- ensembldb::select(orgDb, keys=entrezIds$ENTREZ_ID, keytype="ENTREZID", columns="GENEID")
      mapped$ENTREZID <- as.character(mapped$ENTREZID)
      entrezIds2 <- merge(entrezIds, mapped, by.x="ENTREZ_ID", by.y="ENTREZID", allow.cartesian=TRUE)
      tT[, ID := as.character(ID)]
      entrezIds2[, ID := as.character(ID)]
      tT <- merge(tT, entrezIds2, by.x="ID", by.y="ID", suffixes=c(".x", ""), allow.cartesian=TRUE, all.x=TRUE)
      if (targetGeneIdColName %in% names(tT)) {
        tT[is.na(get(targetGeneIdColName)) & !is.na(GENEID) & GENEID != "", (targetGeneIdColName):=GENEID]
        tT[, GENEID:=NULL]
      } else {
        setnames(tT, "GENEID", targetGeneIdColName)
      }
      mergeIdxWithMissingGeneIds <- getMergeIdxWithMissingGeneId(tT, targetGeneIdColName)
      
      
      print(paste(length(unique(entrezIds$ID)),"probe ids ->", length(unique(entrezIds$ENTREZ_ID)), "entrez gene ids;",
                  length(unique(entrezIds2$ENTREZ_ID)), "entrez gene ids ->", length(unique(entrezIds2$GENEID)) , "ensemble gene ids",
                  paste0("(",countNumberProbesAnnotated(tT, "ID", targetGeneIdColName),
                         "/", length(unique(tT$ID)), " probe ids annotated)")))
    }, error=function(e) {print(e)})}
  }
  return(list(tT=tT, gene.id.col=targetGeneIdColName))
}

getMergeIdxWithMissingGeneId <- function(tT, geneIdCol) {
  if (!is.null(geneIdCol) && geneIdCol %in% names(tT))
    tT[is.na(tT[,get(geneIdCol)]) | tT[,get(geneIdCol)] == "", merge_idx]
  else
    tT$merge_idx
}

ensureGeneIdColumn <- function(tT, currentOrganism, ah, ensembl_release, orgDb=NULL) {
  print("ensureGeneIdColumn")
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(detail='Identifying and mapping gene annotations to ensembl gene ids')
  
  if (is.null(orgDb)) {
    queryResult <- query(x=query(ah, "EnsDb"), pattern=currentOrganism)
    orgDb <- queryResult[[grep(paste("Ensembl", ensembl_release, "EnsDb"), mcols(queryResult)$title, fixed=T)]]
  }
  
  tT[,merge_idx:=1:nrow(tT)]
  
  targetGeneIdColName <- "ENSEMBLGENEID_HOUSEKEEPR"
  tmp <- findEnsemblGeneIdColumn(tT, currentOrganism, orgDb, targetGeneIdColName)
  
  tmp <- mapEnsemblTranscriptToGeneId(
    tmp$tT, 
    currentOrganism, 
    orgDb,
    targetGeneIdColName, 
    getMergeIdxWithMissingGeneId(tmp$tT, tmp$gene.id.col)
  )
  
  tmp <- mapEntrezIdToEnsemblGeneId(
    tmp$tT, 
    currentOrganism, 
    orgDb,
    targetGeneIdColName, 
    getMergeIdxWithMissingGeneId(tmp$tT, tmp$gene.id.col)
  )
  
  # if (is.null(gene.id.col)) {
  #   # ensembl id
  #   # Find Gene ID column name
  #   gene.id.col <- grep("^ensembl_id",tolower(colnames(tT)))
  #   if (length(gene.id.col) == 0)
  #     gene.id.col <- NULL
  #   else
  #     gene.id.col <- colnames(tT)[gene.id.col]
  # }
  
  # # try extracting entrez ID from gene_assignment column / ensembl transcript id
  # if (is.null(gene.id.col) && ("gene_assignment" %in% names(tT))) {
  #   tmp <- copy(tT)
  #   f <- function(x, idx){
  #     if (idx > length(x))
  #       return(list(NA))
  #     list(trimws(x[[idx]]))
  #   }
  #   tmp[, ENSEMBL_ID:=as.character(str_extract(pattern="ENSRNOT[[:digit:]]+", string=tmp$gene_assignment))]
  # 
  #   # verify that we extract indeed gene symbols by trying to get entrez IDs from symbols
  #   symbols <- unique(na.omit(tmp$ENSEMBL_ID))
  #   mapped <- select(orgDb,
  #                    keys=symbols,
  #                    columns=c("ENTREZID"),
  #                    keytype="TXID")
  #   mappedGrouped <- na.omit(mapped) %>% dplyr::group_by(TXID) %>% dplyr::summarise(n=n())
  #   if (nrow(mappedGrouped) >= 250 || (nrow(mappedGrouped)/length(symbols) >= .50)) {
  #     tT[, ENSEMBL_ID:=tmp$ENSEMBL_ID]
  #     tT <- merge(tT, mapped, by.x="ENSEMBL_ID", by.y="TXID")
  #     gene.id.col <- "ENTREZID"
  #   }
  # }
  # replace tmp$gene.id.col by tmp$tT[[tmp$gene.id.col]] 
  if (is.null(tmp$tT[[tmp$gene.id.col]] )) {
    print("Could not find gene IDs of data set. Skipping.")
    return(NULL)
    #return(list(tT=tT, gene.id.col=NULL))
  }
  
  setnames(tmp$tT, tmp$gene.id.col, "gene_id")
  tmp$tT[, gene_id := as.character(tmp$tT$gene_id)]
  
  return(tmp$tT)
}

removeWithMissingEntrezId <- function(tT, colName) {
  tT <- na.omit(tT, cols=colName)
  #tT[, get(colName) != ""]
  tT <- tT[tT[,get(colName) != ""]]
  tT <- cSplit(tT, splitCols=colName,sep="//+", fixed=F,direction="long", type.convert = F)
  tT
}

replaceInf <- function(x) {
  if (is.infinite(x)) return(NA) else return(x)
}

extractParalogs <- function(tT, org_mart, all_paralogs, paralogs_groups, geneIdColumn) {
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(detail='Aggregating identical gene copies')
  print(paste("Identifying paralogue genes using column", geneIdColumn))
  attrs <- listAttributes(org_mart)
  
  attrTargetGene <- attrs[grep("_paralog_ensembl_gene$", attrs$name),"name"]
  attrTargetGeneName <- attrs[grep("_paralog_associated_gene_name$", attrs$name),"name"]
  attrSourcePerc <- attrs[grep("_paralog_perc_id$", attrs$name),"name"]
  attrTargetPerc <- attrs[grep("_paralog_perc_id_r1$", attrs$name),"name"]
  
  uniqueGeneIds <- unique(tT[,get(geneIdColumn)])
  known_paralogs <- all_paralogs[gene_id %in% uniqueGeneIds]
  setnames(known_paralogs, 1:2, c("ensembl_gene_id", "paralog_gene_id"))
  
  unknownGeneIds <- setdiff(uniqueGeneIds, known_paralogs$ensembl_gene_id)
  if (length(unknownGeneIds) > 0) {
    new_paralogs <- na.omit(data.table(getBM(attributes = c("ensembl_gene_id", 
                                                            "external_gene_name",
                                                            attrTargetGene,
                                                            attrTargetGeneName,
                                                            attrSourcePerc,
                                                            attrTargetPerc),
                                             filters = "ensembl_gene_id",
                                             values = unknownGeneIds,
                                             mart = org_mart)))
    new_paralogs <- new_paralogs[
      get(attrSourcePerc)==100 & get(attrTargetPerc)==100, (c("ensembl_gene_id", attrTargetGene)), with=F]
    setnames(new_paralogs, 1:2, c("ensembl_gene_id", "paralog_gene_id"))
    idsWithParalogs <- unique(new_paralogs$ensembl_gene_id)
    idsWithoutParalogs <- setdiff(unknownGeneIds, idsWithParalogs)
    new_trivial_paralogs <- data.table(ensembl_gene_id=unknownGeneIds, paralog_gene_id=unknownGeneIds)
    all_paralogs <- unique(rbindlist(list(all_paralogs, new_paralogs, new_trivial_paralogs), use.names = F))
    #all_paralogs <- unique(rbindlist(list(all_paralogs, new_paralogs, new_paralogs[, c(2,1)], new_trivial_paralogs)))
    
    print(paste("Annotating probes with paralogue groups using column", geneIdColumn))
    currentMaxGroup <- if(nrow(na.omit(paralogs_groups)) == 0) 0 else max(paralogs_groups$paralog_group, na.rm=T)
    new_paralog_groups <- merge(new_paralogs, paralogs_groups, by.x="paralog_gene_id", by.y="gene_id", all.x=T)
    new_paralog_groups[, paralog_group:=(if (any(!is.na(paralog_group))) min(paralog_group, na.rm=T) else Inf), by=ensembl_gene_id]
    todo <- unique(new_paralog_groups[is.infinite(paralog_group), ensembl_gene_id])
    while (length(todo) > 0) {
      gene_with_paralogs <- todo[1]
      paralogs <- new_paralog_groups[ensembl_gene_id == gene_with_paralogs, paralog_gene_id]
      
      group <- currentMaxGroup+1
      new_paralog_groups[ensembl_gene_id %in% c(gene_with_paralogs,paralogs), paralog_group:=group]
      currentMaxGroup <- group
      
      todo <- unique(new_paralog_groups[is.infinite(paralog_group), ensembl_gene_id])
    }
    
    paralogs_groups <- unique(rbindlist(
      list(
        subset(paralogs_groups, !(gene_id %in% new_paralog_groups$paralog_gene_id) & !(gene_id %in% new_paralog_groups$ensembl_gene_id)),
        new_paralog_groups[, .(paralog_gene_id, paralog_group)],
        new_paralog_groups[, .(ensembl_gene_id, paralog_group)]), use.names = F
    ))
    
    if (length(idsWithoutParalogs) > 0) {
      
      paralogs_groups <- unique(rbindlist(
        list(
          paralogs_groups,
          data.table(gene_id=idsWithoutParalogs, 
                     paralog_group = (1:length(idsWithoutParalogs))+currentMaxGroup))))
    }
  } else {
    new_paralogs <- NULL
    new_trivial_paralogs <- NULL
  }
  
  paralogs <- unique(rbind(new_paralogs,new_trivial_paralogs,
                           known_paralogs))
  if('paralog_group' %in% names(tT))
    tT[,paralog_group:=NULL]
  tT.paralogs <- merge(tT, paralogs, by.x=geneIdColumn, by.y="ensembl_gene_id", all.x=T, allow.cartesian=T)
  setnames(tT.paralogs, c(geneIdColumn, "paralog_gene_id"), c("paralog_of_gene_id", geneIdColumn))
  
  tT.paralogs <- merge(tT.paralogs, paralogs_groups, by.x="paralog_of_gene_id", by.y="gene_id", all.x=T, allow.cartesian=T)
  tT.paralogs[, paralog_of_gene_id:=NULL]
  
  list(tT.paralogs=tT.paralogs, all_paralogs=all_paralogs, paralogs_groups=paralogs_groups)
}

mapToTargetOrganismEntrez <- function(tT, org, targetOrg, target_orgDb, source_org_mart, target_org_mart) {
  # retrieve the entrez ID for the target organism
  if (org != targetOrg) {
    if(!is.null(getDefaultReactiveDomain()))
      setProgress(detail=paste('Mapping ensembl gene ids of', org, 'to', targetOrg))
    raw_mapping <-getLDS(
      attributes=c("ensembl_gene_id"), 
      filters="ensembl_gene_id", 
      values= tT$gene_id,
      mart=source_org_mart, 
      attributesL=c("ensembl_gene_id"), 
      martL=target_org_mart)
    
    raw_mapping <- as.data.table(raw_mapping)
    
    setnames(raw_mapping, 
             c("Gene.stable.ID", "Gene.stable.ID.1"), 
             c("source_entrez", "target_entrez"))
    
    raw_mapping <- na.omit(raw_mapping)
    
    raw_mapping[,source_entrez := as.character(source_entrez)]
    raw_mapping[,target_entrez := as.character(target_entrez)]
    
    tT <- merge(tT, raw_mapping, by.x="gene_id", by.y="source_entrez", all.x=T, allow.cartesian=T)
    tT <- na.omit(tT, cols="target_entrez")
    tT <- tT[target_entrez != ""]
    tT <- unique(tT)
    setnames(tT, "gene_id", "source_entrez")
    tT <- cSplit(tT, splitCols="target_entrez",sep="///", fixed=T,direction="long", type.convert = F)
  }
  # if the source is already the target organism
  else {
    # just copy over the entrez IDs
    tT[,target_entrez:=as.character(tT$gene_id)]  
    setnames(tT, "gene_id", "source_entrez")
  }
  return(tT)
}

aggregateByParalogGroup <- function(tT) {
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(detail='Aggregating gene expression of paralogue genes')
  print("aggregateByParalogGroup")
  tT_aggr <- tT[, .(target_entrez, source_entrez, paralog_group, logAveExpr, logFC, variances)][, c("target_entrez",
                      "source_entrez",
                      "logAveExpr",
                      # "adj.P.Val",
                      # "P.Value",
                      # "t",
                      # "B",
                      "logFC",
                      "variances"
                      ) := 
                    list(paste(unique(target_entrez), collapse=", "),
                         paste(unique(source_entrez), collapse=", "),
                         mean(logAveExpr),
                         # mean(adj.P.Val),
                         # mean(P.Value),
                         # mean(t),
                         # mean(B),
                         mean(logFC),
                         mean(variances)),
                  by=paralog_group][!duplicated(paralog_group), 
                                    .(paralog_group,
                                      target_entrez, 
                                      source_entrez, 
                                      logAveExpr, 
                                      # adj.P.Val, 
                                      # P.Value, 
                                      # t, 
                                      # B, 
                                      logFC, 
                                      variances)]
  # tT_aggr[,AveExpr:=2^logAveExpr]
  # tT_aggr <- data.table(tT[, lapply(.SD, mean, na.rm=TRUE), by=c("gene_symbol"),
  #                          .SDcol=c("AveExpr","adj.P.Val","P.Value",
  #                                   "t","B","logFC","variances")])
  
  # Calculate the Fold Change from the log value
  tT_aggr[, FC:= 2^logFC]
  # The Fold Change is treated with traingluar function to be able to
  # use it later in the scoring function
  tT_aggr[, FC:= ifelse(FC>1,FC,1/FC)]
  # Normalizing the AveExpr, FC and variances values to calculate scores
  # tT_aggr[, AveExpr := scale(tT_aggr$AveExpr)]
  #tT_aggr$FC <- scale(tT_aggr$FC)
  # tT_aggr[, variances := scale(tT_aggr$variances)]
  
  tT_aggr$target_entrez <- as.character(tT_aggr$target_entrez)
  setkey(tT_aggr, "target_entrez")
  return(tT_aggr)
}

calculateGeneScoresAndRanks <- function(tT_aggr) {
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(detail='Calculating and ranking gene house-keeping scores')
  print("calculateGeneScoresAndRanks")
  # Calculating scores (AveExpr/(FC*variances))
  # add pseudo counts
  #tT_aggr[, score:=AveExpr/((FC+.Machine$double.eps)*(variances+.Machine$double.eps))]
  #tT_aggr[, score:=FC*sqrt(variances)/AveExpr]
  #tT_aggr$score <- scale(tT_aggr$score)
  # Calculating p-values for each gene based on the normalized scores
  tT_aggr[, gene_rank:=rank(as.double(rank(sqrt(variances)))*as.double(rank(FC)))]
  tT_aggr[, pvalue:=1-ecdf(gene_rank)(gene_rank)]
  setorder(tT_aggr, gene_rank)
}

constructRankingMatrices <- function(tT_list) {
  # adding the ranks from each dataset
  ranking_matrix <- data.table(paralog_group=integer(0))
  ranking_matrix_long <- data.table(paralog_group=integer(0), gene_rank=double(0), dataset=character(0))
  m <- 0
  for (dataset in tT_list) {
    if (!is.null(getDefaultReactiveDomain()))
      setProgress(message=paste('Processing data set', names(tT_list)[m+1]), detail='Ranking genes by house-keepingness')
    
    # g <- 1
    # for (dataset in dataset_gpls) {
    #   if (length(dataset_gpls) > 1) {
    #     ds_name <- paste(names(tT_list)[m+1],g,sep=".")
    #   } else {
    #     ds_name <- names(tT_list)[m+1]
    #   }
    ds_name <- names(tT_list)[m+1]
      
    # use data.tables setter to order by gene_rank
    setorder(dataset, gene_rank, na.last = TRUE)
    ds_sub <- dataset[,.(source_entrez,target_entrez, paralog_group, gene_rank)]
    
    ranking_matrix <- merge(ranking_matrix, ds_sub[, .(paralog_group, gene_rank)],
                            all.x=T,
                            all.y=T,
                            by.x=c("paralog_group"), by.y=c("paralog_group"))
    setnames(ranking_matrix, c(tail(names(ranking_matrix), 1)),
             paste0(c(tail(names(ranking_matrix), 1)), ".", ds_name))
    
    ds_sub2 <- data.table::copy(ds_sub)
    ds_sub2[, dataset:=ds_name]
    ranking_matrix_long <- rbindlist(list(ranking_matrix_long, 
                                          ds_sub2[, .(paralog_group, gene_rank, dataset)]))
    #   g <- g + 1
    # }
    # 
    m <- m+1
    if (!is.null(getDefaultReactiveDomain()))
      setProgress(value=0.2+(0.2/length(tT_list)*m))
  }
  
  ranking_matrix[,na_ratio:= (Reduce(`+`, lapply(.SD,function(x) is.na(x)))/
                                (length(grep("gene_rank", names(ranking_matrix))))),
                 .SDcols=grep("gene_rank", names(ranking_matrix))]
  
  
  ranking_matrix[, num_na_factor:= Reduce(`+`, lapply(.SD,function(x) is.na(x))), .SDcols=grep("gene_rank", names(ranking_matrix))]
  # Calculating rank mean and variance of all dataset rankings
  ranking_matrix[,rank_mean:=apply(.SD,1,mean,na.rm = TRUE), .SDcols=grep("gene_rank",names(ranking_matrix))]
  #ranking_matrix[,rank_mean:=rank_mean*(num_na_factor+1)]
  ranking_matrix[,rank_var:=apply(.SD,1,sd,na.rm = TRUE), .SDcols=grep("gene_rank",names(ranking_matrix))]
  
  # Calculating total rank among all datasets (highest mean, lowest var)
  ranking_matrix[, total_rank:=order(order(rank_mean, rank_var))]
  ranking_matrix[, total_rank:=total_rank*(num_na_factor+1)]
  #setorder(ranking_matrix, "rank_mean", "rank_var", na.last=T)
  setorder(ranking_matrix, total_rank, na.last=T)
  
  
  return(list(ranking_matrix=ranking_matrix, ranking_matrix_long=ranking_matrix_long))
}

calculateBootstrapMatrix <- function(ranking_matrix, bootstrap_replications, bootstrap_sample_size, num_hkg) {
  ### Apply random sampling with replacement (bootstrapping), to reduce the
  # variance and bias produced by arrangement or single datasets having higher
  # influence on ranking
  ranking_matrix_raw <- ranking_matrix[,1:(length(ranking_matrix)-1)]
  bootstrap_matrix <- ranking_matrix[,.(paralog_group, na_ratio, num_na_factor)]
  bootstrap_matrices <- list(data.table(paralog_group=integer(0), gene_rank=double(0), bootstrap_sample=integer(0)))
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(message='Calculating robustness of house-keeping gene ranks', detail='Calculating bootstrap samples and calculating robustness of house-keeping gene ranks')
  #TODO: I don't like that the bootstrap samples are added and deleted from the bootstrap_matrix
  for (i in 1:bootstrap_replications) {
    if(!is.null(getDefaultReactiveDomain()))
      setProgress(detail=paste('Bootstrap sample ',paste(i,'/', bootstrap_replications, sep="")))
    sample_cols <- sample(grep("gene_rank",
                               names(ranking_matrix_raw)), 
                          bootstrap_sample_size, 
                          replace = TRUE)
    sample_ranking <- subset(ranking_matrix_raw, select=sample_cols)
    sample_ranking[,c("rank_mean", "rank_var"):=
                     list(rowMeans(.SD,na.rm = TRUE), 
                          rowSds(.SD, na.rm=TRUE)),
                   .SDcols=grep("gene_rank",names(sample_ranking))]
    sample_ranking[,total_rank:=order(order(rank_mean, rank_var, na.last=TRUE))]
    bootstrap_matrix[,sample:=sample_ranking$total_rank]
    setnames(bootstrap_matrix, "sample", paste0("sample_", i))
    
    bootstrap_matrices[[length(bootstrap_matrices)+1]] <- data.table(
      paralog_group=bootstrap_matrix$paralog_group,
      gene_rank=sample_ranking$total_rank,
      bootstrap_sample=i)
    
    if(!is.null(getDefaultReactiveDomain()))
      setProgress(value = 0.4+(0.3/bootstrap_replications*i))
    else
      print(sprintf("%.0f%%", i*100 / bootstrap_replications))
  }
  
  if(!is.null(getDefaultReactiveDomain()))
    setProgress(message='Gathering bootstrap samples', detail = '', value = 0.7)
  
  # Calculate the mean and SD of the ranking at each time
  # and then get the global ranking based on lowest mean rank and lowest SD
  bootstrap_matrix[,rank_mean:=apply(.SD,1,mean,na.rm = TRUE)*(1+na_ratio),
                   .SDcols=grep("sample_",names(bootstrap_matrix))]
  bootstrap_matrix[,rank_var:=apply(.SD,1,sd,na.rm = TRUE),
                   .SDcols=grep("sample_",names(bootstrap_matrix))]
  bootstrap_matrix[,total_rank:=order(order(rank_mean, rank_var))]
  setorder(bootstrap_matrix, total_rank, na.last = T)
  bootstrap_matrix[,grep("sample_",names(bootstrap_matrix)):=NULL]
  
  bootstrap_matrix_long <- rbindlist(
    lapply(bootstrap_matrices, function(x) x[paralog_group %in% bootstrap_matrix[1:num_hkg, paralog_group]]))
  rm(bootstrap_matrices)
  gc()
  
  list(bootstrap_matrix=bootstrap_matrix[1:num_hkg], bootstrap_matrix_long=bootstrap_matrix_long)
}

annotateTargetParalogGroupsWithSymbols <- function(m, paralogue_groups, geneIdColumn, target_orgDb, 
                                                   paralogue_groups_to_be_annotated=paralogue_groups[,get("paralog_group")]) {
  # get symbols for paralog groups
  mapped <- na.omit(ensembldb::select(target_orgDb, keys=paralogue_groups[paralog_group %in% paralogue_groups_to_be_annotated,
                                                               get(geneIdColumn)], keytype="GENEID", column="SYMBOL"))
  mapped <- data.table(mapped %>%
                         dplyr::group_by(GENEID) %>%
                         dplyr::summarise(symbols = paste(unique(SYMBOL), collapse=", ")), key="GENEID")
  paralogue_groups[match(mapped$GENEID, get(geneIdColumn)),symbols:=mapped$symbols]
  
  paralogs_groups_aggr <- na.omit(paralogue_groups)[paralog_group %in% paralogue_groups_to_be_annotated, lapply(.SD, function(x) {
    paste(unique(x), collapse=", ")
  }), by=paralog_group]
  setkey(paralogs_groups_aggr, paralog_group)
  m[paralog_group %in% paralogue_groups_to_be_annotated, 
           c("symbols", "gene_ids"):=list(
             paralogs_groups_aggr[J(.BY$paralog_group), symbols],
             paralogs_groups_aggr[J(.BY$paralog_group), gene_id]
             ), 
           by=paralog_group]
}


getGeneSymbolsWithoutEnsembleId <- function(dt) {
  dt %>%
    dplyr::group_by(GENE_SYMBOL) %>%
    dplyr::summarise(n=n_distinct(ENSEMBLGENEID, na.rm=TRUE)) %>%
    dplyr::filter(n==0) %>%
    dplyr::pull(GENE_SYMBOL)
}



getUniqueSessionId <- function(pool) {
  sessions <- dbGetQuery(pool, 
                         "SELECT session_uuid FROM sessions;")$session_uuid
  
  new_session_id <- sprintf("%.0f", floor(runif(1)*1e20))
  while (new_session_id %in% sessions) {
    new_session_id <- sprintf("%.0f", floor(runif(1)*1e20))
  }
  new_session_id
}

getSessionStatus <- function(pool, session_id) {
  if (is.null(session_id))
    return("invalid-id")
  res <- dbGetQuery(pool,
                    "SELECT finish_ts FROM sessions WHERE session_uuid = $1;",
                    params=list(session_id))
  if (nrow(res)==0)
    "unknown-id"
  else if (is.na(res$finish_ts))
    "started"
  else
    "finished"
}

getFinishedSessions <- function(pool) {
  res <- dbGetQuery(pool,
                    "SELECT session_uuid FROM sessions WHERE finish_ts is not null;")
  res
}

getSessionStarttime <- function(pool, session_id) {
  if (is.null(session_id))
    return(NA)
  res <- dbGetQuery(pool,
                    "SELECT create_ts FROM sessions WHERE session_uuid = $1;",
                    params=list(session_id))
  if (nrow(res)==0)
    NA
  else
    as.POSIXct(res$create_ts, origin="1970-01-01")
}

getSessionFinishtime <- function(pool, session_id) {
  if (is.null(session_id))
    return(NA)
  res <- dbGetQuery(pool,
                    "SELECT finish_ts FROM sessions WHERE session_uuid = $1;",
                    params=list(session_id))
  if (nrow(res)==0)
    NA
  else
    as.POSIXct(res$finish_ts, origin="1970-01-01")
}

#### NormFinder ####
# Version history: see bottom of file
# 
# 
# Normfinder=function(dat0,Groups=TRUE,ctVal=TRUE,pStabLim=0.25){
#   #
#   # If Groups is TRUE the last row contains the group identifier,
#   # and the last row must be started by a name for that row.
#   # No spaces are allowed in the gene names, sample names and group identifier.
#   #
#   #dat0=read.table(filename,header=TRUE,row.names=1,colClasses="character")
#   #
#   ntotal=dim(dat0)[2] # number of samples
#   k0=dim(dat0)[1] # number of rows
#   #
#   if (Groups){
#     ngenes=k0-1 # number of genes
#     genenames=rownames(dat0)[-k0]
#     grId=dat0[k0,]
#     dat0=dat0[-k0,]
#   } else {
#     ngenes=k0 # number of genes
#     genenames=rownames(dat0)
#     grId=rep(1,ntotal)
#   }
#   #
#   dat=matrix(as.numeric(unlist(dat0)),ngenes,ntotal) # matrix instead of list
#   #
#   if (!ctVal){dat=log2(dat)} # transform to log2 values
#   #
#   samplenames=colnames(dat0)
#   grId=factor(unlist(grId))  # group identifier
#   groupnames=levels(grId)  # group names
#   ngr=length(levels(grId)) # number of groups
#   # Number of samples in each group:
#   nsamples=rep(0,ngr)
#   for (group in 1:ngr){nsamples[group]=sum(grId==groupnames[group])}
#   #
#   #
#   MakeStab=function(da){
#     ngenes=dim(da)[1]
#     # Sample averages
#     sampleavg=apply(da,2,mean)
#     # Gene averages within group
#     genegroupavg=matrix(0,ngenes,ngr)
#     for (group in 1:ngr){
#       genegroupavg[,group]=apply(da[,grId==groupnames[group]],1,mean)}
#     # Group averages
#     groupavg=rep(0,ngr)
#     for (group in 1:ngr){groupavg[group]=mean(da[,grId==groupnames[group]])}
#     
#     # Variances 
#     GGvar=matrix(0,ngenes,ngr)
#     for (group in 1:ngr){
#       grset=(grId==groupnames[group])
#       a=rep(0,ngenes)
#       for (gene in 1:ngenes){
#         a[gene]=sum((da[gene,grset]-genegroupavg[gene,group]-
#                        sampleavg[grset]+groupavg[group])^2)/(nsamples[group]-1)
#       }
#       GGvar[,group]=(a-sum(a)/(ngenes*ngenes-ngenes))/(1-2/ngenes)
#     }
#     #
#     # Change possible negative values
#     genegroupMinvar=matrix(0,ngenes,ngr)
#     for (group in 1:ngr){
#       grset=(grId==groupnames[group])
#       z=da[,grset]
#       for (gene in 1:ngenes){
#         varpair=rep(0,ngenes)
#         for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
#         genegroupMinvar[gene,group]=min(varpair[-gene])/4
#       }
#     }
#     #
#     # Final variances
#     GGvar=ifelse(GGvar<0,genegroupMinvar,GGvar)
#     #
#     # Old stability measure for each gene is calculated:
#     #
#     dif=genegroupavg
#     difgeneavg=apply(dif,1,mean)
#     difgroupavg=apply(dif,2,mean)
#     difavg=mean(dif)
#     for (gene in 1:ngenes){
#       for (group in 1:ngr){
#         dif[gene,group]=dif[gene,group]-difgeneavg[gene]-difgroupavg[group]+difavg
#       }
#     }
#     #
#     nsampMatrix=matrix(rep(nsamples,ngenes),ngenes,ngr,byrow=T)
#     vardif=GGvar/nsampMatrix
#     gamma=sum(dif*dif)/((ngr-1)*(ngenes-1))-sum(vardif)/(ngenes*ngr)
#     gamma=ifelse(gamma<0,0,gamma)
#     #
#     difnew=dif*gamma/(gamma+vardif)
#     varnew=vardif+gamma*vardif/(gamma+vardif)
#     Ostab0=abs(difnew)+sqrt(varnew)
#     Ostab=apply(Ostab0,1,mean)
#     #
#     # Measure of group differences:
#     mud=rep(0,ngenes)
#     for (gene in 1:ngenes){
#       mud[gene]=2*max(abs(dif[gene,]))
#     }
#     # Common variance:
#     genevar=rep(0,ngenes)
#     for (gene in 1:ngenes){
#       genevar[gene]=sum((nsamples-1)*GGvar[gene,])/(sum(nsamples)-ngr)
#     }
#     Gsd=sqrt(genevar)
#     #
#     # Return results:
#     #
#     return(cbind(mud,Gsd,Ostab,rep(gamma,ngenes),GGvar,dif))
#   }    # End of function MakeStab
#   #
#   #
#   MakeComb2=function(g1,g2,res){
#     gam=res[1,4]
#     d1=res[g1,(4+ngr+1):(4+ngr+ngr)]; d2=res[g2,(4+ngr+1):(4+ngr+ngr)]
#     s1=res[g1,(4+1):(4+ngr)]; s2=res[g2,(4+1):(4+ngr)]
#     rho=abs(gam*d1/(gam+s1/nsamples)+gam*d2/(gam+s2/nsamples))*
#       sqrt(ngenes/(ngenes-2))/2
#     rho=rho+sqrt(s1/nsamples+gam*s1/(nsamples*gam+s1)+
#                    s2/nsamples+gam*s2/(nsamples*gam+s2))/2
#     return(mean(rho))
#   }
#   #
#   #
#   MakeStabOne=function(da){
#     ngenes=dim(da)[1]
#     # Sample averages
#     sampleavg=apply(da,2,mean)
#     # Gene averages
#     geneavg=apply(da,1,mean)
#     totalavg=mean(da)
#     #
#     # Variances 
#     genevar0=rep(0,ngenes)
#     for (gene in 1:ngenes){
#       genevar0[gene]=
#         sum((dat[gene,]-geneavg[gene]-sampleavg+totalavg)^2)/
#         ((ntotal-1)*(1-2/ngenes))
#     }
#     genevar=genevar0-sum(genevar0)/(ngenes*ngenes-ngenes)
#     #
#     # Change possible negative values
#     geneMinvar=rep(0,ngenes)
#     z=da
#     for (gene in 1:ngenes){
#       varpair=rep(0,ngenes)
#       for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
#       geneMinvar[gene]=min(varpair[-gene])/4
#     }
#     # Final variances
#     genevar=ifelse(genevar<0,geneMinvar,genevar)
#     #
#     return(genevar)
#   }
#   #     End of function MakeStabOne
#   #
#   #################################################
#   #
#   # Main part
#   #
#   if (ngr>1){   # More than one group.
#     #
#     res=MakeStab(dat)
#     #
#     gcand=c(1:ngenes)[res[,3]<pStabLim]
#     ncand=length(gcand)
#     if (ncand<4){
#       if (ngenes>3){
#         li=sort(res[,3])[4]
#         gcand=c(1:ngenes)[res[,3]<=li]
#         ncand=length(gcand)
#       } else {
#         gcand=c(1:ngenes)
#         ncand=length(gcand)
#       }
#     }
#     #
#     vv2=c()
#     #
#     for (g1 in 1:(ncand-1)){
#       for (g2 in (g1+1):ncand){
#         qmeas=MakeComb2(gcand[g1],gcand[g2],res)
#         vv2=rbind(vv2,c(gcand[g1],gcand[g2],qmeas))
#       }}
#     #
#     ord=order(res[,3])
#     FinalRes=list(Ordered=
#                     data.frame("GroupDif"=round(res[ord,1],2),"GroupSD"=round(res[ord,2],2),
#                                "Stability"=round(res[ord,3],2),row.names=genenames[ord]),
#                   UnOrdered=
#                     data.frame("GroupDif"=round(res[,1],2),"GroupSD"=round(res[,2],2),
#                                "Stability"=round(res[,3],2),
#                                "IGroupSD"=round(sqrt(res[,(4+1):(4+ngr)]),2),
#                                "IGroupDif"=round(res[,(4+ngr+1):(4+ngr+ngr)],2),
#                                row.names=genenames),
#                   PairOfGenes=
#                     data.frame("Gene1"=genenames[vv2[,1]],"Gene2"=genenames[vv2[,2]],
#                                "Stability"=round(vv2[,3],2)))
#     #
#     return(FinalRes)
#     #
#   } else {    # End of more than one group: next is for one group only.
#     #
#     #
#     sigma=sqrt(MakeStabOne(dat))
#     #
#     siglim=(min(sigma)+0.1)
#     gcand=c(1:ngenes)[sigma<siglim]
#     ncand=length(gcand)
#     #
#     if ((ncand>=2)&(ngenes>3)){
#       #
#       vv2=c()
#       #
#       for (g1 in 1:(ncand-1)){
#         for (g2 in (g1+1):ncand){
#           dat1=rbind(dat[-c(gcand[g1],gcand[g2]),],
#                      apply(dat[c(gcand[g1],gcand[g2]),],2,mean))
#           qmeas=sqrt(MakeStabOne(dat1))
#           vv2=rbind(vv2,c(gcand[g1],gcand[g2],qmeas[ngenes-1]))
#         }}
#       ord=order(sigma)
#       FinalRes=list(Ordered=
#                       data.frame("GroupSD"=round(sigma[ord],2),row.names=genenames[ord]),
#                     PairOfGenes=
#                       data.frame("Gene1"=genenames[vv2[,1]],"Gene2"=genenames[vv2[,2]],
#                                  "GroupSD"=round(vv2[,3],2)))
#     } else { # No combined genes to consider
#       ord=order(sigma)
#       FinalRes=list(Ordered=
#                       data.frame("GroupSD"=round(sigma[ord],2),row.names=genenames[ord]))
#     } # End ncand<2 or ngenes<=3
#     #
#     return(FinalRes)
#     #
#   }  # End one group only
#   #
# } # End of main function
# 
# # Version history:
# # Version dated 17/12-2013.
# # Version dated 18/06-2014: Correction in line 121: s2 is squared.
# # Version dated 23/08-2014: Correction to above correction in line 121: 
# #   s1 and s2 are not squared. 
# #   We thank John R. Young for pointing to the error in line 121.
# # Version dated 05/01-2015: sum()/2 changed to mean() in line 126.
# #   This misprint appears in the Supplementary Information as well.
# #   We thank Mark Birtwistle for pointing out this error.
