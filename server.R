PIPELINE_STATUS <- list(IDLE = 1, STARTED = 2, DOWNLOAD_DATASETS = 3, DOWNLOAD_ANNOTATIONS = 4, FIND_HOUSEKEEPING_GENES = 5, FINISHED = 6)

server <- function(input, output, session) {
  enableBookmarking("url")
  source("housekeepr-init.R", local=T, echo=T)
  source("housekeepr-functions.R", local=T, echo=T)
  #source("housekeepr-mesh.R", local = T, echo = T)
  load("init_data.RData") #contains allOrganisms, allTisues and allConditions
  waiter_hide()

  disablePrettyCheckbox <- function(p) {
    p$children[[1]]$children[[1]]$attribs <- append(p$children[[1]]$children[[1]]$attribs, list(disabled="disabled"))
    p
  }
  
  make_color <- function(name, color = 'SteelBlue'){
    return(shiny::HTML(sprintf("<font style='color: %s'>%s</font>", color, name)))
  }
  
  extract_hkg <- function(
    num_hkg = 10,
    GSE_list = list(),
    GSE_samples_annotations = list(),
    bootstrap_sample_size = round(length(GSE_list)/2),
    bootstrap_replications = 50,
    hkgDataDir = "hkg-data") {
    
    print("extract_hkg")
    setProgress(message='Analysis started',
                detail='We are now processing your analysis. This may take a while ...')
    
    # STATIC VARIABLES
    # top tables of all data sets
    tT_list <- list()
    
    # keep track of organisms of each data set
    dataset_organisms <- list()
    
    # mapping of gene id -> paralog gene id
    all_paralogs <- data.table(gene_id=character(0), paralog_gene_id=character(0))
    
    # mapping of gene id -> "paralog group"
    paralogs_groups <- data.table(gene_id=character(0), paralog_group=numeric(0), key="gene_id")
    
    # USER INPUT VALIDATION
    
    if(length(GSE_list) != length(GSE_samples_annotations)){
      return("Please provide annotations as many as datasets")
    }
    
    setProgress(value = 0.0)
    for (i in 1:length(GSE_list)) {
      tryCatch({
        gse_gpl <- GSE_list[[i]]
        gse <- strsplit(gse_gpl, ".", fixed=T)[[1]][1]
        gpl <- strsplit(gse_gpl, ".", fixed=T)[[1]][2]
        setProgress(message=paste('Processing data set', gse, '(', gpl, ')'), detail='Loading ...')
        
        # load series and platform data from GEO
        gset_raw <- reactive_vals$gsets[[gse]]
        
        gsetEx <- getExpressions(gset_raw, gpl=gpl)[[1]]
        
        if (is.null(gsetEx) | is.null(gsetEx$ex)){
          showNotification(paste0("Could not load expression values for ", gse, ". Maybe the expression data is not available in the series matrix?"), type = "error", duration = NULL)
          next()
        }
        gset <- gsetEx$gset
        
        currentOrganism <- tolower(identifyOrganism(gset))
        if (is.null(currentOrganism)){
          showNotification(paste0("Could not find organism for ", gse, "."), type = "error", duration = NULL)
          next()
        }
        dataset_organisms[gse] <- currentOrganism
        
        # make proper column names to match toptable
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # if first dataset cannot be processed these variables were not initialized with i == 1
        if (is.null(reactive_vals$template_organism) | i == 1) {
          reactive_vals$template_organism <- tolower(dataset_organisms[[1]])
          targetOrg <- isolate(reactive_vals$template_organism)
          target_orgDb <- isolate(reactive_vals$orgAnno[[targetOrg]])
          target_org_mart <- isolate(reactive_vals$orgMart[[targetOrg]])
        }
        
        if (currentOrganism != targetOrg) {
          source_orgDb <- isolate(reactive_vals$orgAnno[[currentOrganism]])
          source_org_mart <- isolate(reactive_vals$orgMart[[currentOrganism]])
        } else {
          source_orgDb <- target_orgDb
          source_org_mart <- target_org_mart
        }
        
        tT <- calculateTopTable(GSE_samples_annotations, gse_gpl, gset)
        tT <- ensureGeneIdColumn(tT, currentOrganism, ah, isolate(reactive_vals$ensembl_release), orgDb=source_orgDb)
        if (is.null(tT)){
          showNotification(paste("It seems like no Ensembl or Entrez gene id could be found. Skipping", gse), type = "error", duration = NULL)
          next()
        }
        tT <- removeWithMissingEntrezId(tT, "gene_id")
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value", "logAveExpr","AveExpr", "variances","t","B","logFC", "gene_id"))
        
        tT <- unique(tT)
        
        
        for (i in 1:5){
          res <- try(extractParalogs(tT, source_org_mart, all_paralogs, paralogs_groups, "gene_id"))
          if (!inherits(res, 'try-error')) break
          if (grepl("biomaRt has encountered an unexpected server error|curl::curl_fetch_memory", res)) {
            showNotification(paste0("Attempt ",i,": Querying Ensembl failed. Retrying in 2 secs... ", type = "error", duration = 20))
            Sys.sleep(2)
          } else { 
            return("Querying Ensembl failed 5 times. You may try again later.")
          }
        }
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
        
        # save all datasets to same list
        tT_list[[gse_gpl]] <- tT_aggr
        
        # end GPL
        setProgress(value = (0.2/length(GSE_list)*i))
      }, error=function(e) {
        
        print(paste("Parsing of data set", gse, '(', gpl, ')', "failed:", e))
        showNotification(paste("Parsing of data set", gse, '(', gpl, ')', "failed"), type = "error", duration = NULL)
      })
    }
    if (length(tT_list) == 0) {
      
      #showNotification("None of the selected data sets / samples could be processed.", type = "error", duration = NULL)
      return('None of the selected data sets / samples could be processed')
    } else if (length(tT_list) == 1) {
      
      #showNotification("Only one of the selected data sets / samples could be processed.", type = "error", duration = NULL)
      return('Only one of the selected data sets / samples could be processed')
    }
    
    setProgress(value = 0.2)
    if (length(tT_list) < bootstrap_sample_size) {
      #output$bootstrapSampleSize <- reactive({length(tT_list)})
      updateNumericInput(session, 'bootstrapSampleSize', value=max(length(tT_list),1))
      bootstrap_sample_size <- max(length(tT_list),1)
    }
    
    print("### Calculating the ranking of genes in each individual dataset")
    res <- constructRankingMatrices(tT_list)
    ranking_matrix <- res$ranking_matrix
    ranking_matrix_long <- res$ranking_matrix_long
    
    print("### Running bootstrapping with the following parameters:")
    print(paste("- Number of HKG genes:", num_hkg))
    print(paste("- Bootstrapping sample size:", bootstrap_sample_size))
    print(paste("- Bootstrapping replications:", bootstrap_replications))
    
    setProgress(value = 0.4)
    
    bootstrapRes <- calculateBootstrapMatrix(ranking_matrix, bootstrap_replications, bootstrap_sample_size, num_hkg)
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
    
    #check that number of hkg given as input is not above the total number of genes
    num_hkg <- ifelse(num_hkg>nrow(candidate_hkg_ranking),
                      nrow(candidate_hkg_ranking), num_hkg)
    
    # data preparations for plotting
    top_hkg <- candidate_hkg_ranking$paralog_group[1:num_hkg]
    
    res <- list(bootstrap_ranking=bootstrap_matrix[1:num_hkg], 
                bootstrap_ranking_long=bootstrap_matrix_long[paralog_group%in%top_hkg], 
                ranking_long=ranking_matrix_long[paralog_group%in%top_hkg])
    setProgress(value = 1.0)
    res
  }
  
  ### DATABASE STUFF ----
  # development
  if (houseekpr_env %in% "development") {
    print("Starting HouseKeepR in development mode")
    # in development we use sqlite
    library(RSQLite)
    pool <- dbPool(
      drv = RSQLite::SQLite(),
      dbname = "hkg.sqlite"
    )
  } else if (houseekpr_env == "production") {
    print("Starting HouseKeepR in production mode")
    # in production we use postgresql
    library(RPostgreSQL)
    pool <- dbPool(
      drv = dbDriver("PostgreSQL"),
      dbname = "housekeepr",
      host="db",
      #host="localhost",
      port=5432,
      user="housekeepr",
      password="vz5ENsAd"
    )
  } else
    stop("Invalid mode specified in HOUSEKEEPR_ENV environmental variable. Valid values are \"development\" and \"production\"")
  
  
  #DBI::dbExecute(pool, "DROP TABLE sessions;")
  if (!dbExistsTable(pool, "sessions")) {
    sessionsTable <- data.frame(session_uuid=character(0), 
                                create_ts=as.POSIXct(character()),
                                finish_ts=as.POSIXct(character()))
    dbWriteTable(pool, "sessions", sessionsTable)
  }
  
  if (!dbExistsTable(pool, "datasets")) {
    datasetTable <- data.frame(dataset_id=character(0), status=character(0))
    dbWriteTable(pool, "datasets", datasetTable)
  }
  
  if (!dbExistsTable(pool, "dataset_locks")) {
    dsLocksTable <- data.frame(
      dataset_id=character(0), 
      session_uuid=character(0))
    dbWriteTable(pool, "dataset_locks", dsLocksTable)
  }
  
  #DBI::dbExecute(pool, "DROP TABLE gene_bootstrap_total_rank_scores;")
  if (!dbExistsTable(pool, "gene_bootstrap_total_rank_scores")) {
    geneBootstrapTotalRankScoreTable <- data.frame(
      session_uuid=character(0), 
      paralog_group=numeric(0), 
      symbols=character(0), 
      gene_ids=character(0),
      missing_datasets=integer(0),
      rank_mean=double(0), 
      rank_var=double(0), 
      rank=double(0))
    dbWriteTable(pool, "gene_bootstrap_total_rank_scores", geneBootstrapTotalRankScoreTable)
    conn <- poolCheckout(pool);
    dbSendQuery(conn, "CREATE INDEX idx1 ON gene_bootstrap_total_rank_scores (session_uuid);")
    poolReturn(conn)
  }
  
  #DBI::dbExecute(pool, "DROP TABLE gene_bootstrap_sample_rank_scores;")
  if (!dbExistsTable(pool, "gene_bootstrap_sample_rank_scores")) {
    geneBootstrapSampleRankScoreTable <- data.frame(
      session_uuid=character(0), 
      paralog_group=numeric(0), 
      symbols=character(0),
      gene_ids=character(0),
      bootstrap_sample=integer(0), 
      rank=double(0))
    dbWriteTable(pool, "gene_bootstrap_sample_rank_scores", geneBootstrapSampleRankScoreTable)
    conn <- poolCheckout(pool);
    dbSendQuery(conn, "CREATE INDEX idx2 ON gene_bootstrap_sample_rank_scores (session_uuid);")
    poolReturn(conn)
  }
  
  #DBI::dbExecute(pool, "DROP TABLE gene_dataset_ranks;")
  if (!dbExistsTable(pool, "gene_dataset_ranks")) {
    geneDatasetRankTable <- data.frame(
      session_uuid=character(0),
      paralog_group=numeric(0), 
      symbols=character(0),
      gene_ids=character(0),
      dataset=character(0), 
      rank=double(0))
    dbWriteTable(pool, "gene_dataset_ranks", geneDatasetRankTable)
    conn <- poolCheckout(pool);
    dbSendQuery(conn, "CREATE INDEX idx3 ON gene_dataset_ranks (session_uuid);")
    poolReturn(conn)
  }
  
  #DBI::dbExecute(pool, "DROP TABLE analysis_parameters;")
  if (!dbExistsTable(pool, "analysis_parameters")) {
    analysisParameterTable <- data.frame(session_uuid=character(0),
                                         tissue_type=character(0),
                                         condition=character(0),
                                         organism=character(0),
                                         num_datasets=integer(0),
                                         number_hkg=integer(0),
                                         bootstrap_sample_size=integer(0),
                                         bootstrap_replications=integer(0),
                                         ensembl_release=character(0),
                                         seed=integer(0))
    dbWriteTable(pool, "analysis_parameters", analysisParameterTable)
    conn <- poolCheckout(pool);
    dbSendQuery(conn, "CREATE INDEX idx4 ON analysis_parameters (session_uuid);")
    poolReturn(conn)
  }
  
  onStop(function() {
    poolClose(pool)
  })
  
  
  reactive_vals <- reactiveValues(session_id = NULL,
                                  invalidate_datasets_table=0, 
                                  analysis_started=FALSE,
                                  analysis_running = FALSE, 
                                  status=PIPELINE_STATUS$IDLE, 
                                  organisms=c(),
                                  conditions=c(),
                                  tissues=c(),
                                  tissue_selected=FALSE,
                                  condition_selected=FALSE,
                                  organism_selected=FALSE,
                                  rentrez_search_result=NULL,
                                  analysis_error=NULL,
                                  gsets=list(),
                                  queryingGEOdatasets=FALSE, available_datasets=NULL, selected_datasets=NULL, num_selected_datasets=0,
                                  queryingGEOsamples=FALSE, available_samples=NULL, 
                                  num_available_samples=c(),
                                  selected_samples=NULL,  
                                  num_selected_samples=c(),
                                  valid_sample_selection=FALSE,
                                  # when loading by bookmarked URL we store values until renderUI retrieved them 
                                  restoreState=NULL,
                                  restore_dataset_selection=NULL,
                                  restore_bootstrap_sample_size=NULL,
                                  restore_conditions_all=NULL,
                                  restore_tissues_all=NULL,
                                  restore_organisms_all=NULL,
                                  restore_conditions_selected=NULL,
                                  restore_tissues_selected=NULL,
                                  restore_organisms_selected=NULL,
                                  # use_normfinder=FALSE,
                                  numberGenes=0,
                                  bootstrap_sample_size=0,
                                  bootstrap_replications=0,
                                  bootstrappingRankTableData=data.table(
                                    paralog_group=integer(0), 
                                    symbols=character(0),
                                    gene_ids=character(0),
                                    missing_datasets=integer(0),
                                    rank_mean=double(0), 
                                    rank_var=double(0), 
                                    rank=double(0)),
                                  bootstrappingLongRankTableData=data.table(),
                                  rankLongTableData=data.table(
                                    dataset=character(0), 
                                    paralog_group=integer(0), 
                                    symbols=character(0), 
                                    gene_ids=character(0),
                                    rank=double(0)),
                                  # trigger
                                  delayCounterRankingTableUpdate=T,
                                  bookmark_exclude=c("use_old_session_button",
                                                     "old_session",
                                                     "waiter_hidden",
                                                     "select_all_datasets",
                                                     "select_datasets_table_columns_selected",
                                                     "select_datasets_table_cells_selected",
                                                     "select_datasets_table_cell_clicked",
                                                     "select_datasets_table_rows_all",
                                                     "select_datasets_table_rows_current",
                                                     "select_datasets_table_rows_selected",
                                                     "select_datasets_table_state",
                                                     "select_datasets_table_search",
                                                     "houseKeepingGenesTable_cell_clicked",
                                                     "houseKeepingGenesTable_rows_all",
                                                     "houseKeepingGenesTable_rows_current",
                                                     "houseKeepingGenesTable_search",
                                                     "houseKeepingGenesTable_state",
                                                     "select_none_datasets",
                                                     "select_datasets_table_row_last_clicked",
                                                     "sidebarCollapsed",
                                                     "sidebarItemExpanded",
                                                     "start_analysis",
                                                     "pretty_1",
                                                     "pretty_2",
                                                     "pretty_3",
                                                     "pretty_4",
                                                     "pretty_5",
                                                     "pretty_6",
                                                     "customEnsemblRelease",
                                                     "seed",
                                                     "menu",
                                                     ".clientValue-default-plotlyCrosstalkOpts",
                                                     ".clientValue-plotly_relayout-A",
                                                     "houseKeepingGenesTable_row_last_clicked"),
                                  ensembl_release=NULL,
                                  available_ensembl_releases=NULL,
                                  rankPlotlyHideOutliers=T,
                                  plotsNumHKG=15,
                                  seed=NULL
  )
  
  
  # add custom values to the bookmarking state
  onBookmark(function(state) {
    
    # we are in the results view
    if (getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")) {
      state$values$session <- reactive_vals$session_id
    }
    else {
      if (!is.null(reactive_vals$selected_datasets$accession))
        state$values$ds <- paste(reactive_vals$selected_datasets$accession,reactive_vals$selected_datasets$gpl,sep=".")
      if (!is.null(reactive_vals$selected_samples))
        state$values$samples <- reactive_vals$selected_samples
      if (!is.null(reactive_vals$ensembl_release))
        state$values$ensembl <- reactive_vals$ensembl_release
      if (!is.null(reactive_vals$seed))
        state$values$seed <- reactive_vals$seed
    }
  })
  
  onBookmarked(function(url) {
    updateQueryString(url)
  })
  
  observeEvent(input$menu, 
               # allow UI to update first
               priority=-100, {
                 
                 # render after menu
                 if (is.null(input$menu))
                   return()
                 
                 
                 if (!is.null(reactive_vals$restoreState)) {
                   print("doRestoreState")
                   state <- isolate(reactive_vals$restoreState)
                   
                   if (!is.null(state$input$select_tissues_all)) {
                     userEntries <- state$input$select_tissues_all[!(state$input$select_tissues_all %in% allTissues)]
                     
                     reactive_vals$restore_tissues_all <- c(unlist(allTissues), userEntries)
                     reactive_vals$restore_tissues_selected <- state$input$select_tissues_all
                   }
                   if (!is.null(state$input$select_conditions_all)) {
                     userEntries <- state$input$select_conditions_all[!(state$input$select_conditions_all %in% allConditions)]
                     
                     reactive_vals$restore_conditions_all <- c(unlist(allConditions), userEntries)
                     reactive_vals$restore_conditions_selected <- state$input$select_conditions_all
                   }
                   if (!is.null(state$input$select_organisms_all)) {
                     userEntries <- state$input$select_organisms_all[!(state$input$select_organisms_all %in% allOrganisms)]
                     
                     reactive_vals$restore_organisms_all <- c(unlist(allOrganisms), userEntries)
                     reactive_vals$restore_organisms_selected <- state$input$select_organisms_all
                   }
                   if (!is.null(state$values$ds)) {
                     reactive_vals$restore_dataset_selection <- state$values$ds
                     
                     if (!is.null(state$values$samples)) {
                       reactive_vals$selected_samples <- state$values$samples[names(state$values$samples) %in% reactive_vals$restore_dataset_selection]
                     }
                   }
                   
                   if (!is.null(state$input$numberGenes)){
                     updateNumericInput(session, 'numberGenes', value=state$input$numberGenes)
                     reactive_vals$restore_numberGenes <- state$input$numberGenes
                   }
                   if (!is.null(state$input$bootstrapRepetitions)){
                     updateNumericInput(session, 'bootstrapRepetitions', value=state$input$bootstrapRepetitions)
                     reactive_vals$restore_bootstrapRepetitions <- state$input$bootstrapRepetitions
                   }
                   if (!is.null(state$input$bootstrapSampleSize))
                     reactive_vals$restore_bootstrap_sample_size <- state$input$bootstrapSampleSize
                   if (!is.null(state$values$ensembl)) {
                     updateSelectizeInput(session, "customEnsemblRelease",
                                          selected=state$values$ensembl)
                     reactive_vals$ensembl_release <- state$values$ensembl
                   }
                   if (!is.null(state$values$seed)) {
                     updateSelectizeInput(session, "seed", 
                                          selected=state$values$seed)
                     reactive_vals$restore_seed <- state$values$seed
                   }
                   
                   
                   if (getSessionStatus(pool, state$values$session) == "finished") {
                     reactive_vals$bootstrappingRankTableData <- as.data.table(dbGetQuery(pool, 
                                                                                          "SELECT paralog_group,symbols,gene_ids,missing_datasets, rank_mean,rank_var,rank FROM gene_bootstrap_total_rank_scores WHERE session_uuid = $1;", 
                                                                                          param=list(reactive_vals$session_id)))
                     reactive_vals$bootstrappingLongRankTableData <- as.data.table(dbGetQuery(pool, 
                                                                                              "SELECT paralog_group,symbols,gene_ids,bootstrap_sample,rank FROM gene_bootstrap_sample_rank_scores WHERE session_uuid = $1;", 
                                                                                              param=list(reactive_vals$session_id)))
                     reactive_vals$rankLongTableData <- as.data.table(dbGetQuery(pool, 
                                                                                 "SELECT paralog_group,symbols,gene_ids,dataset,rank FROM gene_dataset_ranks WHERE session_uuid = $1;", 
                                                                                 param=list(reactive_vals$session_id)))
                     
                     analysisParameters <- dbGetQuery(pool, 
                                                      "SELECT tissue_type,condition,organism,num_datasets,number_hkg,bootstrap_sample_size,bootstrap_replications,ensembl_release,seed FROM analysis_parameters WHERE session_uuid = $1;", 
                                                      param=list(reactive_vals$session_id))
                     
                     reactive_vals$tissues <-strsplit(analysisParameters$tissue_type, ";")
                     reactive_vals$conditions <-strsplit(analysisParameters$condition, ";")
                     reactive_vals$organisms <-strsplit(analysisParameters$organism, ";")
                     reactive_vals$num_selected_datasets <- analysisParameters$num_datasets
                     reactive_vals$numberGenes <- analysisParameters$number_hkg
                     reactive_vals$bootstrap_sample_size <- analysisParameters$bootstrap_sample_size
                     reactive_vals$bootstrap_replications <- analysisParameters$bootstrap_replications
                     if (is.na(analysisParameters$ensembl_release))
                       reactive_vals$ensembl_release <- NULL
                     else
                       reactive_vals$ensembl_release <- analysisParameters$ensembl_release
                     reactive_vals$seed <- analysisParameters$seed
                   }
                   
                   if (!is.null(state$input$plotsNumHKG))
                     reactive_vals$restorePlotsNumHKG <- state$input$plotsNumHKG
                   
                   if (!is.null(state$input$rankPlotlyHideOutliers))
                     updateCheckboxInput(session, "rankPlotlyHideOutliers", value=state$input$rankPlotlyHideOutliers)
                   
                   if(!is.null(state$input$houseKeepingGenesTable_rows_selected))
                     reactive_vals$restore_hkg_table_selection <- state$input$houseKeepingGenesTable_rows_selected
                   
                   reactive_vals$restoreState <- NULL
                   reactive_vals$doRestoreState <- F
                 } else print('noRestoreState')
               })
  
  onRestore(function(state) {
    print("onRestore")
    
    if (!is.null(state$values$session)) {
      status <- getSessionStatus(pool, state$values$session)
      if (status == "finished")
        updateTabItems(session, "menu", selected="results")
      else if (status == "started")
        updateTabItems(session, "menu", selected="status")
      else {
        print("invalid session id")
        return()
      }
      reactive_vals$session_id <- state$values$session
    }
    
    reactive_vals$restoreState <- state
  })
  
  ### RENDERER ----
  
  output$menu <- renderMenu({
    print("output$menu")
    
    session_status <- getSessionStatus(pool, reactive_vals$session_id)
    
    # Status
    if (session_status == "started")
      menu <- sidebarMenu(id = "menu",
                          menuItem("Status", selected = T, tabName = "status"))
    # Results
    else if (session_status == "finished")
      menu <- sidebarMenu(
        id = "menu",
        menuItem(
          "Start new Analysis",
          icon = icon("home"),
          href = housekeepr_external_url,
          newtab = F
        ),
        menuItem("Results", selected = T, tabName = "results")
      )
    # Start new
    else if (housekeepr_mode == "computation"){
      example_session_id <- dbGetQuery(pool, 
                                       "SELECT session_uuid FROM analysis_parameters NATURAL JOIN sessions WHERE tissue_type= 'D001921' AND condition = 'D002545;D007511;ischemic brain' AND organism = 'Mus musculus;Rattus norvegicus' AND num_datasets = '12' AND number_hkg = '50' AND bootstrap_sample_size = '12' AND bootstrap_sample_size = '12' AND bootstrap_replications = '10000' AND ensembl_release = '92' AND finish_ts IS NOT NULL LIMIT 1;")
      if (nrow(example_session_id) == 0) 
        old_session_input <- textInput(
          inputId = "old_session",
          label = make_color("Session ID to browse previous results"),
          placeholder = "86667541647329918976"
        )
      else 
        old_session_input <- textInput(
          inputId = "old_session",
          label = make_color("Session ID to browse previous results"),
          value = example_session_id[1,1]
        )
      menu <- sidebarMenu(id = "menu",
                          menuItem(tags$b("Browse old results"),
                                   old_session_input,
                                   uiOutput("use_old_session"),
                                   icon = icon("globe"),
                                   tabName = "browse_sessions",
                                   startExpanded = T),
                          menuItem(tags$b("Use example data"),
                                   icon = icon("table"),
                                   href = sprintf("%s/?_inputs_&select_conditions_all=[\"D002545\",\"D007511\",\"ischemic brain\"]&select_tissues_all=\"D001921\"&bootstrapRepetitions=10000&numberGenes=50&select_organisms_all=[\"Mus musculus\",\"Rattus norvegicus\"]&bootstrapSampleSize=12&_values_&ds=[\"GSE100235.GPL17117\",\"GSE78731.GPL15084\",\"GSE93376.GPL1261\",\"GSE52001.GPL14746\",\"GSE61616.GPL1355\",\"GSE58720.GPL10787\",\"GSE21136.GPL1355\",\"GSE33725.GPL7294\",\"GSE17929.GPL85\",\"GSE17929.GPL341\",\"GSE4206.GPL341\",\"GSE5315.GPL85\"]&ensembl=\"92\"&samples={\"GSE5315.GPL85\":[\"GSM120462\",\"GSM120461\"],\"GSE52001.GPL14746\":[\"GSM1257055\",\"GSM1257056\",\"GSM1257057\",\"GSM1257058\",\"GSM1257059\",\"GSM1257060\"],\"GSE58720.GPL10787\":[\"GSM1418664\",\"GSM1418665\",\"GSM1418666\",\"GSM1418667\",\"GSM1418668\",\"GSM1418669\"],\"GSE61616.GPL1355\":[\"GSM1509427\",\"GSM1509428\",\"GSM1509429\",\"GSM1509430\",\"GSM1509432\",\"GSM1509431\",\"GSM1509433\",\"GSM1509434\",\"GSM1509436\",\"GSM1509435\"],\"GSE78731.GPL15084\":[\"GSM2074821\",\"GSM2074822\",\"GSM2074823\",\"GSM2074824\",\"GSM2074825\",\"GSM2074826\",\"GSM2074827\",\"GSM2074828\",\"GSM2074829\",\"GSM2074830\",\"GSM2074831\"],\"GSE93376.GPL1261\":[\"GSM2452175\",\"GSM2452176\",\"GSM2452177\"],\"GSE100235.GPL17117\":[\"GSM2675546\",\"GSM2675547\",\"GSM2675548\",\"GSM2675549\",\"GSM2675551\",\"GSM2675550\"],\"GSE17929.GPL85\":[\"GSM447894\",\"GSM447896\",\"GSM447900\"],\"GSE17929.GPL341\":[\"GSM447899\",\"GSM447897\",\"GSM447895\"],\"GSE21136.GPL1355\":[\"GSM528795\",\"GSM528818\",\"GSM528817\",\"GSM528816\",\"GSM528815\",\"GSM528814\",\"GSM528813\",\"GSM528812\",\"GSM528811\",\"GSM528810\",\"GSM528809\",\"GSM528808\",\"GSM528807\",\"GSM528806\",\"GSM528805\",\"GSM528804\",\"GSM528803\",\"GSM528802\",\"GSM528801\",\"GSM528800\",\"GSM528799\",\"GSM528798\",\"GSM528797\",\"GSM528796\"],\"GSE33725.GPL7294\":[\"GSM833452\",\"GSM833453\",\"GSM833454\",\"GSM833456\",\"GSM833457\",\"GSM833459\",\"GSM833460\",\"GSM833462\",\"GSM833463\",\"GSM833467\",\"GSM833465\",\"GSM833466\"],\"GSE4206.GPL341\":[\"GSM95954\",\"GSM95955\",\"GSM95958\",\"GSM95959\"]}",
                                                  housekeepr_external_url),
                                   newtab = F),
                          menuItem(
                            tags$b("Start new Analysis"),
                            selected = T,
                            icon = icon("home"),
                            tabName = "start"
                          ))
    }
      
    menu
  })
  
  output$uiSidebar <- renderUI({
    # render after menu
    if (is.null(input$menu))
      return()
    print("renderUiSidebar")
    session_id <- isolate(reactive_vals$session_id)
    if (is.null(session_id)) {
      
      res <- div(
        
        uiOutput("select_tissues_all") %>% withSpinner(size=.25, proxy.height="30px"),
        uiOutput("select_tissues_all_alert"),
        uiOutput("select_conditions_all") %>% withSpinner(size=.25, proxy.height="30px"),
        uiOutput("select_conditions_all_alert"),
        uiOutput("select_organisms_all") %>% withSpinner(size=.25, proxy.height="30px"),
        uiOutput("select_organisms_all_alert"),
        # tags$div(style="margin-top: 15px; margin-left: 15px; margin-right: 15px;", 
        #          tags$b(style = "margin-bottom: 5px;", make_color("Normalization Algorithm")),
        #          radioButtons("algorithm", label = "",
        #                       choiceNames = list(HTML(make_color("HouseKeepR")), HTML(make_color("Normfinder"))),
        #                       choiceValues = c("HouseKeepR", "Normfinder"))),
        # #TODO: here make conditional panel for Normfinder or HouseKeepR
        {
          iDef <- 10
          iMin <- 1
          iMax <- 50
          iTitle <- "Number of candidate genes"
          i <-
            numericInput(
              "numberGenes",
              make_color(iTitle),
              ifelse(
                is.null(isolate(reactive_vals$restore_numberGenes)),
                iDef,
                isolate(reactive_vals$restore_numberGenes)
              ),
              min = iMin,
              max = iMax,
              step = 1
            ) %>%
            shinyInput_label_embed(
              shiny_iconlink() %>%
                bs_embed_popover(
                  title = iTitle, content = div(
                    p("Specify here the number of candidate house-keeping genes you are interested in."),
                    p("This setting is mandatory for controlling the required computational resources of your analysis. Increasing this number will increase the required computational resources."),
                    p(strong("Allowed values"),sprintf(" lie within [%d,%d]", iMin, iMax))
                  ),
                  placement = "right",
                  container = "body",
                  html = TRUE
                )
            )
          i
        },
        uiOutput("bootstrapSampleSize"),
        numericInput(
          "bootstrapRepetitions",
          make_color("Bootstrapping replications"),
          ifelse(
            is.null(isolate(
              reactive_vals$restore_bootstrapRepetitions
            )),
            100,
            isolate(reactive_vals$restore_bootstrapRepetitions)
          )
        ) %>% 
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_attach_modal(
                id_modal='modal_bootstrapping'
              )
          ),
        selectizeInput("customEnsemblRelease", make_color("Ensembl release"), 
                       choices=isolate(reactive_vals$available_ensembl_releases), selected=isolate(reactive_vals$ensembl_release)) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_attach_modal(
                id_modal='modal_ensembl'
              )
          ),
        numericInput("seed",
                     make_color("Random seed"),
                     ifelse(
                       is.null(isolate(reactive_vals$restore_seed)),
                       as.integer(floor(runif(1, max = .Machine$integer.max))),
                       isolate(reactive_vals$restore_seed)
                     )) %>% 
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_attach_modal(
                id_modal='modal_seed'
              )
          ),
        actionButton("start_analysis", "Start analysis"),
        # textInput("old_session", "Session ID to browse old results", value = "", width = NULL,
        #           placeholder = "39977781916968525824"),
        # TODO: remove debugg button
        # actionButton("browser", "browser"),
        use_bs_tooltip(),
        use_bs_popover()
      )
    } else if (getSessionStatus(pool, session_id) == "finished") {
      # showing results,
      res <- div(
        uiOutput("plotsNumHKG") #%>% withSpinner()
      )
    } else {
      res <- NULL
    }
    
    res
  })
  
  output$use_old_session <- renderUI({
    print("output$use_old_session")
    if (is.null(reactive_vals$old_session) ||
        reactive_vals$old_session == "")
      res <-
        shiny::HTML(
          "<div class=\"alert alert-warning\" style=\"padding: 5px; margin-bottom: 0px; margin-top: 0; margin-left: 15px; margin-right: 15px;\">Please insert a session id<br>corresponding to a finished<br>analysis.</div>"
        )
    else
      res <-
        actionButton(
          "use_old_session_button",
          " Browse previous result",
          #icon("th"),
          onclick = paste0(
            "location.href='",
            # session$clientData$url_protocol,
            # "//", session$clientData$url_hostname, ":",
            # session$clientData$url_port, session$clientData$url_pathname,
            housekeepr_external_url,
            "/?_values_&session=\"",
            reactive_vals$old_session,
            "\"';"
          )
        )
    res
  })
  
  output$plotsNumHKG <- renderUI({
    value <- isolate(reactive_vals$plotsNumHKG)
    
    if (!is.null(isolate(reactive_vals$restorePlotsNumHKG))) {
      value <- isolate(reactive_vals$restorePlotsNumHKG)
      reactive_vals$plotsNumHKG <- value
      reactive_vals$restorePlotsNumHKG <- NULL
    }
    sliderInput("plotsNumHKG", "Number of top genes to show", value=value, min=1, step=1, max=reactive_vals$numberGenes)
  })
  
  output$select_conditions_all <- renderUI({
    print("output$select_conditions_all")
    choices <- reactive_vals$restore_conditions_all
    if (!is.null(choices)) {
      selection <- isolate(reactive_vals$restore_conditions_selected)
      #reactive_vals$restore_conditions_all <- NULL
      # reactive_vals$restore_conditions_selected <- NULL
    } else {
      choices <- unlist(allConditions)
      selection <- NULL
    }
    
    res <- selectizeInput(inputId = "select_conditions_all",
                          label = make_color("Condition"), choices = choices, multiple = TRUE, selected=selection,
                          options = list(selectOnTab = F, maxOptions = 1000, openOnFocus = T,
                                         closeAfterSelect = F, hideSelected = T, create = F))
    res
  })
  
  output$select_organisms_all <- renderUI({
    print("output$select_organisms_all")
    choices <- reactive_vals$restore_organisms_all
    if (!is.null(choices)) {
      selection <- isolate(reactive_vals$restore_organisms_selected)
      #reactive_vals$restore_organisms_all <- NULL
      # reactive_vals$restore_organisms_selected <- NULL
    } else {
      choices <- unlist(allOrganisms)
      selection <- NULL
    }
    
    res <- selectizeInput(inputId = "select_organisms_all",
                          label = make_color("Organisms"), choices = choices, multiple = TRUE, selected=selection,
                          options = list(selectOnTab = F, maxOptions = 1000, openOnFocus = T,
                                         closeAfterSelect = F, hideSelected = T, create = F))
    res
  })
  
  output$select_tissues_all <- renderUI({
    print("output$select_tissues_all")
    choices <- reactive_vals$restore_tissues_all
    if (!is.null(choices)) {
      selection <- isolate(reactive_vals$restore_tissues_selected)
      #reactive_vals$restore_tissues_all <- NULL
      # reactive_vals$restore_tissues_selected <- NULL
    } else {
      choices <- unlist(allTissues)
      selection <- NULL
    }
    res <- selectizeInput(inputId = "select_tissues_all",
                          label = make_color("Tissue type"), choices = choices, multiple = TRUE, selected=selection,
                          options = list(selectOnTab = F, maxOptions = 1000, openOnFocus = T,
                                         closeAfterSelect = F, hideSelected = T, create = F))
    res
  })
  
  # output$old_session_id_alert <- renderUI({
  #   
  #   if (is.null(reactive_vals$old_session) || reactive_vals$old_session == "")
  #     shiny::HTML("<div class=\"alert alert-warning\" style=\"padding: 5px; margin-bottom: 0px; margin-top: 0; margin-left: 15px; margin-right: 15px;\">Please insert a session id<br>
  #                                                                                                                                                       corresponding to a finished<br>analysis.</div>")
  # })
  
  output$select_tissues_all_alert <- renderUI({
    if (length(reactive_vals$tissues) == 0)
      tags$div(class="alert alert-warning", style="padding: 5px; margin-bottom: 0px; margin-top: -15px; margin-left: 15px; margin-right: 15px;", "Please select a tissue")
  })
  
  output$select_conditions_all_alert <- renderUI({
    if (length(reactive_vals$conditions) == 0)
      tags$div(class="alert alert-warning", style="padding: 5px; margin-bottom: 0px; margin-top: -15px; margin-left: 15px; margin-right: 15px;", "Please select a condition")
  })
  
  output$select_organisms_all_alert <- renderUI({
    if (length(reactive_vals$organisms) == 0)
      tags$div(class="alert alert-warning", style="padding: 5px; margin-bottom: 0px; margin-top: -15px; margin-left: 15px; margin-right: 15px;", "Please select an organism")
  })
  # 
  # output$uiCustomEnsemblRelease <- renderUI({
  #   selectizeInput("customEnsemblRelease", "Use a specific ensembl release:", 
  #                  choices=NULL, selected=NULL)
  # })

  
  output$uiAnalysisError <- renderUI({
    if (!is.null(reactive_vals$analysis_error)) {
      fluidRow(
        box(status="danger",
            title="Analysis failed",
            solidHeader = T,
            width=12,
            p(class="lead", "HouseKeepR failed with the following error message:"),
            p(class="lead", paste(reactive_vals$analysis_error))
        )
      )
    }
  })
  
  output$houseKeepingGenesTable <- DT::renderDataTable({
    # render after sidebar
    if (is.null(input$plotsNumHKG))
      return()
    
    print("renderHouseKeepingGenesTable")
    
    plotsNumHKG <- reactive_vals$plotsNumHKG
    num_selected_datasets <- isolate(reactive_vals$num_selected_datasets)
    
    d <- subset(reactive_vals$bootstrappingRankTableData,
                select=-1,
                subset = rank <= plotsNumHKG)
    d[, missing_datasets := paste0(
      num_selected_datasets-missing_datasets,
      "/",
      num_selected_datasets)]
    d[, rank_var := round(rank_var, 4)]
    d[, rank_mean := round(rank_mean, 4)]
    
    d$symbols <- sapply(d$symbols, FUN=function(x) {
      paste0(strsplit(x, "###",fixed=T)[[1]],collapse="\n")
    })
    
    setnames(d, c("Gene Symbol(s)", "Ensembl ID(s)", "#Data sets", "Rank Mean", "Rank Variance", "Rank"))
    dt <- datatable(d, selection=list("row", selected=isolate(reactive_vals$restore_hkg_table_selection)), rownames=F, options = list(
      order = list(list(5, 'asc'))
    ))
    
    reactive_vals$restore_hkg_table_selection <- NULL
    reactive_vals$doRenderTable <- F
    reactive_vals$delayCounterRenderResultDatasetChart <- 0
    dt
  },
  server=T)
  
  output$bootstrapSampleSize <- renderUI({
    div(
      numericInput("bootstrapSampleSize", 
                   make_color("Bootstrapping sample size"), 
                   value = if(is.null(reactive_vals$restore_bootstrap_sample_size)) {
                     reactive_vals$num_selected_datasets
                   } else {
                     min(reactive_vals$restore_bootstrap_sample_size,
                         nrow(reactive_vals$available_datasets))
                   },#max(1, min(10, reactive_vals$num_selected_datasets-1)), 
                   min = 1, 
                   max = reactive_vals$num_selected_datasets#reactive_vals$num_selected_datasets-1
      ) %>%
        shinyInput_label_embed(
          shiny_iconlink() %>%
            bs_attach_modal(
              id_modal='modal_bootstrapping'
            )
        ),
      use_bs_tooltip(),
      use_bs_popover()
    )
  })
  
  output$queryingGEOdatasets <- renderUI({
    if (!reactive_vals$queryingGEOdatasets)
      return(NULL)
    
    div(style="text-align: center;",
        shiny::img(src = "ajax-loader.gif"),
        p(class="lead", "Querying GEO for data sets"))
  })
  
  output$queryingGEOsamples <- renderUI({
    if (!reactive_vals$queryingGEOsamples)
      return(NULL)
    
    div(style="text-align: center;",
        shiny::img(src = "ajax-loader.gif"),
        p(class="lead", "Querying GEO for data set samples"))
  })
  
  output$dataSetSelection <- renderUI({
    if (!reactive_vals$tissue_selected || !reactive_vals$condition_selected || !reactive_vals$organism_selected){
      inner <- p(class="lead", "Please make all selections in the sidebar.")
      #inner <- NULL
    }
    else if (reactive_vals$queryingGEOdatasets) {
      #TODO: uiOutput didnt work
      #inner <- uiOutput("queryingGEOdatasets")
      inner <- div(style="text-align: center;",
                   shiny::img(src = "ajax-loader.gif"),
                   p(class="lead", "Querying GEO for data sets"))
    }
    # else if (is.null(reactive_vals$available_datasets))
    #   inner <- NULL
    # added !is.null(reactive_vals$available_datasets) because it is null if nothing is found
    else if (!is.null(reactive_vals$available_datasets)) {# & nrow(reactive_vals$available_datasets) > 0) {
      inner <- div(
        p(class="lead", "We found the following data sets having the conditions, tissue types and organisms you specified. Which data sets should be used to identify housekeeping genes?"),
        DT::dataTableOutput('select_datasets_table'),
        actionButton("select_all_datasets", "Select all"),
        actionButton("select_none_datasets", "Select none")
      )
    } else {
      inner <- p(class="lead", "We could not find any data sets for your criteria.")
    }
    inner
  })
  
  output$sampleSelection <- renderUI({
    if(reactive_vals$num_selected_datasets < 2 || !reactive_vals$tissue_selected || !reactive_vals$condition_selected || !reactive_vals$organism_selected)
      inner <- p(class="lead", "Please make all selections in the sidebar and select at least two data sets.")
    else if (reactive_vals$queryingGEOsamples) {
      #TODO: uiOutput does not work 
      #inner <- uiOutput("queryingGEOsamples")
      inner <-  div(style="text-align: center;",
                    shiny::img(src = "ajax-loader.gif"),
                    p(class="lead", "Querying GEO for data set samples"))
    } else {
      all_sel_samples <- isolate(reactive_vals$selected_samples)
      
      inner <- div(p(class="lead", "Please specify condition and control samples."),
                   {
                     res <- list()
                     for (gse in unique(reactive_vals$available_samples$accession)) {
                       for (GPL in unique(reactive_vals$available_samples[accession==gse]$gpl)) {
                         gse_gpl <- paste(gse, GPL, sep=".")
                         ds <- reactive_vals$available_samples[accession==gse & gpl==GPL]
                         if (is.null(all_sel_samples) || !(gse_gpl %in% names(all_sel_samples)))
                           ds_sel_samples <- c()
                         else
                           ds_sel_samples <- all_sel_samples[[gse_gpl]]
                         res[[length(res)+1]] <- div(
                           h4(paste("Data set", gse, "(", GPL , ")")),
                           wellPanel(
                             apply(ds, MARGIN=1, function(s) {
                               if (!(gse_gpl %in% reactive_vals$selected_datasets$gse_gpl))
                                 return()
                               
                               inputId <- paste("sample_radio", s["samples.accession"], sep="_")
                               reactive_vals$bookmark_exclude <- unique(c(reactive_vals$bookmark_exclude, inputId))
                               rbs <- radioButtons(inputId = inputId,
                                                   paste(s["samples.title"],paste("(",s["samples.accession"],")",sep="")),
                                                   c("Condition" = "disease_sample", "Control" = "control_sample"),
                                                   selected=if (s["samples.accession"] %in% ds_sel_samples) "disease_sample" else "control_sample",
                                                   inline=T)
                               observeEvent(input[[inputId]], {
                                 all_sel_samples <- isolate(reactive_vals$selected_samples)
                                 
                                 if (is.null(all_sel_samples))
                                   all_sel_samples <- list()
                                 
                                 sId <- unname(s["samples.accession"])
                                 
                                 gse_gpl <- paste(reactive_vals$available_samples[samples.accession==sId,.(accession,gpl)],collapse=".")
                                 
                                 if (input[[inputId]] == "disease_sample")
                                   all_sel_samples[[gse_gpl]] <- unique(c(all_sel_samples[[gse_gpl]], sId))
                                 else if (length(all_sel_samples[[gse_gpl]]) > 1)
                                   all_sel_samples[[gse_gpl]] <- setdiff(all_sel_samples[[gse_gpl]], sId)
                                 else
                                   all_sel_samples[[gse_gpl]] <- NULL
                                 reactive_vals$selected_samples <- all_sel_samples
                                 session$doBookmark()
                               }, ignoreInit=T)
                               rbs
                             })
                           )
                         )
                       }
                     }
                     res
                   })
    }
    inner
  })
  
  output$uiStartAnalysis <- renderUI({
    p("In order to start a new analysis you need to",
      {p <- prettyCheckbox(
        inputId = "pretty_1", 
        label = "Select a tissue type, condition, and organism in the sidebar",
        icon = icon("check"),
        status = "success", 
        outline = TRUE,
        value = reactive_vals$tissue_selected && reactive_vals$organism_selected && reactive_vals$condition_selected
      )
      disablePrettyCheckbox(p)},
      {p <- prettyCheckbox(
        inputId = "pretty_2", 
        label = "Specify the desired number of house-keeping genes",
        icon = icon("check"),
        status = "success", 
        outline = TRUE,
        value = !is.null(reactive_vals$numberGenes) && 
          !is.na(reactive_vals$numberGenes) && 
          reactive_vals$numberGenes != ""
      )
      disablePrettyCheckbox(p)},
      {p <- prettyCheckbox(
        inputId = "pretty_3", 
        label = "Configure the bootstrap sample size and the number of bootstrap replications",
        icon = icon("check"),
        status = "success", 
        outline = TRUE,
        value = (!is.null(reactive_vals$bootstrap_sample_size) && 
                   !is.na(reactive_vals$bootstrap_sample_size) && 
                   reactive_vals$bootstrap_sample_size != "") && 
          (!is.null(reactive_vals$bootstrap_replications) &&
             !is.na(reactive_vals$bootstrap_replications) && 
             reactive_vals$bootstrap_replications != "")
      )
      disablePrettyCheckbox(p)},
      {p <- prettyCheckbox(
        inputId = "pretty_4", 
        label = "Select at least two data sets",
        icon = icon("check"),
        status = "success", 
        outline = TRUE,
        value = reactive_vals$num_selected_datasets > 1
      )
      disablePrettyCheckbox(p)},
      {p <- prettyCheckbox(
        inputId = "pretty_5", 
        label = "Select at least one disease and one control sample for each selected data set",
        icon = icon("check"),
        status = "success", 
        outline = TRUE,
        value = reactive_vals$num_selected_datasets > 1 && reactive_vals$valid_sample_selection
      )
      disablePrettyCheckbox(p)},
      {p <- prettyCheckbox(
        inputId = "pretty_6",
        label = "Hit 'Start analysis'",
        icon = icon("check"),
        status = "success", 
        outline = TRUE,
        value = reactive_vals$analysis_started
      )
      disablePrettyCheckbox(p)}
    )
  })
  
  output$uiAnalysisRunning <- renderUI({
    fluidRow(
      box(title="Analysis running", width = 12,
          solidHeader = T,
          collapsible = T,
          status="primary",
          p(class="lead", "We are currently processing your analysis. This may take a while ..."),
          # p("Once it is ready your result will be available ",
          #   a(href=paste0(housekeepr_external_url,"/?_values_&session=\"", reactive_vals$session_id, "\""), "here")
          # )
          p(class="lead", strong("Please do not close this tab or start another HouseKeepR analysis until this is finished."))
      )
    )
  })
  
  output$uiRankPlotly <- renderUI({
    # render after table
    if (is.null(input$houseKeepingGenesTable_search))
      return()
    
    print("renderBootstrapPlot")
    
    if (length(reactive_vals$selected_hkgs) > 0)
      plotsNumHKG <- length(reactive_vals$selected_hkgs)
    else 
      plotsNumHKG <- reactive_vals$plotsNumHKG
    numberDS <- isolate(reactive_vals$num_selected_datasets)
    
    div(
      jqui_resizable(
        plotOutput("rankPlotly", height=paste(as.character(max(220, 150+25*plotsNumHKG)),"px",sep="")) %>% withSpinner(proxy.height="50px")
      ),
      div(style="display: inline-block;",
          checkboxInput("rankPlotlyHideOutliers", "Hide outliers in box plot", value=T)
      ),
      div(style="display: inline-block;",
          downloadButton("downloadBootstrapRanks", "Download Bootstrap ranks")
      )
    )
  })
  
  output$rankPlotly <- renderPlot({
    print("renderRankPlotly")
    if (nrow(reactive_vals$bootstrappingLongRankTableData) > 0) {
      if (length(reactive_vals$selected_hkgs) > 0) {
        tab <- reactive_vals$bootstrappingLongRankTableData[paralog_group %in% reactive_vals$selected_hkgs]
      } else
        tab <- reactive_vals$bootstrappingLongRankTableData
      
      if(!is.null(reactive_vals$bootstrappingRankTableData))
        tab <- tab[paralog_group %in% reactive_vals$bootstrappingRankTableData[rank <= reactive_vals$plotsNumHKG, paralog_group]]
      
      tab <- tab[!is.null(symbols) & !is.na(symbols)]
      setorder(tab, rank)
      tab$symbols <- factor(tab$symbols, levels=
                              rev(unique(reactive_vals$bootstrappingRankTableData[order(rank), symbols])))
      
      levels(tab$symbols) <- sapply(
        levels(tab$symbols), 
        function(x) {
          # first element: symbols
          # second element: ensembl ids
          spl0 <- str_split(x, "###")[[1]]
          for (i in 1:length(spl0)) {
            spl <- str_split(spl0[i], ", ")[[1]]
            all <- c(spl[!startsWith(spl, "LOC")], spl[startsWith(spl, "LOC")])
            if (length(all) > 1) {
              if (i == 1)
                spl0[i] <- (paste(all[1], "...", sep=", ")) 
              else if (i == 2)
                spl0[i] <- (paste(all[1], "...)", sep=", "))
            }
            else 
              spl0[i] <- (all[1])
          }
          paste0(spl0, collapse="###")
        }
      )
      
      if (reactive_vals$rankPlotlyHideOutliers)
        yMax <- max(tab[, yMax:=boxplot.stats(rank)$stats[5], by=symbols]$yMax,na.rm=T)
      else
        yMax <- max(tab$rank)
      
      levels(tab$symbols) <- sapply(levels(tab$symbols), FUN=function(x) {
        paste0(strsplit(x, "###",fixed=T)[[1]],collapse="\n")
      })
      
      ggplot(tab, aes(x=symbols, y=rank)) + 
        geom_boxplot(outlier.shape=if (reactive_vals$rankPlotlyHideOutliers) NA else (19)) + 
        xlab("Gene") + 
        ylab("Ranks across Bootstrap Samples") + 
        scale_fill_ptol() +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylim(1,yMax) + coord_flip(expand=F)
    }
    else {
      NULL
    }
  }, res=96)
  
  output$uiRankLongPlotly <- renderUI({
    # render after table
    if (is.null(input$houseKeepingGenesTable_search))
      return()
    
    print("renderDatasetPlot")
    
    if (length(reactive_vals$selected_hkgs) > 0)
      plotsNumHKG <- length(reactive_vals$selected_hkgs)
    else 
      plotsNumHKG <- reactive_vals$plotsNumHKG
    numberDS <- isolate(reactive_vals$num_selected_datasets)
    div(
      jqui_resizable(
        plotlyOutput("rankLongPlotly", height=paste(as.character(max(220, 150+25*plotsNumHKG)),"px",sep="")) %>% withSpinner(proxy.height="50px")
      ),
      downloadButton("downloadDatasetRanks", "Download Dataset ranks")
    )
  })
  
  output$rankLongPlotly <- renderPlotly({
    print("renderRankLongPlotly")
    if (nrow(reactive_vals$rankLongTableData) > 0) {
      if (length(reactive_vals$selected_hkgs) > 0) {
        tab <- reactive_vals$rankLongTableData[paralog_group %in% reactive_vals$selected_hkgs]
      } else
        tab <- reactive_vals$rankLongTableData
      
      if(!is.null(reactive_vals$bootstrappingRankTableData))
        tab <- tab[paralog_group %in% reactive_vals$bootstrappingRankTableData[rank <= reactive_vals$plotsNumHKG, paralog_group]]
      
      tab <- tab[!is.null(symbols) & !is.na(symbols)]
      setorder(tab, rank)
      tab$symbols <- factor(tab$symbols, levels=
                              rev(unique(reactive_vals$bootstrappingRankTableData[order(rank), symbols])))
      
      levels(tab$symbols) <- sapply(
        levels(tab$symbols), 
        function(x) {
          # first element: symbols
          # second element: ensembl ids
          spl0 <- str_split(x, "###")[[1]]
          for (i in 1:length(spl0)) {
            spl <- str_split(spl0[i], ", ")[[1]]
            all <- c(spl[!startsWith(spl, "LOC")], spl[startsWith(spl, "LOC")])
            if (length(all) > 1) {
              if (i == 1)
                spl0[i] <- (paste(all[1], "...", sep=", ")) 
              else if (i == 2)
                spl0[i] <- (paste(all[1], "...)", sep=", "))
            }
            else 
              spl0[i] <- (all[1])
          }
          paste0(spl0, collapse="###")
        }
      )
      
      tab[,negLogRank := -log2(rank)]
      tab[rank <= 20, c("rankStr", "transRank") := list("<= 20", 1)]
      tab[rank > 20 & rank <= 50, c("rankStr", "transRank") := list("<= 50", 0.85)]
      tab[rank > 50 & rank <= 100, c("rankStr", "transRank") := list("<= 100", 0.7)]
      tab[rank > 100, c("rankStr", "transRank") := list("> 100", 0.5)]
      tab[, rankStr := factor(rankStr, levels = c("<= 20","<= 50","<= 100","> 100"))]
      
      dsWithMultiGPL <- unlist(lapply(strsplit(unique(tab$dataset), ".", fixed=T),function(x) x[1]))
      dsWithMultiGPL <- dsWithMultiGPL[duplicated(dsWithMultiGPL)]
      
      matchesMultiGPL <- function(col) {
        sapply(col, function(dataset) any(stri_startswith(dataset, fixed=dsWithMultiGPL)))
      }
      
      extractGse <- function(col) {
        sapply(col, function(gse_gpl) strsplit(gse_gpl, ".", fixed=T)[[1]][1])
      }
      
      tab[!matchesMultiGPL(dataset), dataset:=extractGse(dataset)]
      
      levels(tab$symbols) <- sapply(levels(tab$symbols), FUN=function(x) {
        paste0(strsplit(x, "###",fixed=T)[[1]],collapse="\n")
      })
      
      ggplotly(
        ggplot(tab, aes(y=symbols,x=dataset, text=paste("Rank:", rank))) + 
          ylab("Gene") + 
          xlab("Dataset") + 
          geom_point(aes(colour=rankStr, size=transRank*16), shape=16) +
          geom_text(aes(label=sapply(rank, function(x) if (x <= 100) paste(x) else ""), size=transRank*6), colour="white", hjust=0, vjust=0) +
          theme_minimal() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_color_manual(name="Rank", 
                             values = c("<= 20"="#005073", 
                                        "<= 50"="#107dac", 
                                        "<= 100"="#189ad3", 
                                        "> 100"="#71c7ec"))
        ,
        dynamicTicks=T,
        tooltip=c("text", "x", "y")
      )
    }
    else {
      NULL
    }
  })
  
  #### boxes ----
  output$boxes <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    
    div(
      shiny::HTML(sprintf(
        "<p class=\"lead\">After the analysis is done, you can come back to the results anytime using the session-id %s. You may want to write it down for later.</p>",
        make_color(reactive_vals$session_id, "red")
      )),
      fluidRow(
        uiOutput("selectedTissueBox"),
        uiOutput("selectedConditionBox"),
        uiOutput("selectedOrganismBox")
      ),
      fluidRow(
        uiOutput("numberDatasetsBox"),
        uiOutput("numberGenesBox"),
        uiOutput("bootstrapSampleSizeBox"),
        uiOutput("bootstrapReplicationsBox"),
        uiOutput("ensemblReleaseBox"),
        uiOutput("resultTimestampBox")
      )
    )
  })
  
  output$selectedTissueBox <-renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    
    m <- match(reactive_vals$tissues, allTissues)
    
    infoBox(
      value=paste0(c(names(allTissues[na.omit(m)]), reactive_vals$tissues[is.na(m)]), collapse=", "), 
      title="Selected Tissue Type(s)", 
      icon = icon("tint"),
      color = "orange"
    )
  })
  
  output$selectedConditionBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    
    m <- match(unlist(reactive_vals$conditions), allConditions)
    
    infoBox(
      value=paste0(
        c(
          names(allConditions[na.omit(m)]), 
          unlist(reactive_vals$conditions)[is.na(m)]), collapse=", "), 
      title="Selected Condition(s)", 
      icon = icon("ambulance"),
      color = "orange"
    )
  })
  
  output$selectedOrganismBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    infoBox(
      value=paste0(unlist(reactive_vals$organisms), collapse=", "), 
      title="Selected Organism(s)", 
      icon = icon("male"),
      color = "orange"
    )
  })
  
  output$numberDatasetsBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    valueBox(
      paste0(reactive_vals$num_selected_datasets), "Number of Selected Datasets", 
      icon = icon("table"),
      color = "purple", width=2
    )
  })
  
  output$numberGenesBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    valueBox(
      paste0(reactive_vals$numberGenes), "Number of House-keeping Genes", 
      icon = icon("key"),
      color = "purple", width=2
    )
  })
  
  output$bootstrapSampleSizeBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    valueBox(
      paste0(reactive_vals$bootstrap_sample_size), "Bootstrap Sample Size", 
      color = "purple", width=2, icon=icon("arrows-alt-v")
    )
  })
  
  output$bootstrapReplicationsBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    valueBox(
      paste0(reactive_vals$bootstrap_replications), "Bootstrap Replications",
      color = "purple", width=2, icon=icon("redo")
    )
  })
  
  output$ensemblReleaseBox <- renderUI({
    if (!(getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished")))
      return(NULL)
    infoBox(
      "Ensembl Release",
      if (is.null(reactive_vals$ensembl_release))
        "Most Recent"
      else
        paste0(reactive_vals$ensembl_release)
      ,
      color = "purple", width=2
    )
  })
  
  output$resultTimestampBox <- renderUI({
    input$menu
    finishTime <- getSessionFinishtime(pool, reactive_vals$session_id)
    if (!is.na(finishTime)) {
      infoBox(
        "Analysis finished", paste0(finishTime),
        icon = icon("clock-o"),
        color = "purple", width=2
      )
    } else {
      startTime <- getSessionStarttime(pool, reactive_vals$session_id)
      if (is.na(startTime))
        return(NULL)
      
      infoBox(
        "Analysis started", paste0(startTime),
        icon = icon("clock-o"),
        color = "purple", width=2
      )
    }
  })
  
  ### /RENDERER
  
  ### OBSERVER ----
  
  observe({
    setBookmarkExclude(reactive_vals$bookmark_exclude)
  })
  
  observeEvent(input$old_session, {
    print("observeEvent(input$old_session")
    old_session_id <- isolate(input$old_session)
    if (!is.null(old_session_id) &&
        (old_session_id %in% getFinishedSessions(pool)$session_uuid))
      reactive_vals$old_session <- old_session_id
    else reactive_vals$old_session <- NULL
  })
  
  # observeEvent(input$algorithm, {
  #   print("observeEvent(input$old_session")
  #   reactive_vals$use_normfinder <- input$algorithm == "Normfinder"
  # })

  
  #TODO: remove debugg button
  # observeEvent(input$browser, {
  #   browser()
  # })
  
  observeEvent(input$select_conditions_all, {
    print("observeEvent(input$select_conditions_all, {")
    session$doBookmark()
    # val <- input$select_conditions_all
    # if (!is.null(val) && length(val) > 0) {
      reactive_vals$conditions <- isolate(input$select_conditions_all) #val
    # }
  }, ignoreInit = T, ignoreNULL=FALSE)
  
  observeEvent(input$select_organisms_all, {
    print("observeEvent(input$select_organisms_all, {")
    session$doBookmark()
    # val <- input$select_organisms_all
    # if (!is.null(val) && length(val) > 0) {
      reactive_vals$organisms <- isolate(input$select_organisms_all) #val
    # }
  }, ignoreInit = T, ignoreNULL=FALSE)
  
  observeEvent(input$select_tissues_all, {
    print("observeEvent(input$select_tissues_all, {")
    session$doBookmark()
    # val <- input$select_tissues_all
    # if (!is.null(val) && length(val) > 0) {
      reactive_vals$tissues <- isolate(input$select_tissues_all) #val
    # }
  }, ignoreInit = T, ignoreNULL=FALSE)
  
  observe({
    # do not query GEO if we are only here to show results
    if (getSessionStatus(pool, isolate(reactive_vals$session_id)) %in% c("started", "finished"))
      return()
    
    mesh_terms <- c('(("gds"[Entry Type] OR "gse"[Entry Type]))')
    if (length(reactive_vals$conditions) > 0 && reactive_vals$conditions != "C") {
      condition <- names(allConditions[which(allConditions%in%reactive_vals$conditions)])
      # user added ones do not have a corresponding entry
      condition <- unique(
        c(condition, 
          reactive_vals$conditions[!(reactive_vals$conditions %in% allConditions)]
        )
      )
      mesh_terms <- c(mesh_terms, paste(
        "(",
        paste(
          c(
            paste(condition, "[MeSH]", sep=""),
            paste(condition, "[Description]", sep="")
          ),
          collapse=" OR "), 
        ")")
      )
    } else {
      reactive_vals$rentrez_search_result <- NULL
      return()
    }
    if (length(reactive_vals$tissues) > 0 && reactive_vals$tissues != "A") {
      tissue <- names(allTissues[which(allTissues%in%reactive_vals$tissues)])
      mesh_terms <- c(mesh_terms, paste(
        "(",
        paste(
          paste(tissue, "[MeSH]", sep=""), 
          collapse=" OR "), 
        ")")
      )
    } else {
      reactive_vals$rentrez_search_result <- NULL
      return()
    }
    if (length(reactive_vals$organisms) > 0 && reactive_vals$organisms != "B") {
      organism <- reactive_vals$organisms
      mesh_terms <- c(mesh_terms, paste(
        "(",
        paste(
          paste(organism, "[Organism]", sep=""), 
          collapse=" OR "), 
        ")")
      )
    } else {
      reactive_vals$rentrez_search_result <- NULL
      return()
    }
    
    if (length(mesh_terms) > 0) {
      if (!isolate(reactive_vals$queryingGEOdatasets)) {
        reactive_vals$queryingGEOdatasets <- TRUE
        # allow the GUI to update
        invalidateLater(1)
        return()
      }
      
      # only find data sets where there is a series matrix available
      # we cross-checked NCBI and this seems to be exclusively the case only for array technologies;
      # -> excluding all high-throughput data sets
      # mesh_terms <- c(mesh_terms, 'NOT "expression profiling by high throughput sequencing"[DataSet Type]')
      # try to only allow expression arrays instead, to exclude Methylation
      # minimum 3 samples per set, otherwise eBays does not work
      mesh_terms <- c(mesh_terms, 'Expression profiling by array [DataSet Type]', '"3"[n_samples] : "1000"[n_samples])')
      print(mesh_terms)
      search_result <- entrez_search(db="gds", term=paste(mesh_terms, collapse = " "), retmax=500)$ids
      if (length(search_result) == 0) {
        reactive_vals$available_datasets <- NULL
        # otherwise error if nothing found
        reactive_vals$queryingGEOdatasets <- FALSE
        print("processed result from NCBI: nothing found")
        return()
      }
      result_list <- list()
      batch_size <- 100
      n_batches <- ceiling(length(search_result)/batch_size)
      for (batch in 1:n_batches) {
        offset <- batch_size*(batch-1)
        # last batch, take the rest
        if (batch == n_batches)
          batch_size <- ifelse(length(search_result)%%batch_size == 0, batch_size, length(search_result)%%batch_size)
        result_list <- c(result_list, entrez_summary(db="gds", id=search_result[(offset+1):(offset+batch_size)], always_return_list = T)) # entrez_fetch
      }
      class(result_list) <- c("esummary_list", "list")
      # same datasets found as before:
      if (!is.null(isolate(reactive_vals$rentrez_search_result))) {
        if (isTRUE(all.equal.character(names(isolate(reactive_vals$rentrez_search_result)), names(result_list))))
          reactive_vals$queryingGEOdatasets <- FALSE
      }
      reactive_vals$rentrez_search_result <- result_list
      print("processed result from NCBI")
    } else {
      reactive_vals$rentrez_search_result <- NULL
      return()
    }
  })
  
  observe({
    reactive_vals$tissue_selected <- !is.null(input$select_tissues_all) && length(input$select_tissues_all) > 0
  })
  
  observe({
    reactive_vals$organism_selected <- !is.null(input$select_organisms_all) && input$select_organisms_all != ""
  })
  
  observe({
    reactive_vals$condition_selected <- !is.null(input$select_conditions_all) && length(input$select_conditions_all) > 0
  })
  
  # available data sets
  observeEvent(reactive_vals$rentrez_search_result, {
    
    search_result_summary <- reactive_vals$rentrez_search_result
    tryCatch({
      if (!is.null(search_result_summary)) {
        attrs <- c("uid", "accession", "pdat", "title", #"summary", 
                   "taxon", "n_samples", 
                   "gpl")
        
        # one esummary for all data sets
        if ("esummary" %in% class(search_result_summary))
          meta <- extract_from_esummary(search_result_summary,attrs, simplify = F)
        # or one esummary per data set batch?
        else if ("esummary_list" %in% class(search_result_summary)) {
          meta <- extract_from_esummary(search_result_summary,attrs, simplify = F)
        }
        else {
          return()
        }
        if (!is.list(meta) || !is.list(meta[[1]]))
          meta <- list(meta)
        
        names(meta) <- lapply(meta, function(x) x$accession)
        meta <- matrix(unlist(meta), ncol = length(attrs), byrow = T)
        colnames(meta) <- attrs
        #meta <- as.data.frame(meta)
        #meta$n_samples <- as.numeric(as.character(meta$n_samples))
        
        meta <- as.data.table(meta)
        meta$n_samples[grep(";", meta$gpl, fixed = T)] <- sprintf("<=%s", meta$n_samples[grep(";", meta$gpl, fixed = T)])
        meta <- cSplit(meta, "gpl", sep=";", direction="long")
        meta[,gpl:=paste0("GPL",gpl)]
        
        reactive_vals$available_datasets <- meta
      } else {
        reactive_vals$available_datasets <- NULL
      }
    }, finally = {
      reactive_vals$queryingGEOdatasets <- FALSE
    })
  }, ignoreInit = T, ignoreNULL = F)
  
  observe({
    if (!is.null(reactive_vals$available_datasets)) {
      newData <- reactive_vals$available_datasets
      selected_datasets <- isolate(reactive_vals$restore_dataset_selection)
      
      dt <- copy(reactive_vals$available_datasets)
      
      #TODO: error in setnames instead of no datasets found but available_datasets not null because of former requests
      
      dt[,uid:=NULL]
      setnames(dt, 
               c("accession",
                 "pdat",
                 "title",
                 "taxon",
                 "n_samples",
                 "gpl"), 
               c("Accession",
                 "Date",
                 "Title",
                 "Taxon",
                 "#Samples",
                 "Platform"))
      
      output$select_datasets_table <- DT::renderDataTable(
        datatable(dt,
                  selection=list(mode='multiple', 
                                 selected = which(paste(newData$accession,newData$gpl,sep=".") %in% selected_datasets)), 
                  rownames=F)
      )
    }
  })
  
  # store selected data sets
  observe({
    ds <- reactive_vals$available_datasets
    if (length(input$select_datasets_table_rows_selected) > 0)
      res <- ds[input$select_datasets_table_rows_selected,]
    else
      res <- NULL
    if (!is.null(res))
      res[,gse_gpl:=paste(accession, gpl, sep=".")]
    reactive_vals$selected_datasets <- res
    # check if samples are valid when there are actually datasets available, otherwise selected samples get deleted when using bookmark
    if (!is.null(res))
      reactive_vals$selected_samples <- reactive_vals$selected_samples[names(reactive_vals$selected_samples) %in% res$gse_gpl]
    
    session$doBookmark()
  })
  
  # output$selected_datasets <- eventReactive(input$select_datasets_table_rows_selected, {
  #   if (length(input$select_datasets_table_rows_selected) == 0)
  #     return("")
  #   
  #   ds <- reactive_vals$available_datasets
  #   res <- ds[input$select_datasets_table_rows_selected,]
  #   paste(as.vector(res$accession),sep=",")
  # }, ignoreNULL = F)
  
  observe({
    if (getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished"))
      return()
    reactive_vals$num_selected_datasets <- length(input$select_datasets_table_rows_selected)
  })
  
  # output$num_selected_datasets <- reactive({
  #   reactive_vals$num_selected_datasets
  # })
  
  observeEvent(input$numberGenes, {
    reactive_vals$numberGenes <- input$numberGenes
  })
  
  observeEvent(input$bootstrapSampleSize, {
    reactive_vals$bootstrap_sample_size <- input$bootstrapSampleSize
  })
  
  observeEvent(input$bootstrapRepetitions, {
    reactive_vals$bootstrap_replications <- input$bootstrapRepetitions
  })
  
  observeEvent(input$select_all_datasets, {
    # select all
    selectRows(DT::dataTableProxy('select_datasets_table'), 1:nrow(isolate(reactive_vals$available_datasets)))  
  })
  
  observeEvent(input$select_none_datasets, {
    # select none
    selectRows(DT::dataTableProxy('select_datasets_table'), NULL)
  })
  
  # available dataset samples
  observe({
    # do not query GEO if we are only here to show results
    if (getSessionStatus(pool, isolate(reactive_vals$session_id)) %in% c("started", "finished"))
      return()
    
    sel_ds <- reactive_vals$selected_datasets
    if(is.null(sel_ds)){
      reactive_vals$available_samples <- NULL
      return()
    }
    search_result_summary <- isolate(reactive_vals$rentrez_search_result)
    if (!is.null(search_result_summary)) {
      sample_summary <- NULL
      if (nrow(sel_ds) > 1) {
        if (!isolate(reactive_vals$queryingGEOsamples)) {
          reactive_vals$queryingGEOsamples <- TRUE
          # allow the GUI to update
          invalidateLater(1)
          return()
        }
        tryCatch({
          
          attrs <- c("gse", "gpl", "accession", "title")
          
          sample_summary <- list()
          n_samples_num <- as.numeric(sel_ds[, gsub('(^<=)?', '', n_samples)])
          max_n_samples <- 100
          curr_n_samples <- 0
          last_idx <- 0
          for (ds_idx in 1:nrow(sel_ds)) {
            # offset <- batch_size*(batch-1)
            curr_n_samples <- curr_n_samples + n_samples_num[ds_idx]
            if (curr_n_samples >= max_n_samples | ds_idx == nrow(sel_ds)) {
              sample_summary <- c(sample_summary, entrez_summary(db="gds", 
                                                                 id=entrez_search(db="gds",
                                                                                  term=paste(sprintf('(%s[Accession] AND %s[Accession] AND "gsm"[Filter])', 
                                                                                                     sel_ds$accession[(last_idx+1):(ds_idx)], 
                                                                                                     sel_ds$gpl[(last_idx+1):(ds_idx)]), 
                                                                                             collapse = " OR "),
                                                                                  retmax=500)$ids, always_return_list = T))
              last_idx <- ds_idx
              curr_n_samples <- 0
            }
          }
          class(sample_summary) <- c("esummary_list", "list")
          sample_summary <- cSplit(rbindlist(lapply(extract_from_esummary(sample_summary, attrs, F), as.data.table)), 
                                   splitCols = 'gse', sep = ";", fixed = T, direction = 'long')
          sample_summary[, gse := sprintf("GSE%s", sample_summary$gse)]
          sample_summary <- sample_summary[gse %in% sel_ds$accession]
          sample_summary[, gpl := sprintf("GPL%s", sample_summary$gpl)]
          sample_summary[, gse_gpl:=sprintf("%s.%s", gse, gpl)]
          setnames(sample_summary, "accession", "samples.accession")
          setnames(sample_summary, "title", "samples.title")
          setnames(sample_summary, "gse", "accession")
          sample_summary[unique(sel_ds[, c("title", "accession")]), on = "accession"]
          setorder(sample_summary, accession, samples.accession)
          
          print("samples done")
        }, finally = {
          reactive_vals$queryingGEOsamples <- FALSE
        })
      }
      
      reactive_vals$available_samples <- sample_summary
      
      } else {
      reactive_vals$available_samples <- NULL
    }
  })
  
  observe({
    available <- reactive_vals$available_samples
    if (is.null(available) || nrow(available) == 0)
      reactive_vals$num_available_samples <- c()
    else {
      reactive_vals$num_available_samples <- as.data.table(available %>%
                                                             dplyr::group_by(accession, gpl) %>%
                                                             dplyr::summarise(n=n()))
      reactive_vals$num_available_samples[, gse_gpl:=paste(accession,gpl,sep=".")]
    }
    #unlist(lapply(available, function(x) nrow(x$samples)))
  })
  
  observe({
    sel <- reactive_vals$selected_samples
    if (length(sel) == 0)
      reactive_vals$num_selected_samples <- c()
    else
      reactive_vals$num_selected_samples <- unlist(lapply(sel, function(x) length(x)))
  })
  
  # output$num_selected_samples <- reactive({
  #   reactive_vals$num_selected_samples
  # })
  
  observe({
    reactive_vals$valid_sample_selection <- length(reactive_vals$num_selected_samples) == reactive_vals$num_selected_datasets && 
      all(reactive_vals$num_selected_samples > 0) &&
      setequal(names(reactive_vals$num_selected_samples), reactive_vals$num_available_samples$gse_gpl) &&
      all(reactive_vals$num_selected_samples < reactive_vals$num_available_samples[match(names(reactive_vals$num_selected_samples), reactive_vals$num_available_samples$gse_gpl)])
  })
  
  # reactive expression for the selected house-keeping genes
  observe({
    reactive_vals$selected_hkgs <- subset(reactive_vals$bootstrappingRankTableData,
                                          subset = rank <= reactive_vals$plotsNumHKG)[input$houseKeepingGenesTable_rows_selected,]$paralog_group
    session$doBookmark()
  })
  
  observeEvent(input$start_analysis, {
    # assign a new session id
    reactive_vals$session_id <- getUniqueSessionId(pool)
    
    set.seed(reactive_vals$seed)
    pool::dbWriteTable(pool, 
                       "sessions",
                       data.frame(session_uuid=reactive_vals$session_id, create_ts=Sys.time()),
                       append=T)
    
    reactive_vals$analysis_started <- TRUE
    reactive_vals$analysis_running <- TRUE
    reactive_vals$analysis_error <- NULL
    session$doBookmark()
  }, ignoreInit = TRUE)
  
  # this has to be in its own observer, such that we can use invalidateLater to allow GUI update before running the script
  observe({
    if (!reactive_vals$analysis_started)
      return()
    # allow the GUI to update before we start our long-running script
    if (is.null(isolate(reactive_vals$delayCounterStartAnalysis))) {
      reactive_vals$delayCounterStartAnalysis <- 1
      invalidateLater(1000)
      return()
    }
    reactive_vals$delayCounterStartAnalysis <- NULL
    
    reactive_vals$status <- PIPELINE_STATUS$STARTED
  })
  
  observeEvent(reactive_vals$status, ignoreInit = T, ignoreNULL = T, handlerExpr = {
    if (reactive_vals$status != PIPELINE_STATUS$STARTED)
      return()
    reactive_vals$status <- PIPELINE_STATUS$DOWNLOAD_DATASETS
  })
  
  observeEvent(reactive_vals$status, ignoreInit = T, ignoreNULL = T, handlerExpr = {
    if (reactive_vals$status != PIPELINE_STATUS$DOWNLOAD_DATASETS)
      return()
    
    withProgress({
      library(GEOquery)
      sel_datasets <- reactive_vals$selected_datasets
      validate(need(!is.null(sel_datasets) && nrow(sel_datasets) > 0, 'Please select at least 1 data set'))
      
      setProgress(message='Downloading data sets', value=0.0)
      i <- 0
      for(dataset_id in sel_datasets$accession) {
        setProgress(value=i/length(sel_datasets$accession))
        poolWithTransaction(pool, function(conn) { ### begin DB transaction
          print(paste("transaction begin", dataset_id))
          
          # get exclusive access to the requested data set
          repeat {
            print(paste("requesting exclusive access to dataset", dataset_id, "in DB"))
            rv <- try(
              dbExecute(
                conn,
                "INSERT INTO dataset_locks (dataset_id, session_uuid) VALUES ($1, $2);",
                params = list( dataset_id, reactive_vals$session_id))
            )
            if(!is(rv, "try-error")) break
            print("failed. retrying in 5secs ...")
            Sys.sleep(5)
          }
          
          # check whether the data set needs downloading / installing
          dsId <- dataset_id
          dataset_status <- dbGetQuery(pool, 
                                       "SELECT status FROM datasets WHERE dataset_id = $1;", 
                                       param=list(dsId))$status
          
          if (length(dataset_status) == 0)
            dataset_status = ""
          
          #print("going to sleep ...")
          #Sys.sleep(300)
          
          is_dataset_available <- dataset_status == "available"
          is_dataset_downloading <- dataset_status == "downloading"
          is_dataset_failed <- dataset_status == "failed"
          
          is_dataset_to_be_downloaded <- !(is_dataset_available || is_dataset_downloading || is_dataset_failed)
          
          #TODO: it seems like sometimes it says dataset available although the file is not in hkg-data
          if (is_dataset_to_be_downloaded) {
            setProgress(detail=paste('Downloading data set', dataset_id), value=0.0)
            print(paste("Downloading data set from GEO: ", dataset_id))
            
            # update status in database
            affectedNRows <- dbExecute(
              conn, 
              "UPDATE datasets SET status=$2 WHERE dataset_id = $1;",
              params = list(dataset_id, "downloading"))
            if (affectedNRows == 0) {
              dbWriteTable(conn, 
                           "datasets",
                           data.frame(dataset_id=dataset_id, status="downloading"),
                           append=T)
            }
            
            print(paste("Downloading ", dataset_id))
            error_thrown <- FALSE
            tryCatch({
              #gset_raw <- getGEO(GEO=dataset_id, GSEMatrix=TRUE, AnnotGPL=T, destdir = hkgDataDir)
              res <- downloadFromGEO(gse=dataset_id, GSEMatrix=TRUE, AnnotGPL=T, destdir = hkgDataDir)
              gset_raw <- res$gset_raw
              if (res$timed_out)
                error_thrown <- TRUE
              else {
                print("success")
                dbClearResult(
                  dbSendStatement(conn, 'update datasets set status=$1 where dataset_id=$2;', params=list("available", dataset_id)))
                names(gset_raw) <- lapply(gset_raw, function(x) x@annotation)
                reactive_vals$gsets[[dataset_id]] <- gset_raw
              }
            }, error=function(e) {
              print(e)
              error_thrown <- TRUE
            })
            if (error_thrown) {
              dbClearResult(
                dbSendStatement(conn, 'update datasets set status=$1 where dataset_id=$2;', params=list("failed", dataset_id))
              )
            }
          } else {
            setProgress(detail=paste('Parsing data set already present locally:', dataset_id), value=0.0)
            print("Data set not to be downloaded.")
            tryCatch({
              fns <- list.files(path = hkgDataDir, 
                                pattern = paste0("^", dataset_id,".*_series_matrix.txt.gz$"),
                                full.names = T)
              gsets <- list()
              for (fn in fns) {
                gpl <- stri_extract(fn, regex="GPL\\d+")
                gset_raw <- getGEO(filename = fn, 
                                   GSEMatrix=TRUE, AnnotGPL=T, getGPL = T, destdir = hkgDataDir)
                if (is.na(gpl))
                  gpl <- gset_raw@annotation
                gsets[[gpl]] <- gset_raw
              }
              reactive_vals$gsets[[dataset_id]] <- gsets
            }, error=function(e) {
              print(e)
            })
          }
          
          dbExecute(
            conn, 
            "DELETE FROM dataset_locks WHERE dataset_id = $1;",
            params = list(dataset_id))
          print(paste("transaction end", dataset_id))
        }) ### end DB transaction
        i <- i+1
      }
      setProgress(value=1.0)
      reactive_vals$status <- PIPELINE_STATUS$DOWNLOAD_ANNOTATIONS
    })
  })
  
  observe({
    if (reactive_vals$status != PIPELINE_STATUS$DOWNLOAD_ANNOTATIONS)
      return()
    
    withProgress({
      setProgress(message='Downloading annotation data', detail='Downloading annotation data for selected organisms', value=0.0)
      
      orgs <- unique(unlist(strsplit(isolate(reactive_vals$selected_datasets)$taxon, "; ")))
      orgs <- setNames(as.list(orgs), tolower(orgs))
      
      reactive_vals$orgAnno <- lapply(
        1:length(orgs),
        FUN=function(i) {
          o <- orgs[[i]]
          setProgress(message='Downloading annotation data', detail=sprintf('Downloading annotation data for %s', o), value=i/length(orgs))
          queryResult <- query(x=query(ah, "EnsDb"), pattern=o)
          ensembl_release <- isolate(reactive_vals$ensembl_release)
          print(paste(o, ": Using Ensembl", ensembl_release ,"EnsDb"))
          queryResult[[grep(paste("Ensembl",ensembl_release,"EnsDb"), mcols(queryResult)$title, fixed=T)]]
        }
      )
      names(reactive_vals$orgAnno) <- names(orgs)
      
      reactive_vals$orgMart <- lapply(
        1:length(orgs),
        FUN=function(i) {
          org <- orgs[[i]]
          setProgress(message='Downloading annotation data', detail=sprintf('Downloading ensembl mart for %s', org), value=i/length(orgs))
          orgAbbr <- strsplit(tolower(org), " ")[[1]]
          orgAbbr[1] <- substr(orgAbbr[1],1,1)
          orgAbbr <- paste0(c(orgAbbr, "_gene_ensembl"), collapse="")
          
          ensembl_release <- isolate(reactive_vals$ensembl_release)
          x <- 1
          max_x <- 5
          mart <- NULL
          while (is.numeric(x)) {
            
            x <- tryCatch({
              if (is.null(ensembl_release)) {
                # take the most recent one for this organism
                ver <- max(as.numeric(listEnsemblArchives()$version), na.rm=T)
                print(paste(org, ": Using ensembl mart", ver))
                mart <- useEnsembl("ensembl", dataset = orgAbbr, mirror="useast")
              }
              else {
                print(paste(org, ": Using ensembl mart", ensembl_release))
                mart <- useEnsembl("ensembl", dataset = orgAbbr, version = ensembl_release)
              }
              mart
            }, error = function(e) {
              if (x > max_x) {
                print(paste('The analysis failed:', e$message))
                setProgress(message = 'The analysis failed', detail =  e$message)
                
                #showNotification('The analysis failed', type = "error", duration = NULL)
                reactive_vals$analysis_error <- sprintf('Using Ensembl mart failed %d times. Error: %s', max_x, e$message)
                showNotification(sprintf("%s: Using Ensembl mart failed %d times", org, max_x), type = "error", duration = NULL)
                reactive_vals$analysis_started <- FALSE
                reactive_vals$analysis_running <- FALSE
                
                
                reactive_vals$session_id <- NULL
                
                # remove session id from URL
                session$doBookmark()
              } else {
                showNotification(sprintf("%s Attempt %d: Using Ensembl mart failed. Retrying in 2 secs... ", org, x), type = "error", duration = 30)
                Sys.sleep(2)
                x + 1
              }
            })
          }
          mart
        })
      names(reactive_vals$orgMart) <- names(orgs)
      
      setProgress(value=1.0)
      
      reactive_vals$status <- PIPELINE_STATUS$FIND_HOUSEKEEPING_GENES
      if (any(sapply(reactive_vals$orgMart, is.null))){
        reactive_vals$status <- PIPELINE_STATUS$IDLE
      }
    })
  })
  
  observeEvent(reactive_vals$status, {
    if (reactive_vals$status != PIPELINE_STATUS$FIND_HOUSEKEEPING_GENES)
      return()
    
    withProgress({
      sel_datasets <- reactive_vals$selected_datasets
      sel_samples <- reactive_vals$selected_samples
      validate(need(length(sel_samples) > 0, 'Please select at least 1 sample'))
      
      annotated_samples <- list()
      for (gse_GPL in unique(reactive_vals$available_samples$gse_gpl)) {
        #gse <- ds_available_samples$accession
        samples <- list()
        for (gsm in reactive_vals$available_samples[gse_gpl==gse_GPL]$samples.accession) {
          if (gsm %in% sel_samples[[gse_GPL]])
            samples[[gsm]] <- "1"
          else
            samples[[gsm]] <- "0"
        }
        annotated_samples[[gse_GPL]] <- samples
      }
      
      tryCatch({
        res <- extract_hkg(num_hkg=reactive_vals$numberGenes,
                           GSE_list = as.character(sel_datasets$gse_gpl),
                           GSE_samples_annotations = annotated_samples,
                           bootstrap_sample_size = reactive_vals$bootstrap_sample_size,
                           bootstrap_replications = reactive_vals$bootstrap_replications,
                           hkgDataDir = hkgDataDir)
        if(is.character(res))
          stop(res)
        bootstrap_ranking <- res$bootstrap_ranking[, c("paralog_group", "symbols", "gene_ids", "num_na_factor", "rank_mean", "rank_var", "total_rank"), with=FALSE]
        setnames(bootstrap_ranking, c("num_na_factor","total_rank"), c("missing_datasets", "rank"))
        setcolorder(bootstrap_ranking, c("paralog_group", "symbols", "gene_ids", "missing_datasets", "rank_mean", "rank_var", "rank"))
        
        # ensure unique symbols by appending ensembl ids
        symbolIsDuplicated <- duplicated(bootstrap_ranking$symbols) | duplicated(bootstrap_ranking$symbols, fromLast=T)
        bootstrap_ranking[symbolIsDuplicated, symbols :=
                            paste0(symbols, "###(", gene_ids, ")")]
        
        mappedUniqueSymbols <- data.table(paralog_group=bootstrap_ranking$paralog_group[symbolIsDuplicated],
                                          newSymbols=bootstrap_ranking$symbols[symbolIsDuplicated])
        
        orderedGeneIds <- bootstrap_ranking$symbols[order(bootstrap_ranking$rank)]
        
        bootstrap_ranking[, symbols := factor(symbols, levels=orderedGeneIds)]
        
        bootstrap_ranking_long <- res$bootstrap_ranking_long[, c("paralog_group", "symbols", "gene_ids", "gene_rank", "bootstrap_sample"), with=FALSE]
        setnames(bootstrap_ranking_long, "gene_rank", "rank")
        setcolorder(bootstrap_ranking_long, c("paralog_group", "symbols", "gene_ids", "bootstrap_sample", "rank"))
        
        # ensure unique symbols by appending ensembl ids
        bootstrap_ranking_long <- merge(bootstrap_ranking_long, mappedUniqueSymbols, all.x=T)
        bootstrap_ranking_long[!is.na(newSymbols), symbols:=newSymbols]
        bootstrap_ranking_long[, symbols := factor(symbols, levels=orderedGeneIds)]
        bootstrap_ranking_long[, newSymbols := NULL]
        
        # ensure unique symbols by appending ensembl ids
        res$ranking_long <- merge(res$ranking_long, mappedUniqueSymbols, all.x=T)
        res$ranking_long[!is.na(newSymbols), symbols:=newSymbols]
        res$ranking_long[, symbols := factor(symbols, levels=orderedGeneIds)]
        res$ranking_long[, newSymbols := NULL]
        setnames(res$ranking_long,c("paralog_group", "rank", "dataset", "symbols", "gene_ids"))
        setcolorder(res$ranking_long, c("paralog_group", "symbols", "gene_ids", "dataset", "rank"))
        
        reactive_vals$analysis_started <- FALSE
        reactive_vals$analysis_running <- FALSE
        reactive_vals$status <- PIPELINE_STATUS$FINISHED
        reactive_vals$bootstrappingRankTableData <- bootstrap_ranking
        reactive_vals$bootstrappingLongRankTableData <- bootstrap_ranking_long
        reactive_vals$rankLongTableData <- res$ranking_long
        reactive_vals$analysis_error <- NULL
        
        # persist results in DB
        dbWriteTable(pool, "gene_bootstrap_total_rank_scores", cbind(session_uuid=reactive_vals$session_id, bootstrap_ranking), append=T)
        dbWriteTable(pool, "gene_bootstrap_sample_rank_scores", cbind(session_uuid=reactive_vals$session_id, bootstrap_ranking_long), append=T)
        dbWriteTable(pool, "gene_dataset_ranks", cbind(session_uuid=reactive_vals$session_id, res$ranking_long), append=T)
        
        dbWriteTable(pool, "analysis_parameters", cbind(session_uuid=reactive_vals$session_id, 
                                                        data.frame(
                                                          tissue_type=paste(reactive_vals$tissues, collapse=";"),
                                                          condition=paste(reactive_vals$conditions, collapse=";"),
                                                          organism=paste(reactive_vals$organisms, collapse=";"),
                                                          num_datasets=reactive_vals$num_selected_datasets,
                                                          number_hkg=reactive_vals$numberGenes,
                                                          bootstrap_sample_size=reactive_vals$bootstrap_sample_size,
                                                          bootstrap_replications=reactive_vals$bootstrap_replications,
                                                          ensembl_release=(
                                                            if(is.null(reactive_vals$ensembl_release)) 
                                                              NA 
                                                            else 
                                                              reactive_vals$ensembl_release),
                                                          seed=reactive_vals$seed)
        ), append=T)
        
        conn <- poolCheckout(pool);
        affectedNRows <- dbExecute(conn, 
                                   if (houseekpr_env == "development") "UPDATE sessions SET finish_ts=$1 WHERE session_uuid = $2;" 
                                   else "UPDATE sessions SET finish_ts=to_timestamp($1) WHERE session_uuid = $2;",
                                   params = list(Sys.time(), reactive_vals$session_id))
        if (affectedNRows == 0) {
          print("Failed to update session finish timestamp in DB")
        }
        poolReturn(conn)
        
        print(paste("Analysis finished (session ", reactive_vals$session_id, ")"))
        
        # make sure the session id is persisted in the URL
        session$doBookmark()
        # workaround for production setup, where we need to make sure 
        # that we send a new request to nginx with the 
        # session id inthe query such that we are forwarded to the presentation app instances
        shinyjs::js$refresh()
      },
      error=function(e) {
        
        print(paste('The analysis failed:', e$message))
        setProgress(message = 'The analysis failed', detail =  e$message)
        
        reactive_vals$analysis_error <- e$message
        reactive_vals$analysis_started <- FALSE
        reactive_vals$analysis_running <- FALSE
        
        reactive_vals$session_id <- NULL
        
        # manually store tissue, organism and condition:
        reactive_vals$restoreState <- list(
          input = list(
            select_tissues_all = isolate(input$select_tissues_all),
            select_conditions_all = isolate(input$select_conditions_all),
            select_organisms_all = isolate(input$select_organisms_all)
          ),
          values = list(
              ds = sprintf("%s.%s",isolate(reactive_vals$selected_datasets$accession), isolate(reactive_vals$selected_datasets$gpl)),
              samples = isolate(reactive_vals$selected_samples)
          )
        )
        
        
        # remove session id from URL
        session$doBookmark()
      })
    })
    # reset seed value
    set.seed(seed=NULL)
  }, ignoreInit = TRUE)
  
  observeEvent(reactive_vals$available_ensembl_releases, {
    updateSelectizeInput(session, "customEnsemblRelease", choices=reactive_vals$available_ensembl_releases)
  })
  
  # stuff depending on selected tab
  observeEvent(input$menu, {
    if (input$menu == "start") {
      if (!is.null(reactive_vals$session_id))
        reactive_vals$session_id <- NULL
      
      reactive_vals$bootstrappingRankTableData <- data.table(
        paralog_group=integer(0), 
        symbols=character(0),
        gene_ids=character(0),
        missing_datasets=integer(0),
        rank_mean=double(0), 
        rank_var=double(0), 
        rank=double(0))
      reactive_vals$bootstrappingLongRankTableData <- data.table()
      reactive_vals$rankLongTableData <- data.table(
        dataset=character(0), 
        paralog_group=integer(0), 
        symbols=character(0),
        gene_ids=character(0),
        rank=double(0))
      
      reactive_vals$available_ensembl_releases <- str_sort(unique(gsub(".*?Ensembl version ([0-9]+).*", "\\1", unlist(mcols(query(ah, "EnsDb"))$description))), numeric = T, decreasing = T)
      
      # does not work if annotationhub does not have current ensembl version
      tryCatch({
        archives <- listEnsemblArchives()
        na.omit(as.numeric(data.table(archives)[grep("^\\d+$", version), version]))
      }, error = function(e){
          reactive_vals$analysis_error <- sprintf("Please try HousKeepR again later. We cannot reach Ensembl: %s", e$message)
      })
      
      gc()
    }
  })
  
  observe({
    if (getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished"))
      return(NULL)
    
    reactive_vals$ensembl_release <- input$customEnsemblRelease
  })
  
  observe({
    if (getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished"))
      return(NULL)
    reactive_vals$seed <- input$seed
  })
  
  # observe({
  #   if (getSessionStatus(pool, reactive_vals$session_id) %in% c("started", "finished"))
  #     return(NULL)
  #   
  #   reactive_vals$ensembl_release
  #   session$doBookmark()
  # })
  
  # deactivate inputs in the sidebar while busy with processing
  # start_analysis button is handled separately below
  observe({
    shouldBeDisabled <- reactive_vals$analysis_running || reactive_vals$queryingGEOdatasets || reactive_vals$queryingGEOsamples
    
    if (shouldBeDisabled) {
      #disable("use_example_data")
      disable("select_tissues_all")
      disable("select_organisms_all")
      disable("select_conditions_all")
      # disable("algorithm")
      disable("numberGenes")
      disable("bootstrapSampleSize")
      disable("bootstrapRepetitions")
      disable("customEnsemblRelease")
      disable("seed")
      disable("old_session")
    } else {
      #enable("use_example_data")
      enable("select_tissues_all")
      enable("select_organisms_all")
      enable("select_conditions_all")
      # enable("algorithm")
      enable("numberGenes")
      enable("bootstrapSampleSize")
      enable("bootstrapRepetitions")
      enable("customEnsemblRelease")
      enable("seed")
      enable("old_session")
    }
    
  })
  
  # deactivate start analysis button if appropriate
  observe({
    # trigger this one, when the button is initialized (necessary, because we use renderUI)
    input$start_analysis
    shouldBeDisabled <- reactive_vals$analysis_running || reactive_vals$queryingGEOdatasets || reactive_vals$queryingGEOsamples ||
      (reactive_vals$num_selected_datasets < 2) ||
      (!reactive_vals$valid_sample_selection);
    
    if (shouldBeDisabled) {
      disable("start_analysis")
    } else {
      enable("start_analysis")
    }
  })
  
  # hide and show sidebar when appropriate
  observe({
    reactive_vals$analysis_running
    reactive_vals$session_id
    
    if (getSessionStatus(pool, reactive_vals$session_id) %in% c("started"))
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    else
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
  })
  
  observeEvent(input$rankPlotlyHideOutliers, {
    reactive_vals$rankPlotlyHideOutliers <- input$rankPlotlyHideOutliers
    session$doBookmark()
  })
  
  observeEvent(input$plotsNumHKG, {
    reactive_vals$plotsNumHKG <- input$plotsNumHKG
    session$doBookmark()
  })
  
  ### /OBSERVER
  
  ### DOWNLOAD HANDLER ----
  
  output$downloadRankingTable <- downloadHandler(
    filename=function() {
      paste0(isolate(reactive_vals$session_id), "_ranking-matrix", ".zip")
    },
    content=function(file) {
      tmpdir <- tempdir()
      cwd <- getwd()
      setwd(tempdir())
      tsv <- paste0(isolate(reactive_vals$session_id), "_ranking-matrix", ".tsv")
      write("# HouseKeepR ranking matrix", tsv)
      write(paste("# Session", reactive_vals$session_id), tsv, append=T)
      write(paste("# Analysis start time", getSessionStarttime(pool, reactive_vals$session_id)), tsv, append=T)
      write(paste("# Analysis finish time", getSessionFinishtime(pool, reactive_vals$session_id)), tsv, append=T)
      write(paste("# File downloaded", Sys.time()), tsv, append=T)
      write.table(reactive_vals$bootstrappingRankTableData[order(rank),-1], tsv, row.names=F,append=T,sep="\t")
      res <- zip(zipfile=file, files=tsv)
      setwd(cwd)
      res
    },
    contentType = "application/zip"
  )
  
  output$downloadDatasetRanks <- downloadHandler(
    filename=function() {
      paste0(isolate(reactive_vals$session_id), "_dataset-ranks", ".zip")
    },
    content=function(file) {
      tmpdir <- tempdir()
      cwd <- getwd()
      setwd(tempdir())
      tsv <- paste0(isolate(reactive_vals$session_id), "_dataset-ranks", ".tsv")
      write("# HouseKeepR dataset ranks", tsv)
      write(paste("# Session", reactive_vals$session_id), tsv, append=T)
      write(paste("# Analysis start time", getSessionStarttime(pool, reactive_vals$session_id)), tsv, append=T)
      write(paste("# Analysis finish time", getSessionFinishtime(pool, reactive_vals$session_id)), tsv, append=T)
      write(paste("# File downloaded", Sys.time()), tsv, append=T)
      write.table(reactive_vals$rankLongTableData[,-1], tsv, row.names=F,append=T,sep="\t")
      res <- zip(zipfile=file, files=tsv)
      setwd(cwd)
      res
    },
    contentType = "application/zip"
  )
  
  output$downloadBootstrapRanks <- downloadHandler(
    filename=function() {
      paste0(isolate(reactive_vals$session_id), "_bootstrap-ranks", ".zip")
    },
    content=function(file) {
      tmpdir <- tempdir()
      cwd <- getwd()
      setwd(tempdir())
      tsv <- paste0(isolate(reactive_vals$session_id), "_bootstrap-ranks", ".tsv")
      write("# HouseKeepR bootstrap ranks", tsv)
      write(paste("# Session", reactive_vals$session_id), tsv, append=T)
      write(paste("# Analysis start time", getSessionStarttime(pool, reactive_vals$session_id)), tsv, append=T)
      write(paste("# Analysis finish time", getSessionFinishtime(pool, reactive_vals$session_id)), tsv, append=T)
      write(paste("# File downloaded", Sys.time()), tsv, append=T)
      write.table(reactive_vals$bootstrappingLongRankTableData[,-1], tsv, row.names=F,append=T,sep="\t")
      res <- zip(zipfile=file, files=tsv)
      setwd(cwd)
      res
    },
    contentType = "application/zip"
  )
  
  outputOptions(output, "queryingGEOdatasets", suspendWhenHidden=T)
  outputOptions(output, "queryingGEOsamples", suspendWhenHidden=T)
  outputOptions(output, "menu", priority=110, suspendWhenHidden=FALSE)
  outputOptions(output, "uiSidebar", priority=100)
}

