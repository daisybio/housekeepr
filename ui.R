source("modal.R", local = TRUE)
library(waiter)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)

ui <- function(request) {
  tagList(
    use_waiter(), # include dependencies
    waiter_show_on_load(tagList(
      spin_fading_circles(),
      tags$div("All required packages and essentials to present the HouseKeepR web application are currently being loaded.",
               tags$br(),
               "This may take a minute. Please be patient...")
    )),
    dashboardPage(title="HouseKeepR",
                  dashboardHeader(
                    title = 
                      tags$div(tags$span(tags$img(src = "HouseKeepR_cropped.png", style = "height: 40px; width: auto;")),
                               tags$span(paste0(HOUSEKEEPR_VER), style = "color: #ddd; font-size: 70%; vertical-align: sub;")
                               )
                  ),
                  dashboardSidebar(
                    sidebarMenuOutput("menu") %>% withSpinner(),
                    uiOutput("uiSidebar") %>% withSpinner()
                  ),
                  dashboardBody(
                    tags$head(
                      tags$link(rel = "stylesheet", type = "text/css", href = "housekeepr.css"),
                      tags$style(".skin-blue .sidebar  a.shiny-download-link { color: #444; }"),
                      HTML("<link rel='shortcut icon' href='favicon.ico' type='image/x-icon' />"),
                      HTML('<link rel="stylesheet" type="text/css" href="cookieconsent.min.css"/><script src="cookieconsent.min.js"></script><script>window.addEventListener("load", function(){window.wpcc.init({"colors":{"popup":{"background":"#cff5ff","text":"#000000","border":"#5e99c2"},"button":{"background":"#5e99c2","text":"#ffffff"}}, "padding":"none","margin":"none","fontsize":"tiny","content":{"href":"https://www.learn-about-cookies.com/"},"position":"top-right"})});</script>')
                    ),
                    useShinyjs(),
                    shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }", functions = c("refresh")),
                    withMathJax(),
                    modal_bootstrapping,
                    modal_ensembl,
                    modal_seed,
                    tabItems(
                      tabItem(tabName = "start",
                              uiOutput("uiAnalysisError"),
                              fluidRow(
                                box(width=12, status="primary",
                                    title="Set up a new HouseKeepR analysis",
                                    solidHeader=T,
                                    p(class="lead", 
                                      "HouseKeepR identifies promising candidate house-keeping genes from gene expression data by array of case and control samples of a set of GEO data sets."),
                                    uiOutput("uiStartAnalysis") %>% withSpinner(proxy.height="200px")
                                )
                              ),
                              fluidRow(
                                box(title="Data set selection", width = 6, collapsible = T,
                                    uiOutput("dataSetSelection") %>% withSpinner(proxy.height="50px")
                                ),
                                box(title="Condition & control samples", width = 6,
                                    collapsible = T,
                                    uiOutput("sampleSelection") %>% withSpinner(proxy.height="50px")
                                )
                              )
                      ),
                      tabItem(tabName = "status",
                              uiOutput("uiAnalysisRunning")
                      ),
                      tabItem(tabName = "results",
                              fluidRow(
                                box(status="primary",
                                    title="Best House-keeping Gene Candidates", 
                                    width = 12,
                                    solidHeader = T,
                                    collapsible = T,
                                    p(class="lead", paste("The top house-keeping gene candidates ranked by total rank score. 
                                              The total rank score is an aggregated score of a gene's ranks over all bootstrap samples.")),
                                    "Select one or multiple genes to update the ranking charts",
                                    DT::dataTableOutput('houseKeepingGenesTable') %>% withSpinner(),
                                    downloadButton("downloadRankingTable", "Download Ranking Table")
                                )
                              ),
                              fluidRow(
                                box(title="Ranks of Candidates on Bootstrap Samples",
                                    collapsible = T,
                                    uiOutput("uiRankPlotly") %>% withSpinner(proxy.height="50px")
                                ),
                                
                                box(title="Ranks of Candidates on Selected Datasets",
                                    collapsible = T,
                                    uiOutput("uiRankLongPlotly") %>% withSpinner(proxy.height="50px")
                                )
                              )
                      ),
                      uiOutput("boxes")
                    )
                  )
    )
  )
}
