library(plotly)
library(shiny)
library(DT)
library(stringr)
require(GSA) 

load("gmtGeneLists.RData")
initxt=paste(readLines("example.tab"), collapse="\n")

shinyUI(fluidPage(theme="bootstrap.min.css",
                  titlePanel("Gene Enrichment"),
                  sidebarLayout(
                    sidebarPanel(width=3,
                                 h4("Paste gene list matrix"),
                                 tags$textarea(id = 'lstGeneLists', placeholder = 'paste Gene List Matrix', initxt ,  style = "width: 100%; height: 300px"),
                                 #h4("Paste background gene list"),
                                 #tags$textarea(id = 'lstBGGeneLists', placeholder = 'Background Gene List',paste(allGenes,collapse="\n") ,  style = "width: 100%; height: 200px"),
                                 selectInput("inpEnrichmentCategory","Select Enrichment Search:", choices=names(gmtGeneLists), selected = "MSigDB_Hallmarks" ),
                                 fileInput('fileBGGeneList', 'Choose background gene list', accept = c('.lst') ),
                                 downloadButton('downloadWB', 'Save as Excel')
                    ),
                    mainPanel(
                      fluidRow(
                      plotlyOutput("volcSelected", height = "600px")
                      ),
                      fluidRow(
                        titlePanel("List Category Enrichment"),
                        dataTableOutput("DE_Enrich")
                      ),
                      fluidRow(
                        titlePanel("Selected Enrichment Detail"),
                        dataTableOutput("DE_Params")
                      ),
                      fluidRow(
                        titlePanel("Genes in Selected List"),
                        dataTableOutput("DE_EnrichSel")
                      )
                    )
                  )
))
