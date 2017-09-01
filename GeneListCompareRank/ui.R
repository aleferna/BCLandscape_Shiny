library(plotly)
library(shiny)
library(DT)
library(stringr)
require(GSA) 

#gmts=str_replace(list.files("../DATA/GeneLists",pattern = ".gmt"),".gmt","")
load("gmtGeneLists.RData")


shinyUI(fluidPage(theme="bootstrap.min.css",
                  titlePanel("Gene Enrichment Rank"),
                  sidebarLayout(
                    sidebarPanel(width=3,
                                 fileInput('fileGeneMtx', 'Choose gene matrix file', accept = c('.tab') ),
                                 selectInput("inpEnrichmentCategory","Select Enrichment Search:", choices=names(gmtGeneLists), selected = "MSigDB_Hallmarks" ),
                                 downloadButton('downloadWB', 'Save as Excel')
                    ),
                    mainPanel(
                      fluidRow(
                      plotlyOutput("volcSelected")
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
