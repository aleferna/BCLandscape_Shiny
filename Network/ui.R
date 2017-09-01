
library(DT)
library(shiny)
library(stringr)
library(ggplot2)
library(plotly)
library(reshape)
library(ComplexHeatmap)
library(pcaMethods)
require(Rmisc)
library(shiny)
#load("../DATA/dtPos.Small.Rdata")
#dt=dtPosSmall

load("DATA/geneListsHJ.RData")
#load("../DATA/dtGeneLevels.Rdata")
load("DATA/dtGeneCategories.Rdata")
load("DATA/gmtGeneLists.RData")


# Define UI for application that draws a histogram

shinyUI(fluidPage(theme="bootstrap.min.css",
                  
                  
                  
                  # Application title
                  titlePanel("Breast Cancer Landscape Gene Correlation Network"),
                  
                  # Sidebar with a slider input for number of bins 
                  sidebarLayout(
                    sidebarPanel(width=3,
                                 #radioButtons("inpPlotType", "Color by", c("Gene List","Gene - Category", "Protein", "mRNA"), selected = "Gene List", inline = T, width = NULL),
                                 radioButtons("inpNetSource","Network source:",choices=c("RNA","Protein"), selected = "Protein"  ),
                                 radioButtons("inpPlotType", "Color by", c("Gene List","Gene Category","Gene Level", "Gene Neighborhood"), selected = "Gene Category", inline = F),
                                 checkboxInput("inpShowEdges","Show edges"),
                                 conditionalPanel(condition="input.inpPlotType == \"Gene Neighborhood\"",      
                                                  selectInput("inpSelectedGene","Select Gene", choices=c() )
                                 ),
                                 
                                 conditionalPanel(condition="input.inpPlotType == \"Gene List\" || input.inpPlotType == \"Gene Neighborhood\" ",
                                                  selectInput("inpEnrichmentCategory","Select Enrichment Search:", choices=names(gmtGeneLists), selected="MSigDB_Hallmarks" )
                                 ),
                                 
                                 conditionalPanel(condition="input.inpPlotType == \"Gene List\"",
                                                  selectInput("inpGeneLists","Select Gene List",choices = c() )
                                 ),
                                 
                                 conditionalPanel(condition="input.inpPlotType == \"Gene List\"  ",
                                                  tags$textarea(id = 'lstGenes', placeholder = 'Gene List', c() ,  style = "width: 100%; height: 300px")
                                 ),
                                 
                                 
                                 conditionalPanel(condition="input.inpPlotType == \"Gene Category\"",      
                                                  selectInput("inpGeneClassLists","Select Gene Classification ", choices=colnames(dtGeneCategories)[2:ncol(dtGeneCategories)] ) ,
                                                  tags$textarea(id = 'lstGeneClasses', placeholder = 'Gene List', c() ,  style = "width: 100%; height: 300px")
                                                  
                                 ),
                                 conditionalPanel(condition="input.inpPlotType == \"Gene Level\"",      
                                                  selectInput("inpGeneLevelLists","Select Gene Level ", choices=c() ) ,
                                                  tags$textarea(id = 'lstGeneLevels', placeholder = 'Gene List', c() ,  style = "width: 100%; height: 300px")
                                 ),
                                 downloadButton('downloadPlot', 'Save as pdf')
                                 
                                 
                                 
                    ),
                    
                    # Show a plot of the generated distribution
                    mainPanel(width=9,
                              plotOutput("NetPlot", width = "100%", height = "700px",click = "plot_click",dblclick = "plot_dblclick",
                                         brush = brushOpts(id = "plot_brush",resetOnNew = TRUE)
                              ),
                              # conditionalPanel(condition="input.inpPlotType == \"Gene Neighborhood\"",      
                              #                  plotOutput("barplot", height="200px")
                              # ),
                              DT::dataTableOutput("GeneNeighborhoodEnrich")
                              #visNetworkOutput("networkPlot", width = "800px", height = "600px")
                    )
                  )
))













