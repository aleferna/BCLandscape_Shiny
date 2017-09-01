
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)
library(d3heatmap)
library(stringr)

load("../DATA/RNA.Rdata")
load("../DATA/Prot.Rdata")
load("../DATA/geneListsHJ.RData")
geneLists[is.na(geneLists)] = F
geneLists$GeneName[ geneLists$FDATargets ]
gns=paste(geneLists$GeneName[ geneLists$FDATargets ],collapse="\n")
lstGeneListNames = colnames(geneLists)[2:ncol(geneLists)]



lstSamps=sort(colnames(dtProt)[2:ncol(dtProt)])
dtMeta= read.table("../DATA/meta.tab",sep="\t",header = T, stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
lstSamps = paste(dtMeta[lstSamps,"ProteinSubtype"],str_split_fixed( lstSamps, "\\.",2)[,2],sep=".")
lstSamps = sort(lstSamps)




shinyUI(fluidPage(theme="../www/bootstrap.min.css",

  # Application title
  titlePanel("Discordance Enrichment"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      
      
      selectInput("inpSample","Select sample",choices = lstSamps, selected = lstSamps[1]),
      selectInput("inpGeneLists","Load Gene List:",choices = lstGeneListNames, selected = "FDATargets"),
      selectInput("inpGeneListsBioGrid","or Load BioGrid interactors:",choices = "" , selected = ""),
      tags$textarea(id = 'lstGenes', placeholder = 'Gene List', gns ,  style = "width: 100%; height: 300px")
      
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("scatterPlot"),
      plotOutput("discRankPlot", height="150px"),
      d3heatmapOutput(outputId = "heatPlot", height="600px")
    )
    
    
    
    
  )
))
