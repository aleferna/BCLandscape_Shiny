#.libPaths(c( "/home/aleferna/R/x86_64-pc-linux-gnu-library/3.2/",.libPaths())  )
#install.packages("rhandsontable")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
#biocLite("pcaMethods")

#install.packages("htmlwidgets")
#install.packages("rpivotTable")
#install.packages("d3heatmap")
#install.packages("plotly")
#install.packages("mixtools")
#citation("mixtools")
#require(mixtools)
library(shiny)
library(stringr)
library(plotly)
library(ggplot2)
library(ComplexHeatmap)

load("DATA/geneListsHJ.RData")
lstGeneListNames = sort(colnames(geneLists)[2:ncol(geneLists)])

#load("DATA/miRNA.Rdata")
#lstmiRNA=dtmiRNA$geneName


busyIndicator2 <- function(text = "Calculation in progress..",img = "busyIndicator.gif", wait=1000) {
  tagList(
    singleton(tags$head(
      tags$link(rel="stylesheet", type="text/css",href="busyIndicator.css")
    ))
    ,div(class="shinysky-busy-indicator",p(text),img(src=img))
    ,tags$script(sprintf(
      "	setInterval(function(){
  		 	 if ($('html').hasClass('shiny-busy')) {
  		    setTimeout(function() {
  		      if ($('html').hasClass('shiny-busy')) {
  		        $('div.shinysky-busy-indicator').show()
  		      }
  		    }, %d)
  		  } else {
  		    $('div.shinysky-busy-indicator').hide()
  		  }
  		},100)
  		",wait)
    )
  )
}





shinyUI(fluidPage(theme="bootstrap.min.css",
                  busyIndicator2(text = "Loading, please wait...",img = "busyIndicator.gif", wait=1000),
                  fluidRow(
                    column(6, align="center", offset = 3,
                           titlePanel("Breast Cancer Landscape"),
                           tags$style(type='text/css', "#button { vertical-align: middle; height: 50px; width: 100%; font-size: 30px;}")
                    )
                  ),
                  sidebarLayout(
                    sidebarPanel(width=3,
                                 radioButtons("inpPlotType",label = "Select Plot Type:",choices = c("Protein","mRNA"),selected="Protein", inline = T),
                                 conditionalPanel("input.inpPlotType != \"miRNA\"",
                                                  selectInput("inpGeneLists","Select Gene List",choices = lstGeneListNames, selected = "AFW_Basal"),
                                                  tags$textarea(id = 'lstGenes', placeholder = 'Gene List', c() ,  style = "width: 100%; height: 300px")
                                 ),
                                 #conditionalPanel("input.inpPlotType == \"miRNA\"",
                                 #                selectInput("inpmiRNA","Select miRNA",choices =lstmiRNA  , selected = sample(lstmiRNA,100), multiple=T, selectize = F)
                                 
                                 #),
                                 downloadButton('downloadPlot', 'Save as pdf')
                    ),
                    mainPanel(
                      shiny::plotOutput("pltHeatmap",height = "800px")
                      
                      #conditionalPanel("input.inpPlotType == \"Protein\"",
                      #                 headerPanel("Protein"),
                      #                 d3heatmapOutput("htProt", height="800px")
                      #),
                      #conditionalPanel("input.inpPlotType == \"mRNA\"",
                      #                 headerPanel("mRNA"),
                      #                 d3heatmapOutput("htmRNA", height="800px")
                      #),
                      #conditionalPanel("input.inpPlotType == \"miRNA\"",
                      #                 headerPanel("miRNA"),
                      #                 d3heatmapOutput("htmiRNA", height="800px")
                      #),
                      #conditionalPanel("input.inpPlotType == \"Reference\"",
                      #                 headerPanel("Reference Plot"),
                      #                 plotOutput("htref")
                      #)
                    )
                    
                  )                  
)
)
