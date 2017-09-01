#.libPaths(c( "/home/aleferna/R/x86_64-pc-linux-gnu-library/3.2/",.libPaths())  )
#install.packages("rhandsontable")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
#biocLite("pcaMethods")

#install.packages("htmlwidgets")
#install.packages("rpivotTable")
#install.packages("d3heatmap")
#install.packages("plotly")
#install.packages("Rmisc")
library(DT)
library(shiny)
library(stringr)
library(plotly)
library(ggplot2)
require(shinysky)


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


load("DATA/DT45.Rdata")


uiScatterPlot = function() {
  sidebarLayout(
    sidebarPanel(width=3,
                 #radioButtons("inpPlotType", "Plot Type", c("Protein_vs_mRNA", "Protein", "mRNA", "Protein_vs_CNA", "mRNA_vs_CNA"), selected = "Protein_vs_mRNA", inline = F, width = NULL),
                 fluidRow(
                   column(6, selectInput("cmbFieldTypeX", "X Axis Field",choices = c("Protein", "mRNA","CNA","miRNA", "Metabolites"), selected = "mRNA",  multiple=F)
                   ),
                   column(6, selectInput("cmbFieldX"    , "X Axis Field",choices =c(), selected = "",  multiple=F)
                   )
                 ),
                 fluidRow(
                   column(6,selectInput("cmbFieldTypeY", "Y Axis Field",choices = c("Protein", "mRNA","CNA","miRNA", "Metabolites"), selected = "Protein",  multiple=F)
                   ),
                   column(6,selectInput("cmbFieldY", "Y Axis Field",choices = c(), selected = "",  multiple=F )
                   )
                 ),
                 
                 selectInput("cmbColorBy",label="Color Field:",choices=c("PAM50","CoTC","ER","PR","HER.2","Metabolic.cluster", "CAAI", "RPPA.subtype"  ), selected = "PAM50" ),
                 dataTableOutput("HighCor"),
                 downloadButton('downloadPlot', 'Publication Ready Figure')
    ),
    mainPanel(
      busyIndicator2(text = "Loading, please wait...",img = "busyIndicator.gif", wait=1000),
      fluidRow(
             plotlyOutput("XYPlot", height="600px")
      )
    )
  )
}



shinyUI(fluidPage(theme="bootstrap.min.css",
                  fluidRow(
                    column(6, align="center", offset = 3,
                           titlePanel("Breast Cancer Landscape"),
                           tags$style(type='text/css', "#button { vertical-align: middle; height: 50px; width: 100%; font-size: 30px;}")
                    )
                  ),
                  tabPanel("Scatter", uiScatterPlot() )
                  
)
)
