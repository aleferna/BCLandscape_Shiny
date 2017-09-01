library(shiny)
library(ggplot2)
library(reshape)
library(plotly)
library(DT)
library(ggrepel)
library(stringr)
require(shinysky)
library(shinyjs)

dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
dtMeta = dtMeta[dtMeta$ProteinSubtype != "",]
dtMeta$TMT.tag = NULL
dtMeta$X = NULL



xcols=c("PAM50","CoTC","ER","PR","HER.2","Metabolic.cluster", "CAAI", "RPPA.subtype"  )

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
                  titlePanel("BC Landscape BoxPlot"),
                  sidebarLayout(
                    sidebarPanel(
                      selectInput("cmbColorBy",label="Group by:",choices=xcols, selected = "subtype" ),
                      tags$div(title="Note: Type gene name slowly",
                               selectInput("lstSelectedGene",label="Select gene:",choices=c("MET"), selected = "MET" , selectize = T)
                      ),
                      downloadButton('downloadPlot', 'Publication Ready Figure')
                      
                    ),
                    mainPanel(
                      plotlyOutput("BoxPlot",height = "600px"),
                      DT::dataTableOutput("tblDiff")
                    )
                  )
))
