library(plotly)
library(shiny)
library(DT)
library(stringr)




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

#dtKeggPaths=read.table("DATA/kegg_entries.tab",sep="\t",quote = "#",stringsAsFactors = F,header = T)
load("DATA/gmtGeneLists.RData")
dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
dtMeta = dtMeta[dtMeta$ProteinSubtype != "",]
dtMeta$TMT.tag = NULL
dtMeta$X = NULL

xcols=c("PAM50","CoTC","ER","PR","HER.2","Metabolic.cluster", "CAAI", "RPPA.subtype"  )

shinyUI(fluidPage(theme="bootstrap.min.css",
                  titlePanel("Differential Expression"),
                  sidebarLayout(
                    sidebarPanel(width=3,
                                 radioButtons('lstDataType',label="Source data:",choices=c("Protein","RNA"),selected = "Protein"),
                                 selectInput("lstCats",label="Select classification:",choices=xcols, selected = "PAM50" ),
                                 selectInput("lstA",label="Group A:",choices=xcols ,multiple = T, selectize = F),
                                 selectInput("lstB",label="Group B:",choices=xcols ,multiple = T, selectize = F),
                                 conditionalPanel("input.tbSet == \"Volcano\"",
                                                  selectInput("inpEnrichmentCategory","Select Enrichment Search:", choices=names(gmtGeneLists), selected="MSigDB_Hallmarks")
                                 ),
                                 shiny::actionButton(inputId="btRefresh","Refresh" ,  icon = icon("refresh")),
                                 downloadButton('downloadPlot', 'Save as pdf'),
                                 downloadButton('downloadWB', 'Save to Excel')
                    ),
                    mainPanel(
                      fluidRow(
                        busyIndicator2(text = "Loading, please wait...",img = "busyIndicator.gif", wait=1000),
                        column(12,
                               plotOutput("volcanoRNA",width = "100%")
                        )
                      ),
                      fluidRow(
                        column(6,dataTableOutput("DiffExpression")),
                        column(6,dataTableOutput("DE_Enrich"))
                      )
                    )
                  )
))
