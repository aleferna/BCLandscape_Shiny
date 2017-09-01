library(stringr)
library(shiny)
library(plotly)
library(DT)
library(networkD3)
load("DATA/geneListsHJ.RData")
load("DATA/gmtGeneLists.RData")
geneLists$All = F
geneLists$Custom = F
#install.packages("amap")



dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
dtMeta = dtMeta[dtMeta$ProteinSubtype != "",]
dtMeta$TMT.tag = NULL
dtMeta$X = NULL
dtMeta$Nstatus = NULL
dtMeta$Grade = NULL
dtMeta$T.status = NULL
dtMeta$CAAI = NULL



shinyUI(fluidPage(theme="bootstrap.min.css",
               
                  # Application title
                  titlePanel("BC Landscape PCA Viewer"),
                  
                  # Sidebar with a slider input for number of bins
                  sidebarLayout(
                    sidebarPanel(width=2,
                                 radioButtons("inpPlotType", "Data source:", c("Protein","mRNA" ), selected = "Protein", inline = T),
                                 selectInput("inpGeneLists","Select Gene List",choices = colnames(geneLists)[2:ncol(geneLists)], selected = "All"),
                                 conditionalPanel("input.inpGeneLists == \"Custom\" ",
                                                  tags$textarea(id = 'lstGenes', placeholder = 'Enter genes here...', c() ,  style = "width: 100%; height: 300px")
                                 ), 
                                 selectInput("cmbColorBy",label="Color by:",choices=colnames(dtMeta), selected = "ProtSub9506" ),
                                 plotOutput("legend", height="300px",),
                                 checkboxInput("inpShowCompAnalysis",label = "Show loadings analysis",value = F),
                                 
                                 conditionalPanel("input.inpShowCompAnalysis",
                                                  selectInput("inpEnrichmentCategory","Select Enrichment Search:", choices=names(gmtGeneLists),selected = "MSigDB_Hallmarks" )
                                 )
                    ),
                    
                    # Show a plot of the generated distribution
                    mainPanel(width=10,
                              fluidRow(
                                column(12,
                                       h3(textOutput("txtStatus"))
                                )
                                
                              ),
                              fluidRow(
                                column(4, plotly::plotlyOutput("pcaPlot12", height="500px", width = "100%")),
                                column(4, plotly::plotlyOutput("pcaPlot34", height="500px", width = "100%")),
                                column(4, forceNetworkOutput ("clusterView", height="500px", width = "100%"))
                              ),
                              DT::dataTableOutput("SampInfo"),
                              conditionalPanel("input.inpShowCompAnalysis",
                                               fluidRow(
                                                 column(4,titlePanel("PCA Loadings Analysis"),offset = 4)
                                               ),
                                               fluidRow(
                                                 column(6, DT::dataTableOutput("PCLoadings")),
                                                 column(6, DT::dataTableOutput("PCEnrich"))
                                               )
                              )
                    )
                  )
)
)
