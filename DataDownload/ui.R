library(shiny)
library(ggplot2)
library(reshape)
library(plotly)
library(DT)
library(ggrepel)
library(stringr)


dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
dtMeta = dtMeta[dtMeta$ProteinSubtype != "",]

xcols=c()
for (i in colnames(dtMeta)){
  ct = dtMeta[,i]
  ct = str_trim(ct)
  ct = ct[ct!=""]
  ct = unique(ct)
  ct = length(ct)
  if (ct < 10 & ct > 1){
    xcols = c(xcols,i)
  }
}
xcols=sort(xcols)



shinyUI(fluidPage(theme="bootstrap.min.css",
  titlePanel("BC Landscape Omics Data Download"),
  sidebarLayout(
    sidebarPanel(
      selectInput("lstCats",label="Select classification:",choices=xcols, selected = c("CoTC","PAM50") , multiple = T),
      
      tags$textarea(id = 'lstGenes', placeholder = 'Gene List', "SEMA3G\nSHANK3\nSKI\nSLIT2\nSLIT3\nSOD3\nSORBS1" ,  style = "width: 100%; height: 300px"),
      downloadButton('downloadPlot', 'Download Excel')
    ),
    mainPanel(
      DT::dataTableOutput("tblExport")
    )
  )
))
