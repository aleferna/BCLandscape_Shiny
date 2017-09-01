library(shiny)
library(doMC)
library(ggplot2)
library(reshape) 
library(plotly)
library(DT)
library(ggrepel)
library(stringr)
require(xlsx)

load("DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL

load("DATA/RNA.Rdata")
rownames(dtRNA) = dtRNA$geneName
dtRNA$geneName = NULL

xcols=intersect(colnames(dtProt),colnames(dtRNA))
dtRNA=dtRNA[,xcols]
dtProt=dtProt[,xcols]

load("DATA/CNA.Rdata")
rownames(dtCNA) = dtCNA$geneName
xcolsMissing = xcols[!xcols %in% colnames(dtCNA)]
dtCNA[,xcolsMissing] = NA
dtCNA = dtCNA[,xcols]


dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
for (i in colnames(dtMeta)){
  dtMeta[,i] = str_trim(dtMeta[,i])
}




shinyServer(function(session,input, output) {
  
  observe({
    updateSelectizeInput(session = session, inputId = "lstSelectedGene",label="Select gene:", choices = rownames(dtRNA), selected = "MET",server = T )
  })
  
  
  getData = reactive({
    xcols=intersect(colnames(dtProt),colnames(dtRNA))
    
    genes = unlist(str_split( input$lstGenes,pattern = "\n"))
    #genes=c("SDPR","SELP","SEMA3G")
    
    idx = genes[genes %in% rownames(dtProt)]
    X=t(dtRNA[idx,])
    colnames(X) = paste0(colnames(X),"_RNA")
    
    idx = genes[genes %in% rownames(dtRNA)]
    Y = t(dtProt[idx,])
    colnames(Y) = paste0(colnames(Y),"_Prot")
    
    lstMetas = input$lstCats
    #lstMetas =sample(colnames(dtMeta),4)
    Z=dtMeta[xcols,lstMetas]
    
    idx = genes[genes %in% rownames(dtCNA)]
    XX=t(dtCNA[idx,])
    colnames(XX) = paste0(colnames(XX),"_CNA")
    
    xres = cbind(Z,X)
    xres = cbind(xres,Y)
    xres = cbind(xres,XX)
    xres
  })
  
  output$tblExport = DT::renderDataTable({
   
    xres=getData()
    datatable(xres, style = 'bootstrap', rownames= T,
              caption = htmltools::tags$caption(
                style = 'caption-side: bottom; text-align: center;',
                paste(input$inpEnrichmentCategory, 'Enrichment:'), htmltools::em('Enrichment by contrast measure.')
              )) 
  })
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 
      "BCLandExtract.xlsx"
    },
    content = function(file) {
      xres=getData()
      wb=createWorkbook()
      ws=createSheet(wb,"BCLandOmics")
      addDataFrame(xres,ws)
      saveWorkbook(wb,file)
    }
  )
  
  
})
