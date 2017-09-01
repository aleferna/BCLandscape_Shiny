
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(stringr)
library(shiny)
library(d3heatmap)


dtMeta= read.table("../DATA/meta.tab",sep="\t",header = T, stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
PAM590cols=c("Basal"="#E31A1C", "Her2"="#FB9A99","LumA"="#1F78B4","LumB"="#A6CEE3","Normal"="#33A02C")
dtMeta$Color = PAM590cols[dtMeta$subtype]


shinyServer(function(input, output) {
  
  getData = function(){
    fl=paste0("Results/",input$lstListSelection,".",input$inpSearchType,".tab")  
    if (file.exists(fl)){   
      dt=read.table(fl,sep="\t", stringsAsFactors = F, col.names = c("Source","Category","ListName", "Sample","hits", "nonhits","pvalue","tstat","direction") , comment.char = "#", quote = "#" )
      dt$Sample = paste(dtMeta[dt$Sample,"ProteinSubtype"],str_split_fixed( dt$Sample, "\\.",2)[,2],sep=".")
      dt$pvalue[dt$direction == "Less"] = dt$pvalue[dt$direction == "Less"]*-1
      nn = cast(dt[, c("ListName","Sample","pvalue")], ListName~Sample,value = "pvalue", mean)
      nn = data.frame(nn, stringsAsFactors = F) 
      for (i in 2:ncol(nn)){
        nn[is.nan(nn[,i]),i] = 0
      }
      rownames(nn) = nn$ListName
      nn$ListName = NULL
      idx = rowSums( abs(nn) > 30 ) > 3
      nn = nn[idx,]
      if (nrow(nn) > 100){
        nn = nn[1:100,]
      }
      nn
    }else{
      data.frame()
    }
  }
  

  output$distPlot <- renderD3heatmap({
    withProgress(message = 'Refreshing plot', value = 0, {
      
    nn=getData()

    
    if (nrow(nn) > 1){
      qq=nn
      nn[nn > 100] = 100
      nn[nn < -100] = -100
      d3heatmap(nn, scale = "none",yaxis_width = 300, xaxis_height = 150, cellnote = qq, na.rm = T, theme = "dark")   
    }else{
      
    }
    
    })
    
  })
  
  output$distBarPlot <- renderPlot({
    
  })
  
  
  
})



