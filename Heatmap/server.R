library(shiny)
library(stringr)
library(ggplot2)
library(plotly)
library(reshape)
library(dendextend)
library(png)
library(ComplexHeatmap)
library(data.table)


#############################################################################
###################  CSS CUSTOM MODIFICATION FOR DARK THEME ###################
####################################################################################
## Modified css for dark
## www/bootstrap.min.css
## Added: 
## .theme-dark .d3heatmap-tip {
## background-color: black;
## color: white;
## }



load("DATA/geneListsHJ.RData")
load("DATA/DT45.Rdata")
rownames(dtRNA) = dtRNA$ID
rownames(dtProt) = dtProt$ID


gnsClass =fread("DATA/geneCats.tab")
gnsClass = data.frame(gnsClass[!is.na(BCLandscape_Protein_Clustering),], stringsAsFactors = F)
gnsClass$Fredlund_Classifier = NULL
gnsClass$NetModK3Prot = NULL
rownames(gnsClass) = gnsClass$Id
gnsClass$Id = NULL



#order by protein clustering
colOrder=c("OSL2U.0351T1","OSL2U.0381T1","OSL2U.0525T1","OSL2U.0524T1","OSL2U.0512T1","OSL2U.0541T1","OSL2U.0439T1","OSL2U.0536T1","OSL2U.0542T1","OSL2U.0334T1","OSL2U.0219T8","OSL2U.0364T1",
           "OSL2U.0030T1","OSL2U.0075T1","OSL2U.0329T1","OSL2U.0194T1","OSL2U.0382T1","OSL2U.0218T1","OSL2U.0407T1","OSL2U.0484T1","OSL2U.0289T1","OSL2U.0352T1","OSL2U.0537T1","OSL2U.0037T1",
           "OSL2U.0085T5","OSL2U.0299T1","OSL2U.0045T1","OSL2U.0260T1","OSL2U.0236T1","OSL2R.3037T4","OSL2U.0201T1","OSL2U.0213T1","OSL2R.3013T1","OSL2U.0392T1","OSL2U.0129T1","OSL2U.0555T1",
           "OSL2U.0050T1","OSL2U.0135T1","OSL2U.0523T1","OSL2U.0429T1","OSL2R.3059T2","OSL2U.0368T1","OSL2R.3045T5","OSL2U.0214T1","OSL2U.0383T1")

dtRefHeatmapGenePos = read.table(file="DATA/refHeatMap.txt",stringsAsFactors = F)


setMethod(f = "draw_dend",
          signature = "Heatmap",
          definition = function(object,
                                which = c("row", "column"), k = 1, max_height = NULL, ...) {
            
            which = match.arg(which)[1]
            
            side = switch(which,
                          "row" = object@row_dend_param$side,
                          "column" = object@column_dend_param$side)
            
            hc = switch(which,
                        "row" = object@row_dend_list[[k]],
                        "column" = object@column_dend)
            
            gp = switch(which,
                        "row" = object@row_dend_param$gp,
                        "column" = object@column_dend_param$gp)
            
            if(length(hc) == 0) {
              return(invisible(NULL))
            }
            
            if(is.null(hc)) return(invisible(NULL))
            
            dend = as.dendrogram(hc)
            n = length(labels(dend))
            if(nnodes(dend) <= 1) {
              return(invisible(NULL))
            }
            
            dend_padding = unit(1, "mm")
            pushViewport(viewport(name = paste(object@name, which, "cluster", k, sep = "_"), gp=gp,...))
            
            if(side == "left") {
              grid.dendrogram(dend, name = paste(object@name, "dend_row", k, sep = "_"), max_height = max_height, facing = "right", order = "reverse", x = dend_padding, width = unit(1, "npc") - dend_padding*2, just = "left")
            } else if(side == "right") {
              grid.dendrogram(dend, name = paste(object@name, "dend_row", k, sep = "_"), max_height = max_height, facing = "left", order = "reverse", x = unit(0, "mm"), width = unit(1, "npc") - dend_padding*2, just = "left")
            } else if(side == "top") {
              grid.dendrogram(dend, name = paste(object@name, "dend_column", sep = "_"), max_height = max_height, facing = "bottom", y = dend_padding, height = unit(1, "npc") - dend_padding*2, just = "bottom")
            } else if(side == "bottom") {
              grid.dendrogram(dend, name = paste(object@name, "dend_column", sep = "_"), max_height = max_height, facing = "top", y = dend_padding, height = unit(1, "npc") - dend_padding*2, just = "bottom")
            }
            
            upViewport()
            
          })



shinyServer(function(input, output,session) {
  
  # drawProtHeatmap = function(){
  #   gns = unlist(str_split( input$lstGenes, pattern = "\n"))
  #   idx = dtProt$Gene.Symbol %in% gns
  #   if (sum(idx) > 1 & paste0(gns, collapse="") != "" & input$inpPlotType == "Protein" ){
  #     dtx=dtProt[idx,colOrder]
  #     colnames(dtx) = paste(dtMeta[colOrder,"subtype"] , str_split_fixed(colOrder,"\\.",2)[,2], sep=".")
  #     rownames(dtx) = dtProt$Gene.Symbol[idx]
  #     dtx=log2(dtx)
  #     qMx=quantile( data.matrix(dtx),0.99)
  #     qMn=quantile( data.matrix(dtx),0.01)
  #     dtx[dtx > qMx ] = qMx
  #     dtx[dtx < qMn ] = qMn
  #     res=d3heatmap(dtx,scale="none", dendrogram = "row",na.rm = T , show_grid = F, anim_duration = 0, colors = c("#0000FF","#000000","#FF0000") , height=600, theme = "dark")
  #   }
  #   else{
  #     res=d3heatmap(data.frame(a=c(1,2),b=c(3,4) ))  
  #   }
  #     
  # }
  # drawmRNAHeatmap = function(){
  #   gns = unlist(str_split( input$lstGenes, pattern = "\n"))
  #   idx = dtRNA$geneName %in% gns
  #   if (sum(idx) > 1 &  paste0(gns, collapse="") != "" & input$inpPlotType == "mRNA"  ){
  #     dtx=dtRNA[idx,colOrder]
  #     colnames(dtx) = paste(dtMeta[colOrder,"subtype"] , str_split_fixed(colOrder,"\\.",2)[,2], sep=".")
  #     rownames(dtx) = dtRNA$geneName[idx]
  #     qMx=quantile( data.matrix(dtx),0.99)
  #     qMn=quantile( data.matrix(dtx),0.01)
  #     dtx[dtx > qMx ] = qMx
  #     dtx[dtx < qMn ] = qMn
  #     d3heatmap(dtx,scale="none", dendrogram = "row",na.rm = T , show_grid = F, anim_duration = 0, colors = c("#0000FF","#000000","#FF0000") , height=600, theme = "dark")
  #   }
  #   else{
  #     d3heatmap(data.frame(a=c(1,2),b=c(3,4) ))  
  #   }
  # }
  # 
  # drawmiRNAHeatmap = function(){
  #   gns = input$inpmiRNA
  #   idx=dtmiRNA$geneName %in% gns
  #   if (sum(idx) > 1  & input$inpPlotType == "miRNA" ){
  #     dtx=dtmiRNA[idx,colOrder]
  #     colnames(dtx) = paste(dtMeta[colOrder,"subtype"] , str_split_fixed(colOrder,"\\.",2)[,2], sep=".")
  #     rownames(dtx) = dtmiRNA$geneName[idx]
  #     qMx=quantile( data.matrix(dtx),0.99)
  #     qMn=quantile( data.matrix(dtx),0.01)
  #     dtx[dtx > qMx ] = qMx
  #     dtx[dtx < qMn ] = qMn
  #     d3heatmap(dtx,scale="row", dendrogram = "row",na.rm = T , show_grid = F, anim_duration = 0, colors = c("#0000FF","#000000","#FF0000"), height=600, theme = "dark" )
  #   }
  #   else{
  #     d3heatmap(data.frame(a=c(1,2),b=c(3,4) ), theme = "dark")  
  #   }
  # }
  # 
  # drawRefHeatmap = function(){
  #   if( input$inpPlotType == "Reference"){
  #     mypic <- readPNG("DATA/refHeatMap.png")
  #     
  #     plot(1,1, type='n', main="Publication Image", xlab="", ylab="", xaxt='n', yaxt='n', xlim=c(0,1467), ylim=c(0,1772) , yaxs='i', xaxs='i')
  #     
  #     lim <- par()
  #     rasterImage(mypic, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
  #     gns = unlist(str_split( input$lstGenes, pattern = "\n"))
  #     #gns=sample( dtRefHeatmapGenePos$V1,6)
  #     #gns=dtRefHeatmapGenePos$V1[1:10]
  #     y=which( dtRefHeatmapGenePos$V1 %in% gns)
  #     y=length(dtRefHeatmapGenePos$V1) - y
  #     y=1010.0*y/length(dtRefHeatmapGenePos$V1)
  #     y=343 + y
  #     
  #     if (length(y) > 1000){
  #       a=paste0("#FFA50005")   
  #     }else{
  #       a=paste0("#FFA50091")   
  #     }
  #     
  #     
  #     
  #     
  #     for (x in y){
  #       lines(c(150,700),c(x,x), col=a) 
  #     }
  #   }
  # 
  #   
  # }
  # 
  #output$htProt <- renderD3heatmap({ drawProtHeatmap() })
  #output$htmRNA <- renderD3heatmap({ drawmRNAHeatmap() })
  #output$htmiRNA <- renderD3heatmap({ drawmiRNAHeatmap() })
  
  doHeatmap = function(xgp){
    gns = toupper(unlist(str_split( input$lstGenes, pattern = "\n")))
    dtx=switch (input$inpPlotType,
                mRNA=dtRNA,
                Protein=dtProt
    )
    idx = dtx$ID %in% gns
    
    if (sum(idx) > 1 &  paste0(gns, collapse="") != "" ){
      dtx=dtx[idx,]
      dtx$ID = NULL
      
      qMx=quantile( data.matrix(dtx),0.99)
      qMn=quantile( data.matrix(dtx),0.01)
      dtx[dtx > qMx ] = qMx
      dtx[dtx < qMn ] = qMn
      
      xcolorRanges$PR = xcolorRanges$PR[c("Unknown","pos","neg")]
      xcolorRanges$ER = xcolorRanges$ER[c("Unknown","pos","neg")]
      ha = HeatmapAnnotation(df = dtMeta[colnames(dtx),c("PAM50","ER","PR")], 
                             col = list(ER=xcolorRanges$ER,
                                        PR=xcolorRanges$PR,
                                        PAM50=xcolorRanges$PAM50 )
      )
      
      ra = rowAnnotation(df=data.frame(ProteinGroup=gnsClass[rownames(dtx),]), 
                         col=list(
                           ProteinGroup=xcolorRanges$ProtColRange
                         ))
      
      #colnames(dtx) = paste(dtMeta[,"PAM50"] , str_split_fixed(colnames(dtx),"\\.",2)[,2], sep=".")
      
      
      ht_global_opt(heatmap_column_names_gp=xgp) 
      ht_global_opt(heatmap_row_names_gp=xgp) 
      ht_global_opt(heatmap_legend_labels_gp=xgp) 
      ht_global_opt(heatmap_legend_title_gp=xgp)
      ht_global_opt(annotation_legend_labels_gp=xgp)
      ht_global_opt(annotation_legend_title_gp=xgp)
      
      
      p= ra + ComplexHeatmap::Heatmap(dtx,row_dend_gp = xgp,
                                      show_row_names = nrow(dtx )< 50,
                                      column_dend_gp = xgp,
                                      top_annotation = ha,
                                      name = input$inpPlotType ) 
      
      p
    }else{
      NULL  
    }
    
  }
  
  
  output$pltHeatmap <- renderPlot({
    par(fontsize = 12, col="white")
    xgp=gpar(fontsize = 12, col="white", fill="white");
    p = doHeatmap(xgp)
    if (is.null(p)){
      p
    }else{
      
      #save(p,file="p.Rdata")
      #load("p.Rdata")
      #ht_global_opt()
      
      
      
      draw(p )
      
      # ,
      #      row_dend_gp = gpar(fontsize = 12, col="white", fill="white"),
      #      column_dend_gp = gpar(fontsize = 12, col="white", fill="white"),
      #      c
      #      row_names_gp= gpar(fontsize = 12, col="white", fill="white"))
    }
    
    #d3heatmap(dtx,scale="none", dendrogram = "row",na.rm = T , show_grid = F, anim_duration = 0, colors = c("#0000FF","#000000","#FF0000") , height=600, theme = "dark")
    
  }, bg = "transparent")
  
  
  #output$htref <- renderPlot({ drawRefHeatmap() },height = 1200)
  observe({
    if  (input$inpGeneLists != ""){
      idx=geneLists[, input$inpGeneLists]  
      gns = sort( unique(geneLists$GeneName[idx])) 
      updateTextInput(session, "lstGenes", label = "Gene list:", value = paste(gns,collapse = "\n")  )
    }
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 
      paste0("BCLandscape.",input$inpPlotType,".pdf")  
    },
    content = function(file) {
      par( col="black", fg="black")
      xgp=gpar(fontsize = 12, col="black", fill="black");  
      pdf(file,12,12)
      p=doHeatmap(xgp)
      draw(p)
      dev.off()
    }
  )
  
})



