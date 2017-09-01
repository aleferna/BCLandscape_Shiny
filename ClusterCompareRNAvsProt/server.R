library(amap)
library(doMC)
library(shiny)
library(plotly)
library(DT)
library(GSA)
library(reshape2)
library(ggplot2)
library(networkD3)
library(stringr)
library(reshape)

# registerDoMC(20)
# 
# load("../DATA/Prot.Rdata")
# rownames(dtProt) = dtProt$Gene.Symbol
# dtProt$Gene.Symbol = NULL
# load("../DATA/dtRNA_Ratio.Rdata")
# dtRNA=dtRNA_Ratio
# 
# dtRNA$geneName = NULL
# xrows = intersect(rownames(dtRNA),rownames(dtProt))
# xcols = intersect(colnames(dtRNA),colnames(dtProt))
# dtxA = dtRNA[xrows, xcols]
# dtxB = dtProt[xrows, xcols]
# dtxA=hclust( Dist(dtxA,method = "correlation",nbproc = 20), method = "average")
# dtxB=hclust( Dist(dtxB,method = "correlation",nbproc = 20), method = "average")
# save(dtxA,file="dtxA.hclust.Rdata")
# save(dtxB,file="dtxB.hclust.Rdata")
# 




load("../DATA/dtGeneCategories.Rdata")
dtGeneCategories=dtGeneCategories[!duplicated(dtGeneCategories$Gene),]
dtGeneCategories=dtGeneCategories[!is.na(dtGeneCategories$Gene),]
rownames(dtGeneCategories) = dtGeneCategories$Gene



shinyServer(function(input, output) {

  calcClusterCompare = function(cluA,cluB,cluC){
     
     # cluA=dtxA
     #  cluB=dtxB
     #  cluC=dtxC

    grps=list()
    for (grp in unique(cluA)  ){
      grps[[grp]] = names(cluA)[cluA == grp]
    }
    
    for (grp in unique(cluB)  ){
      grps[[grp]] = names(cluB)[cluB == grp]
    }
    
    # for (grp in unique(cluC)  ){
    #   grps[[grp]] = names(cluC)[cluC == grp]
    # }
    # 
        
    nodes = data.frame(name=names(grps),idx=1:length(grps) )
    nodes$idx = nodes$idx -1
    rownames(nodes) = nodes$name
    
    grpAll = intersect(names(cluA),names(cluB))
    
    edges = foreach (a = unique(cluA), .combine = rbind) %:%
      foreach (b =  unique(cluB), .combine = rbind) %do% {
        if (a > b){
          grpA = grps[[a]]
          grpB = grps[[b]]
          grpA = intersect(grpA,grpAll)
          grpB = intersect(grpB,grpAll)
          weight = sum(grpA %in% grpB)
          gns=paste(intersect(grpB,grpA),collapse=",")
          c(nodes[a,"idx"],nodes[b,"idx"],weight)
        }
      }

        
    edges=data.frame(edges, stringsAsFactors = F, row.names = 1:nrow(edges))
    colnames(edges) = c("source","target","weight")
    edges=edges[edges$weight > 0,]
    list(nodes=nodes,edges=edges)
  }
  
  
  output$sankPlot <- renderSankeyNetwork({
    k=10
    k=input$ctClus
    
    load("dtxA.hclust.Rdata")
    load("dtxB.hclust.Rdata")
    xrows=dtxA$labels
    
    dtxA=cutree(dtxA,k=k)
    dtxA = paste0("RNA." , dtxA)
    names(dtxA) = xrows
    
    
    dtxB=cutree(dtxB,k=k)
    dtxB = paste0("PROT." , dtxB)
    names(dtxB) = xrows
    
    dtxC=dtGeneCategories$BCLandscape_Protein_Clustering
    names(dtxC) = dtGeneCategories$Gene
    dtxC[is.na(dtxC)] = "Other"
    
    
    res=calcClusterCompare(dtxA,dtxB,dtxC)
    sankeyNetwork(Links = res$edges, Nodes = res$nodes, Source = "source",
                  Target = "target", Value = "weight", NodeID = "name",fontSize = 14)

  })
  
})
