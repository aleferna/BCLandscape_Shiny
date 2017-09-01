
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
#install.packages("ggrepel")
library(ggrepel)
library(shiny)
library(stringr)
library(ggplot2)
library(doMC)
library(stringr)
library(RCurl)
library(d3heatmap)

load("../DATA/Prot.Rdata")
load("../DATA/RNA.Rdata")
load("../DATA/geneListsHJ.RData")
load("../DATA/dtRNA_Ratio.Rdata")
load("../DATA/Disc.Rdata")

rownames(dtProt) = dtProt$Gene.Symbol
xcols=intersect(colnames(dtProt), colnames(dtRNA_Ratio))
xrows=intersect(rownames(dtProt), rownames(dtRNA_Ratio))
dtProt = dtProt[xrows,xcols]
dtRNA = dtRNA_Ratio[xrows,xcols]
dtDisc = dtDisc[xrows,xcols]
dtDiscRank = dtDisc
for (s in xcols){
  dtDiscRank[,s] = rank(dtDiscRank[,s])
}

dtMeta= read.table("../DATA/meta.tab",sep="\t",header = T, stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
xcols = paste(dtMeta[xcols,"ProteinSubtype"],str_split_fixed( xcols, "\\.",2)[,2],sep=".")
colnames(dtRNA) = xcols
colnames(dtProt) = xcols
colnames(dtDisc) = xcols
colnames(dtDiscRank) = xcols



getInteractors = function (gene){
  #gene="UBE2C"
  url=paste0("http://webservice.thebiogrid.org/interactions?searchNames=true&geneList=",gene,"&includeInteractors=true&includeInteractorInteractions=false&taxId=9606&accesskey=415c368b99fd6a7b25901d156c0061cf")
  
  res=getURI(url)
  res=unlist(str_split(res,"\n"))
  res=str_split_fixed(res,"\t",12)[,8]
  unique(toupper( res))
}

shinyServer(function(session , input, output) {
  updateSelectizeInput(session, 'inpGeneListsBioGrid', choices = rownames(dtRNA),selected = "PTEN", server = TRUE)
  
  
  observe({
    if  (input$inpGeneLists != ""){
      idx=geneLists[, input$inpGeneLists]  
      gns = sort( unique(geneLists$GeneName[idx])) 
      updateTextInput(session, "lstGenes", label = "Gene list:", value = paste(gns,collapse = "\n")  )
    }
  })
  
  observe({
    if  (input$inpGeneListsBioGrid != ""){
      gns = getInteractors(input$inpGeneListsBioGrid)
      updateTextInput(session, "lstGenes", label = "Gene list:", value = paste(gns,collapse = "\n")  )
    }
  })
  
  
  
  
  output$scatterPlot <- renderPlot({
    s = input$inpSample
    
    dtx=data.frame(Prot=log2(dtProt[,s]),RNA=dtRNA[,s], gene=rownames(dtRNA))
    lst=unlist(str_split( input$lstGenes,"\n"))
    dtxS=dtx[dtx$gene %in% lst,]
    d <- ggplot(dtx, aes(Prot,RNA))
    d = d + stat_binhex() + geom_smooth() + geom_point(data=dtxS, aes(Prot,RNA), col="#00FF0090" ) + xlim(-5,5) + ylim(-5,5)
    d = d + theme_dark()  + theme(plot.background = element_rect(fill = "transparent",colour = NA), 
                           text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
                           legend.text=element_text(size=16, color="white"),
                           legend.background=element_rect(fill = "transparent",colour = NA),
                           legend.key = element_rect(fill = "transparent", colour=NA),
                           legend.key.size = unit(3,"lines"),
                           legend.position = "right"
    )
    
    d
  }, bg="transparent")
  
  output$discRankPlot2 <- renderPlot({
    s = input$inpSample
    #s = colnames(dtProt)[3]
    #lst = rownames(dtProt)[1:100]
    
    d = log2(dtProt[,s])  - dtRNA[,s]  
    dtx=data.frame(ProteinAccumulation=d,Rank=rank(d), gene=rownames(dtRNA))
    lst=unlist(str_split( input$lstGenes,"\n"))
    dtx$Z = "All" 
    dtx$Z[dtx$gene %in% lst] = "Selected" #input$inpGeneLists
    
    d <- 
      ggplot(dtx, aes(Rank,ProteinAccumulation)) + geom_line()
    d = d + geom_point(data=dtxS, aes(Rank,ProteinAccumulation), col="#00FF0090", size=2 ) 
    
    d = d + theme_dark()  + theme(plot.background = element_rect(fill = "transparent",colour = NA), 
                                  text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
                                  legend.text=element_text(size=16, color="white"),
                                  legend.background=element_rect(fill = "transparent",colour = NA),
                                  legend.key = element_rect(fill = "transparent", colour=NA),
                                  legend.key.size = unit(3,"lines"),
                                  legend.position = "right"
    )
    
    d
  }, bg="transparent")
  
  
  output$discRankPlot <- renderPlot({
    s = input$inpSample
    
    d = log2(dtProt[,s])  - dtRNA[,s]  
    
    dtx=data.frame(ProteinAccumulation=d,Rank=rank(d), gene=rownames(dtRNA))
    lst=unlist(str_split( input$lstGenes,"\n"))
    dtx$Z = "All" 
    dtx$Z[dtx$gene %in% lst] = "Selected" #input$inpGeneLists
    
    # d <- ggplot(dtx, aes(Rank,Disc))
    # d = d + geom_point(data=dtxS, aes(Rank,Disc), col="#00FF0090", size=2 ) 
    # 
    # print(d)
    
    d <-ggplot(dtx, aes(x=ProteinAccumulation, fill=Z ))+ xlim(-5,5) + 
      geom_density(alpha=.8)+ scale_fill_manual(values = c("#142B5C80", "#00FF0090"))
    
    d = d + theme_dark()  + theme(plot.background = element_rect(fill = "transparent",colour = NA), 
                                  text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
                                  legend.text=element_text(size=16, color="white"),
                                  legend.background=element_rect(fill = "transparent",colour = NA),
                                  legend.key = element_rect(fill = "transparent", colour=NA),
                                  legend.key.size = unit(3,"lines"),
                                  legend.position = "right"
    )
    
    d
  }, bg="transparent")
  
  
  
  
  output$heatPlot = renderD3heatmap({
    #dtDisc[1:3,1:4]
    #lst=sample(rownames(dtDiscRank),100)
    lst=unlist(str_split( input$lstGenes,"\n"))
    res=NULL
    
    
    for (s in colnames(dtDiscRank))  {
      tb=dtDiscRank[rownames(dtDiscRank) %in% lst,s]
      tb=floor(10*tb/nrow(dtDiscRank))
      tb = paste0(10*tb,"%")
      tb=data.frame(table(tb) , stringsAsFactors = F )
      colnames(tb) = c("grp",s)
      
      if (is.null(res)){
        res=tb
      }else{
        res=merge(tb,res,all=T,by="grp")
      }
    }
    
    rownames(res) = res$grp
    
    res = res[paste0(0:9*10,"%"),]
    res[is.na(res)] = 0
    res$grp=NULL
    rownames(res) = NULL
    res = t(res) 
    resLab = res
    res = res/sum(rownames(dtDiscRank) %in% lst)
    
    d3heatmap(Colv = NULL, res, cellnote = resLab  ,scale = "none",labCol = c("RNA","----","---","--","-","+","++","+++","++++","Prot"), theme = "dark")
    
  })
  
  
})
