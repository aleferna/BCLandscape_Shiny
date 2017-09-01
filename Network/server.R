library(shiny)
library(stringr)
library(ggplot2)
library(plotly)
library(reshape)
library(ComplexHeatmap)
library(pcaMethods)
library(Rmisc)
library(shiny)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(d3Network)
library(GSA)
#install.packages("GSA")
library(grid)
require(reshape)
require(doMC)
library(DT)
registerDoMC(10)

load("DATA/gmtGeneLists.RData")


load("DATA/RNA.Pearson7.D1.Rdata")
edgesRNA=edges
nodesRNA=nodes

load("DATA/Prot.Pearson7.D1.Rdata")
edgesProt=edges
nodesProt=nodes


load("DATA/Prot.Rdata")
dtMeta= read.table("DATA/meta.tab",sep="\t",header = T, stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
PAM590cols=c("Basal"="#E31A1C", "Her2"="#FB9A99","LumA"="#1F78B4","LumB"="#A6CEE3","Normal"="#33A02C")


load("DATA/geneListsHJ.RData")
load("DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol


#load("../DATA/dtNetEdges.ProtFull.Rdata")
#load("../DATA/dtSubNetworkModularity.Rdata")
load("DATA/dtGeneCategories.Rdata")
#load("../DATA/dtGeneLevels.Rdata")



# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  
  updateSelectizeInput(session, 'inpSelectedGene', choices = rownames(dtProt), selected = "MET", server = TRUE)
  
  theme_black_ggplot = theme(axis.ticks = element_blank(), 
                             axis.text.x = element_blank(), 
                             axis.text.y = element_blank(),
                             panel.background = element_rect(fill = "black",colour = NA),
                             plot.background = element_rect(fill = "black",colour = NA),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             legend.position = "right",
                             legend.text=element_text(size=16, color="white"),
                             legend.background=element_rect(fill = "black",colour = NA),
                             legend.key = element_rect(fill = "black", colour=NA),
                             legend.key.size = unit(3,"lines"))
  
  updateSelected = function(gn){
    if (input$inpNetSource == "RNA"){
      edges=edgesRNA
      nodes=nodesRNA
    }else{
      edges=edgesProt
      nodes=nodesProt
    }
    
    i=1
    nodes$Selected[nodes$name == gn] = gn
    if (gn %in% SubNetWorkModularity$index[,1] & gn != ""){
      lst = paste(SubNetWorkModularity$index[gn,2])
      while(lst %in% names(SubNetWorkModularity)){
        print(lst)
        gns = SubNetWorkModularity[[lst]]
        idx = nodes$name %in% gns & is.na(nodes$Selected )
        nodes$Selected[idx] = lst
        i = i + 1
        lst = str_split(lst,"\\.")[[1]]
        lst = paste(lst[1:length(lst)-1],collapse=".")
      }
    }
    nodes$Selected[is.na( nodes$Selected )] = "Other"
    nodes
  }
  
  
  drawNetwork = function(){
    
    
    if (input$inpNetSource == "RNA"){
      edges=edgesRNA
      nodes=nodesRNA
    }else{
      edges=edgesProt
      nodes=nodesProt
    }
    
    nodes$Selected = "NA"
    nodes$labeltext = nodes$name
    
    if (input$inpPlotType == "Gene List"){
      genes = unlist(str_split( input$lstGenes,pattern = "\n"))
      genes = unlist(str_split_fixed(genes,"\t",2)[,1])
      nodes$Selected = paste(nodes$name %in% genes)
      myColors <- c("#80808080","yellow")
      names(myColors) <- c("FALSE","TRUE")
      colScale <- scale_colour_manual(name = "Selected",values = myColors)
      xguide = guides(colour = guide_legend(override.aes = list(size=12)))
    }
    
    if (input$inpPlotType == "Gene Category"){
      dtCat = unlist(str_split( input$lstGeneClasses,pattern = "\n"))
      dtCat = str_split_fixed(dtCat,"\t",2)
      dtCat = data.frame(dtCat,stringsAsFactors = F )
      colnames(dtCat) = c("Gene","Class")      
      dtCat = dtCat[!duplicated( dtCat$Gene) ,] 
      rownames(dtCat) = dtCat$Gene
      cats=c(unique(dtCat$Class),"Unknown")
      nodes$Selected = "Unknown"
      idx = nodes$name %in% dtCat$Gene
      nodes$Selected[idx] =   dtCat[nodes$name,"Class"][idx]
      myColors <- c(rainbow( length(cats)-1),"#80808080")
      names(myColors) <- cats
      colScale <- scale_colour_manual(name = "Selected",values = myColors)
      xguide = guides(colour = guide_legend(override.aes = list(size=8)))
    }
    
    if (input$inpPlotType == "Gene Level"){
      dtCat = unlist(str_split( input$lstGeneLevels, pattern = "\n"))
      dtCat = str_split_fixed(dtCat,"\t",2)
      dtCat = data.frame(dtCat,stringsAsFactors = F )
      colnames(dtCat) = c("Gene","Level")      
      dtCat = dtCat[!duplicated( dtCat$Gene) ,] 
      rownames(dtCat) = dtCat$Gene
      cats=c(unique(dtCat$Class),"Unknown")
      nodes$Selected = 0
      idx = nodes$name %in% dtCat$Gene
      nodes$Selected[idx] =   as.double(dtCat[nodes$name,"Level"][idx])
      Qs=quantile(nodes$Selected, c(0.05,0.95))
      nodes$Selected[nodes$Selected < Qs[1]] = Qs[1]
      nodes$Selected[nodes$Selected > Qs[2]] = Qs[2]
      colScale <- scale_color_gradient(low = "darkred", high="yellow" )
      xguide = guides( guide_colorbar( title=input$inpGeneLevelLists,barwidth = unit(2,"cm"), barheight=unit(4,"cm"),title.theme = element_text(color = "white", angle = 0) ))
    }
    
    if (input$inpPlotType == "Gene Neighborhood"){
      gn=isolate(input$inpSelectedGene)
      nodes=updateSelected(gn)
      nodes$labeltext = nodes$name
      
      i=length(unique(nodes$Selected))
      if (i > 2){
        xcol=colorRampPalette(colors= c("orange","darkred")  )(i-2)
        myColorsB <- c("green", xcol,"#80808080")  
      }else{
        myColorsB <- c("green","#80808080")  
      }
      
      clus=sort(setdiff(unique(nodes$Selected),c(gn,"Other")))
      names(myColorsB) <- c(gn,clus,"Other")
      colScale <- scale_colour_manual(name = "Selected",values = myColorsB)
      xguide = guides(colour = guide_legend(override.aes = list(size=12)))
      
      
    }
    if (length(unique(nodes$Selected)) == 1){
      p=ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + theme_black_ggplot
      return(p)
    }
    
    
    
    p = ggplot(data = nodes, aes(x = x, y = y)) 
    if (input$inpShowEdges){
      p = p + 
        geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges[edges$pearson < 0 & !edges$BioGrid,], size = 0.3, color="#80000080") +
        geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges[edges$pearson > 0 & !edges$BioGrid,], size = 0.3, color="#00008080") +
        geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges[edges$pearson > 0 &  edges$BioGrid,], size = 0.5, color="#0080FF80") 
    }
    
    p = p + 
      labs(x=NULL , y=NULL) + 
      theme_black_ggplot
    
    
    szPoint = 1.5
    
    
    if (!is.null(ranges$x)){
      szPoint = abs(ranges$x[1] - ranges$x[2])
      szPoint = 1.5 * (max(nodes$x) - min(nodes$x))  / szPoint
    }
    
    p = p + geom_point(aes(text = name, color=Selected), size = szPoint) + colScale + xguide     
    
    if (!is.null(ranges$x)){
      p = p + xlim(ranges$x[1], ranges$x[2]) + ylim(ranges$y[1],ranges$y[2])
      #p = p + geom_text_repel(aes(text = labeltext, label=labeltext), fontface = 'bold', color="white" ,max.iter = 100)
      p = p + geom_text(aes(text = labeltext, label=labeltext), fontface = 'bold', color="white")
    }
    p
  }
  
  getGMT = reactive({
    gmt= gmtGeneLists[[input$inpEnrichmentCategory]]
    return(gmt)
  })
  
  observe({
    if  (input$inpGeneLists != "" & input$inpEnrichmentCategory != ""){
      #idx=geneLists[, input$inpGeneLists]  
      gmt = getGMT()
      idx=gmt$geneset.names == input$inpGeneLists
      
      gns = unlist(gmt$genesets[idx])
      gns = gns[gns != ""]
      gns = str_split_fixed(gns,",",2)[,1]
      updateTextInput(session, "lstGenes", label = "Gene list:", value = paste(gns,collapse = "\n")  )
    }
  })
  
  observe({
    if (input$inpEnrichmentCategory != ""){
      
      gmt = getGMT()
      #gmt=GSA.read.gmt(paste0("../DATA/GeneLists/CORUM.gmt"))
      updateSelectizeInput(session, "inpGeneLists", label = "Select Gene List:", choices=gmt$geneset.names  , server = T)
      
      #idx = !is.na(dtGeneCategories[,input$inpGeneClassLists])
      #gnCats=paste(paste(dtGeneCategories$Gene[idx],dtGeneCategories[idx,input$inpGeneClassLists],sep="\t") ,collapse = "\n")  
      #updateTextInput(session, "lstGeneClasses", label = "Gene classification:", value = paste(gnCats,collapse = "\n")  )
      
    }
  })
  
  
  observe({
    if (input$inpGeneClassLists != ""){
      idx = !is.na(dtGeneCategories[,input$inpGeneClassLists])
      gnCats=paste(paste(dtGeneCategories$Gene[idx],dtGeneCategories[idx,input$inpGeneClassLists],sep="\t") ,collapse = "\n")  
      updateTextInput(session, "lstGeneClasses", label = "Gene classification:", value = paste(gnCats,collapse = "\n")  )
    }
  })

  observe({
    if (input$inpNetSource == "RNA"){
      edges=edgesRNA
      nodes=nodesRNA
    }else{
      edges=edgesProt
      nodes=nodesProt
    }
    
    updateSelectInput(session, inputId = "inpGeneLevelLists", label = "Select Gene Level ", choices = colnames(dtGeneLevels), selected = colnames(dtGeneLevels)[1]  )
  })
  
  observe({
    if (input$inpNetSource == "RNA"){
      edges=edgesRNA
      nodes=nodesRNA
    }else{
      edges=edgesProt
      nodes=nodesProt
    }
    
    
    
    
    dtGeneLevels$kcore
    if (input$inpGeneLevelLists != ""){
      xcat = input$inpGeneLevelLists
      idx = !is.na(dtGeneLevels[,xcat])
      gns=rownames(dtGeneLevels)[idx]
      xvals=dtGeneLevels[idx,xcat]
      xvals=as.integer(1000*xvals)/1000  
      gnLeves=paste(paste(gns,xvals,sep="\t") ,collapse = "\n")  
      updateTextInput(session, "lstGeneLevels", label = "Gene level:", value = paste(gnLeves,collapse = "\n")  )
    }
  })
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  
  
  
  
  observeEvent(input$plot_brush, {
    brush = input$plot_brush
    if (!is.null(brush) && brush$xmax - brush$xmin > 50 &&  brush$ymax - brush$ymin > 50 ) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    }else{
      ranges$x <- NULL
      ranges$y <- NULL  
    }
  })
  
  
  
  observeEvent(input$plot_click, {
    if (input$inpNetSource == "RNA"){
      edges=edgesRNA
      nodes=nodesRNA
    }else{
      edges=edgesProt
      nodes=nodesProt
    }
    
    dtx=nearPoints(nodes, input$plot_click,"x","y", threshold = 300, maxpoints = 1)
    if (nrow(dtx) > 0){
      
      if (input$inpPlotType != "Gene Neighborhood"){
        updateRadioButtons(session, "inpPlotType", selected="Gene Neighborhood")  
      }
      updateSelectizeInput(session, 'inpSelectedGene', choices = nodes$name, selected = dtx$name, server = TRUE)
      
    }
    
    print(paste("Selected" ,dtx$name))
  })
  
  observeEvent(input$plot_dblclick, {
    ranges$x <- NULL
    ranges$y <- NULL
  })
  
  forceRedraw <- reactiveValues(ok = 1)
  
  observe({
    
    if (input$inpSelectedGene != ""){
      forceRedraw$ok = runif(1)
    }
    
  })
  
  output$NetPlot <- renderPlot( { 
    
    ok=forceRedraw$ok
    drawNetwork()
    
    
    
  })
  
  
  
  getEnrichment = function(dtx,lvl){
    
    if (input$inpNetSource == "RNA"){
      edges=edgesRNA
      nodes=nodesRNA
    }else{
      edges=edgesProt
      nodes=nodesProt
    }
    
    
    
    
    gnsClu = SubNetWorkModularity[[lvl]]
    flgmt = input$inpEnrichmentCategory
    #flgmt="MSigDB_Hallmarks"
    
    if (!dir.exists(paste0("Cache/GeneEnrichment/",flgmt))){
      dir.create(paste0("Cache/GeneEnrichment/",flgmt), recursive = T)
    }
    
    flCache=paste0("Cache/GeneEnrichment/",flgmt,"/",input$inpNetSource,".",lvl,".Rdata")
    if (!file.exists(flCache)){
      gmt=getGMT()
      withProgress(message=paste("Updating enrichment for:", flgmt),value=0,{
        res=foreach (i = 1:length(gmt$genesets), .combine=rbind) %do% {
          incProgress(1/length(gmt$genesets), detail = gmt$geneset.names[i]  )
          gns = gmt$genesets[[i]]
          gns = unlist(str_split_fixed(gns,",",2)[,1])
          gns = gns[gns %in% dtx$name]
          ctHits=sum(gnsClu %in% gns)
          ctSel=length(gnsClu)
          ctTot=nrow(dtx)
          ctLst=length(gns)
          if (ctHits > 3){
            pval = 1.0 - phyper(ctHits, ctLst, ctTot - ctLst, ctSel)
          }else{
            pval = 1
          }
          if (pval == 0){
            pval = 300
          }else{
            pval=-10*log10(pval)  
          }
          incProgress(1/length(gmt$genesets), detail = gmt$geneset.names[i]  )
          c(lvl[length(lvl)], gmt$geneset.names[i], pval, ctHits, ctLst, ctSel)
        }
      }
      )
      save(res,file=flCache)
    }else{
      load(flCache)
    }
    res
  }
  
  
  EnrichDataTable =  reactive({
    
    gn=input$inpSelectedGene
    
    if (is.null(gn) | gn == ""){
      return(NULL);
    }
    dtx=updateSelected(gn)
    clus=setdiff(unique(dtx$Selected),c(gn,"Other"))
    if (length(clus)== 0)
      return(NULL)
    
    xa=foreach (clu = clus, .combine=rbind) %do% {
      
      getEnrichment(dtx,clu)
    }
    
    
    xa=data.frame(xa, stringsAsFactors = F)
    colnames(xa)=c("Neighborhood","List","pval", "ctHits","ctLst","ctSel" ) 
    xa$pval = as.integer( as.double(xa$pval)*10)/10
    xa=cast(xa,List~Neighborhood, value="pval")
    
    
    idx = rowMax(data.matrix(xa[,2:ncol(xa)])) > 30
    if (sum(idx)>1){
      xa=xa[idx,]
      idx=order(rowMax(data.matrix(xa[,2:ncol(xa)])), decreasing = T)
      rownames(xa) = xa$List
      xa$List=NULL
      if (length(colnames(xa)) > 1){
        return (xa[idx, ])  
      }
    }
    return (NULL)
  })
  
  
  
  observe({
    
    rows = input$GeneNeighborhoodEnrich_rows_selected
    if (! is.null(rows)){
      isolate({dtx=EnrichDataTable()})
      
      rows = rownames(dtx[rows,])
      if (!is.null(rows)){
        isolate({flGmt=input$inpEnrichmentCategory})
        gmt=gmtGeneLists[[flGmt]]
        gns=c()
        for (r in rows){
          idx=which(r==gmt$geneset.names)
          gns = sort(unique(c(gns,gmt$genesets[[idx]])))
          gns = unlist(str_split_fixed(gns,",",2)[,1])
        }
        gns = gns[gns != ""]
        updateTextInput(session, "lstGenes", label = "Gene list:", value = paste(gns,collapse = "\n")  )
        updateRadioButtons(session, "inpPlotType", selected="Gene List")  
      }else{
        updateTextInput(session, "lstGenes", label = "Gene list:", value = "" )
      }
    }
    
  })
  
  
  
  output$GeneNeighborhoodEnrich = DT::renderDataTable({
    xa=EnrichDataTable()
    xa = cbind(ListName=rownames(xa), xa)
    xa=datatable(xa,style = 'bootstrap',rownames = F, selection="single",autoHideNavigation = T,
                 caption = htmltools::tags$caption(
                   style = 'caption-side: bottom; text-align: center;',
                   'Enrichment tool: ', htmltools::em('Highlight genes in plot.')
                 )
    )
    
    xa
  })
  
  
  output$barplot = renderPlot({
    gn=input$inpSelectedGene
    if (!is.null(gn) & gn != "" & gn %in% dtProt$Gene.Symbol){
      idx = dtProt$Gene.Symbol == gn
      
      dtx=log2(dtProt[idx,2:ncol(dtProt)])
      
      dtx=data.frame(t(dtx))
      dtx$Protein = dtx[,1]
      dtx$samp = rownames(dtx)
      dtx$Type = dtMeta[dtx$samp,"subtype"]
      for (xtype in unique(dtx$Type)){
        idx = dtx$Type == xtype
        dtx$xType[idx] = 1:sum(idx)
      }
      
      #dtx$ProteinSubtype = dtMeta[dtx$samp,"ProteinSubtype"]
      
      ggplot(dtx, aes(x=xType, y=Protein,fill=Type) )+
        geom_bar(stat = "identity", position = "dodge")+
        ylim(-1,2)+
        ggtitle(gn)+
        facet_wrap(~Type,nrow=1)+
        theme_dark()+
        theme(axis.ticks = element_blank(), 
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(),
              panel.background = element_rect(fill = "black",colour = NA),
              plot.background = element_rect(fill = "black",colour = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
              legend.position="none"
        )+ 
        scale_fill_manual(name = "Type",values = PAM590cols)
      
    }else{
      par(bg = 'black')
      plot(0)
    }
    #
    
    
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste('NetworkPlot.pdf', sep='') },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = 14, height = 10)
      ggsave(file, plot = drawNetwork(), device = device)
    }
  )
  
  
  
  
  
  # #  Too slow!!!!
  # output$networkPlot <- renderVisNetwork({
  #   gns=geneLists$GeneName[geneLists$FDATargets]
  #   gns=gns[gns %in% edges$node1 | gns %in% edges$node2]
  #   nodes=data.frame(name=gns, stringsAsFactors = F)
  #   rownames(nodes) = nodes$name
  #   nodes$id = 1:nrow(nodes)
  #   nodes$group = sample(c("A","B"),nrow(nodes), replace = T)
  #   idx = edges$node1 %in% gns |  edges$node2 %in% gns 
  #   
  #   edges=dtE[idx,]
  #   edges$to = nodes[edges$node1,"id"]
  #   edges$from = nodes[edges$node2,"id"]
  #   
  #   
  #   visNetwork(nodes, edges, width = "100%")
  #
  #
  #  #
  #  #
  #  #
  #  # d3ForceNetwork(Nodes = MisNodes,
  #  #                Links = MisLinks,
  #  #                Source = "source", Target = "target",
  #  #                Value = "value", NodeID = "name",
  #  #                Group = "group", width = 400, height = 500,
  #  #                 standAlone = FALSE,
  #  #                parentElement = '#networkPlot')
  # #
  # 
  # 
  #  d3ForceNetwork(Links = edges,
  #                 Nodes = nodes,
  #                 Source = "to",
  #                 Target = "from",
  #                 NodeID = "id",
  #                 width = 550,
  #                 height = 400,
  #                 opacity = 0.9,
  #                 zoom = TRUE
  #  )
  # # #
  # #
  # })
  
})




