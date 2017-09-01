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
library(ggplot2)
library(gridExtra)
library(grid)

registerDoMC(20)
load("DATA/colorRanges.Rdata")

xcolorRanges$subType = xcolorRanges[["subtype"]]

#dtMeta$Nstatus = stringr::str_trim(dtMeta$Nstatus,side = "both")
#dtMeta$T.status = stringr::str_trim(dtMeta$T.status,side = "both")
#unique(dtMeta$T.status) %in% names(xcolorRanges[["T.status"]])
#unique(dtMeta$Nstatus) %in% names(xcolorRanges[["Nstatus"]])

load("DATA/gmtGeneLists.RData")
load("DATA/Prot.Rdata")
rownames(dtProt) = toupper( dtProt$Gene.Symbol)
dtProt$Gene.Symbol = NULL
load("DATA/RNA.Rdata")
rownames(dtRNA) = toupper( dtRNA$geneName)
dtRNA$geneName = NULL
load("DATA/dtGeneCategories.Rdata")
dtGeneCategories$BCLandscape_Protein_Clustering
dtGeneCategories=dtGeneCategories[!duplicated(dtGeneCategories$Gene),]
dtGeneCategories=dtGeneCategories[!is.na(dtGeneCategories$Gene),]
rownames(dtGeneCategories) = dtGeneCategories$Gene

load("DATA/geneListsHJ.RData")
geneLists$All = F
geneLists$Custom = F

dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id

shinyServer(function(input, output, session) {
  
  output$clusterView <- renderForceNetwork({
    genes=rownames(data())
    
    dtx=data()
    smps = colnames(dtx)
    
    dtx=hclust(Dist(t(dtx), method = "correlation", nbproc = 10), method = "average")
    dtx=as.dendrogram(dtx)
    
    ul <- function(x, parent ) {
      if (is.list(x)) {
        i=0
        foreach (y = x, .combine=rbind) %do% {
          i = i + 1 
          child=paste0(parent,".",i)
          rbind(ul(y, child ),c(parent,child,attr(y,"height")))
        }
      }else{
        name = attr(x, "label")
        c(parent,name,0)
      }
      
    }
    dtx=ul(dtx,"root")
    
    edges=data.frame(dtx,stringsAsFactors = F)
    colnames(edges)=c("source","target","weight")
    idxBad = grepl(edges$target, pattern = "OSL")
    Bad=edges[idxBad,]
    edges=edges[!idxBad,]
    for (i in 1:nrow(Bad)){
      idx = edges$target == Bad$source[i] 
      edges$target[idx] = Bad$target[i] 
    } 
    
    #xscale = 'd3.scale.ordinal().domain(["Basal","Her2","LumA","LumB","Normal","Group"]).range(["#E31A1C", "#FB9A99","#1F78B4","#A6CEE3","#33A02C","#808080"])'
    #xvals='"Basal","Her2","LumA","LumB","Normal","Group"'
    cls=input$cmbColorBy
    if (cls %in% labels(xcolorRanges)){
      cls=xcolorRanges[[cls]]
      xvals=paste0(paste0('"',c(labels(cls),"Group"),'"'), collapse=',')
      xcols=paste0(paste0('"',c(cls,"#808080"),'"'), collapse=',')
      xscale = paste0('d3.scale.ordinal().domain([',xvals,']).range([',xcols,"])")
    }else{
      x=dtMeta[smps,input$cmbColorBy]
      x=x[!is.na(x)]
      a=min(x)
      b=max(x)
      xscale = paste0('d3.scale.linear().domain([',a,",",b,']).range(["#d73027","#1a9850"])')
    }
    
    
    
    nodes=data.frame(name=sort(unique(c(edges$source,edges$target))),subtype="" , stringsAsFactors = F)
    nodes$subtype = dtMeta[nodes$name,input$cmbColorBy]
    idx=!nodes$name %in% smps
    nodes$subtype[idx] = "Group"
    nodes[!idx,"sz"] = 100
    nodes[idx,"sz"] = 2
    
    nodes$id = 0:(nrow(nodes)-1)
    rownames(nodes) = nodes$name
    edges$source = nodes[edges$source,"id"]
    edges$target = nodes[edges$target,"id"]
    edges$weight=10*as.double( edges$weight)+0.1
    
    
    
    
    
    forceNetwork(legend = F, zoom = T, linkDistance = 1, Links = edges, Nodes = nodes, Source = "source",
                 Target = "target", Value = "weight", NodeID = "name", Nodesize = "sz",fontSize = 20,
                 Group = "subtype", 
                 opacity = 1 , 
                 colourScale  = xscale )
    
  })
  
  
  data = reactive({
    res=NULL
    genes=c()
    if (input$inpGeneLists != "All"){
      genes = unlist(str_split( input$lstGenes,pattern = "\n"))
      genes = unlist(str_split_fixed(genes,"\t",2)[,1])
      genes = toupper(genes)
      genes = str_trim(genes,"both")
      genes = genes[genes != ""]
    }else{
      if (input$inpPlotType == "Protein"){
        genes = rownames(dtProt)
      }else{
        genes = rownames(dtRNA)
      }
    }
    genes=unique(genes)
    
    if (input$inpPlotType == "Protein"){
      idx = genes %in% rownames(dtProt)
      genes = genes[ idx]
      if (length(genes) > 0){
        res=dtProt[genes,]    
        res=log2(res)
      }
    }else{
      idx=genes %in% rownames(dtRNA) 
      genes = genes[idx]
      
      xcols=intersect(colnames(dtRNA),colnames(dtProt))
      if (length(genes) > 0){
        res=dtRNA[genes,xcols]
      }
    }
    res
  })
  
  
  output$txtStatus = reactive({
    res=NULL
    if (input$inpGeneLists != "All"){
      genes = unlist(str_split( input$lstGenes,pattern = "\n"))
      genes = unlist(str_split_fixed(genes,"\t",2)[,1])
      genes = toupper(genes)
      genes = str_trim(genes,"both")
      genes = genes[genes != ""]
      genes=unique(genes)
      
      if (input$inpPlotType == "Protein"){
        idx = genes %in% rownames(dtProt)
      }else{
        idx = genes %in% rownames(dtRNA)
      }
      
      if (sum(!idx) > 0 ){
        missing=paste0(genes[!idx], collapse=", ")
        if (nchar(missing) > 50) missing = paste0(substr(missing,1,50),"...")
        hdr=paste0("Analyzing ",  sum(idx)," Genes, Omiting (",sum(!idx),"): ",missing)
      }else{
        hdr=paste0("Analyzing ", sum(idx)," Genes")
      }
    }else{
      if (input$inpPlotType == "Protein"){
        hdr=paste0("Analyzing ", nrow(dtProt)," Genes")
      }else{
        hdr=paste0("Analyzing ", nrow(dtRNA)," Genes")
      }
    }
    
  })
  
  dopca = reactive({
    dtx = data()
    if (!is.null(dtx)){
      prcomp(t(dtx), scale=F)  
    }else{
      NULL
    }
  })
  
  scores = reactive({
    pca = dopca()
    if (!is.null(pca)){
      subType=dtMeta[rownames(pca$x),input$cmbColorBy]
      #subTypeProt=dtMeta[rownames(pca$x),"ProteinSubtype"]
      data.frame(subType, pca$x[,1:4])
    }else{
      NULL  
    }
  })
  
  loadings = reactive({
    pca <- dopca()
    if (!is.null(pca)){
      idx = pca$rotation[,1] == -1000
      for (i in 1:4){
        sdx=sd(pca$rotation[,i])
        mnx=mean(pca$rotation[,i])
        idx=idx | pca$rotation[,i] > mnx + 3*sdx | pca$rotation[,i] < mnx - 3*sdx
      }
      if (sum(idx) > 1) 
        data.frame(pca$rotation[idx,1:4])
      else 
        data.frame(PC1=0,PC2=0,PC3=0,PC4=0)
    }else{
      NULL
    }
  })
  
  getColorScale = reactive({
    cls = input$cmbColorBy 
    
    if (cls %in% names(xcolorRanges)){
      scale_color_manual(name = "subType", values = xcolorRanges[[cls]]  ,  na.value = "black" ) 
    }else{
      scale_colour_gradient(name = "subType", low ="#d73027", high="#1a9850" ,  na.value = "black")
    }
    
  })
  
  output$legend = renderPlot({
    
    dtx = scores()
    if (is.null(dtx)){
      return(NULL)
    }
    dtx$ID = rownames(dtx)
    gx=ggplot(data=dtx) +
      geom_point(aes(x=PC1, y=PC2,  colour=subType, text=ID))+
      getColorScale() + 
      theme_dark() + 
      theme(legend.background=element_rect(fill = "transparent",colour = NA), 
            legend.position="left",
            legend.text=element_text(size=24),
            legend.key.size=unit(1,"cm"),
            text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
            legend.title =element_blank(), 
            panel.background = element_rect(fill = "transparent",colour = NA), 
            plot.background = element_rect(fill = "transparent",colour = NA)
      )+
      guides(colour = guide_legend(override.aes = list(size=10, name='', title='')))
    
    
    
    
    #Extract Legend 
    tmp <- ggplot_gtable(ggplot_build(gx)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    
    grid.draw(legend)
    
  }, bg="transparent")
  
  output$pcaPlot12 = renderPlotly({
    dtx = scores()
    pca = dopca()
    dtx$ID = rownames(dtx)
    #dtx$subType[is.na(dtx$subType)] = "NA"
    
    ggplot(data=dtx) +
      geom_point(aes(x=PC1, y=PC2,  colour=subType, text=ID, size =3 ) )+
      getColorScale() +
      theme_dark() + 
      theme(legend.position="none", 
            plot.background = element_rect(fill = "transparent",colour = NA), 
            text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ))+
      xlab(paste0("PC1: ",as.integer(1000*pca$sdev[1]^2/sum(pca$sdev^2))/10,"%"))+
      ylab(paste0("PC2: ",as.integer(1000*pca$sdev[2]^2/sum(pca$sdev^2))/10,"%"))
  })
  
  output$pcaPlot34 = renderPlotly({
    dtx=scores()
    pca=dopca()
    dtx$ID= rownames(dtx)
    #dtx$subType[is.na(dtx$subType)] = "NA"
    ggplot(data=dtx) +
      geom_point(aes(x=PC3, y=PC4, colour=subType, text=ID, size =3 ) )+
      getColorScale() +
      theme_dark() + 
      theme(legend.position="none", 
            plot.background = element_rect(fill = "transparent",colour = NA), 
            text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ))+
      xlab(paste0("PC3: ",as.integer(1000*pca$sdev[3]^2/sum(pca$sdev^2))/10,"%"))+
      ylab(paste0("PC4: ",as.integer(1000*pca$sdev[4]^2/sum(pca$sdev^2))/10,"%"))
  })
  
  DoEnrichment = function(dtx){
    if (!input$inpShowCompAnalysis) {
      return(NULL)
    }
    gmt = gmtGeneLists[[input$inpEnrichmentCategory]]
    Category=input$inpEnrichmentCategory
    Source=input$inpEnrichmentCategory
    
    flCacheEnrich = paste0("Cache/ByGene.",input$inpPlotType,".",input$inpGeneLists,".",input$inpEnrichmentCategory,".RData")
    if (file.exists(flCacheEnrich)){
      load(flCacheEnrich)
    }else{
      withProgress(message = "Calculating enrichment", value=0, expr = {
        result=foreach (i = 1:length(gmt$geneset.names), .combine=rbind, .inorder = F) %dopar% {
          incProgress(1/length(gmt$genesets), detail = gmt$geneset.names[i]  )
          
          gsName = gmt$geneset.names[i]
          gns = gmt$genesets[[i]]
          gns = str_split_fixed(gns,",",2)[,1]
          idx = toupper(rownames(dtx)) %in% gns
          xres = NULL
          if (sum(idx) > 5 && sum(idx) < 2000){
            xres=foreach (s = colnames(dtx), .combine=rbind, .errorhandling = "remove", .inorder = F) %do% {
              res = NULL
              X = dtx[idx, s ]
              Y = dtx[!idx, s]
              pGT = wilcox.test( X, Y , alternative =  "greater")
              pLT = wilcox.test( X, Y , alternative =  "less")
              if (pGT$p.value < 0.001){
                res = c(Source,Category,gsName, s, sum(idx), sum(!idx), -10*log10(pGT$p.value), pGT$statistic,"Greater" )
              }
              if (pLT$p.value < 0.001){
                res = c(Source,Category,gsName, s, sum(idx), sum(!idx), -10*log10(pLT$p.value), pLT$statistic,"Less" )
              }
              res
            }
            return(data.frame(xres,stringsAsFactors = F))
          }else{
            return(NULL)
          }
        }
        save(result, file=flCacheEnrich,compress = "gzip")
      })
    }
    return(result)
  }
  

  output$SampInfo = DT::renderDataTable({

    evt = input$`.clientValue-plotly_hover-A`


    if (!is.null(evt)){
      evt = jsonlite::fromJSON( evt)
      dtx=scores()
      idx = as.integer(100*evt$x) == as.integer(100*dtx$PC1) & as.integer(100*evt$y) == as.integer(100*dtx$PC2) | 
        as.integer(100*evt$x) == as.integer(100*dtx$PC3) & as.integer(100*evt$y) == as.integer(100*dtx$PC4) 
      xid = rownames(dtx[idx,])
      dtx=dtMeta
      dtx = datatable(dtx[xid,2:15], style = 'bootstrap',
                      rownames = FALSE,
                      options = list(dom = 't')
      )
      dtx
    }else{
      # dtx=dtMeta
      # dtx = datatable(dtx[1,2:15], style = 'bootstrap',
      #                 rownames = FALSE,
      #                 options = list(dom = 't')
      # )
      # dtx
      NULL
    }
    
    
  })
  
  
  output$PCEnrich = DT::renderDataTable({
    if (!input$inpShowCompAnalysis) {
      return(NULL)
    }
    dtx= dopca()
    dtx=dtx$rotation[,1:4]
    dt=DoEnrichment(dtx)
    
    
    
    if (is.null(dt) | nrow(dt) < 1){
      dtx = datatable(data.frame(PC1=0,PC2=0,PC3=0,PC4=0), style = 'bootstrap',
                      caption = htmltools::tags$caption(
                        style = 'caption-side: bottom; text-align: center;',
                        paste(input$inpEnrichmentCategory, 'No significant enrichments found...'), htmltools::em('')
                      ))
      return(dtx)
    }
    
    colnames(dt) = c("Source","Category","ListName", "Sample","hits", "nonhits","pvalue","tstat","direction")
    dt$pvalue[dt$direction == "Less"] = as.double( dt$pvalue[dt$direction == "Less"])*-1
    dt$pvalue = as.double(dt$pvalue)
    nn = cast(dt[, c("ListName","Sample","pvalue")], ListName~Sample,value = "pvalue", mean)
    nn = data.frame(nn, stringsAsFactors = F) 
    
    
    
    for (x in paste0("PC",1:4)){
      if (x %in% colnames(nn)){
        nn[is.nan(nn[,x]),x] = 0  
      }else{
        nn[,x] = 0
      }
    }
    # nn=nn[,paste0("PC",1:4)]
    # rownames(nn)
    
    dtx = datatable(nn, style = 'bootstrap',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;',
                      paste(input$inpEnrichmentCategory, 'Enrichment:'), htmltools::em('Enrichment by component loading.')
                    ))
    
    dtx %>% 
      formatRound('PC1', 3) %>%
      formatRound('PC2', 3) %>% 
      formatRound('PC3', 3) %>% 
      formatRound('PC4', 3) 
    
    
  })
  
  output$PCLoadings = DT::renderDataTable({
    if (!input$inpShowCompAnalysis) {
      return(NULL)
    }
    dtx=loadings()
    X=c(dtx$PC1,dtx$PC2,dtx$PC3,dtx$PC4)
    X=quantile(X,c(0.01,0.5,0.99) )
    
    
    #dtx = as.integer(dtx*1000)/1000
    dtx= datatable(dtx, style = 'bootstrap',
                   caption = htmltools::tags$caption(
                     style = 'caption-side: bottom; text-align: center;',
                     'Top PC Loadings: ', htmltools::em('Gene weights with statistically significant discrimitating power.')
                   ))
    
    
    
    dtx %>% 
      formatRound('PC1', 4) %>%
      formatRound('PC2', 4) %>% 
      formatRound('PC3', 4) %>% 
      formatRound('PC4', 4)
  })  
  
  
  observe({
    if  (input$inpGeneLists != ""){
      idx=geneLists[, input$inpGeneLists]  
      gns = geneLists$GeneName[idx]
      updateTextInput(session, "lstGenes", label = "Gene list:", value = paste(gns,collapse = "\n")  )
    }
  })
  
})
