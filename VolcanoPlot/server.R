require(doMC)
require(shiny)
require(ggplot2)
require(reshape)
require(plotly) 
require(DT)
require(ggrepel)
require(stringr)
#require(xlsx)
#require(rJava)
require(openxlsx)
#install.packages("openxlsx")
require(digest)
#library(gage)
#library(pathview)
#library(gridExtra)
require(png)
require(grid) 
require(gridExtra)
require(gtable)
#library(rasterGrob)
#install.packages("rasterGrob")
#registerDoMC(detectCores())

load("DATA/gmtGeneLists.RData")

# dtKeggPaths=read.table("DATA/kegg_entries.tab",sep="\t",quote = "#",stringsAsFactors = F,header = T)
# kegglst=list()
# for (i in 1:nrow(dtKeggPaths)){
#   kegglst[[dtKeggPaths$Name[i]]] = unlist(str_split(dtKeggPaths$Genes[i],","))
# }


load("DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL
load("DATA/RNA.Rdata")
gns = dtRNA$geneName
dtRNA$geneName = NULL
dtRNA = data.frame(2^dtRNA)
rownames(dtRNA) = gns

xcols=intersect(colnames(dtProt),colnames(dtRNA))
dtRNA=dtRNA[,xcols]
dtProt=dtProt[,xcols]

dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
for (i in colnames(dtMeta)){
  dtMeta[,i] = str_trim(dtMeta[,i])
}



PAM590cols=c("Basal"="#E31A1C", "Her2"="#FB9A99","LumA"="#1F78B4","Lum"="#A6CEE3","LumB"="#A6CEE3","Normal"="#33A02C","Normal-like"="#33A02C")

#pca=prcomp(t(dtProt), scale=F)
#save(pca,file="DATA/PCA.RData")
#load("DATA/PCA.RData")

shinyServer(function(input, output,session) {
  
  values <- reactiveValues(doRefresh = 0)
  
  
  
  
  observe ({
    if  (input$btRefresh ){
      isolate({
        values$doRefresh <- runif(1)  
      })
    }
    
  })
  
  observe ({
    xcol=input$lstCats
    x=dtMeta[,xcol]
    x=x[x!=""]
    x=str_trim(x,"both")
    x=sort(unique(x))
    updateSelectInput(session, "lstA",choices=x, selected = x[1])
    updateSelectInput(session, "lstB",choices=x, selected = x[2])
    
    
  })
  
  # output$pcaPlot <- renderPlotly({
  #   if (FALSE){ 
  #     subType=dtMeta[rownames(pca$x),"subtype"]
  #     subTypeProt=dtMeta[rownames(pca$x),"ProteinSubtype"]
  #     dtx = data.frame(subType,subTypeProt, pca$x, group="x",selected="", stringsAsFactors = F)
  #     dtx$ID= rownames(dtx)
  #     group = dtMeta[rownames(pca$x),input$lstCats]
  #     #dtx$group[group %in% input$lstA] = "A" 
  #     #dtx$group[group %in% input$lstB] = "B" 
  #     dtx$group = group
  #     dtx$selected[group %in% input$lstA] = "A"
  #     dtx$selected[group %in% input$lstB] = "B"
  #     
  #     
  #     p = ggplot(data=dtx)
  #     #p=p+geom_point(aes(x=PC1, y=PC2,  colour=selected,   size=75 ) )
  #     #p=p+scale_colour_manual(name = "group",values = PAM590cols)  
  #     p = p + geom_point(aes(x=PC1, y=PC2,  colour=group, text=ID ), size=6 )  +  ggtitle("PCA Analysis")
  #     if (input$lstCats == "subtype"){
  #       p = p + scale_colour_manual(name = "Subtype",values = PAM590cols)  
  #     }
  #     
  #     p = p + geom_text(aes(x=PC1, y=PC2,text = selected, label=selected), fontface = 'bold', color="white")
  #     p = p + guides(colour = guide_legend(override.aes = list(size=10)))
  #     
  #     
  #     p=p+theme_dark()  + theme(plot.background = element_rect(fill = "transparent",colour = NA), 
  #                               text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
  #                               legend.text=element_text(size=16, color="white"),
  #                               legend.background=element_rect(fill = "transparent",colour = NA),
  #                               legend.key = element_rect(fill = "transparent", colour=NA),
  #                               legend.key.size = unit(3,"lines"),
  #                               legend.position = "right"
  #     )
  #     save( p, file= "DATA/pca.RData")
  #   }
  #   
  #   load("DATA/pca.RData")
  #   p
  # })
  # 
  # 
  
  
  getDiffExpr = function(dtx,lstA, lstB){
    
    
    
    
    
    isolate({
      withProgress(message = "Calculating differential expression...",expr = {
        gns=toupper( rownames(dtx))
        dtx= foreach (i = 1:nrow(dtx),.combine = rbind) %do% {
          incProgress(1/nrow(dtx), message=paste("Processing: ",gns[i])   )
          A=unlist(dtx[i,lstA])
          B=unlist(dtx[i,lstB])
          wL=wilcox.test(A,B, paired=F,alternative = "less",exact = T)
          wG=wilcox.test(A,B, paired=F,alternative = "greater",exact = T)
          pL=wL$p.value
          pG=wG$p.value
          if (pL < pG){
            p=wL
            alternative = "Less"
          }else{
            p=wG
            alternative = "Greater"
          }
          rnk=rank(c(A,B))
          t.stat= mean(rnk[lstA]) - mean(rnk[lstB])
          FC=log2( (1+mean(A))/(1+mean(B)) )
          p.val=-10*log10(p$p.value)
          if (is.infinite(p.val)){
            p.val=500
          }
          if (p.val > 10){
            data.frame(gene=gns[i], p.val=p.val, t.stat, log2FC=FC, stringsAsFactors = F, alternative)  
          }
        }
      })
      #dtx=data.frame(dtx, stringsAsFactors = F)
      #colnames(dtx) = c("p.val","alternative","t.stat","log2FC")
      dtx[order(dtx$p.val,decreasing = T), ]
      
    })
  }
  
  # gageDTX  = reactive({
  #   a=input$btRefresh
  #   isolate({
  #     if (input$lstDataType == "RNA"){
  #       dtx =dtRNA  
  #     }else{
  #       dtx =dtProt  
  #     }
  #     
  #     idxA = dtMeta[,input$lstCats] %in% input$lstA
  #     idxB = dtMeta[,input$lstCats] %in% input$lstB
  #     if (sum(idxA) < 3){return(NULL)}
  #     if (sum(idxB) < 3){return(NULL)}
  #     idxA =  dtMeta$id[idxA]
  #     idxB =  dtMeta$id[idxB]
  #     idxA = idxA[idxA %in% colnames(dtx)]
  #     idxB = idxB[idxB %in% colnames(dtx)]
  #     flCache = paste0( input$lstDataType,".",input$lstCats,".",paste0(input$lstA,collapse=","),".vs.",paste(input$lstB,collapse=","))
  #     flCache = digest(flCache,algo="md5")
  #     flCache = paste0("Cache/GAGE.",flCache,".Rdata")
  #     if (file.exists(flCache)){
  #       load(flCache)
  #       return(dtx)
  #     }else{
  #       x=data.matrix(dtx)
  #       idxA = which(colnames(x) %in% idxA)
  #       idxB = which(colnames(x) %in% idxB)
  #       
  #       #idxA=c(2:8)
  #       #idxB=c(9:15)
  #       dtx=gage(x,gsets = kegglst,ref=idxA,samp=idxB, rank.test = T,compare = "unpaired" )
  #       dtxG=data.frame(Name=rownames(dtx$greater),type="Greater",dtx$greater[,1:5], stringsAsFactors = F)
  #       dtxL=data.frame(Name=rownames(dtx$greater),type="Less",dtx$less[,1:5], stringsAsFactors = F)
  #       dtx=rbind(dtxG,dtxL)
  #       rownames(dtx)= NULL
  #       save(dtx,file=flCache)
  #     }
  #     
  #     dtx=dtx[order(dtx$p.val,decreasing = T),]
  #     
  #     dtx
  #   })
  # })
  
  
  diffDtx = reactive({
    if (values$doRefresh == 0){
      return(NULL)
    }
    
    
    isolate({
      if (input$lstDataType == "RNA"){
        dtx =dtRNA  
      }else{
        dtx =dtProt  
      }
      idxA = dtMeta[,input$lstCats] %in% input$lstA
      idxB = dtMeta[,input$lstCats] %in% input$lstB
      if (sum(idxA) < 3){return(NULL)}
      if (sum(idxB) < 3){return(NULL)}
      idxA =  dtMeta$id[idxA]
      idxB =  dtMeta$id[idxB]
      idxA = idxA[idxA %in% colnames(dtx)]
      idxB = idxB[idxB %in% colnames(dtx)]
      flCache = paste0(input$lstDataType,".",input$lstCats,".",paste0(input$lstA,collapse=","),".vs.",paste(input$lstB,collapse=","))
      flCache = digest(flCache,algo="md5")
      flCache = paste0("Cache/",flCache,".Rdata")
      if (file.exists(flCache)){
        load(flCache)
      }else{
        dtx=getDiffExpr(dtx,idxA,idxB)  
        save(dtx,file=flCache)
      }
      
      dtx=dtx[order(dtx$p.val,decreasing = T),]
      dtx
    })
  })
  
  getGmtGenes = function(lsts){
    
    gmt=gmtGeneLists[[input$inpEnrichmentCategory]]
    gns=c()
    for (lst in lsts){
      idx = gmt$geneset.names == lsts
      gns=unique(c(gmt$genesets[idx],gns))
    }
    return(gns)
  }
  
  drawVolcanoPlot = function(theme = "black"){
    #load("/Z/ShinyApps/BCLandscape/Cache/VolcanoPlot/b3a5b0afe7fa6f6051ab645c1b6246ca.Rdata")
    
    
    dtx=diffDtx()
    isolate({
      xtitle = paste( paste0(input$lstA,collapse=","),"vs",paste(input$lstB,collapse=",") )
    })
    
      if (is.null(dtx)) return()
      idx=dtx$p.val > 10
      dtx=dtx[idx,]
      #dtx$alternative[dtx$p.val < 30] =
      #idxBG = !idx & dtx$p.val > 10
      
      
      xmax=max(abs(dtx$log2FC))
      if (is.infinite(xmax) ) xmax=10
      ymax=max(dtx$p.val)
      if (is.infinite(ymax) ) ymax=300
      
      
      #Plot black dots
      
      p=ggplot(dtx[dtx$p.val < 30,], aes(log2FC, p.val)) +
        xlim(c(-1*xmax,xmax))+  ylim(c(10,ymax))+
        #stat_bin2d() + 
        geom_point(aes(log2FC, p.val), alpha=1)+
        guides(colour = guide_legend(override.aes = list(size=12))) 
      
      
      
      if (theme == "black"){
        p=p+theme_dark() +
          theme(plot.background = element_rect(fill = "transparent",colour = NA), 
                text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
                legend.text=element_text(size=16, color="white"),
                legend.background=element_rect(fill = "transparent",colour = NA),
                legend.key = element_rect(fill = "transparent", colour=NA),
                legend.key.size = unit(3,"lines"),
                legend.position = "right")
      }else{
        p=p+theme_light()
      }
      
      #Plot sig dif genes
      p = p + geom_point(data=dtx[dtx$p.val >= 30,], aes(text=gene, colour=alternative), show.legend=T  )  
      
      
      #Plot selected genes
      rows = input$DE_Enrich_rows_selected
      if (length(rows) > 0){
        dtE = dtEnrichment()
        #p=p+ggtitle(dtE$ListName[rows])
        gns=getGmtGenes(dtE$ListName[rows])
        idx = dtx$gene %in% unlist(gns)
        
        p = p + geom_point(data=dtx[idx,], aes(text=gene, colour="Selected") ,show.legend = T )  
        p = p + scale_colour_manual(name="Diff. Expr",
                                    values=c(Less="#67a9cf", Greater="#ef8a62", Selected="yellow")) 
        xtitle =  paste0(xtitle,"\n", dtE$ListName[rows]) 
        #p = p + ggtitle(dtE$ListName[rows])
      }else{
        
        p = p + scale_colour_manual(name="Diff.Expr",
                                    values=c(Less="#67a9cf", Greater="#ef8a62", Selected="yellow")) 
      }
      
      p=p+ggtitle(xtitle) 
      
      rows = input$DiffExpression_rows_selected
      dtx=diffDtx()
      
      #if (length(rows) > 1E9){
      ######### Not working in shiny cloud only on tamarindo, must fix
      if (length(rows) > 0){
        gns=dtx[rows,"gene"]
        idx = dtx$gene %in% gns                              
        p = p + geom_point(data=dtx[idx,], aes(text=gene), colour="green"  )  
        p = p + geom_text_repel(data=dtx[idx,], aes(log2FC, p.val, label=gene, labeltext=gene, text=gene), color="white")  
      }
      p
    
  }
  
  output$volcanoRNA <- renderPlot({
    drawVolcanoPlot()
  }, bg="transparent")
  
  
  
  
  output$DE_Enrich = DT::renderDataTable({
    dt=dtEnrichment()
    doStop = is.null(dt)
    if (!doStop){
      doStop = nrow(dt) <= 1
    }
    if (doStop){
      return(NULL)
    }
    
    
    dtx = datatable(dt, style = 'bootstrap', rownames= F, selection = "single",
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;',
                      paste(input$inpEnrichmentCategory, 'Enrichment:'), htmltools::em('Enrichment by contrast measure.')
                    ))
    dtx %>% 
      formatRound('pvalue', 1) %>%
      formatRound('tstat', 1)  
    
  })
  
  output$DiffExpression = DT::renderDataTable({
    dtx=diffDtx()
    doStop = is.null(dtx)
    if (!doStop){
      doStop = nrow(dtx) <= 1
    }
    if (doStop){
      return(NULL)
    }
    rownames(dtx) = dtx$gene
    dtx= datatable(dtx, style = 'bootstrap', rownames= F,
                   caption = htmltools::tags$caption(
                     style = 'caption-side: bottom; text-align: center;',
                     'Differentially expressed genes: ', htmltools::em('Gene weights with statistically significant discrimitating power.')
                   ))
    dtx %>% 
      formatRound('p.val', 1) %>%
      formatRound('log2FC', 1) %>%
      formatRound('t.stat', 1)  
  })  
  
  
  drawCover = function(){
    txtCite="Please cite:\n   Breast Cancer Landscape Paper, Henrik et al, Nature 2016"
    txtMatMethods=paste("Materials and Methods:\n ", 
                        "Produced by http://lehtiolab.se/toolsx/44-omics-databases-x/breast-cancer-landscape/XXX??\n",
                        "Generated: ",date())
    img <- readPNG("DATA/lehtiolab.png")
    logo <- rasterGrob(img, interpolate=TRUE)
    
    dtx=data.frame(a=c(1,100),b=c(1,100))
    p=qplot(a, b, data = dtx, geom = "blank")+
      theme_bw() +
      xlab("") +
      ylab("") +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    p + annotate("text", x = 50, y = 90, label = txtCite  ) +
      annotate("text", x = 50, y = 50, label = txtMatMethods  )+
      annotate("text", x = 10, y = 25, label = paste0("ListA:",paste(input$lstA, collapse=",") ) )+
      annotate("text", x = 10, y = 20, label = paste0("ListB:",paste(input$lstB, collapse=",") ) )+
      annotate("text", x = 10, y = 15, label = "Statistic: Wilcox Rank test" )+
      annotate("text", x = 10, y = 10, label = paste0("Src Data: ", input$lstDataType) )+
      annotation_custom(logo, xmin = 70, xmax=100, ymin=1, ymax=10)
    
    
  }
  
  
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$lstCats , 'DiffExprPlot.pdf', sep='') },
    content = function(file) {
      #device <- function(..., width, height) grDevices::pdf(..., width = 16, height = 12,bg = "darkgrey")
      #ggsave(file, plot = drawVolcanoPlot("white"), device = device)
      
      pdf(file, onefile = TRUE, width = 12, height = 8)
      grid.arrange( drawCover() )
      grid.arrange( drawVolcanoPlot("white") )
      
      #df=data.frame(gender=c(1,2),freq=c(3,4), fraud=c(32,23))
      df=diffDtx()[1:10,]
      t1 <- tableGrob(df,rows=NULL)
      title <- textGrob(paste("Top 10 differentially expressed gene") )
      padding <- unit(5,"mm")
      table <- gtable_add_rows(t1,   heights = grobHeight(title) + padding, pos = 0)
      table <- gtable_add_grob(table, title, 1, 1, 1, ncol(table))
      grid.newpage()
      grid.draw(table)
      
      df=dtEnrichment()[1:10,]
      t1 <- tableGrob(df,rows=NULL)
      title <- textGrob(paste("Top 10 differentially expressed gene lists") )
      padding <- unit(5,"mm")
      table <- gtable_add_rows(t1,   heights = grobHeight(title) + padding, pos = 0)
      table <- gtable_add_grob(table, title, 1, 1, 1, ncol(table))
      grid.newpage()
      grid.draw(table)
      
      
      
      
      dev.off()
      
      
      
    }
  )
  
  
  output$downloadWB <- downloadHandler(
    filename = function() { 
      "GeneListEnrichCompare.xlsx"
    },
    content = function(file) {
      
      
      wb=createWorkbook()
      
      ws=addWorksheet(wb,"Parameters")
      
      dtx=data.frame(c("ListA","ListB","DataSrc", "Citation"),
                     c(paste(input$lstA,collapse=","),paste(input$lstB,collapse=","),input$lstDataType,"Breast Cancer Landscape Paper, Henrik et al, Nature 2016")
      )
      #dtx=data.frame(a=1:10,b=11:20)
      colnames(dtx) = c("X","Y")
      #dtx=data.frame(ListA=,ListB=paste(input$lstB,sep=","),Src=input$lstDataType)
      writeData(wb,ws,dtx)
      
      ws=addWorksheet(wb,"Gene Diff. Exp.")
      dtx=diffDtx()
      writeData(wb,ws,dtx)
      
      dtx=dtEnrichment()
      ws=addWorksheet(wb,input$inpEnrichmentCategory)
      writeData(wb,ws,dtx)
      
      
      saveWorkbook(wb,file)
    }
  )
  
  
  # output$Kegg_Enrich  = DT::renderDataTable({
  #   print(input$inpKeggCategory)
  #   dtx=gageDTX()
  #   if (!is.null(dtx) ){
  #     
  #     
  #     dtx$p.val =-10*log( dtx$p.val,10)
  #     dtx$q.val =-10*log( dtx$q.val,10)
  #     dtx= datatable(dtx, style = 'bootstrap', rownames= F, selection = "single",
  #                    caption = htmltools::tags$caption(
  #                      style = 'caption-side: bottom; text-align: center;',
  #                      'Differentially expressed genes KEGG (GAGE): ', htmltools::em('')
  #                    ))
  #     dtx %>% 
  #       formatRound('p.val', 1) %>%
  #       formatRound('q.val', 1) %>%
  #       formatRound('p.geomean', 1) %>%
  #       formatRound('stat.mean', 1)
  #     
  #     
  #   }else{
  #     return(NULL)
  #   }
  #   
  #   
  # })
  # 
  # 
  # output$keggPlot <- renderImage({
  #   a=input$btRefreshc
  #   
  #   tblid=input$Kegg_Enrich_rows_selected
  #   dtDTX = gageDTX()
  #   tblid = dtDTX$Name[tblid]
  #   tblid = dtKeggPaths$ID[ dtKeggPaths$Name == tblid]
  #   tblid = str_pad(paste(tblid),5,"left","0")
  #   
  #   outfl="./Untitled.png"
  #   isolate({
  #     if (input$lstDataType == "RNA"){
  #       dtx =dtRNA
  #     }else{
  #       dtx =dtProt
  #     }
  #     idxA = dtMeta[,input$lstCats] %in% input$lstA
  #     idxB = dtMeta[,input$lstCats] %in% input$lstB
  #     if (sum(idxA) < 3){return(NULL)}
  #     if (sum(idxB) < 3){return(NULL)}
  #     idxA =  dtMeta$id[idxA]
  #     idxB =  dtMeta$id[idxB]
  #     idxA = idxA[idxA %in% colnames(dtx)]
  #     idxB = idxB[idxB %in% colnames(dtx)]
  #     A = rowMeans( dtx[,idxA]) - rowMeans(dtx[,idxB])
  #     x <- org.Hs.egSYMBOL2EG
  #     mapped_genes <- mappedkeys(x)
  #     xx <- as.list(x[mapped_genes])
  #     names(A) = xx[names(A)]
  #     outfl=as.integer(runif(1)*1E7)
  #     
  #     
  #     #pathid="04110"
  #     #pathid=str_pad(,5,"left","0")
  #     pathview(gene.data = A, pathway.id = tblid , species = "hsa", out.suffix = outfl, kegg.native = T, kegg.dir = "Cache/" )
  #     outfl=paste0("hsa",tblid,".",outfl,".png")
  #   })
  #   list(src = outfl, contentType = 'image/png', width = 900, height = 600, alt = "")
  # }, deleteFile = TRUE)
  
  dtEnrichment = reactive({
    
    
    if (values$doRefresh == 0){
      return(NULL)
    }
    
    
    isolate({
      dtx=diffDtx()
      if (is.null(dtx)) return(NULL)
      rownames(dtx) = dtx$gene
      flCache = paste0(input$lstDataType,".",input$inpEnrichmentCategory,'.',input$lstCats,".",paste0(input$lstA,collapse=","),".vs.",paste(input$lstB,collapse=","))
      flCache = digest(flCache,algo="md5")
      flCache = paste0("Cache/",flCache,".Rdata")
      #flCache="../Cache/VolcanoPlot/45ceb53b2058aefa370035c4ff281849.Rdata"
      if (file.exists(flCache)){
        load(flCache)
        return(dt)
      }
      
      dtx=dtx[,c('p.val','log2FC','t.stat')]
      
      if (input$inpEnrichmentCategory == ""){
        return(NULL)
      }
      
      gmt=gmtGeneLists[[input$inpEnrichmentCategory]]
      Category=input$inpEnrichmentCategory
      Source=input$inpEnrichmentCategory
      
      
      
      withProgress(message = "Calculating enrichment", value=0, expr = {
        result=foreach (i = 1:length(gmt$geneset.names), .combine=rbind, .inorder = F) %do% {
          incProgress(1/length(gmt$genesets), detail = gmt$geneset.names[i]  )
          gsName = gmt$geneset.names[i]
          gns = gmt$genesets[[i]]
          gns = str_split_fixed(gns,",",2)[,1]
          idx = rownames(dtx) %in% gns
          xres = NULL
          if (sum(idx) > 5 && sum(idx) < 2000){
            #xres=foreach (s = colnames(dtx), .combine=rbind, .errorhandling = "remove", .inorder = F) %do% {
            s = "t.stat"
            res = NULL
            X = dtx[idx, s ]
            Y = dtx[!idx, s]
            pGT = wilcox.test( X, Y , alternative =  "greater")
            pLT = wilcox.test( X, Y , alternative =  "less")
            if (pGT$p.value < 0.001){
              xres = c(Source,Category,gsName, s, sum(idx), sum(!idx), -10*log10(pGT$p.value), pGT$statistic,"Greater" )
            }
            if (pLT$p.value < 0.001){
              xres = c(Source,Category,gsName, s, sum(idx), sum(!idx), -10*log10(pLT$p.value), pLT$statistic,"Less" )
            }
            
            #}
          }
          xres
        }
      })
      
      dt = data.frame(result,stringsAsFactors = F)
      colnames(dt) = c("Source","Category","ListName", "Sample","szList", "nonhits","pvalue","tstat","direction")
      dt$pvalue = as.double(dt$pvalue)
      dt$szList = as.double(dt$szList)
      dt$tstat = as.double(dt$tstat)
      dt$pvalue[dt$direction == "Less"] = dt$pvalue[dt$direction == "Less"]*-1
      dt = dt[,c("ListName", "szList","pvalue","tstat","direction")]
      dt = dt[order(dt$pvalue, decreasing = T),]
      dt = dt[dt$pvalue > 30,]
      save(dt,file=flCache)
      return(dt)
    })
  })
  
  
  observe({
    
    if (values$doRefresh == 0 ){
      invalidateLater(1000, session)  
    }
    
    
    if (values$doRefresh == 0 &  !is.null(input$lstA) ){
      values$doRefresh <- runif(1)  
     
    }
    
  })
  
})




