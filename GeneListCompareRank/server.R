library(doMC)
library(shiny)
library(ggplot2)
library(reshape)
library(plotly)
library(DT)
library(ggrepel)
require(xlsx)
registerDoMC(20) 
library(stringr)


load("gmtGeneLists.RData")


shinyServer(function(input, output,session) {
  
  dtxGeneLists = reactive({
    inFile <- input$fileGeneMtx
    if (is.null(inFile))
      return(read.table("example.tab",sep = "\t",header = T, stringsAsFactors = F, row.names = 1))
    txt2=readLines(inFile$datapath)
    txt3=str_split(txt2,'\t')
    result=data.frame(do.call(rbind, txt3[2:length(txt3)]), stringsAsFactors = F)
    colnames(result) = txt3[[1]]
    rownames(result) = result[,1]
    result[,2:ncol(result)]
  })
  

  
  dtxGetEnrichment = reactive({
    res=dtxGetEnrichmentFull()
    if (is.null(res)){return(NULL)} 
    if (nrow(res) < 1){return(NULL)}
    idx = res$alternative == "Less"
    res$pval[idx] = res$pval[idx]*-1
    res=cast(res, GmtList~List,  value = "pval", fun.aggregate = mean,fill = 0)
    
  })
  
  
  
  dtxGetEnrichmentFull = reactive({
    dtx=dtxGeneLists()
    if (nrow(dtx) < 3){return(NULL)}
    #dtx=result
    #x=sample(colnames(dtx),1)
    #gmt=GSA.read.gmt(paste0("../DATA/GeneLists/CORUM.gmt"))

    flgmt=input$inpEnrichmentCategory
    gmt=gmtGeneLists[[flgmt]]

    allGenes = rownames(dtx)
    
    withProgress(message=paste("Updating enrichment for:", flgmt),value=0,{
      res=foreach (x = colnames(dtx), .combine = rbind ) %do% { 
        incProgress(1/ncol(dtx),message =  x  )
        foreach (i = 1:length(gmt$genesets), .combine=rbind) %dopar% {
          GmtList = gmt$geneset.names[i]
          gns = gmt$genesets[[i]]
          gns = unlist(str_split_fixed(gns,",",2)[,1])
          if (sum(gns %in% allGenes) > 10){
            gns = gns[gns %in% allGenes]
            idx = allGenes %in% gns
            A = as.double(dtx[idx,x])
            B = as.double(dtx[!idx,x])
          
            wL=wilcox.test(A,B, paired=F,alternative = "less",exact = F)
            wG=wilcox.test(A,B, paired=F,alternative = "greater",exact = F)
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
            t.stat= mean(rnk[1:length(A)]) - mean(rnk[1:length(B)])
            FC=log2( (1+mean(A))/(1+mean(B)) )
            p.val=-10*log10(p$p.value)
            if (is.infinite(p.val)){
              p.val=300
            }
            p.val = as.integer(min(p.val,300)*10)/10.0
            if (p.val > 10){
              List=x
              data.frame(GmtList,List, pval=p.val, t.stat, log2FC=FC, stringsAsFactors = F, alternative)  
            }
          }
        }
      }
      res=data.frame(res, stringsAsFactors = F)
      res$pval = as.double(res$pval)
      res
    })
    
    
  })
  
  output$DE_Enrich = DT::renderDataTable({
    
    
    dt=dtxGetEnrichment()
    doStop = is.null(dt)
    if (!doStop){
      doStop = nrow(dt) <= 1
    }
    if (doStop){
      dtx = datatable(data.frame(p.val=0,log2FC=0), style = 'bootstrap',
                      caption = htmltools::tags$caption(
                        style = 'caption-side: bottom; text-align: center;',
                        paste(input$inpEnrichmentCategory, 'No significant enrichments found...'), htmltools::em('')
                      ))
      return(dtx)
    }
    dtx = datatable(dt, style = 'bootstrap', rownames= F,selection = "single",
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;',
                      paste(input$inpEnrichmentCategory, 'Enrichment:'), htmltools::em('Enrichment by gene list.')
                    ))
    
    dtx
  })
  
  output$DE_EnrichSel = DT::renderDataTable({
    dt=getSelectedEnrichment()
    doStop = is.null(dt)
    if (!doStop){
      doStop = nrow(dt) <= 1
    }
    if (doStop){
      
      dtx = datatable(data.frame(p.val=0,log2FC=0), style = 'bootstrap',
                      caption = htmltools::tags$caption(
                        style = 'caption-side: bottom; text-align: center;',
                        paste(input$inpEnrichmentCategory, 'No significant enrichments found...'), htmltools::em('')
                      ))
      return(dtx)
    }
    dtx = datatable(dt, style = 'bootstrap', rownames= F,selection = "none",
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;',
                      '', htmltools::em('Enrichment detail.')
                    ))
    dtx
  })
  
  output$DE_Params = DT::renderDataTable({
    dt=getSelectedEnrichmentPhyper()
    doStop = is.null(dt)
    if (!doStop){
      doStop = nrow(dt) <= 1
    }
    if (doStop){
      dtx = datatable(data.frame(p.val=0,log2FC=0), style = 'bootstrap',
                      caption = htmltools::tags$caption(
                        style = 'caption-side: bottom; text-align: center;',
                        '','' 
                      ))
      return(dtx)
    }
    dtx = datatable(dt, style = 'bootstrap', rownames= F,selection = "none",
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;',
                      'No selection...','' 
                    ))
    dtx
  })
  
  getSelectedEnrichment = reactive({
    rows = input$DE_Enrich_rows_selected
    if (! is.null(rows)){
      isolate({dtxE=dtxGetEnrichment()})
      List = dtxE[rows[1],"GmtList"]
      #flgmt=paste0("../DATA/GeneLists/",input$inpEnrichmentCategory ,".gmt")
      #gmt=GSA.read.gmt(flgmt)

      flgmt=input$inpEnrichmentCategory
      gmt=gmtGeneLists[[flgmt]]

      idx = which(gmt$geneset.names == List)
      gmtGenes = gmt$genesets[[idx]]
      gmtGenes = unlist(str_split_fixed(gmtGenes,",",2)[,1])
      
      isolate({dtxLists=dtxGeneLists()})
      #gns=unique(dtx[,List])
      #gns=gns[gns %in% allGenes]
      
      
      dtx=data.frame(Genes=gmtGenes, Valid=T , stringsAsFactors = F)
      for (x in colnames(dtxLists)){
        dtx[,x] = dtx$Genes %in% dtxLists[,x]
      }
      dtx
      
    }
    
  })
  
  
  output$volcSelected = renderPlotly({
    dtx=dtxGetEnrichmentFull()
    #load("dtx.Rdata")
    if (is.null(dtx) ){
      return(NULL)
    }
    
    xcols=rainbow( length(unique(dtx$List))+1) 
    p=ggplot(dtx, aes(log2FC, pval, colour=List,label=GmtList)) +
      geom_point(aes(log2FC, pval), alpha=0.75)+
      ylim(c(0,300))+
      scale_colour_manual(values = xcols)+
      theme_dark() +
      theme(plot.background = element_rect(fill = "transparent",colour = NA), 
            text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
            legend.text=element_text(size=16, color="white"),
            legend.background=element_rect(fill = "transparent",colour = NA),
            legend.key = element_rect(fill = "transparent", colour=NA),
            legend.key.size = unit(3,"lines"),
            legend.position = "left"
      )+
      guides(colour = guide_legend(override.aes = list(size=12), title = ""))
    
    
    rows = input$DE_Enrich_rows_selected
    if (! is.null(rows)){
      isolate({dtxE=dtxGetEnrichment()})
      lstName  = dtxE[rows[1],"GmtList"]
      idx = dtx$GmtList == lstName
      p = p + ggtitle(label = lstName)
      p = p + geom_point(data = dtx[idx,], aes(log2FC, pval, colour=List), size=4 ) 
    }
    
    
    p
  })
  
  
  getSelectedEnrichmentPhyper = reactive({
    rows = input$DE_Enrich_rows_selected
    if (! is.null(rows)){
      
      isolate({dtxE=dtxGetEnrichment()})
      isolate({dtxF=dtxGetEnrichmentFull()})
      List = dtxE[rows[1],"GmtList"]
      
      
      flgmt=input$inpEnrichmentCategory
      gmt=gmtGeneLists[[flgmt]]
      
      idx = which(gmt$geneset.names == List)
      idx = dtxF$GmtList == List
      dtxF[idx,]
      
      
      
      
      
    }
    
  })
  
  
  output$downloadWB <- downloadHandler(
    filename = function() { 
      "BCLandEnrichExtract.xlsx"
    },
    content = function(file) {
      
      wb=createWorkbook()
      ws=createSheet(wb,input$inpEnrichmentCategory)
      dtxE=dtxGetEnrichment()
      addDataFrame(dtxE,ws)
      
      rows = input$DE_Enrich_rows_selected
      if (! is.null(rows)){
        List = dtxE[rows[1],"List"]
        ws=createSheet(wb,List)
        addDataFrame(getSelectedEnrichmentPhyper(),ws)
        ws=createSheet(wb,paste0("Genes_",List)  )
        addDataFrame(getSelectedEnrichment(),ws)
      }
      
      
      
      saveWorkbook(wb,file)
    }
  )
  
  
})
