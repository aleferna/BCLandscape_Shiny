library(doMC)
library(shiny)
library(ggplot2)
library(reshape)
library(plotly)
library(DT)
library(ggrepel)
require(xlsx)
#require(GSA)
registerDoMC(20) 
library(stringr)

load("RNA.Rdata")
load("gmtGeneLists.RData")



shinyServer(function(input, output,session) {
  
  allGenes = dtRNA$geneName
  
  #allGenes = unlist(str_split(allGenes,pattern="[ ,\t\n]"))
  
  
  dtxGeneLists = reactive({
    
    txt=input$lstGeneLists
    txt2=unlist(str_split(txt,"\n"))
    #txt2=readLines("/Z/ShinyApps/BCLandscape/DATA/RAW/ExtraCategories.tab")
    
    txt3=str_split(txt2,'\t')
    
    result=data.frame(do.call(rbind, txt3[2:length(txt3)]), stringsAsFactors = F)
    colnames(result) = txt3[[1]]
    result
  })
  
  observe({
    
    inFile <- input$fileBGGeneList
    
    if (is.null(inFile))
      return(NULL)
    
    allGenes <<- readLines(inFile$datapath)
    allGenes <<- unlist(unlist(str_split(allGenes,pattern="[ ,\t\n]")))
  })
  
  
  dtxGetEnrichment = reactive({
    res=dtxGetEnrichmentFull()
    if (is.null(res)){return(NULL)} 
    if (nrow(res) < 1){return(NULL)}
    
    cast(res, GmtList~List,  value = "pval", fun.aggregate = mean,fill = 0)
    
  })
  
  
  
  dtxGetEnrichmentFull = reactive({
    
    dtx=dtxGeneLists()
    
    if (nrow(dtx) < 3){return(NULL)}
    
    #dtx=result
    #x=sample(colnames(dtx),1)
    #gmt=GSA.read.gmt(paste0("../DATA/GeneLists/CORUM.gmt"))

    flgmt=input$inpEnrichmentCategory
    gmt=gmtGeneLists[[flgmt]]
    
    withProgress(message=paste("Updating enrichment for:", flgmt),value=0,{
      res=foreach (x = colnames(dtx), .combine = rbind ) %dopar% { 
        gnsClu = unique(dtx[,x])
        gnsClu = gnsClu[gnsClu %in% allGenes]
        foreach (i = 1:length(gmt$genesets), .combine=rbind) %do% {
          incProgress(1/length(gmt$genesets), message = gmt$geneset.names[i]  )
          
          GmtList = gmt$geneset.names[i]
          gns = gmt$genesets[[i]]
          gns = unlist(str_split_fixed(gns,",",2)[,1])
          gns = gns[gns %in% allGenes]
          ctHits=sum(gnsClu %in% gns)
          ctSel=length(gnsClu)
          ctTot=length(allGenes)
          ctLst=length(gns)
          
          exp = 1.0 * ctLst / ctTot 
          exp = ctSel * exp
          log2FC = log2( (1+ctHits)/(1+exp))
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
          pval = as.integer(pval*10)/10
          log2FC = as.integer(log2FC*10)/10
          exp = as.integer(exp*10)/10
          List = x
          data.frame(List,pval, ctHits, exp, log2FC, ctSel, GmtList, ctLst, ctTot, stringsAsFactors = F)
        }
      }
      res=data.frame(res, stringsAsFactors = F)
      res$pval = as.double(res$pval)
      res
    }
    
    )
    
    
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
      
      flgmt=input$inpEnrichmentCategory
      gmt=gmtGeneLists[[flgmt]]
      
      idx = which(gmt$geneset.names == List)
      gmtGenes = gmt$genesets[[idx]]
      gmtGenes = unlist(str_split_fixed(gmtGenes,",",2)[,1])
      
      isolate({dtxLists=dtxGeneLists()})
      #gns=unique(dtx[,List])
      #gns=gns[gns %in% allGenes]
      
      
      dtx=data.frame(Genes=gmtGenes, Valid=gmtGenes %in% allGenes , stringsAsFactors = F)
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