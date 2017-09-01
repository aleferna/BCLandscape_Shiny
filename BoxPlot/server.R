library(shiny)
library(doMC)
library(ggplot2)
library(reshape)
library(plotly)
library(DT)
library(ggrepel)
library(stringr)
library(gridExtra)
library(png)
library(grid)
library(gridExtra)
library(gtable)

load("DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL
load("DATA/RNA.Rdata")
rownames(dtRNA) = dtRNA$geneName
dtRNA$geneName = NULL
xcols=intersect(colnames(dtProt),colnames(dtRNA))
dtRNA=dtRNA[,xcols]
dtProt=log2(dtProt[,xcols])


dtMeta=read.table(file="DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
for (i in colnames(dtMeta)){
  dtMeta[,i] = str_trim(dtMeta[,i])
}
load("DATA/colorRanges.Rdata")

#load("../DATA/colorRanges.Rdata")
#xcolorRanges$CoTC = xcolorRanges$ProtSub9506
#xcolorRanges$CoCT = NULL
#save(xcolorRanges,file="../DATA/colorRanges.Rdata")

thm=theme_dark()  + 
  theme(plot.background = element_rect(fill = "transparent",colour = NA), 
        text = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0), color="white" ),
        legend.text=element_text(size=16, color="white"),
        legend.background=element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent", colour=NA),
        legend.key.size = unit(3,"lines"),
        legend.position = "none")



shinyServer(function(session,input, output) {

  observe({
    withProgress(message = 'Loading fields', value = 0, {
      gns = unique(c(rownames(dtRNA),rownames(dtProt)))  
      
      
      
      
      
      updateSelectizeInput(session = session, inputId = "lstSelectedGene",label="Select gene:", choices = gns, selected = "MET",server = T  )
      
    })
  })
  
  
  drawBoxPlot = function(isPDF=FALSE){
    gn=input$lstSelectedGene
    cls=input$cmbColorBy
    
    if (gn == "" | cls == "") {
      p = ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlab("") + ylab("") +thm
      return(p)
    } 
    
    dtxA=NULL
    if (gn %in% rownames(dtRNA) ){
      y=t(dtRNA[gn,])
      
      dtxA=data.frame(y=y,x=dtMeta[rownames(y),cls], src="RNA",stringsAsFactors = F)
    }
    
    dtxB=NULL
    if(gn %in% rownames(dtProt)){
      y=t(dtProt[gn,])
      dtxB=data.frame(y=y,x=dtMeta[rownames(y),cls], src="Protein",stringsAsFactors = F)
    }
    dtx=rbind(dtxA,dtxB)
    
    colnames(dtx) = c("y","group","src")
    dtx=dtx[!is.na(dtx$group), ]
    p=ggplot(data = dtx, aes(x=group, y=y)) + 
      ggtitle(paste(gn,"by",cls))+
      geom_boxplot(aes(fill=group)) +
      guides(colour = guide_legend(title=cls, override.aes = list(size=10))) + 
      xlab("") + 
      ylab("")+ facet_grid(src ~ .)+
      scale_fill_manual(name = cls, values = xcolorRanges[[cls]]  ) 
    if (!isPDF){
      p=p+thm
      
    }
    
    p
  }
  
  output$BoxPlot <- renderPlotly({
    drawBoxPlot()
  })
  
  getDiffTable = function(){
    
    gn=input$lstSelectedGene
    cls=input$cmbColorBy
    emptydtx = datatable(data.frame(), style = 'bootstrap',
                         caption = htmltools::tags$caption(
                           style = 'caption-side: bottom; text-align: center;',
                           paste(cls, 'No significant contrasts found...'), htmltools::em('')
                         ))
    if (gn == "" | cls == "") {
      #return(emptydtx)
      return(NULL)
    }
    
    
    xres=foreach (tp = c("RNA","Protein"), .combine=rbind) %do% {
      grps=unique(dtMeta[,cls])
      grps=grps[!is.na(grps)]
      grps=grps[grps!=""]
      foreach (a = grps, .combine=rbind) %do% {
        foreach (b = grps,  .combine=rbind) %do% {
          res=NULL
          if (a > b){
            if(tp == "RNA"){
              dtx=dtRNA  
            }else{
              dtx=dtProt 
            }
            if (gn %in% rownames(dtx)){
              A = dtMeta$id[ dtMeta[,cls] == a]
              A = intersect(colnames(dtx),A)
              A = unlist(dtx[gn,A])
              
              B = dtMeta$id[ dtMeta[,cls] == b]
              B = intersect(colnames(dtx),B)
              B = unlist(dtx[gn,B])
              
              
              p=wilcox.test(A,B)
              if (p$p.value < 0.05){
                res=c(a,b, length(A), length(B), p$p.value,tp)  
              }
            }
            res
          }
        }
      }
    }
    if (is.null(xres))return(NULL)
    if (nrow(xres) < 1)return(NULL)
    
    xres=data.frame(xres,stringsAsFactors = F)
    colnames(xres) = c("A","B","ctA","ctB","pval","src")
    xres$pval=as.double(xres$pval)
    #xres$p.val_code = -1*as.integer(log10(xres$p.val))
    #xres$p.val_code = str_dup("*",xres$pval_code)
    xres$pval <- format(xres$pval, scientific = T,digits = 2)
    xres
  }
  
  output$tblDiff = DT::renderDataTable({
    xres=getDiffTable()
    if (is.null(xres)) return(NULL)
    
    datatable(xres, style = 'bootstrap', rownames= F,
              caption = htmltools::tags$caption(
                style = 'caption-side: bottom; text-align: center;',
                paste(input$cmbColorBy, 'Difference between groups:'), htmltools::em('Statistical significance (two sided Wilcox) of the difference between types')
              )) 
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
      annotation_custom(logo, xmin = 70, xmax=100, ymin=1, ymax=10)
    
    
  }
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 
      paste('BoxPlot',input$lstSelectedGene,'by',input$cmbColorBy,'pdf',sep=".") 
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = 14, height = 10,bg = "darkgrey")
      
      pdf(file, onefile = TRUE)
      grid.arrange( drawCover() )
      grid.arrange( drawBoxPlot(T) )
      
      #df=data.frame(gender=c(1,2),freq=c(3,4), fraud=c(32,23))
      df=getDiffTable()
      t1 <- tableGrob(df,rows=NULL)
      title <- textGrob(paste("Differential expression of \n ",input$lstSelectedGene,'by',input$cmbColorBy) )
      padding <- unit(5,"mm")
      table <- gtable_add_rows(t1,   heights = grobHeight(title) + padding,
                               pos = 0)
      table <- gtable_add_grob(table, title, 1, 1, 1, ncol(table))
      
      grid.newpage()
      grid.draw(table)
      
      dev.off()
      
      
    }
  )
  
  
})
