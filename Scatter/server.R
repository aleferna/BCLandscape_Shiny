library(gridExtra)
library(shiny)
library(stringr) 
library(ggplot2)
library(plotly)
library(reshape)
library(png)
library(grid)
library(RColorBrewer)
library(data.table)
library(DT)
#devtools::install_github('hadley/ggplot2')
load("DATA/DT45.Rdata")
load("DATA/corMtx.Rdata")

#xcolorRanges$CoTC = xcolorRanges$ProtSub9506

shinyServer(function(input, output,session) {
  
  
  getX = function(){
    X=switch( input$cmbFieldTypeX, 
              Protein = dtProt,
              mRNA =  dtRNA,
              CNA = dtCNA,
              miRNA = dtmiRNA,
              Metabolites = dtMetabolites
    )
    
    if(input$cmbFieldX %in% X$ID)
      return(X[X$ID == input$cmbFieldX,])
    else
      return(NULL)
  }
  
  getY = function(){
    Y=switch( input$cmbFieldTypeY, 
              Protein = dtProt,
              mRNA =  dtRNA,
              CNA = dtCNA,
              miRNA = dtmiRNA,
              Metabolites = dtMetabolites
    )
    if(input$cmbFieldY %in% Y$ID)
      return(Y[Y$ID == input$cmbFieldY,])
    else
      return(NULL)
  }
  
  
  
  
  
  drawXY = function(thm="dark"){
    #p=qqplot(0,0)
    X = getX()
    Y = getY()
    
    
    if (!is.null(X) && !is.null(Y)) {
      
      
      
      Xlab = switch( input$cmbFieldTypeX, 
                     Protein = paste0(input$cmbFieldX , " (Protein)"),
                     mRNA =  paste0(input$cmbFieldX , " (RNA)"),
                     miRNA = paste0(input$cmbFieldX , " (miRNA)") ,
                     CNA = paste0(input$cmbFieldX , " (CNA)"),
                     Metabolites = paste0(input$cmbFieldX , " \n(Metabolites)")
      )
      
      
      Ylab =  switch( input$cmbFieldTypeY, 
                      Protein = paste0(input$cmbFieldY, " (Protein)"),
                      mRNA =  paste0(input$cmbFieldY , " (RNA)"),
                      miRNA = paste0(input$cmbFieldY , " (miRNA)"),
                      CNA = paste0(input$cmbFieldY , " (CNA)"),
                      Metabolites = paste0(input$cmbFieldY , " \n(Metabolites)")
      )
      
      ov = intersect(colnames(X) , colnames(Y) ) 
      ov = ov[grep(ov,pattern = "OSL") ]
      
      
      X=t(X[1,ov])
      Y=t(Y[1,ov])
      dtx=data.frame(X,Y)
      colnames(dtx) = c("X","Y")
      
      dtx$subType = dtMeta[ov,input$cmbColorBy]
      dtx$subType[is.na(dtx$subType)] = "Unknown" 
      dtx$Desc = rownames(dtx)
      
      
      p=ggplot(data = dtx, aes(x = X, y = Y)) +
        geom_point(aes(text = Desc, color=subType), size = 4) +
        labs(x=Xlab , y=Ylab) + 
        scale_colour_manual(name = input$cmbColorBy, values = xcolorRanges[[input$cmbColorBy]]) 
      
      if (thm!="white"){
        #theme(plot.margin=unit(c(10,10,10,10),"mm")) +
        p=p+
          theme_dark()+  
          theme(legend.position="bottom", 
                plot.background = element_rect(fill = "transparent",colour = NA), 
                legend.background= element_rect(fill = "transparent",colour = NA), 
                text = element_text(size=20, face="bold", 
                                    #margin = margin(10, 0, 10, 0), 
                                    color="white" ))
        
      }
      p
    }else{
      ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100)
    }
  }
  
  
  
  output$XYPlot <- renderPlotly({
    
    p=drawXY("black")
    p = ggplotly(p, tooltip=c("text")) %>% 
      config(displayModeBar = F)
    p
    
  })
  
  observe({
    gnsX = switch( input$cmbFieldTypeX, 
                   Protein = dtProt$ID,
                   mRNA =  dtRNA$ID,
                   CNA = dtCNA$ID,
                   miRNA = dtmiRNA$ID,
                   Metabolites = dtMetabolites$ID
    )
    if (is.null(gnsX)) return(NULL)
    isolate({
      if (input$cmbFieldX %in% gnsX)
        selx = input$cmbFieldX
      else
        selx = sample(gnsX,1)
      
    })
    
    withProgress(message = 'Loading fields', value = 0, {
      updateSelectizeInput(session,inputId = "cmbFieldX",choices = gnsX, selected=selx, server = T)
    })
  })
  
  
  observe({  
    gnsY = switch( input$cmbFieldTypeY, 
                   Protein = dtProt$ID ,
                   mRNA =  dtRNA$ID,
                   CNA = dtCNA$ID,
                   miRNA = dtmiRNA$ID,
                   Metabolites = dtMetabolites$ID
    )
    
    if (is.null(gnsY)) return(NULL)
    
    isolate({
      if (input$cmbFieldY %in% gnsY)
        sely = input$cmbFieldY
      else
        sely = sample(gnsY,1)
      
      withProgress(message = 'Loading fields', value = 0, {
        updateSelectizeInput(session,inputId = "cmbFieldY",choices = gnsY, selected=sely, server = T)
      })
      
    })
    
  })
  
  drawCover = function(){
    txtCite="Please cite:\n   Breast Cancer Landscape Paper, Henrik et al, Cell 2017"
    txtMatMethods=paste("Materials and Methods:\n Produced by http://www.breastcancerlandscape.org/ \nGenerated: ",date())
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
    filename = function() { paste("BCLandscape.Scatter",input$inpPlotType,input$cmbColorBy,input$cmbFieldX,input$cmbFieldY,"pdf",sep=".")  },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = 10, height = 5)
      #ggsave(file, plot = , device = device)
      
      
      pdf(file, onefile = TRUE,width = 8,height = 4)
      grid.arrange(drawCover())
      grid.arrange( drawXY("white") )
      dev.off()
      
      
    }
  )
  
  observe({
    idx= input$HighCor_rows_selected
    
    
    if (is.null(idx)) return(0)
    
    X=getCor()
    if (nrow(X) < 1) return(0)
    
    isolate({
      gnsY = switch( input$cmbFieldTypeY, 
                     Protein = dtProt$ID ,
                     mRNA =  dtRNA$ID,
                     CNA = dtCNA$ID,
                     miRNA = dtmiRNA$ID,
                     Metabolites = dtMetabolites$ID
      )
      
      
      if (X$Type[idx] != input$cmbFieldTypeY ){
        updateSelectizeInput(session,inputId = "cmbFieldY",     choices = gnsY,  selected=X$Field[idx], server = T)
        updateSelectizeInput(session,inputId = "cmbFieldTypeY", choices = c("Protein", "mRNA","CNA","miRNA", "Metabolites"), selected=X$Type[idx], server = T)  
        updateSelectizeInput(session,inputId = "cmbFieldY",     choices = gnsY,  selected=X$Field[idx], server = T)
      }else{
        updateSelectizeInput(session,inputId = "cmbFieldY",     choices = gnsY,  selected=X$Field[idx], server = T)
      }
      
      
    })
    
  })
  
  getCor = reactive({
    idxA = input$cmbFieldX == corMtx$FieldA & input$cmbFieldTypeX == corMtx$TypeA
    idxB = input$cmbFieldX == corMtx$FieldB & input$cmbFieldTypeX == corMtx$TypeB 
    
    A = corMtx[idxA,c("TypeB","FieldB","value")]
    colnames(A) = c("Type","Field","value")
    
    B = corMtx[idxB,c("TypeA","FieldA","value")]
    colnames(B) = c("Type","Field","value")
    X=rbind(A,B)
    
    X=X[order(-abs(value))]
  
    X
  })
  
  output$HighCor = renderDataTable({
    X = getCor()
    colnames(X) = c("Type","Field","Pearson")
    X = datatable(X, style = 'bootstrap', rownames= F, selection = "single",
                  options = list(searching=FALSE, serverSide=TRUE, info=FALSE, pagingType="simple", 
                                 paging=T, pageLength=7, lengthChange=F),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: center;',
                    "",
                    htmltools::em(paste0("High (Anti)Correlated to ",input$cmbFieldX, " ", input$cmbFieldTypeX  ))
                  ))
    X %>% 
      formatRound('Pearson', 2) 
  })
  
})



