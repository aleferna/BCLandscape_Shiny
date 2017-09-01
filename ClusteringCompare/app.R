library(doMC)
library(shiny)
library(networkD3)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(shiny)
require(clValid)
require(stringr)

dtMeta=read.table(file="../DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id

load("../DATA/Prot.Rdata")
rownames(dtProt) =dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL
dtProt=log2(t(dtProt))


load("../DATA/colorRanges.Rdata")

dtMeta=read.table(file="../DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id



ui <- shinyUI(fluidPage(
  
  # Application title
  titlePanel("Clustering vs Gene "),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId="cmbRefClust",label="Reference Clustering", choices = colnames(dtMeta), selected="PAM50" ),
      selectInput(inputId = "cmbAlgorithm",label = "Algorithm",choices = c("hierarchical","kmeans", "diana", "fanny", "som", "model", "sota", "pam", "clara","agnes"), selected="hierachical"),
      selectInput(inputId = "cmbDistance", label = "Distance", choices = c("euclidean", "correlation", "manhattan"), selected="correlation"),
      sliderInput("cmbNumClu", "Number of clusters:", min = 3, max = 7, value = 5),
      tags$textarea(id = 'lstCatA', placeholder = 'Clustering A', "" ,  style = "width: 100%; height: 100px")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      #plotOutput("distPlot")
      tabsetPanel(
        tabPanel("Consensus", plotOutput("ConsensusPlot")),
        tabPanel("Sankey", sankeyNetworkOutput("sankPlot",height = "800px"))
      )
    )
  )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
  
  
  updateTextInput(session, "lstCatA",label = "Selected Genes", value = paste(readLines("ProtSubTypeGeneList.txt"), collapse="\n"))
  
  getClusterInfo = function(d,m,n,lst){
    if (lst == ""){
      dtx = dtProt
    }else{
      dtx = dtProt[,colnames(dtProt) %in% lst]
    }
    
    x=clValid(dtx, nClust=n, metric=d, clMethods=m, method="average" ,validation=c("internal"))
    obj=x@clusterObjs[[m]]
    if (paste(n) %in% names(obj)) {
      obj = obj[[paste(n)]]
    }
    obj
  }
  
  clusterInfo = reactive({
    d=input$cmbDistance
    m=input$cmbAlgorithm
    n= input$cmbNumClu
    lst=unlist(str_split(input$lstCatA,"\n"))
    getClusterInfo(d,m,n,lst)
  })
  
  getClassification <- function(obj, d,m,n){
    
    switch(m,
           kmeans = {
             obj$cluster
           },
           hierarchical = {
             cutree(obj,k = n)
           },
           diana = {
             lst=cutree(obj,n)
             names(lst) = obj$order.lab
             lst
           },
           fanny = {
             obj$clustering
           },
           som = {
             lst=obj$unit.classif
             names(lst) = rownames(obj$data)
             lst
           },
           model={
             obj$classification
           },
           sota={
             lst=obj$clust 
             names(lst)=rownames(obj$data)
             lst
           },
           pam={
             obj$clustering 
           },
           clara={
             obj$clustering
           },
           agnes={
             lst=cutree(obj,n)
             names(lst) = obj$order.lab
             lst
           }
    )
  }
  
  
  
  
  calcClusterCompare = function(cluA,cluB){
    grps=list()
    for (grp in unique(cluA)  ){
      grps[[grp]] = names(cluA)[cluA == grp]
    }
    
    for (grp in unique(cluB)  ){
      grps[[grp]] = names(cluB)[cluB == grp]
    }
    
    nodes = data.frame(name=names(grps),idx=1:length(grps) )
    nodes$idx = nodes$idx -1
    rownames(nodes) = nodes$name
    
    grpAll = intersect(names(cluA),names(cluB))
    
    edges = foreach (a = unique(cluA), .combine = rbind) %:%
      foreach (b = unique(cluB), .combine = rbind) %do% {
        grpA = grps[[a]]
        grpB = grps[[b]]
        grpA = intersect(grpA,grpAll)
        grpB = intersect(grpB,grpAll)
        weight = sum(grpA %in% grpB)
        c(nodes[a,"idx"],nodes[b,"idx"],weight)
      }
    edges=data.frame(edges, stringsAsFactors = F, row.names = 1:nrow(edges))
    colnames(edges) = c("source","target","weight")
    edges=edges[edges$weight > 0,]
    list(nodes=nodes,edges=edges)
  }
  
  
  output$sankPlot <- renderSankeyNetwork({
    obj=clusterInfo()
    d=input$cmbDistance
    m=input$cmbAlgorithm
    n= input$cmbNumClu
    cA = getClassification(obj, d,m,n)
    cx = paste0(input$cmbAlgorithm, cA)
    names(cx) = names(cA)
    cA=cx
    
    
    cB = dtMeta[,input$cmbRefClust]
    names(cB) = dtMeta$id
    cB = cB[cB != ""]
    res=calcClusterCompare(cA,cB)
    
    sankeyNetwork(Links = res$edges, Nodes = res$nodes, Source = "source",
                  Target = "target", Value = "weight", NodeID = "name",fontSize = 14)
    
    #cB=input$lstRefCat
  })
  output$ConsensusPlot = renderPlot({
    lst = unlist(str_split(input$lstCatA,"\n"))
    dtx = t(dtProt[,colnames(dtProt) %in% lst])
    
    dtx = dtx[,!colnames(dtx) %in% c("OSL2U.0439T1","OSL2U.0555T1","OSL2U.0523T1")]
    
    
    #dtx = t(dtProt)
    
    #dtx = sweep(dtx,1, apply(dtx,1,median,na.rm=T))
    
    results = ConsensusClusterPlus(data.matrix( dtx), maxK=7,reps=50,pItem=0.8,pFeature=1,
                                   title= "./X", clusterAlg="hc",distance="spearman",seed=1262118388.71279, plot="png")
    
    ConsPlot = results[[ input$cmbNumClu ]]
    xcol=rainbow(input$cmbNumClu)
    names(xcol) = paste0("CSNS",1:input$cmbNumClu)
    
    
    haT = HeatmapAnnotation(df = data.frame(ProtSub9506= dtMeta[colnames(dtx),"ProtSub9506" ]),  col=list(ProtSub9506 = xcolorRanges$ProtSub9506))
    haL = rowAnnotation(df = data.frame(PAM50= dtMeta[colnames(dtx),"PAM50" ]),  col=list(PAM50 = xcolorRanges$PAM50))
    haB =  HeatmapAnnotation(df = data.frame(Consensus=paste0("CSNS",ConsPlot$consensusClass)),  col=list(Consensus=xcol))
    
    
    X=results[[input$cmbNumClu]][["consensusMatrix"]]
    rownames(X) = colnames(dtx)
    
    Heatmap(X, bottom_annotation = haB, top_annotation = haT, row_names_side = "left" ) + haL
    
    
  }, width = 1000, height = 600)
  
})

# Run the application 
shinyApp(ui = ui, server = server)

