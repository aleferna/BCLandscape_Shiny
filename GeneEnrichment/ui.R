
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
require(d3heatmap)
require(reshape)
library(shiny)
library(stringr)

lst=list.files("Results/", pattern = "*.tab", full.names = F)
lst=str_replace( lst,".Discordance.tab","")
lst=str_replace( lst,".Protein.tab","")
lst=str_replace( lst,".RNA.tab","")
lst=str_replace( lst,".Random.tab","")
lst=sort(unique(lst))

shinyUI(fluidPage(theme="../www/bootstrap.min.css",
  
  titlePanel("GeneList Enrichment Viewer"),
  fluidRow(
    column(4,selectInput(inputId="lstListSelection", label="Select list category:", choices = lst,selected = "Hallmarks" )),
    column(8,radioButtons(inputId="inpSearchType",label="Discordance Type:",choices=c("Protein","RNA","Discordance","Random") ,selected = "Discordance", inline = T))
  ),
  d3heatmapOutput("distPlot",height = "800px" )
))
