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

shinyUI(fluidPage(
  titlePanel("Cluster Compare"),
  sankeyNetworkOutput ("sankPlot",height = "800px",width="100%"),
  shiny::sliderInput("ctClus","Number of clusters",4,50,10,step=1)
))
