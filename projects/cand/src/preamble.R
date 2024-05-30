rm(list = ls())
options(bookdown.theorem.preamble = FALSE) #https://bookdown.org/yihui/bookdown/markdown-extensions-by-bookdown.html
knitr::opts_chunk$set(packages = c("amsmath", "amsfonts", "amssymb", "dsfont","pstricks"))
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
knitr::opts_chunk$set(fig.pos = "H",fig.align = "center",fig.width = 4,fig.height = 4)
# Set working directory
setwd(dirname(rstudioapi::documentPath()))
# Libraries
library(ggplot2)
library(dplyr)
library(latex2exp)
library(sysfonts)
library(showtext)
library(kableExtra)
library(ggpubr)
library(reshape2)
#remotes::install_github("bvancilku/kubrand")
library(kubrand)
library(fields)
# Python
library(reticulate)
font_add_google("Crimson Pro","crimson")
font_add_google("Tinos","tinos")
showtext_auto() # listen!
# LaTeX native Tikz
# https://matthew-parker.rbind.io/post/2018-08-21-easy-latex-titles/
# Theme
theme_custom <- function(size = 1){ 
  # https://rpubs.com/mclaire19/ggplot2-custom-themes
  font <- "crimson"   #assign font family up front
  font <- "serif"
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      #panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 12*size,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 10*size,
        hjust = 0,                #left align
        vjust = 2),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 10*size,                 #font size
        hjust = 0),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10*size),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 10*size),                #font size
      
      legend.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 10*size),
      legend.title = element_text(              #axis text
        family = font,            #axis famuly
        face = 'bold',
        size = 10*size),
      legend.position = "bottom",
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0),
      legend.box.spacing = unit(0, "pt")
    )
}
# https://en.wikipedia.org/wiki/Stochastic_differential_equation#Linear_SDE:_general_case

#Tables
library(kableExtra)
library(xtable)

xtable2kable <- function(x) {
  out <- capture.output(print(x, table.placement = NULL))[-(1:2)]
  out <- paste(out, collapse = "\n")
  structure(out, format = "latex", class = "knitr_kable")
}

wrapFig_plot <- function(
    plot = NULL, filename, caption = NULL, label = NULL,
    width= 3,height=3,pos="r",
    doc_width = 0.5,res = 3,innerwidth = 0.9
) {
  if (!is.null(plot)) {
    png(filename,width = width*res*100,height = height*res*100, units = "px",res = res*100)
    print(plot)
    dev.off()
  }
  # Remenber result = 'asis' in chunck!
  cat(paste0("
  \\begin{wrapfigure}{",pos,"}{",doc_width,"\\textwidth}
    \\begin{center}
      \\includegraphics[width=",round(doc_width*innerwidth,digits = 2),"\\textwidth]{",filename,"}
    \\end{center}",ifelse(is.null(caption),"",paste0(
      "\\caption{",caption,"}",ifelse(is.null(label),"",paste0("\\label{fig:",label,"}")))),"
  \\end{wrapfigure}
  "))
}

# Get functions
source("src/functions.R")