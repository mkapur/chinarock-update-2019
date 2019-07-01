---
  title: "Status of China rockfish off the U.S. Pacific Coast in 2015"
author: ''
date: ''
output:
  pdf_document:
  fig_caption: yes
highlight: haddock
includes:
  before_body: Titlepage.tex
in_header: header.tex
keep_tex: yes
latex_engine: xelatex
number_sections: yes
toc: yes
toc_depth: 4
html_document:
  toc: yes
word_document: default
documentclass: article
fontsize: 12pt
geometry: margin=1in
csl: CJAFS.csl
bibliography: chinarockfish2015.bib
---
  ```{r, echo=FALSE, message=FALSE, warning=FALSE}
##read in necessary libraries
# library(r4ss)
library(xtable)
library(ggplot2)
library(reshape2)
options(xtable.comment = FALSE)  #turns off xtable comments
options(scipen=999) # turn off scientific notation
#########################################################################
#CHANGE ME
FirstYR=2006  # minimum year for tables
LastYR=2015  # maximum year for tables
#########################################################################

##############################################################################################################
#### loads workspace image of the BASE MODEL(S) from SS_output in r4ss    
load("./r4ss/China_SS_output2019.RData")

#number of workspaces
n_models = 3
#check the structure of the data here to make sure it contains "mod," the SS_output variable
##############################################################################################################
```

