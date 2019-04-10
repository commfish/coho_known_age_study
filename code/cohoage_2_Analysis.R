#IN PROGRESS. This will be removed once code is completely updated. 
# Script 2 - Analysis
#' ---
#' title: "Coho Known Age LDA Analysis"
#' author: "Sara Miller & Justin Priest"
#' date: "`r format(Sys.Date())`"
#' ---

# Email addresses: Sara Miller (sara.miller@alaska.gov) and Justin Priest (justinpriest.ak@gmail.com)


#Run LDA function  
#To perform LDA, one must ensure that the within-group covariance matrices of
#the explanatory variables are homogeneous, a condition that is frequently violated
#with ecological data ;If one the groups defined by the dependent variable has greater 
#dispersion than others, cases will tend to be overclassified in it.

#function betadisper in package vegan
#MASS packake LDA and QDA or cluster analysis?




library(here)
library (MASS)
library (vegan)

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))
system.time(source(here::here("code/OLD_LDA_SEM_Draft.R")))
