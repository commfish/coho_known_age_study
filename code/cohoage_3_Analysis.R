#IN PROGRESS. This will be removed once code is completely updated. 
# Script 3 - Analysis
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
library(MASS)
library(vegan)
library(klaR)

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))



# Following approach from: https://www.statmethods.net/advstats/discriminant.html
fit_lda <- lda(Age ~ Length + Q5, data=coho_scales_aukelake,
           na.action="na.omit", CV=TRUE)
fit_lda 
ct <- table(coho_scales_aukelake$Age[!is.na(coho_scales_aukelake$Length)], fit_lda$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))

fit_qda <- qda(as.factor(Age) ~ Location + Length + Q9, data=coho_scales_fulldata,
           prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)
fit_qda
klaR::partimat(as.factor(Age)~Length+Q9,data=coho_scales_fulldata,method="qda") 
klaR::partimat(as.factor(Age)~Length+Q9,data=coho_scales_berners,method="qda") 
klaR::partimat(as.factor(Age)~Length+Q9,data=coho_scales_aukelake,method="qda") 
klaR::partimat(as.factor(Age)~Length+Q9,data=coho_scales_hughsmith,method="qda") 

coho_scales_fulldata <- coho_scales_fulldata %>% 
  mutate(pred_age = predict(fit_qda, coho_scales_fulldata)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))
coho_scales_fulldata %>% group_by(accuracy) %>% tally()


