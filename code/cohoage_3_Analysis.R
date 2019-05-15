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
#MASS package LDA and QDA or cluster analysis?




library(here)
library(MASS)
library(vegan)
library(klaR)

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))


# Transform Variables - Create New Datasets
coho_scales_fulldata <- coho_scales_fulldata %>%
  mutate(Q2plus_trans = asin(sqrt(Q2plus)),
         Q9abs_trans = asin(sqrt(Q9abs)))

coho_scales_berners <- coho_scales_fulldata %>% filter(Location == "BR")
coho_scales_hughsmith <- coho_scales_fulldata %>% filter(Location == "HS")




# Run QDA
fit_qda_final <- qda(as.factor(Age) ~ Location + Q2plus_trans + Q9abs_trans, data=coho_scales_fulldata,
                prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)

# Create blank predicted values (for plotting)
predictedvals <- expand.grid(Q9abs_trans =seq(min(coho_scales_fulldata$Q9abs_trans), 
                                          max(coho_scales_fulldata$Q9abs_trans), length.out = 100), 
                         Q2plus_trans=seq(min(coho_scales_fulldata$Q2plus_trans), 
                                          max(coho_scales_fulldata$Q2plus_trans), length.out = 100),
                         Location = c("AL", "HS", "BR"))

predictedvals$predval <- predict(fit_qda_final, predictedvals)$class

coho_scales_berners <- coho_scales_berners %>% 
  mutate(pred_age = predict(fit_qda_final, coho_scales_berners)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))

coho_scales_hughsmith <- coho_scales_hughsmith %>% 
  mutate(pred_age = predict(fit_qda_final, coho_scales_hughsmith)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "BR"), aes(x=Q9abs_trans, y=Q2plus_trans, fill=predval)) +
  geom_text(data = coho_scales_berners, aes(x=Q9abs_trans, y=Q2plus_trans, label=Age, color=accuracy)) + 
  scale_color_manual(values = c("black", "red"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("#89b6ff", "#ff9260")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  ggtitle("Berners River")

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs_trans, y=Q2plus_trans, fill=predval)) +
  geom_text(data = coho_scales_hughsmith, aes(x=Q9abs_trans, y=Q2plus_trans, label=Age, color=accuracy)) + 
  scale_color_manual(values = c("black", "red"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("#89b6ff", "#ff9260")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  ggtitle("Hugh Smith")



