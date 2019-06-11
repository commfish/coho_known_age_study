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
#devtools::install_github("ben-williams/FNGr")
library(FNGr)
library(processx)
library(fs)
library(cowplot)
library(extrafont)
windowsFonts(Times=windowsFont("TT Times New Roman")) 
theme_set(theme_sleek())

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))

# Transform Variables - Create New Datasets
coho_scales_bothriv <- coho_scales_bothriv %>%
  mutate(Q2plus_trans = (Q2plus), # Leaving this untransformed for now
         Q9abs_trans = (Q9abs)) # Not transformed

coho_scales_berners <- coho_scales_bothriv %>% filter(Location == "BR")
coho_scales_hughsmith <- coho_scales_bothriv %>% filter(Location == "HS")

# Run QDA
fit_qda_final <- qda(as.factor(Age) ~ Location + Q2plus_trans + Q9abs_trans, data=coho_scales_bothriv,
                prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)

# Create blank predicted values (for plotting)
predictedvals <- expand.grid(Q9abs_trans =seq(min(coho_scales_bothriv$Q9abs_trans), 
                                          max(coho_scales_bothriv$Q9abs_trans), length.out = 100), 
                         Q2plus_trans=seq(min(coho_scales_bothriv$Q2plus_trans), 
                                          max(coho_scales_bothriv$Q2plus_trans), length.out = 100),
                         Location = c("HS", "BR"))

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
  scale_color_manual(values = c("grey60", "black"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  ggtitle("Berners River")+ xlab ("Q9abs (transformed)")+ylab("Q2plus (transformed)") +
  annotate("text", x = 0.44, y=0.40, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.32, y=1.3, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.48, y=1.3, label="A)", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none") -> plot1

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs_trans, y=Q2plus_trans, fill=predval)) +
  geom_text(data = coho_scales_hughsmith, aes(x=Q9abs_trans, y=Q2plus_trans, label=Age, color=accuracy)) + 
  scale_color_manual(values = c("grey60", "black"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  ggtitle("Hugh Smith")+ xlab ("Q9abs (transformed)")+ylab("Q2plus (transformed)") +
  annotate("text", x = 0.44, y=0.40, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.32, y=1.3, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.48, y=1.3, label="B)", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none") -> plot2

cowplot::plot_grid(plot1, plot2,  align = "h", nrow = 1, ncol=2) 
ggsave("figures/model_pred.png", dpi = 500, height = 4 , width =6, units = "in")


coho_scales_hughsmith %>% group_by(accuracy) %>% tally()
coho_scales_berners %>% group_by(accuracy) %>% tally()
1-117/(1850 + 117)
1-62/(790+62)



