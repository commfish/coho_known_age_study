#IN PROGRESS. This will be removed once code is completely updated. 
# Script 3 - Analysis
#' ---
#' title: "Coho Known Age LDA Analysis"
#' author: "Sara Miller & Justin Priest"
#' date: "`r format(Sys.Date())`"
#' ---

# Email addresses: Sara Miller (sara.miller@alaska.gov) and Justin Priest (justinpriest.ak@gmail.com)


# Run QDA function  
# To perform QDA, one must ensure that the within-group covariance matrices of
# the explanatory variables are homogeneous, a condition that is frequently violated
# with ecological data ;If one the groups defined by the dependent variable has greater 
# dispersion than others, cases will tend to be overclassified in it.


library(here)
library(MASS)
library(vegan)
library(klaR)
library(FNGr) #devtools::install_github("ben-williams/FNGr") 
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
                                          max(coho_scales_bothriv$Q9abs_trans)+0.05, length.out = 100), #adding a little bit for plotting later
                         Q2plus_trans=seq(min(coho_scales_bothriv$Q2plus_trans), 
                                          max(coho_scales_bothriv$Q2plus_trans)+0.1, length.out = 100), #adding a little bit for plotting later
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
  #geom_text(data = coho_scales_berners, aes(x=Q9abs_trans, y=Q2plus_trans, label=Age, color=accuracy)) + 
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Correct"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="grey60") + 
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Incorrect"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="black") +
  #scale_color_manual(values = c("grey60", "black"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(0.09, 0.72), breaks = c(seq(0.1, 0.7, by = 0.1))) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.07, 0.23), breaks = c(seq(0.08, 0.2, by = 0.04))) +
  ggtitle("A) Berners River")+ xlab ("Q9abs (transformed)")+ylab("Q2plus (transformed)") +
  annotate("text", x = 0.19, y=0.13, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.1, y=0.68, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  # annotate("text", x = 0.08, y=0.7, label="B)", family="Times New Roman", colour="black", size=5) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> plot1

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs_trans, y=Q2plus_trans, fill=predval)) +
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Correct"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="grey60") + 
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Incorrect"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="black") +
  # I did it this way so that the black text would be on top and not buried
  #geom_text(data = coho_scales_hughsmith, aes(x=Q9abs_trans, y=Q2plus_trans, label=Age, color=accuracy)) + 
  #scale_color_manual(values = c("grey60", "black"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(0.09, 0.72), breaks = c(seq(0.1, 0.7, by = 0.1))) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.07, 0.23), breaks = c(seq(0.08, 0.2, by = 0.04))) +
  ggtitle("B) Hugh Smith")+ xlab ("Q9abs (transformed)") + ylab("Q2plus (transformed)") +
  annotate("text", x = 0.19, y=0.13, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.1, y=0.68, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  # annotate("text", x = 0.08, y=0.7, label="B)", family="Times New Roman", colour="black", size=5) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> plot2

cowplot::plot_grid(plot1, plot2,  align = "h", nrow = 1, ncol=2) 
ggsave("figures/model_pred1.png", dpi = 500, height = 4 , width = 6, units = "in")


coho_scales_hughsmith %>% group_by(accuracy) %>% tally()
coho_scales_berners %>% group_by(accuracy) %>% tally()
paste0("BR SQRT ARCSINE: ", round((1-117/(1850 + 117))*100, 2), "%")
paste0("HS SQRT ARCSINE: ", round((1-62/(790+62))*100, 2), "%")

paste0("BR NO TRANS: ", round((1-118/(1849 + 118))*100, 2), "%")
paste0("HS NO TRANS: ", round((1-60/(792+60))*100, 2), "%")

paste0("BR SQRT: ", round((1-119/(1848 + 119))*100, 2), "%")
paste0("HS SQRT: ", round((1-102/(750+102))*100, 2), "%")








