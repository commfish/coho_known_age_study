#IN PROGRESS. This will be removed once code is completely updated. 
# Script 3 - Analysis
#' ---
#' title: "Coho Known Age LDA Analysis"
#' author: "Sara Miller & Justin Priest"
#' date: "`r format(Sys.Date())`"
#' ---

# Email addresses: Sara Miller (sara.miller@alaska.gov) and Justin Priest (justinpriest.ak@gmail.com)
library(here)
#library(MASS)
library(FNGr) #devtools::install_github("ben-williams/FNGr") 
library(cowplot)
library(extrafont)
library(mgcv)
windowsFonts(Times=windowsFont("TT Times New Roman")) 
theme_set(theme_sleek())
options(scipen=999)

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))
source(here::here("code/functions.R"))




# Run Binomial Regression
binomfit <- glm((Age-1)~Location+Q2+Q9abs, data=coho_scales_bothriv, family="binomial")
summary(binomfit)



coho_scales_berners <- coho_scales_bothriv %>% 
  filter(Location == "BR") %>%
  dplyr::select(-c(Q10:Q44))

coho_scales_berners <- coho_scales_berners %>% 
  mutate(pred_binom = predict(binomfit, coho_scales_berners, type = "response"),
         pred_age = ifelse(pred_binom > 0.5, 2, 1),
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))


coho_scales_hughsmith <- coho_scales_bothriv %>% 
  filter(Location == "HS") %>%
  dplyr::select(-c(Q10:Q44))

coho_scales_hughsmith <- coho_scales_hughsmith %>% 
  mutate(pred_binom = predict(binomfit, coho_scales_hughsmith, type = "response"),
         pred_age = ifelse(pred_binom > 0.5, 2, 1),
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))

# Evaluate accuracy
coho_scales_berners %>% group_by(Age, accuracy) %>% tally()
coho_scales_hughsmith %>% group_by(Age, accuracy) %>% tally()




predictedvals <- expand.grid(Q9abs=seq(min(coho_scales_bothriv$Q9abs), max(coho_scales_bothriv$Q9abs), length.out = 100), 
                             Q2=seq(min(coho_scales_bothriv$Q2), max(coho_scales_bothriv$Q2), length.out = 100),
                             Location=c("BR", "HS"))



predictedvals$predval.binom <- predict(binomfit, predictedvals, type = "response")
predictedvals <- predictedvals %>% 
  mutate(predval = as.factor(ifelse(predval.binom > 0.5, 2, 1)))


ggplot() + 
  #geom_tile(data = predictedvals %>% filter(Location == "BR"), aes(x=Q9abs, y=Q2, fill=predval)) +
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Correct"), aes(x=Q9abs, y=Q2, label=Age), color="grey60") + 
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Incorrect"), aes(x=Q9abs, y=Q2, label=Age), color="black") +
  geom_contour(data = predictedvals %>% filter(Location == "BR"), aes(x=Q9abs, y=Q2, z=predval.binom), color="black", breaks = 0.5) +
  geom_contour(data = predictedvals %>% filter(Location == "BR"), aes(x=Q9abs, y=Q2, z=predval.binom), color="black", lty=2, breaks = c(0.005, 0.995)) +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(0.09, 0.72), breaks = c(seq(0.1, 0.7, by = 0.1))) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.07, 0.23), breaks = c(seq(0.08, 0.2, by = 0.04))) +
  ggtitle("A) Berners River")+ xlab ("Circulus 3–Circulus 8 (mm)")+ylab("Circulus 3–Scale Edge  (mm)") +
  annotate("text", x = 0.19, y=0.13, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.1, y=0.68, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> BRmodelpred

ggplot() + 
  #geom_tile(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs, y=Q2, fill=predval)) +
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Correct"), aes(x=Q9abs, y=Q2, label=Age), color="grey60") + 
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Incorrect"), aes(x=Q9abs, y=Q2, label=Age), color="black") +
  geom_contour(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs, y=Q2, z=predval.binom), color="black", breaks = 0.5) +
  geom_contour(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs, y=Q2, z=predval.binom), color="black", lty=2, breaks = c(0.005, 0.995)) +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(0.09, 0.72), breaks = c(seq(0.1, 0.7, by = 0.1))) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.07, 0.23), breaks = c(seq(0.08, 0.2, by = 0.04))) +
  ggtitle("B) Hugh Smith Lake")+ xlab ("Circulus 3–Circulus 8 (mm)") + ylab("Circulus 3–Scale Edge  (mm)") +
  annotate("text", x = 0.19, y=0.13, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.1, y=0.68, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> HSmodelpred

cowplot::plot_grid(BRmodelpred, HSmodelpred,  align = "h", nrow = 1, ncol=2) 
#ggsave("figures/model_pred_binom.png", dpi = 500, height = 4 , width = 6, units = "in")






cor(coho_scales_berners$dayofyear, coho_scales_berners$Q2plus-coho_scales_berners$Q2 )

