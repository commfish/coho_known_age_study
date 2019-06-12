#IN PROGRESS. This will be removed once code is completely updated. 
# Script 3 - Analysis
#' ---
#' title: "Coho Known Age LDA Analysis"
#' author: "Sara Miller & Justin Priest"
#' date: "`r format(Sys.Date())`"
#' ---

# Email addresses: Sara Miller (sara.miller@alaska.gov) and Justin Priest (justinpriest.ak@gmail.com)
library(here)
library(MASS)
library(vegan)
library(klaR)
library(FNGr) #devtools::install_github("ben-williams/FNGr") 
library(processx)
library(fs)
library(cowplot)
library(extrafont)
library(rcompanion)
library(mgcv)
library(car)
library(modEvA) #pseudo-R squared
library(broom)
library(boot)
windowsFonts(Times=windowsFont("TT Times New Roman")) 
theme_set(theme_sleek())

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))
source(here::here("code/functions.R"))

#Generalized Linear Models Diagnostics https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885900/
binomfit <- glm((Age-1)~Location+Q9abs+Q2plus, data=coho_scales_fulldata, family="binomial")
summary(binomfit)
Anova(binomfit)
RsqGLM(binomfit)#peudo R2 
outlierTest(binomfit) #Bonferroni p-values
residualPlots(binomfit) #lack-of fit curvature test
influenceIndexPlot(binomfit) #these plots are meade below in ggplot
marginalModelPlots(binomfit) #marginal model plots
par(mfrow=c(1,1))
png(filename="figures/glm_diagnostics1.png", width = 6, height = 6, units = 'in', res = 400)
windowsFonts(A = windowsFont("Times New Roman"))
influencePlot(binomfit, xlab="Hat-values", ylab="Studentized residuals", family="A")
dev.off()

coho_scales_fulldata %>% 
  do(A2 = glm((Age-1) ~ Location+Q9abs+Q2plus, data = coho_scales_fulldata, family = binomial(link=logit))) -> lm_out
head(augment(lm_out,coho_scales_fulldata, type.residuals="pearson")) #another way to write same model

lm_out %>% 
  tidy(A2) %>% 
  mutate(model = "binom") -> A2
write_csv(A2, "output/lm.csv")

lm_out %>% 
  glance(A2) %>% 
  mutate(model = "binom") -> A2
write_csv(A2, "output/lm_R2.csv")

lm_out %>% # Pearson residuals against covariate
  augment(A2) %>% 
  mutate(resid = (.resid)) %>% 
  ggplot(aes(x = Q9abs, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3), limits = c(0, 0.3))+
  scale_y_continuous(breaks = c(-3, -2.5,-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3), limits = c(-3,3)) +
  geom_smooth(aes(colour = Q9abs, fill = Q9abs)) +
  labs(y = "Pearson residuals", x =  "Q9abs")-> plot1

lm_out %>% # Pearson residuals against covariate
  augment(A2) %>% 
  mutate(resid = (.resid)) %>% 
  ggplot(aes(x = Q2plus, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7), limits = c(0.1, 0.7))+
  scale_y_continuous(breaks = c(-3.5,-3,-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5,3, 3.5), limits = c(-3.5,3.5)) +
  geom_smooth(aes(colour = Q9abs, fill = Q9abs)) +
  labs(y = "Pearson residuals", x =  "Q2plus")-> plot2


lm_out %>% #Pearson residuals against fitted
  augment(A2) %>% 
  mutate(resid = (.resid),
         fit = (.fitted)) %>% 
  ggplot(aes(x = fit, y = resid)) +
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4,0.5), limits = c(0, 0.5))+
  scale_y_continuous(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit)) +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values")-> plot3

lm_out %>% #Cook's distance plot
  augment(A2) %>% 
  mutate(cooksd = (.cooksd),
         count = 1:3117,
         name= ifelse(cooksd >0.011, count, "")) %>% #this displays cooks distances>0.011 (this can be changed)
  ggplot(aes(x = count, y = cooksd, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1)) + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03), limits = c(0,0.03)) +
  labs(y = "Cook's distance", x =  "Index") -> plot4

lm_out %>% #leverage plot
  augment(A2) %>% 
  mutate(hat= (.hat),
         count = 1:3117,
         name= ifelse(hat >0.01, count, "")) %>% 
  ggplot(aes(x = count, y = hat, label=name)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  scale_y_continuous(breaks = c(0, 0.005, 0.01, 0.015, 0.02), limits = c(0,0.02)) +
  labs(y = "hat-values", x =  "Index") -> plot5

lm_out %>% #Pearson by index
  augment(A2) %>% 
  mutate(resid = (.resid),
         count = 1:3117) %>% 
  ggplot(aes(x = count, y = resid)) +
  geom_bar(stat = "identity", colour = "grey50", 
           fill = "lightgrey",alpha=.7,
           width = 0.8, position = position_dodge(width = 0.2)) + 
  scale_y_continuous(breaks = c(-3, -2.5,-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3), limits = c(-3,3)) +
  labs(y = "Pearson residuals", x =  "Index") -> plot6
cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,  align = "vh", nrow = 3, ncol=2)
ggsave("figures/glm_diagnostics.png", dpi = 500, height = 6, width = 8, units = "in")

#Residual plots as above
binom <- glm.diag(binomfit)
glm.diag.plots(binomfit, binom)

#Boxplot Figures 
coho_scales_fulldata %>% filter(Location == "BR") -> BR
coho_scales_fulldata %>% filter(Location == "HS") -> HS
BR %>% 
  ggplot(aes(x = as.factor(Age), y = Q2plus, color = as.factor(Age))) + labs(x = "Age",
                                                              y = "Q2plus")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +ggtitle("Berners") +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot1

BR %>% 
  ggplot(aes(x = as.factor(Age), y = Q9abs, color = as.factor(Age))) + labs(x = "Age",
                                                                             y = "Q9abs")  +
  geom_boxplot() + theme(legend.position= "none") +ggtitle("Berners") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot2

HS %>% 
  ggplot(aes(x = as.factor(Age), y = Q2plus, color = as.factor(Age))) + labs(x = "Age",
                                                                             y = "Q2plus")  +
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) +ggtitle("Hugh Smith") +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot3

HS %>% 
  ggplot(aes(x = as.factor(Age), y = Q9abs, color = as.factor(Age))) + labs(x = "Age",
                                                                            y = "Q9abs")  +
  geom_boxplot() + theme(legend.position= "none") +ggtitle("Hugh Smith") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" ))-> plot4

cowplot::plot_grid(plot1, plot2, plot3, plot4, align = "vh", nrow = 2, ncol=2)
ggsave("figures/boxplot.png", dpi = 500, height = 6, width = 8, units = "in")

#cool plot if you can get it to work with our model (needs work!)
outer.data <- data.frame(Q2plus = seq(0, 1, 0.1))# Create a temporary data frame of hypothetical values
outer.data2 <- data.frame(Q9abs = seq(0, 1, 0.1))#
data<-cbind(outer.data, outer.data2)
Location<-data.frame(rep("BR", 11))
data<-cbind(data, Location)
data %>% 
  rename(
    Location = rep..BR...11.)-> data1
outer.data <- data.frame(Q2plus = seq(0, 1, 0.1))# Create a temporary data frame of hypothetical values
outer.data2 <- data.frame(Q9abs = seq(0, 1, 0.1))#
data<-cbind(outer.data, outer.data2)
Location<-data.frame(rep("HS", 11))
data<-cbind(data, Location)
data %>% 
  rename(
    Location = rep..HS...11.)-> data2
data<-rbind(data1, data2)

predicted.data <- as.data.frame(predict(binomfit, newdata = data, # Predict the fitted values given the model and hypothetical data
                                        type="link", se=TRUE))
new.data <- cbind(data, predicted.data)# Combine the hypothetical data and predicted values
std <- qnorm(0.95 / 2 + 0.5)# Calculate confidence intervals
new.data$ymin <- binomfit$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- binomfit$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- binomfit$family$linkinv(new.data$fit)  # Rescale to 0-1

coho_scales_fulldata %>% 
  ggplot(aes(x=Q9abs, y=Age-1)) +
  geom_point() + 
  geom_ribbon(data=new.data, aes(y=fit, ymin=ymin, ymax=ymax), alpha=0.5) + 
  geom_line(data=new.data, aes(y=fit)) + 
  labs(x="Q9abs", y="Age") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4,0.5), limits = c(0, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8,1.0), limits = c(0, 1.0)) -> plot3







# Run QDA function  
# To perform QDA, one must ensure that the within-group covariance matrices of
# the explanatory variables are homogeneous, a condition that is frequently violated
# with ecological data ;If one the groups defined by the dependent variable has greater 
# dispersion than others, cases will tend to be overclassified in it.


#Transformed Data; https://rcompanion.org/handbook/I_12.html: Q2plus age1 Berners
Box = boxcox(coho_scales_berners$Q2plus[coho_scales_berners$Age == 1] ~ 1, # Transform Turbidity as a single vector
             lambda = seq(-6,6,0.1)) # Try values -6 to 6 by 0.1)
Cox = data.frame(Box$x, Box$y) # Create a data frame with the results
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
Cox2[1,]# Display the lambda with the greatest log likelihood
lambda1 = Cox2[1, "Box.x"] # Extract that lambda
x<-coho_scales_berners$Q2plus[coho_scales_berners$Age == 1]
T_box = (x^lambda1 - 1)/lambda1# Transform the original data
eda.norm(T_box)
eda.norm(coho_scales_berners$Q2plus[coho_scales_berners$Age == 1]) #compare to untransformed

#Q2plus age2 Berners
Box = boxcox(coho_scales_berners$Q2plus[coho_scales_berners$Age == 2] ~ 1, 
             lambda = seq(-6,6,0.1)) 
Cox = data.frame(Box$x, Box$y) 
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda2 = Cox2[1, "Box.x"] 
x<-coho_scales_berners$Q2plus[coho_scales_berners$Age == 2]
T_box = (x^lambda2 - 1)/lambda2
eda.norm(T_box)
eda.norm(coho_scales_berners$Q2plus[coho_scales_berners$Age == 2]) 

#Q9abs age1 Berners
Box = boxcox(coho_scales_berners$Q9abs[coho_scales_berners$Age == 1] ~ 1,
             lambda = seq(-6,6,0.1)) 
Cox = data.frame(Box$x, Box$y) 
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda3 = Cox2[1, "Box.x"] 
x<-coho_scales_berners$Q9abs[coho_scales_berners$Age == 1]
T_box = (x^lambda3 - 1)/lambda3
eda.norm(T_box)
eda.norm(coho_scales_berners$Q9abs[coho_scales_berners$Age == 1])

#Q9abs age2 Berners
Box = boxcox(coho_scales_berners$Q9abs[coho_scales_berners$Age == 2] ~ 1, 
             lambda = seq(-6,6,0.1))
Cox = data.frame(Box$x, Box$y)
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda4 = Cox2[1, "Box.x"] # Extract that lambda
x<-coho_scales_berners$Q9abs[coho_scales_berners$Age == 2]
T_box = (x^lambda4 - 1)/lambda4
eda.norm(T_box)
eda.norm(coho_scales_berners$Q9abs[coho_scales_berners$Age == 2]) 

#Transformed Data; Q2plus age1 Hughsmith
Box = boxcox(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 1] ~ 1,
             lambda = seq(-6,6,0.1)) 
Cox = data.frame(Box$x, Box$y)
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda5 = Cox2[1, "Box.x"]
x<-coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 1]
T_box = (x^lambda5 - 1)/lambda5
eda.norm(T_box)
eda.norm(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 1]) 

#Q2plus age2 hughsmith
Box = boxcox(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 2] ~ 1, 
             lambda = seq(-6,6,0.1)) 
Cox = data.frame(Box$x, Box$y) 
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda6 = Cox2[1, "Box.x"] 
x<-coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 2]
T_box = (x^lambda6 - 1)/lambda6
eda.norm(T_box)
eda.norm(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 2]) 

#Q9abs age1 hughsmith
Box = boxcox(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 1] ~ 1,
             lambda = seq(-6,6,0.1)) 
Cox = data.frame(Box$x, Box$y) 
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda7 = Cox2[1, "Box.x"] 
x<-coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 1]
T_box = (x^lambda7 - 1)/lambda7
eda.norm(T_box)
eda.norm(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 1])

#Q9abs age2 hughsmith
Box = boxcox(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 2] ~ 1, 
             lambda = seq(-6,6,0.1))
Cox = data.frame(Box$x, Box$y)
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] 
Cox2[1,]
lambda8 = Cox2[1, "Box.x"] 
x<-coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 2]
T_box = (x^lambda8 - 1)/lambda8
eda.norm(T_box)
eda.norm(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 2]) 

# Transformed Variables - Create New Datasets
coho_scales_berners<- coho_scales_bothriv %>% 
  filter(Location == "BR") %>%
  mutate(Q2plus_trans = ifelse(Age==1, ((Q2plus^lambda1-1)/lambda1),
                ifelse(Age==2, ((Q2plus^lambda2-1)/lambda2),NA)),
  Q9abs_trans = ifelse(Age==1, ((Q9abs^lambda3-1)/lambda3),
                      ifelse(Age==2, ((Q9abs^lambda4-1)/lambda4),NA)))

# Transformed Variables - Create New Datasets
coho_scales_hughsmith<- coho_scales_bothriv %>% 
  filter(Location == "HS") %>%
  mutate(Q2plus_trans = ifelse(Age==1, ((Q2plus^lambda5-1)/lambda5),
                               ifelse(Age==2, ((Q2plus^lambda6-1)/lambda6),NA)),
         Q9abs_trans = ifelse(Age==1, ((Q9abs^lambda7-1)/lambda7),
                              ifelse(Age==2, ((Q9abs^lambda8-1)/lambda8),NA)))

coho_scales_bothriv <-rbind(coho_scales_berners, coho_scales_hughsmith)


# Run QDA
fit_qda_final <- qda(as.factor(Age) ~ Location + Q2plus + Q9abs, data=coho_scales_bothriv,
                prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)

# Create blank predicted values (for plotting)
predictedvals <- expand.grid(Q9abs =seq(min(coho_scales_bothriv$Q9abs), 
                                          max(coho_scales_bothriv$Q9abs)+0.05, length.out = 100), #adding a little bit for plotting later
                         Q2plus=seq(min(coho_scales_bothriv$Q2plus), 
                                          max(coho_scales_bothriv$Q2plus)+0.1, length.out = 100), #adding a little bit for plotting later
                         Location = c("HS", "BR"))

predictedvals$predval <- predict(fit_qda_final, predictedvals)$class

coho_scales_berners <- coho_scales_berners %>% 
  mutate(pred_age = predict(fit_qda_final, coho_scales_berners)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))

coho_scales_hughsmith <- coho_scales_hughsmith %>% 
  mutate(pred_age = predict(fit_qda_final, coho_scales_hughsmith)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))
final<-rbind(coho_scales_berners,coho_scales_hughsmith)
write.csv(write.csv(final, "data/check.csv"))

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "BR"), aes(x=Q9abs, y=Q2plus, fill=predval)) +
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Correct"), aes(x=Q9abs, y=Q2plus, label=Age), color="grey60") + 
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Incorrect"), aes(x=Q9abs, y=Q2plus, label=Age), color="black") +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(0.09, 0.72), breaks = c(seq(0.1, 0.7, by = 0.1))) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.07, 0.23), breaks = c(seq(0.08, 0.2, by = 0.04))) +
  ggtitle("A) Berners River")+ xlab ("Q9abs")+ylab("Q2plus") +
  annotate("text", x = 0.19, y=0.13, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.1, y=0.68, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> plot1

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs, y=Q2plus, fill=predval)) +
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Correct"), aes(x=Q9abs, y=Q2plus, label=Age), color="grey60") + 
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Incorrect"), aes(x=Q9abs, y=Q2plus, label=Age), color="black") +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(0.09, 0.72), breaks = c(seq(0.1, 0.7, by = 0.1))) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.07, 0.23), breaks = c(seq(0.08, 0.2, by = 0.04))) +
  ggtitle("B) Hugh Smith")+ xlab ("Q9abs") + ylab("Q2plus") +
  annotate("text", x = 0.19, y=0.13, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = 0.1, y=0.68, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> plot2

cowplot::plot_grid(plot1, plot2,  align = "h", nrow = 1, ncol=2) 
ggsave("figures/model_pred.png", dpi = 500, height = 4 , width = 6, units = "in")

# Run QDA Transformed
fit_qda_final <- qda(as.factor(Age) ~ Location + Q2plus_trans + Q9abs_trans, data=coho_scales_bothriv,
                     prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)


coho_scales_hughsmith %>% group_by(accuracy) %>% tally()
coho_scales_berners %>% group_by(accuracy) %>% tally()
paste0("BR SQRT ARCSINE: ", round((1-117/(1850 + 117))*100, 2), "%")
paste0("HS SQRT ARCSINE: ", round((1-62/(790+62))*100, 2), "%")

paste0("BR NO TRANS: ", round((1-118/(1849 + 118))*100, 2), "%")
paste0("HS NO TRANS: ", round((1-60/(792+60))*100, 2), "%")

paste0("BR SQRT: ", round((1-119/(1848 + 119))*100, 2), "%")
paste0("HS SQRT: ", round((1-102/(750+102))*100, 2), "%")



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
final<-rbind(coho_scales_berners,coho_scales_hughsmith)
write.csv(write.csv(final, "data/check_trans.csv"))

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "BR"), aes(x=Q9abs_trans, y=Q2plus_trans, fill=predval)) +
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Correct"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="grey60") + 
  geom_text(data = coho_scales_berners %>% filter(accuracy == "Incorrect"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="black") +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_y_continuous(expand = c(0,0), limits = c(-3, -0.4), breaks = c(seq(-3,-0.4, by = 0.2))) + 
  scale_x_continuous(expand = c(0,0), limits = c(-1, -0.8), breaks = c(seq(-1,-0.8, by = 0.1))) + 
  ggtitle("A) Berners River")+ xlab ("Q9abs (transformed)")+ylab("Q2plus (transformed)") +
  annotate("text", x = -0.85, y=-2.8, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = -0.84, y=-0.5, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> plot1

ggplot() + 
  geom_tile(data = predictedvals %>% filter(Location == "HS"), aes(x=Q9abs_trans, y=Q2plus_trans, fill=predval)) +
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Correct"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="grey60") + 
  geom_text(data = coho_scales_hughsmith %>% filter(accuracy == "Incorrect"), aes(x=Q9abs_trans, y=Q2plus_trans, label=Age), color="black") +
  scale_fill_manual(name="Predicted\nAge", values = c("grey90", "grey95")) +
  scale_x_continuous(expand = c(0,0), limits = c(-3, -1), breaks = c(seq(-3,-1, by = 0.5))) + 
  scale_y_continuous(expand = c(0,0), limits = c(-1, -0.3), breaks = c(seq(-1,-0.3, by = 0.1))) + 
  ggtitle("B) Hugh Smith")+ xlab ("Q9abs (transformed)") + ylab("Q2plus (transformed)") +
  annotate("text", x = -2.6, y=-0.4, label="predicted age 1", family="Times New Roman", colour="black", size=3) +
  annotate("text", x = -1.3, y=-0.35, label="predicted age 2", family="Times New Roman", colour="black", size=3) +
  theme(legend.position="none", text=element_text(family="Times New Roman", size=12)) -> plot2

cowplot::plot_grid(plot1, plot2,  align = "h", nrow = 1, ncol=2) 
ggsave("figures/model_pred_trans.png", dpi = 500, height = 4 , width = 6, units = "in")