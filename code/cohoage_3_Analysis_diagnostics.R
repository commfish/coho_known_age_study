#Binomial regression diagnostics


# Email addresses: Sara Miller (sara.miller@alaska.gov) and Justin Priest (justinpriest.ak@gmail.com)
library(here)
library(MASS)
library(FNGr) #devtools::install_github("ben-williams/FNGr") 
library(cowplot)
library(extrafont)
library(rcompanion)
library(mgcv)
library(car)
library(modEvA) #pseudo-R squared
library(broom)
library(boot)


# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))


#Generalized Linear Models Diagnostics https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885900/
binomfit <- glm((Age-1)~Location+Q9abs+Q2, data=coho_scales_fulldata, family="binomial")
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
  do(A2 = glm((Age-1) ~ Location+Q9abs+Q2, data = coho_scales_fulldata, family = binomial(link=logit))) -> lm_out
head(augment(lm_out,coho_scales_fulldata, type.residuals="pearson")) #another way to write same model

lm_out %>% 
  tidy(A2) %>% 
  mutate(model = "binom") -> A2


lm_out %>% 
  glance(A2) %>% 
  mutate(model = "binom") -> A2


lm_out %>% # Pearson residuals against covariate
  augment(A2) %>% 
  mutate(resid = (.resid)) %>% 
  ggplot(aes(x = Q9abs, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3), limits = c(0, 0.3))+
  scale_y_continuous(breaks = seq(-3, 3, 0.5), limits = c(-3,3)) +
  geom_smooth(aes(colour = Q9abs, fill = Q9abs)) +
  labs(y = "Pearson residuals", x =  "Q9abs") -> plot1

lm_out %>% # Pearson residuals against covariate
  augment(A2) %>% 
  mutate(resid = (.resid)) %>% 
  ggplot(aes(x = Q2, y = resid)) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7), limits = c(0.1, 0.7))+
  scale_y_continuous(breaks = seq(-3.5, 3.5, 0.5), limits = c(-3.5,3.5)) +
  geom_smooth(aes(colour = Q9abs, fill = Q9abs)) +
  labs(y = "Pearson residuals", x =  "Q2") -> plot2


lm_out %>% #Pearson residuals against fitted
  augment(A2) %>% 
  mutate(resid = (.resid),
         fit = (.fitted)) %>% 
  ggplot(aes(x = fit, y = resid)) +
  geom_point(color ="grey50") + 
  scale_x_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  scale_y_continuous(breaks = seq(-2, 2, 0.5), limits = c(-2,2)) +
  geom_smooth(aes(colour = fit, fill = fit)) +
  geom_hline(yintercept = 0, lty=2) + 
  labs(y = "Pearson residuals", x =  "Fitted values") -> plot3

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
  scale_y_continuous(breaks = seq(-3, 3, 0.5), limits = c(-3,3)) +
  labs(y = "Pearson residuals", x =  "Index") -> plot6
cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6,  align = "vh", nrow = 3, ncol=2)
#ggsave("figures/glm_diagnostics.png", dpi = 500, height = 6, width = 8, units = "in")

#Residual plots as above
binom <- glm.diag(binomfit)
glm.diag.plots(binomfit, binom)

#Boxplot Figures 
coho_scales_fulldata %>% filter(Location == "BR") -> BR
coho_scales_fulldata %>% filter(Location == "HS") -> HS
BR %>% 
  ggplot(aes(x = as.factor(Age), y = Q2, color = as.factor(Age))) + 
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  scale_fill_manual(values=c("#999999", "#E69F00" )) +
  ggtitle("Berners") +labs(x = "Age", y = "Q2") -> plot1

BR %>% 
  ggplot(aes(x = as.factor(Age), y = Q9abs, color = as.factor(Age))) + 
  geom_boxplot() + theme(legend.position= "none") + 
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) +
  ggtitle("Berners") + labs(x = "Age", y = "Q9abs") -> plot2

HS %>% 
  ggplot(aes(x = as.factor(Age), y = Q2, color = as.factor(Age))) + 
  geom_boxplot() + theme(legend.position= "none") +
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  scale_fill_manual(values=c("#999999", "#E69F00" )) +
  ggtitle("Hugh Smith") + labs(x = "Age", y = "Q2") -> plot3

HS %>% 
  ggplot(aes(x = as.factor(Age), y = Q9abs, color = as.factor(Age))) + 
  geom_boxplot() + theme(legend.position= "none") + 
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_fill_manual(values=c("#999999", "#E69F00" )) +
  ggtitle("Hugh Smith") + labs(x = "Age", y = "Q9abs")  -> plot4

cowplot::plot_grid(plot1, plot2, plot3, plot4, align = "vh", nrow = 2, ncol=2)
#ggsave("figures/boxplot.png", dpi = 500, height = 6, width = 8, units = "in")



#cool plot if you can get it to work with our model (needs work!)

outer.data <- data.frame(Q2 = seq(0, 1, 0.1))# Create a temporary data frame of hypothetical values
outer.data2 <- data.frame(Q9abs = seq(0, 1, 0.1))#
data<-cbind(outer.data, outer.data2)
Location<-data.frame(rep("BR", 11))
data<-cbind(data, Location)
data %>% 
  rename(
    Location = rep..BR...11.)-> data1
outer.data <- data.frame(Q2 = seq(0, 1, 0.1))# Create a temporary data frame of hypothetical values
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



# JTP: I think that you can replace all of the above with:
predicted.data <- expand.grid(Q9abs=seq(min(coho_scales_bothriv$Q9abs), max(coho_scales_bothriv$Q9abs), length.out = 20), 
                              Q2=seq(min(coho_scales_bothriv$Q2), max(coho_scales_bothriv$Q2), length.out = 20),
                              Location=c("BR", "HS"))
predicted.data <- cbind(predicted.data, predict(binomfit, newdata = predicted.data, # Predict the fitted values given the model and hypothetical data
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
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8,1.0), limits = c(0, 1.0))

