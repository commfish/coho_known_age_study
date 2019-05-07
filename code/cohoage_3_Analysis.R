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



# Run Linear and Quadratic Discriminant Analysis (LDA & QDA) 
# Following approach from: https://www.statmethods.net/advstats/discriminant.html
fit_lda <- lda(Age ~ Location + Length + Q9, data=coho_scales_fulldata, #Q9 chosen after extensive testing
           na.action="na.omit", CV=TRUE)
ct <- table(coho_scales_fulldata$Age[!is.na(coho_scales_fulldata$Length)], fit_lda$class)
diag(prop.table(table(coho_scales_fulldata$Age[!is.na(coho_scales_fulldata$Length)], fit_lda$class), 1))
# total percent correct
sum(diag(prop.table(table(coho_scales_fulldata$Age[!is.na(coho_scales_fulldata$Length)], fit_lda$class)))) 
#Performs pretty well! 94% accuracy
klaR::partimat(as.factor(Age)~Length+Q9,data=coho_scales_fulldata,method="lda") 


# Perhaps model could be improved by accounting for quadratic curvature
fit_qda <- qda(as.factor(Age) ~ Location + Length + Q9, data=coho_scales_fulldata,
           prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)
# fits better as quadratic. And relaxed assumptions about homoscedascity


klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_fulldata, method="qda", prec=250, main="All Sites") 
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_aukelake, method="qda", prec=250, main="Auke Lake") 
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_berners, method="qda", prec=250, main="Berners River") 
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_hughsmith, method="qda", prec=250, main="Hugh Smith") 

coho_scales_fulldata <- coho_scales_fulldata %>% 
  mutate(pred_age = predict(fit_qda, coho_scales_fulldata)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))
coho_scales_fulldata %>% group_by(accuracy) %>% tally()



#########
# Same analysis but uses test and train subsets

# This function takes a dataframe, splits it into two portions (train and test)
# and then saves those dataframes. It will create two new dataframes at the end

err_rate <- data.frame(err_rate =as.numeric())
for(i in 1:1000){ #lower this to run faster
  traintest_subset <- function(df, trainproportion){ #proportion should be ~60%?
    .randrow <- sample(nrow(df), round(nrow(df) * trainproportion))
    .df_train <- df[.randrow,]
    .df_test <- df[-.randrow,]
    assign(paste0(deparse(substitute(df)), "_train"), get(".df_train"), envir = .GlobalEnv)
    assign(paste0(deparse(substitute(df)), "_test"), get(".df_test"), envir = .GlobalEnv)
  }
  
  traintest_subset(coho_scales_fulldata, 0.6) 
  
  fit_qda_train <- qda(as.factor(Age) ~ Location + Length + Q9, data=coho_scales_fulldata_train,
                       prior=c(0.66,0.34))
  
  err_rate[i,1] <- (coho_scales_fulldata_test %>% 
    mutate(pred_age = predict(fit_qda_train, coho_scales_fulldata_test)$class,
           accuracy = ifelse(pred_age == Age, "Correct", "Incorrect")) %>%
    group_by(accuracy) %>% tally())[3,2]
  
}
hist(err_rate$err_rate, breaks = 25)
err_rate$total <- err_rate$err_rate / nrow(coho_scales_fulldata_test)
quantile(err_rate$total, c(0.025, 0.975)) #95% CI for error rate is 5.5-7.9%

rm(list = c(ls(pattern = "._train"), ls(pattern = "._test")), "err_rate") 

# Summary, similar results as without using train/test datasets


#######
#manually recreate partimat() plot in ggplot

testpred <- expand.grid(Q9=seq(min(coho_scales_fulldata$Q9), max(coho_scales_fulldata$Q9), length.out = 100), 
            Length=seq(min(coho_scales_fulldata$Length, na.rm = TRUE), max(coho_scales_fulldata$Length, na.rm = TRUE), by=1),
            Location = c("AL", "HS", "BR"))

testpred$predval <- predict(fit_qda, testpred)$class

coho_scales_aukelake <- coho_scales_aukelake %>% 
  mutate(pred_age = predict(fit_qda, coho_scales_aukelake)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))
coho_scales_berners <- coho_scales_berners %>% 
  mutate(pred_age = predict(fit_qda, coho_scales_berners)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))
coho_scales_hughsmith <- coho_scales_hughsmith %>% 
  mutate(pred_age = predict(fit_qda, coho_scales_hughsmith)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))

write.csv(coho_scales_berners, "berners.csv")

ggplot() + 
  geom_tile(data = testpred %>% filter(Location == "BR"), aes(x=Q9, y=Length, fill=predval)) +
  geom_text(data = coho_scales_berners, aes(x=Q9, y=Length, label=Age, color=accuracy)) + 
  scale_color_manual(values = c("black", "red"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("#89b6ff", "#ff9260")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  ggtitle("Berners River")

ggplot() + 
  geom_tile(data = testpred %>% filter(Location == "HS"), aes(x=Q9, y=Length, fill=predval)) +
  geom_text(data = coho_scales_hughsmith, aes(x=Q9, y=Length, label=Age, color=accuracy)) + 
  scale_color_manual(values = c("black", "red"), guide=FALSE) +
  scale_fill_manual(name="Predicted\nAge", values = c("#89b6ff", "#ff9260")) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  ggtitle("Hugh Smith")



coho_scales_hughsmith %>% group_by(Age, accuracy) %>% tally()
