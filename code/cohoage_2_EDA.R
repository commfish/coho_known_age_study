# load libraries
library(corrplot)
library(GGally)
library(here)
library(mgcv)
library(mvnormtest)
library(HH)
library(klaR)
library(MASS)
library(mgcv)
library(mvoutlier)
library(prabclus)

# Run previous script to import data
source(here::here("code/functions.R"))
source(here::here("code/cohoage_1_DataImport_JTP.R"))




#----------------------------# 
#### BEGIN EDA ####


# Correlation between variables
for(q in c("HS", "BR", "AL")) { # Print a few correlation plots between variables
  print(corrplot.mixed(cor(coho_scales_fulldata %>%
                     filter(Location == q) %>%
                     dplyr::select("Q32_sum":"Q44"), use="complete.obs"), title=q, mar=c(0,0,1,0) ))
  
  # Now explore how these correlations and "dropoffs" change btwn locations
  .temp4 <- cor(coho_scales_fulldata %>%
                  filter(Location == q) %>% 
                  dplyr::select("Q4":"Q31"), use="complete.obs")
  .temp4[upper.tri(.temp4)] <- NA

  shift <- function(x, n){c(x[-(seq(n))], rep(NA, n))}
  
  .temp4 <- as.data.frame(.temp4)
  for(i in 2:29){.temp4[,i] <- shift(.temp4[,i], i-1)}
  
  .temp4["rownum"] <- seq(1:29)
  .temp4 <- .temp4 %>% gather("Qnum", "corr", Q4:Q31)
  
  print(ggplot(.temp4, aes(x=rownum, y = corr, color=Qnum)) +
    geom_line() + 
    labs(title = q))
} 
#Summary, clear correlation between close circuli variables


#caution, this takes a while to eval! Uncomment to run
#ggpairs(coho_scales_berners %>% dplyr::select("Zone1":"Count_Zone2", "Year", "Age":"Q10"))


# EDA SUMMARY ON VARIABLE SELECTION
# Based on discussion with LW and viewing data, we probably want to explore variables
# that are a few circuli out (circuli 5- 15?), because not all scales will have more 
# than a dozen circuli (esp if Age 1). This approx corresponds to Q6 - Q16. 
# Also, Auke Lake data might not have had true known ages on them.
# Exclude Auke Lake ("AL") for most analyses. But can be useful to see how it performs. 


# --------------------------------------------#
### Exploratory Modeling ###
# Explore a few plots and model types to get a general feel for the data
ggplot(coho_scales_berners, aes(y = Q9, group = Age)) +
  geom_boxplot() #clear difference in Q9 between ages

summary(lm(data=coho_scales_aukelake, Age~Length+Q5+Q12))
summary(gam(data=coho_scales_aukelake, Age ~ Length + Q5))
summary(gam(data=coho_scales_aukelake, Age ~ Length + s(Q5))) # deviance explained rises quite a bit
summary(gam(data=coho_scales_aukelake, (Age-1) ~ Length + Q5, family = "binomial"))
summary(glm(data=coho_scales_aukelake, (Age-1) ~ log(Length) + log(Q5), family = "binomial")) # Top linear model (Len&Q5 only)
summary(gam(data=coho_scales_aukelake, (Age-1) ~ Length + s(Q5), family = "binomial")) # Top non-linear model (Len&Q5 only)

boxplot(Q2~Age, data = coho_scales_berners)
boxplot(Length~Age, data = coho_scales_berners)
summary(lm(Age ~ Length, data = coho_scales_berners))
summary(lm(Age ~ Q2, data = coho_scales_berners))
summary(lm(Age ~ Q9, data = coho_scales_berners))


# --------------------------------------------#
### Correlation ###
test1 <- coho_scales_aukelake %>% filter(!is.na(Length)) #remove NAs from Length to test corr with other vars
cor(test1$Q5, test1$Length) # High corr
cor(test1$Q9, test1$Length) # Top predictor variables are very correlated
cor(test1$Q5, test1$Age) 
cor(test1$Length, test1$Age) # But are also good predictors of Age
rm(test1)

cor(coho_scales_berners$Q9abs, coho_scales_berners$Q2plus)
cor(coho_scales_hughsmith$Q9abs, coho_scales_hughsmith$Q2plus)


# --------------------------------------------#
### Rough Histograms ###
hist(coho_scales_aukelake$Length, breaks=15)
hist((coho_scales_aukelake %>% filter(Age==1))$Length, breaks=15)
hist((coho_scales_aukelake %>% filter(Age==2))$Length, breaks=15)
hist(coho_scales_aukelake$Q9, breaks=15)
hist(log(coho_scales_aukelake$Length), breaks=15)
hist(log((coho_scales_aukelake %>% filter(Age==1))$Length), breaks=15)
hist(log((coho_scales_aukelake %>% filter(Age==2))$Length), breaks=15)
hist(log(coho_scales_aukelake$Q9), breaks=15)


# --------------------------------------------#
#### UNIVARIATE NORMALITY ASSESSMENTS ####
coho_scales_berners2 <- coho_scales_fulldata %>% filter(Location == "BR")
coho_scales_hughsmith2 <- coho_scales_fulldata %>% filter(Location == "HS")
qqnorm((coho_scales_berners2 %>% filter(Age==1))$Q9abs); qqline((coho_scales_berners2 %>% filter(Age==1))$Q9abs)
qqnorm((coho_scales_berners2 %>% filter(Age==2))$Q9abs); qqline((coho_scales_berners2 %>% filter(Age==2))$Q9abs)
qqnorm((coho_scales_berners2 %>% filter(Age==1))$Q2plus); qqline((coho_scales_berners2 %>% filter(Age==1))$Q2plus)
qqnorm((coho_scales_berners2 %>% filter(Age==2))$Q2plus); qqline((coho_scales_berners2 %>% filter(Age==2))$Q2plus)

qqnorm((coho_scales_hughsmith2 %>% filter(Age==1))$Q9abs); qqline((coho_scales_hughsmith2 %>% filter(Age==1))$Q9abs)
qqnorm((coho_scales_hughsmith2 %>% filter(Age==2))$Q9abs); qqline((coho_scales_hughsmith2 %>% filter(Age==2))$Q9abs)
qqnorm((coho_scales_hughsmith2 %>% filter(Age==1))$Q2plus); qqline((coho_scales_hughsmith2 %>% filter(Age==1))$Q2plus)
qqnorm((coho_scales_hughsmith2 %>% filter(Age==2))$Q2plus); qqline((coho_scales_hughsmith2 %>% filter(Age==2))$Q2plus)


# --------------------------------------------#
#### NORMALITY TESTING ####
mvnormtest::mshapiro.test(t(coho_scales_hughsmith %>% dplyr::select(Age, Length, Q2, Q5, Q9) %>% na.omit %>% as.matrix()))

#try centering data
#coho_scales_berners_sub <-coho_scales_berners %>%
#  dplyr::select(Q2plus, Q9abs)
#preproc.param <- coho_scales_berners_sub %>%
#  preProcess(method = c("center", "scale"))
#coho_scales_berners_trans <- preproc.param %>% predict(coho_scales_berners_sub)
#coho_scales_berners2<-coho_scales_berners %>%
#  dplyr::select(Location,Age, Length)
#coho_scales_berners <-cbind(coho_scales_berners2, coho_scales_berners_trans)

eda.norm(coho_scales_aukelake$Q9) # Not normal
eda.norm(log(coho_scales_aukelake$Q9)) # Normal
eda.norm(log(coho_scales_aukelake$Length))
eda.norm((coho_scales_aukelake$Length)^0.25)

eda.norm((coho_scales_aukelake$Length[coho_scales_aukelake$Age == 1])) #Not normal, but slightly better log trans
eda.norm((coho_scales_aukelake$Length[coho_scales_aukelake$Age == 2])) #Not normal, but slightly better log trans
eda.norm((coho_scales_berners$Length[coho_scales_berners$Age == 1])) # VERY Not normal, but slightly better log trans
eda.norm((coho_scales_berners$Length[coho_scales_berners$Age == 2])) #Not normal, but worse with log trans
eda.norm((coho_scales_hughsmith$Length[coho_scales_hughsmith$Age == 1])) # Normal!
eda.norm((coho_scales_hughsmith$Length[coho_scales_hughsmith$Age == 2])) # Normal

eda.norm(coho_scales_berners$Q2plus[coho_scales_berners$Age == 1]) #very not normal
eda.norm(coho_scales_berners$Q2plus[coho_scales_berners$Age == 2]) # not normal
eda.norm(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 1]) #normal
eda.norm(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 2]) #almost normal

eda.norm(coho_scales_berners$Q9abs[coho_scales_berners$Age == 1]) # normal
eda.norm(coho_scales_berners$Q9abs[coho_scales_berners$Age == 2]) # normal
eda.norm(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 1]) # not normal
eda.norm(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 2]) # normal

eda.norm(asin(sqrt(coho_scales_berners$Q2plus[coho_scales_berners$Age == 1]))) #not normal but improvement
eda.norm(asin(sqrt(coho_scales_berners$Q2plus[coho_scales_berners$Age == 2]))) # not normal, about same
eda.norm(asin(sqrt(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 1]))) #normal
eda.norm(asin(sqrt(coho_scales_hughsmith$Q2plus[coho_scales_hughsmith$Age == 2]))) #almost normal, same

eda.norm(asin(sqrt(coho_scales_berners$Q9abs[coho_scales_berners$Age == 1]))) # almost normal, worse
eda.norm(asin(sqrt(coho_scales_berners$Q9abs[coho_scales_berners$Age == 2]))) # normal, worse
eda.norm(asin(sqrt(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 1]))) # not normal, great improvement tho
eda.norm(asin(sqrt(coho_scales_hughsmith$Q9abs[coho_scales_hughsmith$Age == 2]))) # normal

#recommend transforming with arcsine square root. But first explore everything without transformations, save for last.


# Homogeneity of variances Tests
hist((coho_scales_berners$Length[coho_scales_berners$Age == 1]), breaks = 25)

bartlett.test(log(Q5)~Age, data=coho_scales_aukelake) 
fligner.test(Length~Age, data=coho_scales_aukelake) 

boxplot(Length~Age, data = coho_scales_berners) # Clear difference in length by ages
boxplot(Length~Age, data = coho_scales_hughsmith)
cov(coho_scales_aukelake$Q5, coho_scales_aukelake$Age)
cov(coho_scales_aukelake$Length[!is.na(coho_scales_aukelake$Length)], 
    coho_scales_aukelake$Age[!is.na(coho_scales_aukelake$Length)])

HH::hov(Q9~as.factor(Age), data=coho_scales_aukelake)
hov(Length~as.factor(Age), data=coho_scales_aukelake %>% filter(!is.na(Length)))

hovPlot(Q5~as.factor(Age), data=coho_scales_aukelake)
hovPlot(Length~as.factor(Age), data=coho_scales_aukelake %>% filter(!is.na(Length)))
HH::hov(Q9abs~as.factor(Age), data=coho_scales_aukelake) #p >0.05 homogeneity of variances (Brown Forsyth test)
HH::hov(Q9abs~as.factor(Age), data=coho_scales_berners)
HH::hov(Q9abs~as.factor(Age), data=coho_scales_hughsmith)

HH::hov(Q2plus~as.factor(Age), data=coho_scales_aukelake)
HH::hov(Q2plus~as.factor(Age), data=coho_scales_berners)
HH::hov(Q2plus~as.factor(Age), data=coho_scales_hughsmith)

HH::hov(Length~as.factor(Age), data=coho_scales_aukelake %>% filter(!is.na(Length)))

hovPlot(Q2plus~as.factor(Age), data=coho_scales_aukelake)
hovPlot(Q2plus~as.factor(Age), data=coho_scales_berners)
hovPlot(Q2plus~as.factor(Age), data=coho_scales_hughsmith)

hovPlot(Q9abs~as.factor(Age), data=coho_scales_aukelake)
hovPlot(Q9abs~as.factor(Age), data=coho_scales_berners)
hovPlot(Q9abs~as.factor(Age), data=coho_scales_hughsmith)


# --------------------------------------------#
### Visualize Departures from Normality ###

hist.normal(coho_scales_aukelake, "Length", 1)
hist.normal(coho_scales_aukelake, "Length", 2)
hist.normal(coho_scales_aukelake, "Length")
hist.normal(coho_scales_aukelake, "Q9", 1)
hist.normal(coho_scales_aukelake, "Q9", 2)
hist.normal(coho_scales_aukelake, "Q9")

hist.normal(coho_scales_berners, "Length", 1)
hist.normal(coho_scales_berners, "Length", 2)
hist.normal(coho_scales_berners, "Length")
hist.normal(coho_scales_berners, "Q2plus", 1)
hist.normal(coho_scales_berners, "Q2plus", 2)
hist.normal(coho_scales_berners, "Q2plus")
hist.normal(coho_scales_berners, "Q9", 1)
hist.normal(coho_scales_berners, "Q9", 2)
hist.normal(coho_scales_berners, "Q9")
hist.normal(coho_scales_berners, "Q9abs", 1)
hist.normal(coho_scales_berners, "Q9abs", 2)

hist.normal(coho_scales_hughsmith, "Length", 1)
hist.normal(coho_scales_hughsmith, "Length", 2)
hist.normal(coho_scales_hughsmith, "Length")
hist.normal(coho_scales_hughsmith, "Q2plus", 1)
hist.normal(coho_scales_hughsmith, "Q2plus", 2)
hist.normal(coho_scales_hughsmith, "Q9", 1)
hist.normal(coho_scales_hughsmith, "Q9", 2)
hist.normal(coho_scales_hughsmith, "Q9")
hist.normal(coho_scales_hughsmith, "Q9abs", 1)
hist.normal(coho_scales_hughsmith, "Q9abs", 2)



# All in all, most times there is not a large violation of normality, at least visually


# --------------------------------------------#
### Check Multivariate Normality ###
ggplot(coho_scales_fulldata %>% filter(Age == 1), aes(Length, Q9)) + 
  stat_bin2d(bins = 20) +
  coord_cartesian(xlim = c(60, 140), ylim = c(0, 1)) +
  facet_wrap(~Location)

ggplot(coho_scales_fulldata %>% filter(Age == 2), aes(Length, Q9)) + 
  stat_bin2d(bins = 20) +
  coord_cartesian(xlim = c(60, 140), ylim = c(0, 1)) +
  facet_wrap(~Location)

# Large differences by location. Clear correlation. Not all combinations are bad, Age2s are pretty good



#### EDA OVERVIEW SUMMARY ####-----------------
# Very large differences between locations (model will need to account for location)
# Many of the variables are approximately normal but fail the 'strict' conditions of tests
# Variance is not even between Age 1 & 2
# It appears that the higher numbered Q variables are not as good of a fit. 
# Explore further with Q variables <Q12 or so.
# There is clear violation of homogeneity of variance by age.
# Some moderate MVN violations, especially for age-1 fish
# Top models appear to be Age ~ Length + Q9. 
# Downside of this is that we would like to use scale-only data (i.e., not use Smolt Length)
# According to LW, scale distance (including plus growth), is good proxy for length
# I have added this as variable "Q2plus"; to reduce multi-collinearity, I used the 
# absolute distance to Circuli 8 (Q9abs), instead of relative distance to Circuli 8 (Q9)
# Thus, I recommend using top model as Age ~ Q9abs + Q2plus


# --------------------------------------------#



# --------------------------------------------#
#### EDA MODELING #### -----------------------


#### LDA & QDA MODELING ####
# Linear and Quadratic Discriminant Analysis (LDA & QDA) 
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

# Visualize
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_fulldata, method="qda", prec=250, main="All Sites") 
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_aukelake, method="qda", prec=250, main="Auke Lake") 
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_berners, method="qda", prec=250, main="Berners River") 
klaR::partimat(as.factor(Age) ~ Length + Q9, data=coho_scales_hughsmith, method="qda", prec=250, main="Hugh Smith") 

coho_scales_fulldata <- coho_scales_fulldata %>% 
  mutate(pred_age = predict(fit_qda, coho_scales_fulldata)$class,
         accuracy = ifelse(pred_age == Age, "Correct", "Incorrect"))
coho_scales_fulldata %>% group_by(accuracy) %>% tally() #93.5% accuracy, even inc Auke Lake

#### Summary
# Use QDA for analysis. QDA fits better and has better data assumptions. 


### Train / Test QDA Analysis
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
# Our sample size is large enough that we don't necessarily have to use train/test method



# --------------------------------------------#

# Manually recreate partimat() plot in ggplot
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
coho_scales_berners %>% group_by(Age, accuracy) %>% tally()
# Very good accuracy!



######### TRANSFORMED
# Explore using transformed versions of the variables
# To address concerns of normality, we should transform variables, using arcsine sqrt

coho_scales_fulldata <- coho_scales_fulldata %>%
  mutate(Q2plus_trans = asin(sqrt(Q2plus)),
         Q9abs_trans = asin(sqrt(Q9abs)))

coho_scales_berners1 <- coho_scales_fulldata %>% filter(Location == "BR")
coho_scales_hughsmith1 <- coho_scales_fulldata %>% filter(Location == "HS")

#QDA with Transformed Variables
fit_qda1 <- qda(as.factor(Age) ~ Location + Q2plus_trans + Q9abs_trans, data=coho_scales_fulldata,
                prior=c(0.66,0.34)) #prior prob 66% fish are Age1 (leave blank for uninformed)


testpred1 <- expand.grid(Q9abs_trans =seq(min(coho_scales_fulldata$Q9abs_trans), 
                                         max(coho_scales_fulldata$Q9abs_trans), length.out = 100), 
                         Q2plus_trans=seq(min(coho_scales_fulldata$Q2plus_trans), 
                                          max(coho_scales_fulldata$Q2plus_trans), length.out = 100),
                         Location = c("AL", "HS", "BR"))

testpred1$predval <- predict(fit_qda1, testpred1)$class

coho_scales_berners1 <- coho_scales_berners1 %>% 
  mutate(pred_age1 = predict(fit_qda1, coho_scales_berners1)$class,
         accuracy = ifelse(pred_age1 == Age, "Correct", "Incorrect"))

coho_scales_hughsmith1 <- coho_scales_hughsmith1 %>% 
  mutate(pred_age1 = predict(fit_qda1, coho_scales_hughsmith1)$class,
         accuracy = ifelse(pred_age1 == Age, "Correct", "Incorrect"))

coho_scales_berners %>% group_by(accuracy) %>% tally() 
coho_scales_berners1 %>% group_by(accuracy) %>% tally() #Slight improvement in accuracy

coho_scales_hughsmith %>% group_by(Age, accuracy) %>% tally()
coho_scales_hughsmith1 %>% group_by(Age, accuracy) %>% tally() #Slight improvement in accuracy


coho_scales_hughsmith1 %>% group_by(Age, accuracy) %>% tally() %>% spread(accuracy, n) %>%
  mutate(PercentCorrect = Correct / (Correct + Incorrect)) 

coho_scales_berners1 %>% group_by(Age, accuracy) %>% tally() %>% spread(accuracy, n) %>%
  mutate(PercentCorrect = Correct / (Correct + Incorrect)) 


#############
#predict using binomial
# --------------------------------------------#
### Binomial Regression ###


library(mgcv)
binomfit <- glm((Age-1)~Location+Q2+Q9abs, data=coho_scales_bothriv, family="binomial")
summary(binomfit)

testpred5 <- expand.grid(Q9abs=seq(min(coho_scales_bothriv$Q9abs), max(coho_scales_bothriv$Q9abs), length.out = 100), 
                         Q2=seq(min(coho_scales_bothriv$Q2), max(coho_scales_bothriv$Q2, na.rm = TRUE), length.out = 100),
                         Location=c("BR", "HS"))



coho_scales_berners1 <- coho_scales_bothriv %>% filter(Location == "BR")
coho_scales_berners1 <- coho_scales_berners1 %>% 
  mutate(pred_binom = predict(binomfit, coho_scales_berners1, type = "response"),
         pred_age5 = ifelse(pred_binom > 0.5, 2, 1),
         accuracy = ifelse(pred_age5 == Age, "Correct", "Incorrect"))

coho_scales_berners1 %>% group_by(accuracy) %>% tally()


coho_scales_hughsmith1 <- coho_scales_fulldata %>% filter(Location == "HS")
coho_scales_hughsmith1 <- coho_scales_hughsmith1 %>% 
  mutate(pred_binom = predict(binomfit, coho_scales_hughsmith1, type = "response"),
         pred_age5 = ifelse(pred_binom > 0.5, 2, 1),
         accuracy = ifelse(pred_age5 == Age, "Correct", "Incorrect"))

coho_scales_hughsmith1 %>% group_by(Age, accuracy) %>% tally()


plot(fitted(binomfit), residuals(binomfit), xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(binomfit), residuals(binomfit)))





###############
#What sorts of things cause our inaccuracy 
coho_scales_berners1 <- coho_scales_berners1 %>%
  mutate(acc1 = ifelse(accuracy == "Correct", 1, 0))
plot(coho_scales_berners1$acc1, coho_scales_berners1$dayofyear)

summary(lm(acc1 ~ dayofyear, data = coho_scales_berners1)) # signif
summary(lm(acc1 ~ Q1, data = coho_scales_berners1)) # signif
summary(lm(acc1 ~ Q3, data = coho_scales_berners1)) # insig, also multicollinearity
summary(lm(acc1 ~ Year, data = coho_scales_berners1)) #insig

coho_scales_hughsmith1 <- coho_scales_hughsmith1 %>%
  mutate(acc1 = ifelse(accuracy == "Correct", 1, 0))

summary(lm(acc1 ~ dayofyear, data = coho_scales_hughsmith1)) #marginally signif
summary(lm(acc1 ~ Q1, data = coho_scales_hughsmith1)) # very signif
summary(lm(acc1 ~ Q3, data = coho_scales_hughsmith1)) # insig, also multicollinearity
summary(lm(acc1 ~ Year, data = coho_scales_hughsmith1)) # insig



boxcox(coho_scales_hughsmith$Age~coho_scales_hughsmith$Q2plus)
boxcox(lm(Age~Q2plus, data=coho_scales_berners),lambda=seq(-10,1,by=.1))





#### FINAL SUMMARY #### 
# Recommend using variables Q2plus and Q9abs to predict age
# If using discriminate analysis, use a QDA to address homoscedascity 
# Do not transform data 








binomfit.hs <- glm((Age-1)~Q2+Q9abs, data=coho_scales_hughsmith, family="binomial")
summary(binomfit.hs)

coho_scales_berners2 <- coho_scales_bothriv %>% filter(Location == "BR")
coho_scales_berners2 <- coho_scales_berners2 %>% 
  mutate(pred_binom = predict(binomfit.hs, coho_scales_berners2, type = "response"),
         pred_age5 = ifelse(pred_binom > 0.5, 2, 1),
         accuracy = ifelse(pred_age5 == Age, "Correct", "Incorrect"))

coho_scales_berners2 %>% group_by(accuracy) %>% tally()





binomfit.br <- glm((Age-1)~Q2+Q9abs, data=coho_scales_berners, family="binomial")
summary(binomfit.br)

coho_scales_hughsmith2 <- coho_scales_bothriv %>% filter(Location == "HS")
coho_scales_hughsmith2 <- coho_scales_hughsmith2 %>% 
  mutate(pred_binom = predict(binomfit.br, coho_scales_hughsmith2, type = "response"),
         pred_age5 = ifelse(pred_binom > 0.5, 2, 1),
         accuracy = ifelse(pred_age5 == Age, "Correct", "Incorrect"))

coho_scales_hughsmith2 %>% group_by(accuracy) %>% tally()








