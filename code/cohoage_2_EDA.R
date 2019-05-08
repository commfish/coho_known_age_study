library(corrplot)
library(GGally)
library(here)
library(mgcv)
library(mvnormtest)
library(HH)

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))


####
#load some functions for usage

hist.normal <- function(data, histcolumn, agefilter=0){ 
  # This function creates a histogram of data, overlaid with a normal curve (to show departures from normal)
  # Note that histcolumn needs to be in quotes when run. Leave agefilter blank to show both ages
  if (agefilter != 1 & agefilter != 2){
    .cleandata <- data %>% filter(!is.na(Length))
    .agenum <- "Both ages"
  }
  else{
    .cleandata <- data %>% filter(Age == agefilter, !is.na(Length))
    .agenum <- paste0("Age ", agefilter)
  }
  .cleandata <- .cleandata %>% mutate(newcol=(!!as.name(histcolumn))) %>% 
    dplyr::select(newcol)
  h <- hist(.cleandata$newcol, breaks=20, 
            main = paste0(deparse(substitute(data)), "\n", histcolumn, " hist - ", .agenum))
  
  xfit <- seq(min(.cleandata$newcol), max(.cleandata$newcol), length = 40)
  yfit <- dnorm(xfit, mean = mean(.cleandata$newcol), sd = sd(.cleandata$newcol)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(.cleandata$newcol) 
  lines(xfit, yfit, col = "red", lwd = 2)
}


eda.norm <- function(x, ...) #Function from Franz Mueter
{
  par(mfrow=c(2,2))
  if(sum(is.na(x)) > 0)
    warning("NA's were removed before plotting")
  x <- x[!is.na(x)]
  hist(x, main = "Histogram and non-\nparametric density estimate", prob = T)
  iqd <- summary(x)[5] - summary(x)[2]
  lines(density(x, width = 2 * iqd))
  boxplot(x, main = "Boxplot", ...)
  qqnorm(x)
  qqline(x)
  plot.ecdf(x, main="Empirical and normal cdf")
  LIM <- par("usr")
  y <- seq(LIM[1],LIM[2],length=100)
  lines(y, pnorm(y, mean(x), sqrt(var(x))))
  shapiro.test(x)
}


#### 
# BEGIN EDA
####

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
  for(i in 2:28){.temp4[,i] <- shift(.temp4[,i], i-1)}
  
  .temp4["rownum"] <- seq(1:28)
  .temp4 <- .temp4 %>% gather("Qnum", "corr", Q4:Q31)
  
  print(ggplot(.temp4, aes(x=rownum, y = corr, color=Qnum)) +
    geom_line() + 
    labs(title = q))
}

#caution, this takes a while to eval!
#ggpairs(coho_scales_berners %>% dplyr::select("Zone1":"Count_Zone2", "Year", "Age":"Q10"))


# Explore a few plots and model types to get a general feel for the data
ggplot(coho_scales_berners, aes(y = Q9, group = Age)) +
  geom_boxplot() #clear difference in Q9 between ages

summary(lm(data=coho_scales_aukelake, Age~Length+Q5+Q12))
summary(gam(data=coho_scales_aukelake, Age ~ Length + Q5))
summary(gam(data=coho_scales_aukelake, Age ~ Length + s(Q5))) # deviance explained rises quite a bit
summary(gam(data=coho_scales_aukelake, (Age-1) ~ Length + Q5, family = "binomial"))
summary(glm(data=coho_scales_aukelake, (Age-1) ~ log(Length) + log(Q5), family = "binomial")) # Top linear model (Len&Q5 only)
summary(gam(data=coho_scales_aukelake, (Age-1) ~ Length + s(Q5), family = "binomial")) # Top non-linear model (Len&Q5 only)



test1 <- coho_scales_aukelake %>% filter(!is.na(Length)) #remove NAs from Length to test corr with other vars
cor(test1$Q5, test1$Length) # High corr
cor(test1$Q9, test1$Length) # Top predictor variables are very correlated
cor(test1$Q5, test1$Age) 
cor(test1$Length, test1$Age) # But are also good predictors of Age
rm(test1)

hist(coho_scales_aukelake$Length, breaks=15)
hist((coho_scales_aukelake %>% filter(Age==1))$Length, breaks=15)
hist((coho_scales_aukelake %>% filter(Age==2))$Length, breaks=15)
hist(coho_scales_aukelake$Q9, breaks=15)
hist(log(coho_scales_aukelake$Length), breaks=15)
hist(log((coho_scales_aukelake %>% filter(Age==1))$Length), breaks=15)
hist(log((coho_scales_aukelake %>% filter(Age==2))$Length), breaks=15)
hist(log(coho_scales_aukelake$Q9), breaks=15)



mvnormtest::mshapiro.test(t(coho_scales_aukelake %>% dplyr::select(Age, Length, Q2, Q5, Q9) %>% na.omit %>% as.matrix()))

shapiro.test(coho_scales_aukelake$Q9) # Not normal
shapiro.test(log(coho_scales_aukelake$Q9)) # Normal
shapiro.test(log(coho_scales_aukelake$Length))
shapiro.test((coho_scales_aukelake$Length)^0.25)

shapiro.test((coho_scales_aukelake$Length[coho_scales_aukelake$Age == 1])) #Not normal, but slightly better log trans
shapiro.test((coho_scales_aukelake$Length[coho_scales_aukelake$Age == 2])) #Not normal, but slightly better log trans
shapiro.test((coho_scales_berners$Length[coho_scales_berners$Age == 1])) # VERY Not normal, but slightly better log trans
shapiro.test((coho_scales_berners$Length[coho_scales_berners$Age == 2])) #Not normal, but worse with log trans
shapiro.test((coho_scales_hughsmith$Length[coho_scales_hughsmith$Age == 1])) # Normal!
shapiro.test((coho_scales_hughsmith$Length[coho_scales_hughsmith$Age == 2])) # Normal


hist((coho_scales_berners$Length[coho_scales_berners$Age == 1]), breaks = 25)


bartlett.test(log(Q5)~Age, data=coho_scales_aukelake)
fligner.test(Length~Age, data=coho_scales_aukelake) 


boxplot(Length~Age, data = coho_scales_aukelake)
cov(coho_scales_aukelake$Q5, coho_scales_aukelake$Age)
cov(coho_scales_aukelake$Length[!is.na(coho_scales_aukelake$Length)], 
    coho_scales_aukelake$Age[!is.na(coho_scales_aukelake$Length)])

HH::hov(Q9~as.factor(Age), data=coho_scales_aukelake)
hov(Length~as.factor(Age), data=coho_scales_aukelake %>% filter(!is.na(Length)))

hovPlot(Q5~as.factor(Age), data=coho_scales_aukelake)
hovPlot(Length~as.factor(Age), data=coho_scales_aukelake %>% filter(!is.na(Length)))

##################
# explore how not normal it is


hist.normal(coho_scales_aukelake, "Length", 1)
hist.normal(coho_scales_aukelake, "Length", 2)
hist.normal(coho_scales_aukelake, "Length")
hist.normal(coho_scales_aukelake, "Q9", 1)
hist.normal(coho_scales_aukelake, "Q9", 2)
hist.normal(coho_scales_aukelake, "Q9")

hist.normal(coho_scales_berners, "Length", 1)
hist.normal(coho_scales_berners, "Length", 2)
hist.normal(coho_scales_berners, "Length")
hist.normal(coho_scales_berners, "Q9", 1)
hist.normal(coho_scales_berners, "Q9", 2)
hist.normal(coho_scales_berners, "Q9")

hist.normal(coho_scales_hughsmith, "Length", 1)
hist.normal(coho_scales_hughsmith, "Length", 2)
hist.normal(coho_scales_hughsmith, "Length")
hist.normal(coho_scales_hughsmith, "Q9", 1)
hist.normal(coho_scales_hughsmith, "Q9", 2)
hist.normal(coho_scales_hughsmith, "Q9")

# All in all, most times there is not a large violation of normality, at least visually

###### Check multivariate normality
ggplot(coho_scales_fulldata %>% filter(Age == 1), aes(Length, Q9)) + 
  stat_bin2d(bins = 20) +
  coord_cartesian(xlim = c(60, 140), ylim = c(0, 1)) +
  facet_wrap(~Location)

ggplot(coho_scales_fulldata %>% filter(Age == 2), aes(Length, Q9)) + 
  stat_bin2d(bins = 20) +
  coord_cartesian(xlim = c(60, 140), ylim = c(0, 1)) +
  facet_wrap(~Location)

# Large differences by location. Clear correlation. Not all combinations are bad, Age2s are pretty good



# Exploratory
boxplot(Q2~Age, data = coho_scales_berners)
boxplot(Length~Age, data = coho_scales_berners)
summary(lm(Age ~ Length, data = coho_scales_berners))
summary(lm(Age ~ Q2, data = coho_scales_berners))
summary(lm(Age ~ Q9, data = coho_scales_berners))

hist.normal(data = coho_scales_berners, "Q2plus", agefilter = 1)

#################################
# EDA Overview
# Many of the variables are approximately normal but fail the 'strict' conditions of tests
# Variance is not even between Age 1 & 2
# It appears that the higher numbered Q variables are not as good of a fit. 
# Explore further with Q variables <Q12 or so

