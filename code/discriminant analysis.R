#This is Kevin McNeil's example code
##############Discriminant Analysis###########
#data loading

datamap = read.table("clipboard", header=T, sep="\t")##pull in data from clipboard

data<-data1[,-c(1,2,7,11)] ##data without unused columns

head(data)                 ##look at first few specimens to checck

data$SPECIES<-factor(data$SPECIES)

New = read.table("clipboard", header=T, sep="\t") #Copy new data to check

#Your data should have columns named:  SPECIES WEIGHT LENGTH GENDER OTOLENGTH OTOHEIGHT OTOWEIGHT


###### For random Test#################################
###### LINEAR DISCRIM from R BOOK pg 821###############

load("G:\\ZZ McNeel\\R\\6 29 2015 152 151 discriminant")
test152<-New[1:20,]
test151<-New[21:39,]

pairs(data)
library (MASS)
library(klaR)
library(lattice)

##Alternative methods####################

partimat(factor(SPECIES2)~WEIGHT:LENGTH+GENDER:WEIGHT+OTOLENGTH:OTOHEIGHT+OTOWEIGHT:OTOLENGTH
,data=New1,method="lda",plot.matrix = TRUE) #### Visualize grouping

##modelq<-qda(factor(SPECIES)~.,data)  #Didn't improve the error based on partimat
##modelCV<- lda(factor(SPECIES)~.,data,CV=TRUE)   #leave one out cross validation

train<-sort(sample(1:nrow(data),round((nrow(data)/2),0))) #model development data
table(data$SPECIES[train])
unused<-data[-train,]     #withhled data for accuracy test 
table(unused$SPECIES)
plot(data)
model<-lda(factor(SPECIES)~GENDER:WEIGHT+LENGTH:WEIGHT+LENGTH:OTOWEIGHT+WEIGHT:OTOLENGTH+OTOLENGTH:OTOWEIGHT+OTOHEIGHT:OTOWEIGHT+OTOLENGTH:OTOHEIGHT+LENGTH+WEIGHT+OTOLENGTH+OTOHEIGHT+OTOWEIGHT+GENDER,data,subset=train)# model with subset

model2<-lda(factor(SPECIES)~GENDER+LENGTH+WEIGHT+OTOLENGTH+OTOHEIGHT+OTOWEIGHT,data,subset=train)# model with subset
model
model2
predict(model)$class
predict(model,unused)
predict(model,New)
New$LD1<-data.frame(predict(model,New)$x)
New$Species<-data.frame(c(
	"Shortraker","Shortraker","Shortraker","Shortraker","Shortraker",
	"Shortraker","Shortraker","Shortraker","Shortraker","Shortraker",
	"Shortraker","Shortraker","Shortraker","Shortraker","Shortraker",
	"Shortraker","Shortraker","Shortraker","Shortraker","Shortraker",
	"Rougheye","Rougheye","Rougheye","Rougheye","Rougheye",
	"Rougheye","Rougheye","Rougheye","Rougheye","Rougheye",
	"Rougheye","Rougheye","Rougheye","Rougheye","Rougheye",
	"Rougheye","Rougheye","Rougheye","Rougheye"))
ls(New)

# Assess the accuracy of the prediction on the current data
# percent correct for each category of SPECIES
 ct2 <- table(data$SPECIES[train], predict(model2)$class)
 diag(prop.table(ct2, 1))
 # total percent correct
 sum(diag(prop.table(ct2)))

# Assess the accuracy of the prediction
# percent correct for each category of SPECIES
 ct <- table(unused$SPECIES, predict(model,unused)$class)
 diag(prop.table(ct, 1))
# total percent correct
 sum(diag(prop.table(ct)))

plot(model2)


results<- data.frame(predict(model2, data1)$posterior) ###P(Class|X)
class<-data.frame(predict(model2,data1)$class)
names(results)<-c("prob151","prob152")
write.table(results, "clipboard", sep="\t") ### export x into clipboard
write.table(class, "clipboard", sep="\t") ### export x into clipboard

resultsN<- data.frame(predict(model, New)$x) ###P(Class|X)
classN<-data.frame(predict(model2,New)$class)
names(resultsN)<-c("prob151","prob152")
write.table(resultsN, "clipboard", sep="\t") ### export x into clipboard
write.table(classN, "clipboard", sep="\t") ### export x into clipboard


