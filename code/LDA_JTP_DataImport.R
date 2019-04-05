#JTP: Re-working this code to be tidyverse compliant, portable ("here"), and removing plyr/reshape2, reshape
#This is converted SAS code that I worked on awhile back (SEM-2-28-2019)
#Converted SAS code to R code: Coho Smolt Scale Data-LDA
#**************************************************************************************************
#PART I:Import data and create dataset by one row/circuli and one row/zone instead of one row/fish
library(here)
library(tidyverse)
library(lubridate)
library(reshape2)

library(reshape)
library (plyr)
library (dplyr)
library (MASS)
library (vegan)

# Will need this function for converting year correctly
convyear <- function(x, year=2000){ m <- year(x) %% 100
  year(x) <- ifelse(m > year %% 100, 1900+m, 2000+m)
  x
}

# First, read in the data and make sure that the dates import correctly. 
coho_scales <- read.csv(here::here("data/AL_BR_HS.csv"), stringsAsFactors = FALSE) %>%
  mutate(Date = ymd(convyear(strptime(Date, format = "%d-%b-%Y", tz="US/Alaska"))))



A1_JTP <- coho_scales %>% dplyr::select("IMAGENAME":"Comment") %>% # InClude only the first 11 columns
  bind_cols(coho_scales %>% dplyr::select(starts_with("Z"))) %>% # Select only the Z columns, then put back together
  gather(key = "Variable", value = "Zone", "Z1":"Z41") # turn from wide to long, each row is a zone

A2_JTP <- coho_scales %>% dplyr::select("IMAGENAME":"Comment") %>% # InClude only the first 11 columns
  bind_cols(coho_scales %>% dplyr::select(-"Comment") %>% dplyr::select(starts_with("C"))) %>% # Select only the C cols, then bind
  gather(key = "Circulus", value = "Distance2", "C1":"C41") # turn from wide to long, each row is a zone

coho_scales_long <- A2_JTP %>% bind_cols(A1_JTP %>% dplyr::select("Zone", "Variable")) %>%
  filter(Circulus != "C1", Circulus != "C2", Age != 3) %>%
  # Drop C1 (dist from focus to C1), C2 (dist from C1 to C2), and age 3 fish (too few samples)
  mutate(Distance = ifelse(Age==1 & Zone==1, Distance2, # Add a new column "Distance" that excludes some distance/age combos
                           ifelse(Age==2 & Zone<=2, Distance2, NA))) %>% 
  arrange(Sample_ID, Circulus) %>% # Sort (order) the data by Sample_ID then Circulus 
  dplyr::select(-"Distance2") %>% # Distance2 is just the plus growth now, drop it
  drop_na("Distance") # remove rows with NAs for Distance

jj_JTP <- coho_scales_long %>% dcast(Sample_ID ~ Circulus, value.var = "Distance", fun = sum) %>% # Now that data are correct, turn back to wide
  select(Sample_ID, num_range("C", range = 3:41)) %>% # Grab just the Sample_ID and "C" columns
  left_join(coho_scales_long %>% dcast(Sample_ID ~ Variable, value.var = "Zone") %>% # Join these with the correct wide data...
            select(Sample_ID, num_range("Z", range = 3:41)), by = c("Sample_ID" = "Sample_ID")) %>% # ... grabbing just Z cols
  left_join(coho_scales %>% dplyr::select("IMAGENAME":"Comment"), by = c("Sample_ID" = "Sample_ID")) %>% # Now merge back w/ header data
  dplyr::select("Sample_ID", "IMAGENAME":"Comment", everything()) # Arrange col order to be as preferred
#rm(A1_JTP, A2_JTP, coho_scales_long)
#new dataset jj_JTP is without plus groups or distance from focus to C2

# NOTE: The long dataframe is different slightly than the same length dataframe. It appears that melt truncates results to 128000 while gather can handle all 
# Runs on just four packages (tidyverse, lubridate, reshape2, and here). Removing reshape functionality.




#**************************************************************************************************
#PART II: Summarize data by sample ID, zone, and circuli distance and count

temp2 <- coho_scales_long %>% filter(Zone == 1 | Zone ==2) %>% group_by(Sample_ID, Zone) %>% 
  summarise(freq = length(Distance), Distance = sum (Distance)) %>% #does not include plus group or distance from focus to C2
  mutate(Zone = replace(Zone, Zone==1, "Zone1"),
         Zone = replace(Zone, Zone==2, "Zone2")) %>%
  rename(value = Distance) 

j_JTP <- left_join(temp2 %>% dplyr::select(-freq) %>% spread(Zone, value = value), #summary of circuli distance by zone and sample_ID
                   temp2 %>% dplyr::select(-value) %>% rename(value = freq) %>%
                   spread(Zone, value = value) %>% 
                   rename(Count_Zone1 = Zone1, Count_Zone2 = Zone2), #count of circuli distance by zone and fish ID
                   by = c("Sample_ID" = "Sample_ID")) %>%
  left_join(jj_JTP, by = c("Sample_ID" = "Sample_ID"))
#new dataset with zone summaries and C1... and Z1...



#**************************************************************************************************
# PART III: Calculate step#1 for variables Q32 and Q33
# PART IV: Merge full dataset with summarized dataset by Sample_ID

j_JTP <- j_JTP %>% dplyr::select(Sample_ID:C40) %>% gather(key = "Varr", value = "Value", "C3":"C40") %>%
  filter(Value != 0, !is.na(Value)) %>%  # Filter to remove 0 & NAs
  mutate(Circulus = as.numeric(gsub("^C", '', Varr)), Distance = Value) %>% #Make "Circulus" just the numbers (ie drop "C") 
  dplyr::select(-Varr, -Value) %>% #drop these two cols
  mutate(NCFAZ = ifelse(Age == 1, Count_Zone1, # Add column NCFAZ
         ifelse(Age == 2, Count_Zone1 + Count_Zone2, NA)),
         NCFAZ_adj = ifelse(Age == 1, Count_Zone1 + 2, # Adjust NCFAZ for focus to C2; data pairs variable includes plus group
                            ifelse(Age == 2, Count_Zone1 + Count_Zone2 + 2, NA)),
         NCFAZ_6 = NCFAZ_adj-5, 
         NCFAZ_7 = NCFAZ_adj-6) %>% # Minimum circulus to count
  arrange(Sample_ID, Circulus) %>%
  mutate(Q32 = ifelse(Circulus >= NCFAZ_6 & Circulus <= NCFAZ_adj, Distance,0),
         Q33 = ifelse(Circulus >= NCFAZ_7 & Circulus <= NCFAZ_adj, Distance,0)) %>% #count distance if btwn NCFAZ_adj & NCFAZ_6)
  group_by(Sample_ID) %>% summarise(Q32_sum = sum(Q32), Q33_sum = sum(Q33)) %>% #dataset of distance from NCFAZ_adj, NCFAZ_6, NCFAZ_7 
  left_join(j_JTP, by = c("Sample_ID" = "Sample_ID"))
  # New dataset that includes Q32_sum and Q33_sum


#NEWTEST. IN PROG
temp4 <- j_JTP %>% dplyr::select(Sample_ID:C40) %>% gather(key = "Varr", value = "Value", "C3":"C40") %>%
  filter(Value != 0, !is.na(Value)) %>%  # Filter to remove 0 & NAs
  mutate(Circulus = as.numeric(gsub("^C", '', Varr)), Distance = Value) %>% #Make "Circulus" just the numbers (ie drop "C") 
  dplyr::select(-Varr, -Value) %>% #drop these two cols
  mutate(NCFAZ = ifelse(Age == 1, Count_Zone1, # Add column NCFAZ
                        ifelse(Age == 2, Count_Zone1 + Count_Zone2, NA)),
         NCFAZ_adj = ifelse(Age == 1, Count_Zone1 + 2, # Adjust NCFAZ for focus to C2; data pairs variable includes plus group
                            ifelse(Age == 2, Count_Zone1 + Count_Zone2 + 2, NA)),
         NCFAZ_6 = NCFAZ_adj-5, 
         NCFAZ_7 = NCFAZ_adj-6,# Minimum circulus to count
         SFAZ = ifelse(Age == 1, Zone1,
                       ifelse(Age == 2, Zone1 + Zone2, NA)), #total distance of circuli by Sample_ID (no plus group or focus to C2)
         SFAZ_0.5 = SFAZ * 0.5) %>% # Distance of half SFAZ
  arrange(Sample_ID, Circulus) %>% 
  mutate(Cum.sum = cumsum(Distance) )

  mutate(Q32 = ifelse(Circulus >= NCFAZ_6 & Circulus <= NCFAZ_adj, Distance,0),
         Q33 = ifelse(Circulus >= NCFAZ_7 & Circulus <= NCFAZ_adj, Distance,0)) %>% #count distance if btwn NCFAZ_adj & NCFAZ_6)
  group_by(Sample_ID) %>% summarise(Q32_sum = sum(Q32), Q33_sum = sum(Q33)) %>% #dataset of distance from NCFAZ_adj, NCFAZ_6, NCFAZ_7 
  left_join(j_JTP, by = c("Sample_ID" = "Sample_ID"))



#**************************************************************************************************

#PART V: Calculate step#1 for variables Q34 & Q35
#VARIABLE Q34
A1 <- j[, c(1:15, 16:53)] 
A1 <- melt(A1, id=c(1:15)) #circulus by distance
A1$value[which(A1$value==0)] = NA
A1<-subset(A1, !is.na(value)) #delete rows with NA 
A1["Circulus"] <-as.numeric(gsub("^C", '', A1$variable))
A1["Distance"] <-A1$value
A1<- subset(A1, select = -c(variable,value))
A1["NCFAZ"] <-ifelse(A1$Age==1,A1$Count_Zone1,
                     ifelse(A1$Age==2,A1$Count_Zone1+A1$Count_Zone2,NA)) #number of circuli by Sample_ID (no plus group or focus to C2)
A1["SFAZ"] <-ifelse(A1$Age==1,A1$Zone1,
                    ifelse(A1$Age==2,A1$Zone1+A1$Zone2,NA)) #total distance of circuli by Sample_ID (no plus group or focus to C2)
A1["SFAZ_0.5"] <-0.5*A1$SFAZ #distance of half of SFAZ
A1["Circulus"]<-as.numeric(A1$Circulus)
A1 <- A1[order(A1$Sample_ID, A1$Circulus),]
A1<-ddply(A1, .(Sample_ID), transform, Cumulative.Sum = cumsum(Distance)) #cumulative sum of distances
A1["Circuli_SFAZ_0.5"] <-ifelse(A1$Cumulative.Sum<=A1$SFAZ_0.5,1,0) #count if distance is <=SFAZ*0.5
A1 <- aggregate(Circuli_SFAZ_0.5~Sample_ID, A1, FUN=sum, na.action=na.omit) #summarize number of circuli in first half of SFAZ

#VARIABLE Q35
B1 <- j[, c(1:15, 16:53)] 
B1 <- melt(B1, id=c(1:15))
B1$value[which(B1$value==0)] = NA
B1<-subset(B1, !is.na(value)) #delete rows with NA 
B1["Circulus"] <-as.numeric(gsub("^C", '', B1$variable))
B1["Distance"] <-B1$value
B1<- subset(B1, select = -c(variable,value))
B1["NCFAZ"] <-ifelse(B1$Age==1,B1$Count_Zone1,
                     ifelse(B1$Age==2,B1$Count_Zone1+B1$Count_Zone2,NA)) #number of circuli by Sample_ID (no plus group or focus to C2)
B1["SFAZ"] <-ifelse(B1$Age==1,B1$Zone1,
                    ifelse(B1$Age==2,B1$Zone1+B1$Zone2,NA))  #total distance of circuli by Sample_ID (no plus group or focus to C2)
B1["SFAZ_0.75"] <-0.75*(B1$SFAZ)#distance of 3/4 of SFAZ
B1["Circulus"]<-as.numeric(B1$Circulus)
B1 <- B1[order(B1$Sample_ID, -B1$Circulus),]
B1<-ddply(B1, .(Sample_ID), transform, Cumulative.Sum = cumsum(Distance))#cumulative sum of distances (start from end of zone)
B1["Circuli_SFAZ_0.75"] <-ifelse(B1$Cumulative.Sum<=B1$SFAZ_0.75,1,0) #count if distance is <=SFAZ*0.75
B1 <- aggregate(Circuli_SFAZ_0.75~Sample_ID, B1, FUN=sum)#summarize number of circuli in last 3/4 of SFAZ
Merge <- merge(A1,B1, by=c("Sample_ID")) #merge Q34 and Q35 variables by Sample_ID
j <- merge(j,Merge, by=c("Sample_ID")) 
rm(A1,B1) 

#**************************************************************************************************

#PART VI: Data check 
# If age 1 there must be zone 1 measurements also zone 2 is optional if plus is present. 
# There must not be any zones but 1 or 2.
# If age 2 there must be zone 1 and zone 2 measurements also zone 3 is optional if plus is present. 
# There must not be any zones but 1 and 2 or 3.

j["Zone1_Check"] <-ifelse(j$Zone1>0,1,0)
j["Zone2_Check"] <-ifelse(j$Zone2>0,1,0)
j["Sum_Zones"] <-j$Zone1_Check+j$Zone2_Check
j<- subset(j, select = -c(Zone1_Check, Zone2_Check))
j["Check"] <-      ifelse(j$Age==1 & j$Sum_Zones>2,"Need to Check",
                          ifelse(j$Age==2 & j$Sum_Zones>3,"Need to Check","OK"))
j["Check1"] <-      ifelse(j$Age==1 & j$Count_Zone2>0,"Need to Check","OK")
j["Check2"] <-      ifelse(j$Age==1 & j$Zone2>0,"Need to Check","OK")
newdata <- j[ which(j$Check=="Need to Check")] #output rows that have too many zones for the age
newdata1 <- j[ which(j$Check1=="Need to Check")] #output rows that have too many zones for the age
newdata2 <- j[ which(j$Check2=="Need to Check")] #output rows that have too many zones for the age
j<- subset(j, select = -c(Check))

#**************************************************************************************************

#PART VII: Calculate Variables (See Linear Discriminant Analysis-Project Overview document in S:\Region1Shared-DCF\Research\Salmon\Coho\Linear Discriminant Analysis.doc for list of variables)
#Note C5 distance is distance from C4 to C5 circulus, therefore to start at C4 you need to start at C5
options("na.actions"=na.omit)
j["Q1"] <-ifelse(j$Age==1,j$Count_Zone1,
                 ifelse(j$Age==2,j$Count_Zone1+j$Count_Zone2,NA)) #circuli count (no plus group or focus to C2)
j["Q2"] <-ifelse(j$Age==1,j$Zone1,
                 ifelse(j$Age==2,j$Zone1+j$Zone2,NA)) #circuli distance (no plus group or focus to C2)
j["Q3"] <-j$Q2/j$Q1
j["Q4"] <-j$C3/j$Q2
j["Q5"] <-(j$C4+j$C3)/j$Q2 
j["Q6"] <-(j$C5+j$C4+j$C3)/j$Q2     
j["Q7"] <-(j$C6+j$C5+j$C4+j$C3)/j$Q2  
j["Q8"] <-(j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q9"] <-(j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2      
j["Q10"] <-(j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q11"] <-(j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q12"] <-(j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q13"] <-(j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2     
j["Q14"] <-(j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2   
j["Q15"] <-(j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q16"] <-(j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2  
j["Q17"] <-(j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2      
j["Q18"] <-(j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q19"] <-(j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2  
j["Q20"] <-(j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2   
j["Q21"] <-(j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2  
j["Q22"] <-(j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q23"] <-(j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q24"] <-(j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2      
j["Q25"] <-(j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2  
j["Q26"] <-(j$C25+j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q27"] <-(j$C26+j$C25+j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q28"] <-(j$C27+j$C26+j$C25+j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2     
j["Q29"] <-(j$C28+j$C27+j$C26+j$C25+j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2  
j["Q30"] <-(j$C29+j$C28+j$C27+j$C26+j$C25+j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 
j["Q31"] <-(j$C30+j$C29+j$C28+j$C27+j$C26+j$C25+j$C24+j$C23+j$C22+j$C21+j$C20+j$C19+j$C18+j$C17+j$C16+j$C15+j$C14+j$C13+j$C12+j$C11+j$C10+j$C9+j$C8+j$C7+j$C6+j$C5+j$C4+j$C3)/j$Q2 

j["Q32"] <-j$Q32_sum/j$Q2 
j["Q33"] <-j$Q33_sum/j$Q2 
j["Q34"] <-j$Circuli_SFAZ_0.5/(j$Q1-j$Circuli_SFAZ_0.5)
j["Q35"] <-j$Circuli_SFAZ_0.75/(j$Q1-j$Circuli_SFAZ_0.75)

j["Q36"] <-(j$C7+j$C6+j$C5)/j$Q2 
j["Q37"] <-(j$C8+j$C7+j$C6)/j$Q2 
j["Q38"] <-(j$C9+j$C8+j$C7)/j$Q2 
j["Q39"] <-(j$C10+j$C9+j$C8)/j$Q2 
j["Q40"] <-(j$C11+j$C10+j$C9)/j$Q2 
j["Q41"] <-(j$C12+j$C11+j$C10)/j$Q2 
j["Q42"] <-(j$C13+j$C12+j$C11)/j$Q2 
j["Q43"] <-(j$C14+j$C13+j$C12)/j$Q2 
j["Q44"] <-(j$C15+j$C14+j$C13)/j$Q2 
j["Q45"] <-(j$C16+j$C15+j$C14)/j$Q2 
j["Q46"] <-(j$C17+j$C16+j$C15)/j$Q2 
j["Q47"] <-(j$C18+j$C17+j$C16)/j$Q2 
j["Q48"] <-(j$C19+j$C18+j$C17)/j$Q2  
j["Q49"] <-(j$C20+j$C19+j$C18)/j$Q2 
j["Q50"] <-(j$C21+j$C20+j$C19)/j$Q2 
j["Q51"] <-(j$C22+j$C21+j$C20)/j$Q2 
j["Q52"] <-(j$C23+j$C22+j$C21)/j$Q2 
j["Q53"] <-(j$C24+j$C23+j$C22)/j$Q2 
j["Q54"] <-(j$C25+j$C24+j$C23)/j$Q2 
j["Q55"] <-(j$C26+j$C25+j$C24)/j$Q2 
j["Q56"] <-(j$C27+j$C26+j$C25)/j$Q2 
j["Q57"] <-(j$C28+j$C27+j$C26)/j$Q2 
j["Q58"] <-(j$C29+j$C28+j$C27)/j$Q2  
j["Q59"] <-(j$C30+j$C29+j$C28)/j$Q2 
j["Q60"] <-(j$C31+j$C30+j$C29)/j$Q2 
j["Q61"] <-(j$C32+j$C31+j$C30)/j$Q2 
j["Q62"] <-(j$C33+j$C32+j$C31)/j$Q2 
j["Q63"] <-(j$C34+j$C33+j$C32)/j$Q2 
j["Q64"] <-(j$C35+j$C34+j$C33)/j$Q2 
j["Q65"] <-(j$C36+j$C35+j$C34)/j$Q2  
j["Q66"] <-(j$C37+j$C36+j$C35)/j$Q2 
j["Q67"] <-(j$C38+j$C37+j$C36)/j$Q2 
j["Q68"] <-(j$C39+j$C38+j$C37)/j$Q2 
j["Q69"] <-(j$C40+j$C39+j$C38)/j$Q2 
j["Q70"] <-(j$C41+j$C40+j$C39)/j$Q2 
j<-j[,c(1:93, 99:167,94:98)] #reorganize variable order


#PART VIII: Create variables Q71:Q72 
#VARIABLE Q71
A1 <- j[, c(1, 16:54, 96)] 
A1 <- melt(A1, id=c(1,41))
A1["Circulus"] <-as.numeric(gsub("^C", '', A1$variable))#delete C from C3...
A1["Distance"] <-A1$value
A1<- subset(A1, select = -c(variable,value))
A1["Q3_0.75"] <-j$Q3*0.75 #circulus distance of 0.75 * Q3
A1 <- A1[order(A1$Sample_ID, A1$Circulus),]
A1["Count_Q3_0.75"] <-ifelse(A1$Distance<A1$Q3_0.75,1,0) #if distance is <0.75*Q3, then count
A1 <- aggregate(Count_Q3_0.75~Sample_ID, A1, FUN=sum, na.action=na.omit) #count; aggregate by Sample_ID
A1["Q71"] <-A1$Count_Q3_0.75
A1 <- subset(A1, select = -c(Count_Q3_0.75)) #delete unneeded variables
j <- merge(j,A1, by=c("Sample_ID"))
rm(A1)#remove datasets

#VARIABLE Q72
A1 <- j[, c(1, 16:54, 96)] 
A1 <- melt(A1, id=c(1,41))
A1["Circulus"] <-as.numeric(gsub("^C", '', A1$variable))#delete C from C1...
A1["Distance"] <-A1$value
A1<- subset(A1, select = -c(variable,value))
A1["Q3_1.25"] <-j$Q3*1.25
A1 <- A1[order(A1$Sample_ID, A1$Circulus),]
A1["Count_Q3_1.25"] <-ifelse(A1$Distance>A1$Q3_1.25,1,0) #if distance is >1.25*Q3, then count
A1 <- aggregate(Count_Q3_1.25~Sample_ID, A1, FUN=sum, na.action=na.omit) #count; aggregate by Sample_ID
A1["Q72"] <-A1$Count_Q3_1.25
A1 <- subset(A1, select = -c(Count_Q3_1.25)) #delete unneeded variables
j <- merge(j,A1, by=c("Sample_ID"))
rm(A1)#remove datasets
j<-j[,c(1:162, 168, 169,163:167)] #reorganize variable order
j<- subset(j, select = -c(Sum_Zones))
write.csv(j, "H:\\Salmon\\LDA.csv")

#subset datasets for different analyses
AL_BR_HS<- subset(j, select = c(1:8, 92:137))#Full dataset (all areas)
AL<-subset(AL_BR_HS, AL_BR_HS$Location == 'AL') #lake system (Auke Lake)
BR<-subset(AL_BR_HS, AL_BR_HS$Location == 'BR') #river system (Berner's River)
HS<-subset(AL_BR_HS, AL_BR_HS$Location == 'HS') #Lake system (High Smith Lake)
Lake<-subset(AL_BR_HS, AL_BR_HS$Location == 'HS'| AL_BR_HS$Location =='AL') #lake versus river systems

#PART IX: Run LDA function  
#To perform LDA, one must ensure that the within-group covariance matrices of
#the explanatory variables are homogeneous, a condition that is frequently violated
#with ecological data ;If one the groups defined by the dependent variable has greater 
#dispersion than others, cases will tend to be overclassified in it.

#function betadisper in package vegan
#MASS packake LDA and QDA or cluster analysis?

