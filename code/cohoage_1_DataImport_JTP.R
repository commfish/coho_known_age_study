#IN PROGRESS. This will be removed once code is completely updated. 
#' ---
#' title: "Coho Known Age Data Import and Cleanup"
#' author: "Sara Miller & Justin Priest"
#' date: "`r format(Sys.Date())`"
#' ---
#' 
# Script 1 - Data Import and Cleaup
# Email addresses: Sara Miller (sara.miller@alaska.gov) and Justin Priest (justinpriest.ak@gmail.com)
# JTP: Re-worked this code to be tidyverse compliant, portable ("here"), and remove plyr/reshape2, reshape

# SEM: This is converted SAS code that I worked on awhile back (SEM-2-28-2019)
# Converted SAS code to R code: Coho Smolt Scale Data-LDA

# Data info: These data come from Coho Salmon scales which were collected between 1994 and 2005 by ADF&G
# at three locations: Auke Lake (AL), Berner's River (BR), and Hugh Smith Lake (HS)
# Scale data were then read by the ADF&G Mark, Tag, and Age Lab
# Define "Circuli"
# Define "Zones"

#PART I:Import data and create dataset by one row/circuli and one row/zone instead of one row/fish
library(here)
library(tidyverse)
library(lubridate)
source("code/functions.R")

# Import data
coho_scales <- read.csv(here::here("data/AL_BR_HS.csv"), stringsAsFactors = FALSE) %>%
  mutate(Date = ymd(convyear(strptime(Date, format = "%d-%b-%Y", tz="US/Alaska"), 
                             year(strptime("10-Jun-05", format = "%d-%b-%Y")))),
         dayofyear = yday(Date)) %>%
  dplyr::select(IMAGENAME:Date, dayofyear, everything()) %>%
  drop_na(Data_Pairs) # Exclude rows with no distances measured

coho_scales %>% 
  mutate (include = ifelse(Sample_ID == 'BR04-0453', 0, 1))-> coho_scales
coho_scales %>% 
  filter (include == 1) %>%
  dplyr::select(-"include") -> coho_scales

A1_JTP <- coho_scales %>% 
  dplyr::select("IMAGENAME":"Comment") %>% # Include only the first 11 columns
  bind_cols(coho_scales %>% 
              dplyr::select(starts_with("Z"))) %>% # Select only the Z columns, then put back together
  gather(key = "Variable", value = "Zone", "Z1":"Z41") # turn from wide to long, each row is a zone

A2_JTP <- coho_scales %>% 
  dplyr::select("IMAGENAME":"Comment") %>% # Include only the first 11 columns
  bind_cols(coho_scales %>% 
              dplyr::select(-"Comment") %>% 
              dplyr::select(starts_with("C"))) %>% # Select only the C cols, then bind
  gather(key = "Circulus", value = "Distance2", "C1":"C41") # turn from wide to long, each row is a zone

coho_scales_long <- A2_JTP %>% 
  bind_cols(A1_JTP %>% 
              dplyr::select("Zone", "Variable")) %>%
  filter(Circulus != "C1", Circulus != "C2", Age != 3) %>% # Drop C1 (dist from focus to C1), C2 (dist from C1 to C2), and age 3 fish (too few samples)
  arrange(Sample_ID, Circulus)%>% # Sort (order) the data by Sample_ID then Circulus 
  mutate(PlusGrowth = ifelse(Distance2 > 0,
                             ifelse(is.na(Zone) | Zone=="", "plusgrowth", paste0("Zone", Zone)), Distance2))%>%
  #JTP added this section, 5/5 to account for plus growth
  drop_na("Distance2") %>%
  mutate(Distance = Distance2) %>%# remove rows with NAs for Distance
  dplyr::select(-"Distance2") %>%
  dplyr::select("IMAGENAME", "Sample_ID",	"Location",	"Year",	"Date",	"dayofyear"	,"Floy_No",	"Tag_Code",	"Age"	,"Length",
                "Sublocation",	"Comment",	"Circulus",	"Zone","Variable",	"Distance","PlusGrowth")
rm(A1_JTP, A2_JTP)
#write.csv(write.csv(coho_scales_long, "data/check.csv")) #outputs data for manual cehgck against original data

jj_JTP <- coho_scales_long %>% 
  dplyr::select(-Variable, -Zone, -PlusGrowth) %>% # Drop these two cols so that spread works correctly
  spread(key = Circulus, value = Distance) %>%
  left_join(coho_scales_long %>% # Next, we'll merge with the Z cols
              dplyr::select(-c(IMAGENAME, Location:Comment, Circulus, Distance, Zone)) %>% # exclude extra header info
              spread(key = Variable, value = PlusGrowth), # Spread by Z cols
            by = c("Sample_ID" = "Sample_ID")) %>% # Merge with the C cols
  dplyr::select("Sample_ID", "IMAGENAME":"Comment", #Finally, reorder all the columns correctly
                num_range("C", range = 3:41), num_range("Z", range = 3:41))
  
#rm(A1_JTP, A2_JTP, coho_scales_long)
#new dataset jj_JTP is without plus groups or distance from focus to C2
#write.csv(write.csv(jj_JTP, "data/check2.csv"))
# NOTE: The JTP long dataframe is different slightly than the same long SEM dataframe. 
# Runs on just three packages (tidyverse, lubridate, and here). Removed reshape functionality.

#PART II: Summarize data by sample ID, zone, and circuli distance and count
temp2 <- coho_scales_long %>% dplyr::select(-Zone) %>%
  rename(Zone = PlusGrowth) %>%
  filter(Zone == "Zone1" | Zone == "Zone2" | Zone == "plusgrowth" ) %>% 
  group_by(Sample_ID, Zone) %>% 
  summarise(freq = length(Distance), Distance = sum (Distance)) %>% #does not include distance from focus to C2
  rename(value = Distance) 

j_JTP <- left_join(temp2 %>% 
                      dplyr::select(-freq) %>% 
                      spread(Zone, value = value) %>% #summary of circuli distance by zone and sample_ID
                      mutate(Zone2 = replace_na(Zone2, 0),
                             plusgrowth = replace_na(plusgrowth, 0)) %>% #JTP added this line Apr9,2019. Ensure that this is correct.
                      rename(Zoneplus = plusgrowth),
                    temp2 %>% 
                      dplyr::select(-value) %>% 
                      rename(value = freq) %>%
                      spread(Zone, value = value) %>% 
                      mutate(Zone2 = replace_na(Zone2, 0),
                             plusgrowth = replace_na(plusgrowth, 0)) %>% #JTP added this line Apr9,2019. Ensure that this is correct.
                      rename(Count_Zone1 = Zone1, Count_Zone2 = Zone2, 
                             Count_Plus = plusgrowth), #count of circuli distance by zone and sample ID
                    by = c("Sample_ID" = "Sample_ID")) %>%
  left_join(jj_JTP, by = c("Sample_ID" = "Sample_ID")) %>% 
  dplyr::select(Sample_ID, Zone1, Zone2, Zoneplus, Count_Zone1, Count_Zone2, Count_Plus, everything()) 

# PART III: Calculate step#1 for variables Q32 and Q33
# PART IV: Merge full dataset with summarized dataset by Sample_ID
# PART V: Calculate step#1 for variables Q34 & Q35
j_JTP <- j_JTP %>%
  dplyr::select(Sample_ID:C40) %>% 
  gather(key = "Varr", value = "Value", "C3":"C40") %>%
  filter(Value != 0, !is.na(Value)) %>%  # Filter to remove 0 & NAs
  mutate(Circulus = as.numeric(gsub("^C", '', Varr)), # Make "Circulus" just the numbers (ie drop "C")
         Distance = Value) %>%
  dplyr::select(-Varr, -Value) %>% #drop these two cols
  mutate(NCFAZ = ifelse(Age == 1, Count_Zone1, # Add column NCFAZ
                        ifelse(Age == 2, Count_Zone1 + Count_Zone2, NA)),
         NCFAZ_adj = ifelse(Age == 1, Count_Zone1 + 2, # Adjust NCFAZ for focus to C2; data pairs variable includes plus group
                            ifelse(Age == 2, Count_Zone1 + Count_Zone2 + 2, NA)),
         NCFAZ_6 = NCFAZ_adj-5,
         NCFAZ_7 = NCFAZ_adj-6,# Minimum circulus to count
         SFAZ = ifelse(Age == 1, Zone1,
                       ifelse(Age == 2, Zone1 + Zone2, NA)), #total distance of circuli by Sample_ID (no plus group or focus to C2)
         SFAZ_0.5 = SFAZ * 0.5, # Distance of half SFAZ
         SFAZ_0.75 = SFAZ * 0.75) %>% # Distance of 3/4 of SFAZ
  arrange(Sample_ID, Circulus) %>%
  group_by(Sample_ID) %>%
  mutate(Cum.sum_0.5 = cumsum(Distance)) %>% # Take cumsum distance by Sample_ID, inner to outer (used for SFAZ0.5)
  arrange(Sample_ID, -Circulus) %>%
  group_by(Sample_ID) %>%
  mutate(Cum.sum_0.75 = cumsum(Distance)) %>% # Take cum distance by Sample_ID, outer to inner (used for SFAZ0.75)
  mutate(Q32 = ifelse(Circulus >= NCFAZ_6 & Circulus <= NCFAZ_adj, Distance,0),
         Q33 = ifelse(Circulus >= NCFAZ_7 & Circulus <= NCFAZ_adj, Distance,0), #count distance if btwn NCFAZ_adj & NCFAZ_6)
         Circuli_SFAZ_0.5 = ifelse(Cum.sum_0.5 <= SFAZ_0.5, 1, 0), # Count if SFAZ_0.5 is less than cum sum
         Circuli_SFAZ_0.75 = ifelse(Cum.sum_0.75 <= SFAZ_0.75, 1, 0)) %>% # Count if SFAZ_0.75 is less than cum sum (outer to inner)
  group_by(Sample_ID) %>%
  summarise(Q32_sum = sum(Q32), Q33_sum = sum(Q33),
            Circuli_SFAZ_0.5 = sum(Circuli_SFAZ_0.5), # Variable Q34
            Circuli_SFAZ_0.75 = sum(Circuli_SFAZ_0.75)) %>% # Variable Q35 
  left_join(j_JTP, by = c("Sample_ID" = "Sample_ID")) %>%
  dplyr::select(Sample_ID, Zone1:Z40, Q32_sum, Q33_sum, Circuli_SFAZ_0.5, Circuli_SFAZ_0.75) # reorganize order

#PART VI: Data check 
# If age 1 there must be zone 1 measurements also zone 2 is optional if plus is present. 
# There must not be any zones but 1 or 2.
# If age 2 there must be zone 1 and zone 2 measurements also zone 3 is optional if plus is present. 
# There must not be any zones but 1 and 2 or 3.

j_JTP <- j_JTP %>% 
  mutate(Sum_Zones = ifelse(Zone1 > 0, 1, 0) + ifelse(Zone2 > 0, 1, 0),
         Check = ifelse(Age == 1 & Sum_Zones > 2,"Need to Check",
                        ifelse(Age == 2 & Sum_Zones > 3,"Need to Check","OK")),
         Check1 = ifelse(Age == 1 & Count_Zone2 > 0,"Need to Check","OK"),
         Check2 = ifelse(Age == 1 & Zone2 > 0,"Need to Check","OK")) 
  
j_JTP %>% filter(Check == "Need to Check" | Check1 == "Need to Check" | Check2 == "Need to Check") %>% 
   head()
# Empty! S U C C E S S F U L L E X E C U T I O N
# JTP FIX
j_JTP <- j_JTP %>% 
  dplyr::select(-Check)

#PART VII: Calculate Variables (See Linear Discriminant Analysis-Project Overview document
#in S:\Region1Shared-DCF\Research\Salmon\Coho\Coho_Known_Age_study\coho_known_age_study\references\reports\Linear Discriminant Analysis.doc for list of variables)
#Note C5 distance is distance from C4 to C5 circulus, therefore to start at C4 you need to start at C5
options("na.actions"= na.omit)

j_JTP <- j_JTP %>% 
  mutate(Q1 = ifelse(Age == 1, Count_Zone1,
                     ifelse(Age == 2, Count_Zone1 + Count_Zone2, NA)), #circuli count (no plus group or focus to C2),
         Q1plus = Count_Zone1 + Count_Zone2 + Count_Plus,
         Q2 = ifelse(Age == 1, Zone1,
                     ifelse(Age == 2, Zone1 + Zone2, NA)), # circuli distance (no plus group or focus to C2)
         Q2plus = Zone1 + Zone2 + Zoneplus,
         Q3 = Q2 / Q1,
         Q4 = C3 / Q2, 
         Q5 = f_sum(., C3, C4, Q2), # This function sums Columns C3 thru C4 and divides by Q2
         Q6 = f_sum(., C3, C5, Q2), # This function sums Columns C3 thru C5 and divides by Q2
         Q7 = f_sum(., C3, C6, Q2), # etc.
         Q8 = f_sum(., C3, C7, Q2), # The end result is that Q8 is Distance to Circuli 7 (excluding Circ1&2) / tot dist (exc Circ1&2)
         Q9 = f_sum(., C3, C8, Q2),
         Q9abs = Q9 * Q2,
         Q10 = f_sum(., C3, C9, Q2),
         Q11 = f_sum(., C3, C10, Q2),
         Q12 = f_sum(., C3, C11, Q2),
         Q13 = f_sum(., C3, C12, Q2),
         Q14 = f_sum(., C3, C13, Q2),
         Q15 = f_sum(., C3, C14, Q2),
         Q16 = f_sum(., C3, C15, Q2),
         Q17 = f_sum(., C3, C16, Q2),
         Q18 = f_sum(., C3, C17, Q2),
         Q19 = f_sum(., C3, C18, Q2),
         Q20 = f_sum(., C3, C19, Q2),
         Q21 = f_sum(., C3, C20, Q2),
         Q22 = f_sum(., C3, C21, Q2),
         Q23 = f_sum(., C3, C22, Q2),
         Q24 = f_sum(., C3, C23, Q2),
         Q25 = f_sum(., C3, C24, Q2),
         Q26 = f_sum(., C3, C25, Q2),
         Q27 = f_sum(., C3, C26, Q2),
         Q28 = f_sum(., C3, C27, Q2),
         Q29 = f_sum(., C3, C28, Q2),
         Q30 = f_sum(., C3, C29, Q2),
         Q31 = f_sum(., C3, C30, Q2),
         Q32 = Q32_sum / Q2,
         Q33 = Q33_sum / Q2,
         Q34 = Circuli_SFAZ_0.5 / (Q1 - Circuli_SFAZ_0.5),
         Q35 = Circuli_SFAZ_0.75 / (Q1 - Circuli_SFAZ_0.75),
         Q36 = (C7 + C6 + C5) / Q2,
         Q37 = (C8 + C7 + C6) / Q2,
         Q38 = (C9 + C8 + C7) / Q2,
         Q39 = (C10 + C9 + C8) / Q2,
         Q40 = (C11 + C10 + C9) / Q2,
         Q41 = (C12 + C11 + C10) / Q2,
         Q42 = (C13 + C12 + C11) / Q2,
         Q43 = (C14 + C13 + C12) / Q2,
         Q44 = (C15 + C14 + C13) / Q2,
         Q45 = (C16 + C15 + C14) / Q2,
         Q46 = (C17 + C16 + C15) / Q2,
         Q47 = (C18 + C17 + C16) / Q2,
         Q48 = (C19 + C18 + C17) / Q2,
         Q49 = (C20 + C19 + C18) / Q2,
         Q50 = (C21 + C20 + C19) / Q2,
         Q51 = (C22 + C21 + C20) / Q2,
         Q52 = (C23 + C22 + C21) / Q2,
         Q53 = (C24 + C23 + C22) / Q2,
         Q54 = (C25 + C24 + C23) / Q2,
         Q55 = (C26 + C25 + C24) / Q2,
         Q56 = (C27 + C26 + C25) / Q2,
         Q57 = (C28 + C27 + C26) / Q2,
         Q58 = (C29 + C28 + C27) / Q2,
         Q59 = (C30 + C29 + C28) / Q2,
         Q60 = (C31 + C30 + C29) / Q2,
         Q61 = (C32 + C31 + C30) / Q2,
         Q62 = (C33 + C32 + C31) / Q2,
         Q63 = (C34 + C33 + C32) / Q2,
         Q64 = (C35 + C34 + C33) / Q2,
         Q65 = (C36 + C35 + C34) / Q2,
         Q66 = (C37 + C36 + C35) / Q2,
         Q67 = (C38 + C37 + C36) / Q2,
         Q68 = (C39 + C38 + C37) / Q2,
         Q69 = (C40 + C39 + C38) / Q2) 

j_JTP <- j_JTP %>% 
  dplyr::select(Sample_ID:Q33_sum, Q1:Q69, Circuli_SFAZ_0.5:Check2) 

#PART VIII: Create variables Q71:Q72 
coho_scales_fulldata <- j_JTP %>% 
  dplyr::select("Sample_ID", "C3":"C40", "Q3") %>% # educated guess as to using Q3 since it appears to be called later
  gather(key = "C_num", value = "Distance", "C3":"C40") %>% 
  drop_na(Distance) %>% #Note that this is different than how SEM did it
  mutate(Circulus = as.numeric(gsub("^C", '', C_num)), # Make "Circulus" just the numbers (ie drop "C")
         Q3_0.75 = Q3 * 0.75, 
         Count_Q3_0.75 = ifelse(Distance < Q3_0.75, 1, 0), #if distance is <0.75*Q3, then count)
         Q3_1.25 = Q3 * 1.25, 
         Count_Q3_1.25 = ifelse(Distance > Q3_1.25, 1, 0)) %>% #if distance is >1.25*Q3, then count) 
  arrange(Sample_ID, Circulus) %>% # Sort (order) the data by Sample_ID then Circulus 
  group_by(Sample_ID) %>%
  summarise(Q71 = sum(Count_Q3_0.75),
            Q72 = sum(Count_Q3_1.25)) %>%
  left_join(j_JTP, by = c("Sample_ID" = "Sample_ID")) %>% 
  dplyr::select("Sample_ID", "Zone1":"Q69", "Q71":"Q72", "Circuli_SFAZ_0.5":"Check2", -"Sum_Zones")

#subset datasets for different analyses
coho_scales_aukelake <- coho_scales_fulldata %>% 
  dplyr::select("Sample_ID":"Length", "Q32_sum":"Q44") %>%
  filter(Location == "AL") #lake system (Auke Lake)
  
coho_scales_berners <- coho_scales_fulldata %>% 
  dplyr::select("Sample_ID":"Length", "Q32_sum":"Q44") %>%
  filter(Location == "BR") #river system (Berner's River)

coho_scales_hughsmith <- coho_scales_fulldata %>% 
  dplyr::select("Sample_ID":"Length", "Q32_sum":"Q44") %>%
  filter(Location == "HS") #Lake system (Hugh Smith Lake)

coho_scales_lakes <- coho_scales_fulldata %>% 
  dplyr::select("Sample_ID":"Length", "Q32_sum":"Q44") %>%
  filter(Location == "HS" | Location == "AL") # both lake systems only

write.csv(write.csv(coho_scales_fulldata, "data/check.csv")) #outputs data for manual cehgck against original data




