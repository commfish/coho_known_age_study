library(corrplot)
library(GGally)
library(here)

# Run previous script to import data
source(here::here("code/cohoage_1_DataImport_JTP.R"))





for(q in c("HS", "BR", "AL")) {
  print(corrplot.mixed(cor(coho_scales_fulldata %>%
                     filter(Location == q) %>%
                     dplyr::select("Q32_sum":"Q44"), use="complete.obs"), title=q, mar=c(0,0,1,0) ))
  
  # Now explore how these correlations and "dropoffs" change btwn locations
  .temp4 <- cor(coho_scales_fulldata %>%
                  filter(Location == q) %>% 
                  dplyr::select("Q4":"Q31"), use="complete.obs")
  .temp4[upper.tri(.temp4)] <- NA

  shift <- function(x, n){
    c(x[-(seq(n))], rep(NA, n))
  }
  
  .temp4 <- as.data.frame(.temp4)
  for(i in 2:28){
    .temp4[,i] <- shift(.temp4[,i], i-1)
  }
  
  .temp4["rownum"] <- seq(1:28)
  .temp4 <- .temp4 %>% gather("Qnum", "corr", Q4:Q31)
  
  print(ggplot(.temp4, aes(x=rownum, y = corr, color=Qnum)) +
    geom_line() + labs(title = q))
}

#caution, this takes a while to eval!
ggpairs(coho_scales_aukelake %>% dplyr::select("Zone1":"Count_Zone2", "Year", "Age":"Q10"))
