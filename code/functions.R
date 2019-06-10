f_sum <- function(data, col1, col2, div){
  col1 = enquo(col1)
  col2 = enquo(col2)
  div = enquo(div)
  
  data %>% 
    dplyr::select(!!col1:!!col2) %>% 
    apply(1, sum, na.rm=T) %>% 
    as.data.frame %>% 
    bind_cols(data) %>% 
    transmute(temp = . / !!div) %>% 
    .$temp
} # Sum function from Ben Williams

convyear <- function(x, year=2000){ # This function converts year correctly
  m <- year(x) %% 100
  year(x) <- ifelse(m > year %% 100, 1900+m, 2000+m)
  x
}


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

eda.norm <- function(x, ...){ #Function from Franz Mueter
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
