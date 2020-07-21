coho_scales_fulldata <- coho_scales_fulldata %>% filter(Location != "AL")
coho_scales_fulldata <- coho_scales_fulldata %>% filter(Location == "HS")

binomfit <- glm((Age-1)~Q9abs+Q2, data=coho_scales_fulldata, family="binomial")
newdata <- expand.grid(Q9abs=seq(min(coho_scales_bothriv$Q9abs), max(coho_scales_bothriv$Q9abs), length.out = 100), 
                             Q2=seq(min(coho_scales_bothriv$Q2), max(coho_scales_bothriv$Q2), length.out = 100))
                 
                 
newdata$prob = predict(binomfit, newdata, type="response")  

ggplot(newdata, aes(Q9abs, Q2, fill=prob)) +
  geom_tile() +
  scale_fill_gradient2(low="red",mid="yellow",high="blue",midpoint=0.5,
                       limits=c(0,1)) +
  geom_contour(aes(x=Q9abs, y=Q2, z=prob), color="black", breaks = c(0.005, 0.1, 0.5, 0.9, 0.995))

#https://stackoverflow.com/questions/14423325/confidence-intervals-for-predictions-from-logistic-regression

#https://cran.r-project.org/web/packages/sjPlot/vignettes/plot_marginal_effects.html
library(sjPlot)
library(ggplot2)
library(margins)
theme_set(theme_sjplot())

plot_model(binomfit, type = "pred", terms = "Q9abs")
plot_model(binomfit, type = "pred", terms = "Q2")
cplot(binomfit, "Q9abs")
cplot(binomfit, "Q2")

#https://socviz.co/modeling.html
coho_scales_fulldata <- coho_scales_fulldata %>% filter(Location != "AL")
coho_scales_fulldata <- coho_scales_fulldata %>% filter(Location == "HS")
binomfit <- glm((Age-1)~Q9abs+Q2, data=coho_scales_fulldata, family="binomial")
#lm(formula = lifeExp ~ gdpPercap + pop + continent, data = gapminder)
min_Q9abs <- min(coho_scales_fulldata$Q9abs)
max_Q9abs <- max(coho_scales_fulldata$Q9abs)
med_Q2 <- median(coho_scales_fulldata$Q2)

pred_df <- expand.grid(Q9abs = (seq(from = min_Q9abs,
                                        to = max_Q9abs,
                                        length.out = 100)),
                       Q2 = med_Q2)

pred_out <- predict(object = binomfit,
                    newdata = pred_df,
                    interval = "predict")
head(pred_out)
pred_df <- cbind(pred_df, pred_out)
p <- ggplot(data = subset(pred_df, Location %in% c("HS", "BR")),
            aes(x = Q9abs,
                y = fit, ymin = lwr, ymax = upr,
                color = Location,
                fill = Location,
                group = Location))

p + geom_point(data = subset(coho_scales_bothriv,
                             Location %in% c("HS", "BR")),
               aes(x = Q9abs, y = Age-1,
                   color = Location),
               alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_line() +
  geom_ribbon(alpha = 0.2, color = FALSE) +
  scale_x_log10(labels = scales::dollar)


predicted.data <- expand.grid(Q9abs=seq(min(coho_scales_bothriv$Q9abs), max(coho_scales_bothriv$Q9abs), length.out = 20), 
                              Q2=median(coho_scales_bothriv$Q2))
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

visreg(binomfit)