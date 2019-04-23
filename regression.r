###################
#Regression
###################
library(MASS)
fit <- lm(dat$Blastocystis_reads ~ dat$age_years + dat$filter_water + dat$sex)
step <- stepAIC(fit, direction="both")
step$anova
