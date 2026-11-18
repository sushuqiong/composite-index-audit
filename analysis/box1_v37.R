library(CompositeIndexAudit)
set.seed(20261023); n <- 1500
d <- data.frame(sex=factor(sample(c('Male','Female'),n,TRUE),levels=c('Male','Female')),
  age=runif(n,20,80), BMI=rnorm(n,27,5), waist_cm=rnorm(n,95,14),
  tg_mmol_L=rlnorm(n,0,0.5), glucose_mg_dL=rnorm(n,100,22))
d$metab <- log10(d$glucose_mg_dL) + log10(d$tg_mmol_L)
d$log_age <- log10(d$age); d$log_bmi <- log10(d$BMI)
d$GLM7 <- d$metab + d$log_age + d$log_bmi
register_index('DemoGLM7','metab + log10(age) + log10(BMI)',
  units=c(metab='-',age='years',BMI='kg/m^2'),overwrite=TRUE)
keys_min <- c(log_age='age|log10|10|years|-',log_bmi='BMI|log10|10|kg/m^2|-')
a <- audit_registered_index('DemoGLM7',keys_min,data=d)
keys_full <- c(metab='metab|identity|-|-|-',keys_min)
b <- audit_registered_index('DemoGLM7',keys_full,data=d)
stopifnot(a$classification=='rank_adding_exact',b$classification=='rank_preserving_exact')
print(c(minimal_base=a$classification,embedded_base=b$classification))
d$LAP <- ifelse(d$sex=='Male',d$waist_cm-65,d$waist_cm-58)*d$tg_mmol_L
# Exploratory only: built-in expansion and bare names; not exact certification.
e <- audit_registered_index('LAP',c('sex','waist_cm','tg_mmol_L'),data=d,strict=FALSE)
stopifnot(e$classification=='reducible_to_interaction')
print(c(LAP_exploratory=e$classification))
d$dx <- rbinom(n,1,plogis(-3+0.3*d$metab+0.2*d$log_age))
fit <- glm(dx ~ metab + log_age + log_bmi + GLM7,data=d,family=binomial())
stopifnot(is.na(coef(fit)['GLM7']))
print(coef(fit)['GLM7'])
cat('BOX1_V37_COMPLETE\n')
