#########################################################################################
# SARS-CoV-2 household analyses for the Netherlands (FFX study). Code is copyrighted by #
# Michiel van Boven and licensed with BSD3 clause (use but mention). Created: 9/9/2020  #
# Here, the household data are analysed using 1) generalised estimating equations of    #
# the secondary attack rate (SAR) using various logistic regressions, and 2) with       #
# a Bayesian final size analysis of a transmission model using Stan (v2.21).            #
# Original code kindly provided by Chris van Dorp. Stan analyses have been validated    #
# using ML analyses developed earlier (de Greeff et al (2012) Epidemiology 23:852-60).  #
#########################################################################################

# Load packages
library(readxl)
library(geepack)
library(emmeans)
library(tidyverse)
#library(broom) 
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(ggplot2)

# Load data 
ffx.data <- read_excel(file.choose())
view(ffx.data)

### 1. CLASSICAL ESTIMATION OF SAR (e.g., Longini et al (2003) AJE, see paper) ###

# Prepare data for SAR estimation using GEE
ffx.data$age_group <- cut(ffx.data$Leeftijd_d1, breaks = c(0,12,18,100), include.lowest = TRUE, right = FALSE, labels = c('young', 'adolescent', 'adult'))
ffx.data.gee <- ffx.data %>%
  mutate(household = Gezin) %>%
  mutate(sex = Geslacht) %>%
  mutate(index = Index) %>%
  mutate(primary = Primair) %>%
  mutate(infected = Conclusie2) %>%
  mutate(health_worker = Gezondheidszorgmedewerker) %>%
  mutate(ECDCcase = CaseECDC_ooit) %>%
  mutate(severity_class = case_when(Severity_ooit == 0 ~ 'mild', # 'asymptomatic',
                                    Severity_ooit == 1 ~ 'mild',
                                    Severity_ooit == 2 ~ 'severe', # 'moderate'
                                    Severity_ooit == 3 ~ 'severe')) %>%
  group_by(household) %>% mutate(household_size = n()) %>%
  select(household, sex, index, primary, age_group, severity_class, health_worker, household_size, ECDCcase, infected) %>%
  mutate(severity_class = replace_na(severity_class, 'not applicable')) %>%  # depends on assumptions which one you choose
  #mutate(severity_class = replace_na(severity_class, 'mild')) %>% # NA is mild
  mutate(health_worker = replace_na(health_worker, -1)) %>%
  mutate(severity_index = severity_class[index == 1]) %>% 
  mutate(severity_primary = severity_class[primary == 1])
view(ffx.data.gee)

# We need to exclude the primary and/or index case for estimation of the SAR. Removing the index case is troublesome when it is not the 
# primary case. We can either 1) base the analyses on the index case (this is usually done in household studies), 2) remove the index and 
# primary case which may introduce bias as you will now remove 2 infections in households where the index case is not the primary case, or 
# 3) think of some clever analysis conditioning on infected index cases that are not the primary case (outside the current scope). 
# I here choose 1) as the GEE analyses are problematic anyway because they do not take the more intricate dependencies in the household 
# into account. To account for this, I use the transmission model.
dim(ffx.data.gee)
ffx.data.gee <- ffx.data.gee %>% 
  subset(index == 0)
  #  %>% subset(primary == 0) # only include this if you would try to estimate SARs from primary case(s)
dim(ffx.data.gee)

# define the reference level for age
ffx.data.gee <- within(ffx.data.gee, age_group <- relevel(age_group, ref = 'adult'))

# Run GEE using geepack

# Based on our own case definition - a selection of models 
# Choose what you need
regression_equation <- formula(infected ~ 1)
regression_equation <- formula(infected ~ age_group) # best model judged by QICC
regression_equation <- formula(infected ~ age_group + severity_primary) 
regression_equation <- formula(infected ~ sex) 
regression_equation <- formula(infected ~ household_size) 
regression_equation <- formula(infected ~ health_worker) 
regression_equation <- formula(infected ~ education) 
regression_equation <- formula(infected ~ severity_index) 
regression_equation <- formula(infected ~ severity_class) 
regression_equation <- formula(infected ~ age_group + sex) # this used to be the best model
regression_equation <- formula(infected ~ age_group + sex + health_worker)  
regression_equation <- formula(infected ~ sex + age_group) 
regression_equation <- formula(infected ~ age_group + sex + household_size) 
regression_equation <- formula(infected ~ age_group + severity_class) 
regression_equation <- formula(infected ~ sex + severity_class) 
regression_equation <- formula(infected ~ age_group + sex + severity_class) 
regression_equation <- formula(infected ~ age_group + sex + severity_index) 
regression_equation <- formula(infected ~ age_group + sex + severity_index + household_size + education + health_worker) 

# based on ECDC case definition - a selection of models
regression_equation <- formula(ECDCcase ~ 1) 
regression_equation <- formula(ECDCcase ~ age_group) 
regression_equation <- formula(ECDCcase ~ sex) 
regression_equation <- formula(ECDCcase ~ household_size) 
regression_equation <- formula(ECDCcase ~ age_group + sex) 
regression_equation <- formula(ECDCcase ~ sex + age_group) 
regression_equation <- formula(ECDCcase ~ age_group + sex + household_size) 
regression_equation <- formula(ECDCcase ~ age_group + sex + severity_index) 

# run the model
fit.gee <- geeglm(regression_equation, id = household, data = ffx.data.gee, family = binomial(link = logit), corstr = "exchangeable")

# variable and model selection
summary(fit.gee)
anova(fit.gee)
geepack::QIC(fit.gee)

# Estimates and 95% confidence intervals for linear predictors 
# Choose what you need
# CHECK: vcov.method = c("vbeta", "vbeta.naiv", "vbeta.j1s", "vbeta.fij", "robust", "naive") # ok
pred <- as_tibble(emmeans(fit.gee, ~ (age_group + sex + severity_index))) 
pred <- as_tibble(emmeans(fit.gee, ~ (age_group + sex))) 
pred <- as_tibble(emmeans(fit.gee, ~ (age_group)))
pred <- as_tibble(emmeans(fit.gee, ~ (1)))
dim(pred)
view(pred)

# estimates and 95% confidence intervals for adjusted marginal secondary attack rates
as.numeric(pred[1,c('emmean', 'asymp.LCL', 'asymp.UCL')]) %>% plogis # group 1 - reference; notice that this can/should be made more abstract
as.numeric(pred[2,c('emmean', 'asymp.LCL', 'asymp.UCL')]) %>% plogis
as.numeric(pred[3,c('emmean', 'asymp.LCL', 'asymp.UCL')]) %>% plogis
as.numeric(pred[4,c('emmean', 'asymp.LCL', 'asymp.UCL')]) %>% plogis
as.numeric(pred[5,c('emmean', 'asymp.LCL', 'asymp.UCL')]) %>% plogis
as.numeric(pred[6,c('emmean', 'asymp.LCL', 'asymp.UCL')]) %>% plogis

### 2. ESTIMATION OF INFECTION RATES USING FINAL SIZE OF SEIR TRANSMISSION MODEL ###

# Prepare data for susceptibility/infectivity estimation using final size of transmission model
ffx.data.fs <- ffx.data %>%
  mutate(household = Gezin) %>%
  mutate(sex = Geslacht) %>%
  mutate(index = Index) %>%
  mutate(primary = Primair) %>%
  mutate(infected = Conclusie2) %>%
  #mutate(severity = minCt) %>%
  select(household, sex, index, primary, age_group, infected)
view(ffx.data.fs)

# A quick view of the data
ffx.data.fs %>%
group_by(household, age_group, primary, index, infected) %>%
tally() %>%
view()  

# Reshape the data for final size analysis using the well-known j, a, n structure for each household (cf stan file)
# where the vector j contains the household infections,n the vector of non-primary infected persons, 
# and a the vector of primary infections. In addition, the vector c contains the non-primary index cases.
# Also, we introduce four types, you never know. Protocode for any number of types is available on request
household.data <- ffx.data.fs %>%
  group_by(household) %>%
  mutate(ones = 1) %>% 
  mutate(conditioning = if_else(index == 1 & primary ==0, 1, 0)) %>% 
  summarise(j1 = sum(infected[age_group == 'young']) - sum(primary[age_group == 'young']),
            j2 = sum(infected[age_group == 'adolescent']) - sum(primary[age_group == 'adolescent']),
            j3 = sum(infected[age_group == 'adult']) - sum(primary[age_group == 'adult']),
            j4 = as.numeric(0),
            a1 = sum(primary[age_group == 'young']),
            a2 = sum(primary[age_group == 'adolescent']),
            a3 = sum(primary[age_group == 'adult']),
            a4 = 0,
            n1 = sum(ones[age_group == 'young']) - sum(primary[age_group == 'young']),
            n2 = sum(ones[age_group == 'adolescent']) - sum(primary[age_group == 'adolescent']),
            n3 = sum(ones[age_group == 'adult']) - sum(primary[age_group == 'adult']),
            n4 = 0,
            c1 = sum(conditioning[age_group == 'young']),
            c2 = sum(conditioning[age_group == 'adolescent']),
            c3 = sum(conditioning[age_group == 'adult']),
            c4 = 0) %>% 
  mutate(conditioning = case_when(c1+c2+c3+c4==0 ~ as.numeric(0), # works bc there is always at most 1 index that is not a primary case
                                  c1 == 1 ~ as.numeric(1),
                                  c2 == 1 ~ as.numeric(2),
                                  c3 == 1 ~ as.numeric(3),
                                  c4 == 1 ~ as.numeric(4))
         )
view(household.data)

# Write data
write_csv(as.data.frame(yourdata), "data/household_data.csv")

# Prepare data for Stan
household.data.list <- list(
  J = as.matrix(household.data[,c("j1", "j2", "j3", "j4")]),
  A = as.matrix(household.data[,c("a1", "a2", "a3", "a4")]),
  N = as.matrix(household.data[,c("n1", "n2", "n3", "n4")]),
  C = as.matrix(household.data[,c("c1", "c2", "c3", "c4")]), 
  H = nrow(household.data),
  conditioning = as.vector(pull(household.data[,c("conditioning")])),
  susceptibility = c(1.0, 1.0),
  #susceptibility_children = 1,
  infectivity = c(1.0, 1.0),
  external_escape = c(0.9999, 0.9999, 0.9999, 0.9999),
  mode = 0
) 

# Initial values for Stan 
# Different scenarios considered in the paper can be evaluated by modifying the stan file and code below
initials = function(){
  return(list(
    #susceptibility_children = 1,
    #susceptibility = c(1, 1),
    #infectivity = c(1, 1),
    #bb = c(0.99, 0.99, 0.99, 0.99),
    #ext_escape = c(0.99),
    #external_escape = c(0.99, 0.99, 0.99, 0.99),
    beta = 1
  ) 
  )
}

# Run Stan model
# On our servers compiled code runs in under a minute
fit <- stan("scripts/model_households.stan", 
            data=household.data.list, 
            init = initials,
            chains = 10,
            iter = 3000, 
            warmup = 1000,  
            thin = 10, 
            refresh = 100, 
            control = list(adapt_delta = 0.85, max_treedepth = 12))
            
# Print fit
print(fit, digits = 4)

# Traceplots 
traceplot(fit)

# Pairwise plots of parameters
pairs(fit, pars=c("beta", "susceptibility", "external_escape"))

# LOO_IC and some checks
loo_output_unstructured = loo(fit, cores = 10, is_method = "psis")
loo_output_unstructured

loo_output_transmissibility = loo(fit, cores = 10, is_method = "psis")
loo_output_transmissibility

loo_output_susceptibility = loo(fit, cores = 10, is_method = "psis")
loo_output_susceptibility

loo_output_fullmodel = loo(fit, cores = 10, is_method = "psis")
loo_output_fullmodel

loo_output_susceptibility = loo(fit, cores = 10, is_method = "psis")
loo_output_susceptibility

loo_output_suschildren = loo(fit, cores = 10, is_method = "psis")
loo_output_suschildren

loo_output_suschildren_ext = loo(fit, cores = 10, is_method = "psis")
loo_output_suschildren_ext

loo_output_suschildren_ext_flat = loo(fit, cores = 10, is_method = "psis")
loo_output_suschildren_ext_flat

loo_compare(loo_output_unstructured, loo_output_transmissibility, loo_output_susceptibility, loo_output_suschildren, loo_output_suschildren_ext, loo_output_suschildren_ext_flat)

# Compare models with loo_compare - gives better indication of which model(s) perform best
plot(loo_output$pointwise[, "elpd_loo"])
plot(loo_output$pointwise[, "mcse_elpd_loo"])
plot(loo_output$pointwise[,])
plot(loo_output)
pareto_k_table(loo_output)
pareto_k_ids(loo_output, threshold = 0.7)
pareto_k_values(loo_output)
psis_n_eff_values(loo_output)
mcse_loo(loo_output, threshold = 0.5)

# Extract parameters
params <- extract(fit)

### include scenario with female versus male adults ###
### insert this in the appropriate place above      ###
view(ffx.data.fs)
ffx.data.fs3 <- ffx.data.fs %>%
  group_by(household) %>%
  mutate(ones = 1) %>% 
  mutate(conditioning = if_else(index == 1 & primary ==0, 1, 0)) %>% 
  summarise(j1 = sum(infected[age_group == 'young']) - sum(primary[age_group == 'young']),
            j2 = sum(infected[age_group == 'adolescent']) - sum(primary[age_group == 'adolescent']),
            j3 = sum(infected[age_group == 'adult' & sex == "M"]) - sum(primary[age_group == 'adult' & sex == "M"]),
            j4 = sum(infected[age_group == 'adult' & sex == "V"]) - sum(primary[age_group == 'adult' & sex == "V"]),
            a1 = sum(primary[age_group == 'young']),
            a2 = sum(primary[age_group == 'adolescent']),
            a3 = sum(primary[age_group == 'adult' & sex == "M"]),
            a4 = sum(primary[age_group == 'adult' & sex == "V"]),
            n1 = sum(ones[age_group == 'young']) - sum(primary[age_group == 'young']),
            n2 = sum(ones[age_group == 'adolescent']) - sum(primary[age_group == 'adolescent']),
            n3 = sum(ones[age_group == 'adult' & sex == "M"]) - sum(primary[age_group == 'adult' & sex == "M"]),
            n4 = sum(ones[age_group == 'adult' & sex == "V"]) - sum(primary[age_group == 'adult' & sex == "V"]),
            c1 = sum(conditioning[age_group == 'young']),
            c2 = sum(conditioning[age_group == 'adolescent']),
            c3 = sum(conditioning[age_group == 'adult' & sex == "M"]),
            c4 = sum(conditioning[age_group == 'adult' & sex == "V"])
            ) %>% 
  mutate(conditioning = case_when(c1+c2+c3+c4==0 ~ as.numeric(0), # works bc there is always at most 1 index that is not a primary case
                                  c1 == 1 ~ as.numeric(1),
                                  c2 == 1 ~ as.numeric(2),
                                  c3 == 1 ~ as.numeric(3),
                                  c4 == 1 ~ as.numeric(4))
  )
view(ffx.data.fs3)

### include scenario with female versus male adults/adolescents ###
### and female versus male children ###
### insert this in the appropriate place above      ###
view(ffx.data.fs)
ffx.data.fs3 <- ffx.data.fs %>%
  group_by(household) %>%
  mutate(ones = 1) %>% 
  mutate(conditioning = if_else(index == 1 & primary == 0, 1, 0)) %>% 
  summarise(j1 = sum(infected[age_group == 'young' & sex == "M"]) - sum(primary[age_group == 'young' & sex == "M"]),
            j2 = sum(infected[age_group == 'young' & sex == "F"]) - sum(primary[age_group == 'young' & sex == "F"]),
            
            
            j3 = sum(infected[age_group == 'adult' & sex == "M"]) - sum(primary[age_group == 'adult' & sex == "M"]),
            j4 = sum(infected[age_group == 'adult' & sex == "V"]) - sum(primary[age_group == 'adult' & sex == "V"]),
            a1 = sum(primary[age_group == 'young']),
            a2 = sum(primary[age_group == 'adolescent']),
            a3 = sum(primary[age_group == 'adult' & sex == "M"]),
            a4 = sum(primary[age_group == 'adult' & sex == "V"]),
            n1 = sum(ones[age_group == 'young']) - sum(primary[age_group == 'young']),
            n2 = sum(ones[age_group == 'adolescent']) - sum(primary[age_group == 'adolescent']),
            n3 = sum(ones[age_group == 'adult' & sex == "M"]) - sum(primary[age_group == 'adult' & sex == "M"]),
            n4 = sum(ones[age_group == 'adult' & sex == "V"]) - sum(primary[age_group == 'adult' & sex == "V"]),
            c1 = sum(conditioning[age_group == 'young']),
            c2 = sum(conditioning[age_group == 'adolescent']),
            c3 = sum(conditioning[age_group == 'adult' & sex == "M"]),
            c4 = sum(conditioning[age_group == 'adult' & sex == "V"])
  ) %>% 
  mutate(conditioning = case_when(c1+c2+c3+c4==0 ~ as.numeric(0), # works bc there is always at most 1 index that is not a primary case
                                  c1 == 1 ~ as.numeric(1),
                                  c2 == 1 ~ as.numeric(2),
                                  c3 == 1 ~ as.numeric(3),
                                  c4 == 1 ~ as.numeric(4))
  )
view(ffx.data.fs3)
