
## Setup ----

# load and view data
source('code/manage_data.R')

ecoevo_df_allWeeks_MeasurementError <- ecoevo_df_allWeeks_imputed %>% 
  # following Barbour et al. 2022, Science supplement average for 2 aphid species and estimate for D. rapae
  # reduce measurement error for weeks where aphids/parasitoids were not present for each variable
  mutate(se_lnAt = ifelse(week == 1, 0.01, 
                          ifelse(Aphids_sub == 1, 0.26, 0.01)), 
         se_lnPt = ifelse(week < 5, 0.01, 
                          ifelse(Ptoids_sub == 1, 0.13, 0.01)),
         # adding measurement error based on standard error for a sample of binomial random variables: 
         # following: https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/what-is-the-standard-error-of-a-sample/
         se_Rt = sqrt(Rt*(1-Rt)/(expm1(lnAt)-1)),
         se_Yt = sqrt(Yt*(1-Yt)/(expm1(lnAt)-1)),
         ind = 1:nrow(.)) %>%
  mutate(se_Rt = ifelse(se_Rt == 0, 0.001, se_Rt),
         se_Yt = ifelse(se_Yt == 0, 0.001, se_Yt)) %>%
  # transform NA to marginal values to facilitate model fitting
  mutate(se_Rt = ifelse(is.na(se_Rt == T), 0.001, se_Rt),
         se_Yt = ifelse(is.na(se_Yt == T), 0.001, se_Yt)) %>%
  filter(week < 12) # to match experimental data

ecoevo_df_allWeeks_MeasurementError # 330 data points

# calculate number of observations with frequencies outside range of 0.2 to 0.8
filter(ecoevo_df_allWeeks_MeasurementError, 
       Treatment %in% c("RYG","RYGP"),
       Rt1 > 0.8 | Yt1 > 0.8 | Gt1 > 0.8,
       Rt1 < 0.2 | Yt1 < 0.2 | Gt1 < 0.2) %>%
  select(week, ID, Treatment, Rt1, Yt1, Gt1)
nrow(filter(ecoevo_df_allWeeks_MeasurementError, Treatment %in% c("RYG","RYGP")))
2/132 # less than 2% (2 of 132) of observations were outside the range of 0.2 to 0.8, suggesting that our linear modeling approach is likely robust

# load other required libraries
library(brms)
library(cowplot)
library(marginaleffects)

# set plot theme
theme_set(theme_cowplot())

## Check for patterns among treatments in plants with minimal growth (covariate) ----

# view pattern in plants with minimal growth over time
ecoevo_df_allWeeks_MeasurementError %>%
  group_by(week) %>%
  summarise(sum_plants_no_growth = sum(plants_no_growth, na.rm = T),
            n_plants_no_growth = n()) %>%
  ungroup() %>%
  mutate(prop_plants_no_growth = sum_plants_no_growth/n_plants_no_growth) %>%
  ggplot(aes(x = week, y = prop_plants_no_growth)) +
  geom_line()
# strong spike at week 6

ecoevo_df_allWeeks_MeasurementError %>%
  group_by(ptoid_in) %>%
  summarise(sum_plants_no_growth = sum(plants_no_growth, na.rm = T),
            n_plants_no_growth = n()) %>%
  mutate(prop = sum_plants_no_growth/n_plants_no_growth)
# no pattern caused by presence/absence of parasitoids in probability of having plants with no growth

# test for differences among groups within just parasitoid
# treatments and also between parasitoid treatments
# and overall.
ptoid_chisq_table <- ecoevo_df_allWeeks_MeasurementError %>%
  filter(week > 3, Aphids_sub == 1) %>% # all cages with aphids (base for data)
  group_by(ptoid_in) %>%
  summarise(Small = sum(plants_no_growth, na.rm = T),
            Normal = n() - Small) %>%
  column_to_rownames(var = "ptoid_in")
ptoid_chisq_table
chisq.test(ptoid_chisq_table)
chisq.test(ptoid_chisq_table)$expected
fisher.test(ptoid_chisq_table)
# no clear evidence for a bias among treatments, particularly
# those with and without aphids

treatment_chisq_table <- ecoevo_df_allWeeks_MeasurementError %>%
  filter(Aphids_sub == 1, Treatment != "RYG") %>% # all cages with aphids (base for data)
  group_by(Treatment) %>%
  summarise(Small = sum(plants_no_growth, na.rm = T),
            Normal = n() - Small) %>%
  column_to_rownames(var = "Treatment")
treatment_chisq_table
chisq.test(treatment_chisq_table)
chisq.test(treatment_chisq_table)$expected
fisher.test(treatment_chisq_table)
# no clear evidence for bias among parasitoid treatments

## Assess collinearity among predictors ----
# ptoid_in == 0 (all below 2)
car::vif(lm(lnAt1 ~ plants_no_growth + lnAt + Rt + Yt,# + lnPt, 
            data = filter(ecoevo_df_allWeeks_MeasurementError, Aphids_sub == 1, ptoid_in == 0)))
car::vif(lm(Rt1 ~ plants_no_growth + lnAt + Rt + Yt,# + lnPt, 
            data = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1, ptoid_in == 0)))
# note that doing a model with Yt1 would be redundant, since collinearity assess predictors and the data is the same as Rt1
# ptoid_in == 1 (all below 1.4)
car::vif(lm(lnAt1 ~ plants_no_growth + lnAt + Rt + Yt + lnPt, 
            data = filter(ecoevo_df_allWeeks_MeasurementError, Aphids_sub == 1, ptoid_in == 1)), type = "predictor")
car::vif(lm(lnPt1 ~ plants_no_growth + lnAt + Rt + Yt + lnPt, 
            data = filter(ecoevo_df_allWeeks_MeasurementError, Ptoids_sub == 1, ptoid_in == 1)), type = "predictor")
car::vif(lm(Rt1 ~ plants_no_growth + lnAt + Rt + Yt + lnPt, 
            data = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1, ptoid_in == 1)), type = "predictor")
# note: Rt1 or Yt1 give same results since the same data is used

## Fit Bayesian multivariate autoregressive model ----

# multivariate formula
MAR1_model <- 
  # aphid ecological dynamics
  bf(lnAt1|subset(Aphids_sub) ~ 
       0 + Intercept + mi(lnAt, idx=ind) + mi(Rt, idx=ind) + mi(Yt, idx=ind) +
       ptoid_in + ptoid_in:(mi(lnAt, idx=ind) + mi(lnPt, idx=ind) + mi(Rt, idx=ind) + mi(Yt, idx=ind)) + 
       plants_no_growth + 
       (1|i|ID)) +
  # parasitoid ecological dynamics
  bf(lnPt1|subset(Ptoids_sub) ~   
       0 + Intercept + mi(lnAt, idx=ind) + mi(lnPt, idx=ind) + mi(Rt, idx=ind) + mi(Yt, idx=ind) + 
       plants_no_growth + 
       (1|i|ID)) +
  # change in red morph frequency (aphid evo dynamics)
  bf(Rt1|subset(Freq_sub) ~   
       0 + Intercept + mi(Rt, idx=ind) + mi(Yt, idx=ind) + mi(lnAt, idx=ind) +  
       ptoid_in + ptoid_in:(mi(Rt, idx=ind) + mi(Yt, idx=ind) + mi(lnAt, idx=ind) + mi(lnPt, idx=ind)) +
       plants_no_growth + 
       (1|i|ID)) +
  # change in yellow morph frequency (aphid evo dynamics)
  bf(Yt1|subset(Freq_sub) ~  
       0 + Intercept + mi(Yt, idx=ind) + mi(Rt, idx=ind) + mi(lnAt, idx=ind) +  
       ptoid_in + ptoid_in:(mi(Yt, idx=ind) + mi(Rt, idx=ind) + mi(lnAt, idx=ind) + mi(lnPt, idx=ind)) +
       plants_no_growth + 
       (1|i|ID)) +
  # measurement error
  bf(lnAt|mi(se_lnAt) + index(ind) ~ 1 + (1|i|ID)) +  
  bf(lnPt|mi(se_lnPt) + index(ind) ~ 1 + (1|i|ID)) +
  bf(Rt|mi(se_Rt) + index(ind) ~ 1 + (1|i|ID)) +
  bf(Yt|mi(se_Yt) + index(ind) ~ 1 + (1|i|ID)) +
  set_rescor(rescor = FALSE)

# visualize priors
get_prior(formula = MAR1_model, data = ecoevo_df_allWeeks_MeasurementError)
hist(rnorm(100, 0.5*7, 3)) # aphid intrinsic growth rates, based on previous experiments
hist(rnorm(100, -2, 2)) # parasitoid intrinsic growth rates, based on Barbour et al. 2022, Science Fig. S5
hist(rnorm(100, 0, 1)) # density and frequency dependence, regularizing toward zero. 
hist(rexp(100, 1)) # random effects

# fit model
ecoevo_dynamics_brm <- brm(formula = MAR1_model, # rescor = F because of different data subsets
                           data = ecoevo_df_allWeeks_MeasurementError,
                           cores = 4,
                           iter = 4000, # increased bulk effective sample size
                           prior = c(
                             # intercepts
                             set_prior("normal(3.5,3)", class = "b", coef = "Intercept", resp = "lnAt1"),
                             set_prior("normal(-2,2)", class = "b", coef = "Intercept", resp = "lnPt1"),
                             set_prior("normal(0,1)", class = "b", coef = "Intercept", resp = "Rt1"),
                             set_prior("normal(0,1)", class = "b", coef = "Intercept", resp = "Yt1"),
                             # intra-specific and intra-trait coefficients
                             set_prior("normal(1,1)", class = "b", coef = "milnAtidxEQind", resp = "lnAt1"),
                             set_prior("normal(1,1)", class = "b", coef = "milnPtidxEQind", resp = "lnPt1"),
                             set_prior("normal(1,1)", class = "b", coef = "miRtidxEQind", resp = "Rt1"),
                             set_prior("normal(1,1)", class = "b", coef = "miYtidxEQind", resp = "Yt1"),
                             # additional coefficients
                             set_prior("normal(0,1)", class = "b", resp = "lnAt1"),
                             set_prior("normal(0,1)", class = "b", resp = "lnPt1"),
                             set_prior("normal(0,1)", class = "b", resp = "Rt1"),
                             set_prior("normal(0,1)", class = "b", resp = "Yt1"),
                             # random effect terms
                             set_prior("exponential(1)", class = "sd", resp = "lnAt1"),
                             set_prior("exponential(1)", class = "sd", resp = "lnPt1"),
                             set_prior("exponential(1)", class = "sd", resp = "Rt1"),
                             set_prior("exponential(1)", class = "sd", resp = "Yt1")),
                           control = list(adapt_delta = 0.95, max_treedepth = 20),
                           file = "analyses/ecoevo_dynamics_brm") 

summary(ecoevo_dynamics_brm)

# calculate Bayesian R2
bayes_R2(ecoevo_dynamics_brm, resp = "lnAt1", newdata = filter(ecoevo_df_allWeeks_MeasurementError, Aphids_sub == 1))
bayes_R2(ecoevo_dynamics_brm, resp = "lnPt1", newdata = filter(ecoevo_df_allWeeks_MeasurementError, Ptoids_sub == 1, is.na(Rt) != T))  
bayes_R2(ecoevo_dynamics_brm, resp = "Rt1", newdata = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1)) 
bayes_R2(ecoevo_dynamics_brm, resp = "Yt1", newdata = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1)) 


## Check posterior predictions for individual cages ----

# aphids
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df_allWeeks_MeasurementError, Aphids_sub == 1), 
            resp = "lnAt1", 
            re_formula = NULL, 
            type = "prediction") %>%
  ggplot(aes(x = week, y = estimate, group = ID)) +
  geom_line(aes(color = Treatment)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), 
              alpha = 0.2) +
  geom_point(aes(y = lnAt1, color = Treatment)) +
  facet_wrap(~ID) +
  geom_hline(yintercept = 0, linetype = "dotted")
# difficulty predicting final week of extinction in some cages

# parasitoids
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df_allWeeks_MeasurementError, Ptoids_sub == 1, is.na(Rt) != T), 
            resp = "lnPt1", 
            re_formula = NULL, 
            type = "prediction") %>%
  ggplot(aes(x = week, y = estimate, group = ID)) +
  geom_line(aes(color = Treatment)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), 
              alpha = 0.2) +
  geom_point(aes(y = lnPt1, color = Treatment)) +
  facet_wrap(~ID) +
  geom_hline(yintercept = 0, linetype = "dotted")
# difficulty predicting weeks 6 and 7 (trends seem incorrect)

# Red morph frequency
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1), 
            resp = "Rt1", 
            re_formula = NULL, 
            type = "prediction") %>%
  ggplot(aes(x = week, y = estimate, group = ID)) +
  geom_line(aes(color = Treatment)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), alpha = 0.2) +
  geom_point(aes(y = Rt1, color = Treatment)) +
  facet_wrap(~ID) +
  geom_hline(yintercept = 0, linetype = "dotted")


# Yellow morph frequency
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1), 
            resp = "Yt1", 
            re_formula = NULL, 
            type = "prediction") %>%
  ggplot(aes(x = week, y = estimate, group = ID)) +
  geom_line(aes(color = Treatment)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), alpha = 0.2) +
  geom_point(aes(y = Yt1, color = Treatment)) +
  facet_wrap(~ID) +
  geom_hline(yintercept = 0, linetype = "dotted")

