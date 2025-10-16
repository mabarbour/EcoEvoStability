## Setup ----

# load and view data
source('code/manage_data.R')
ecoevo_df 

# load other required libraries
library(brms)
library(cowplot)
library(marginaleffects)

# set plot theme
theme_set(theme_cowplot())


## Fit Bayesian multivariate autoregressive model ----

# formula for aphid ecological dynamics
A_bf <- bf(lnAt1|subset(Aphids_sub) ~ offset(lnAt) + 
             0 + Intercept + lnAt + Rt + Yt +
             ptoid_in + ptoid_in:(lnAt + lnPt + Rt + Yt) + 
             (1|i|ID)) 

# formula for parasitoid ecological dynamics
P_bf <- bf(lnPt1|subset(Ptoids_sub) ~ offset(lnPt) + 
             0 + Intercept + lnAt + lnPt + Rt + Yt + 
             (1|i|ID)) 

# formula for change in red morph frequency (aphid evo dynamics)
R_bf <- bf(Rt1|subset(Freq_sub) ~ offset(Rt) +
             0 + Intercept + Rt + Yt + lnAt + 
             ptoid_in + ptoid_in:(Rt + Yt + lnAt + lnPt) +
             (1|i|ID)) 

# formula for change in yellow morph frequency (aphid evo dynamics)
Y_bf <- bf(Yt1|subset(Freq_sub) ~ offset(Yt) +
             0 + Intercept + Yt + Rt + lnAt +
             ptoid_in + ptoid_in:(Yt + Rt + lnAt + lnPt) + 
             (1|i|ID)) 

# visualize priors
hist(rnorm(100, 0.5*7, 3)) # aphid intrinsic growth rates, based on previous experiments
hist(rnorm(100, -2, 2)) # parasitoid intrinsic growth rates, based on Barbour et al. 2022, Science Fig. S5
hist(rnorm(100, 0, 1)) # density and frequency dependence, regularizing toward zero. 
hist(rexp(100, 1)) # random effects

# fit model
ecoevo_dynamics_brm <- brm(formula = mvbf(A_bf, P_bf, R_bf, Y_bf, rescor = F), # rescor = F because of different data subsets
                           data = ecoevo_df,
                           cores = 4,
                           prior = c(
                             # intercepts
                             set_prior("normal(3.5,3)", class = "b", coef = "Intercept", resp = "lnAt1"),
                             set_prior("normal(-2,2)", class = "b", coef = "Intercept", resp = "lnPt1"),
                             set_prior("normal(0,1)", class = "b", coef = "Intercept", resp = "Rt1"),
                             set_prior("normal(0,1)", class = "b", coef = "Intercept", resp = "Yt1"),
                             # coefficients
                             set_prior("normal(0,1)", class = "b", resp = "lnAt1"),
                             set_prior("normal(0,1)", class = "b", resp = "lnPt1"),
                             set_prior("normal(0,1)", class = "b", resp = "Rt1"),
                             set_prior("normal(0,1)", class = "b", resp = "Yt1"),
                             # random effect terms
                             set_prior("exponential(1)", class = "sd", resp = "lnAt1"),
                             set_prior("exponential(1)", class = "sd", resp = "lnPt1"),
                             set_prior("exponential(1)", class = "sd", resp = "Rt1"),
                             set_prior("exponential(1)", class = "sd", resp = "Yt1")),
                           control = list(adapt_delta = 0.95),
                           file = "analyses/ecoevo_dynamics_brm") 

summary(ecoevo_dynamics_brm)
bayes_R2(ecoevo_dynamics_brm)

## Check posterior predictions for individual cages ----

# aphids
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df, Aphids_sub == 1), 
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

# parasitoids
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df, Ptoids_sub == 1), 
            resp = "lnPt1", 
            re_formula = NULL, 
            type = "prediction") %>%
  filter(is.na(estimate) == F) %>% # limit to data where we can make predictions. 
  # this excludes a few data points where aphids were not detected, therefore, we cannot measure morph frequency,
  # which is part of the model for parasitoid abundance.
  ggplot(aes(x = week, y = estimate, group = ID)) +
  geom_line(aes(color = Treatment)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), 
              alpha = 0.2) +
  geom_point(aes(y = lnPt1, color = Treatment)) +
  facet_wrap(~ID) +
  geom_hline(yintercept = 0, linetype = "dotted")


# Red morph frequency
predictions(ecoevo_dynamics_brm, 
            newdata = filter(ecoevo_df, Freq_sub == 1), 
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
            newdata = filter(ecoevo_df, Freq_sub == 1), 
            resp = "Yt1", 
            re_formula = NULL, 
            type = "prediction") %>%
  ggplot(aes(x = week, y = estimate, group = ID)) +
  geom_line(aes(color = Treatment)) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), alpha = 0.2) +
  geom_point(aes(y = Yt1, color = Treatment)) +
  facet_wrap(~ID) +
  geom_hline(yintercept = 0, linetype = "dotted")


