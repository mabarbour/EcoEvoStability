
## Setup ----

# load posterior predictions (and bayesian MAR1 model of eco-evo dynamics)
load("analyses/posteriors_ecoevo_model.RData")

# load and view data
source('code/manage_data.R')

# copy-pasted from fit_ecoevo_model_MeasurementError.R
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
  # testing, transform NA to marginal values to facilitate model fitting
  mutate(se_Rt = ifelse(is.na(se_Rt == T), 0.001, se_Rt),
         se_Yt = ifelse(is.na(se_Yt == T), 0.001, se_Yt)) %>%
  filter(week < 12) # to match experimental data

# load other required libraries
library(cowplot)
library(marginaleffects)

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Plot of eco-evolutionary dynamics in aphid and aphid-parasitoid systems (Fig. 3) ----

# get expected value of aphid population dynamics in polycultures
epred_lnAt1 <- predictions(ecoevo_dynamics_brm, 
                           by = c("week","Treatment"),
                           newdata = filter(ecoevo_df_allWeeks_MeasurementError, Aphids_sub == 1), 
                           resp = "lnAt1", 
                           re_formula = NULL, 
                           type = "response") %>%
  mutate(resp = "lnAt1")

# get expected value of parasitoid population dynamics in polycultures
epred_lnPt1 <- predictions(ecoevo_dynamics_brm, 
                           by = c("week","Treatment"),
                           newdata = filter(ecoevo_df_allWeeks_MeasurementError, Ptoids_sub == 1, is.na(Rt) != T), 
                           resp = "lnPt1", 
                           re_formula = NULL, 
                           type = "response") %>%
  mutate(resp = "lnPt1")

# get raw aphid population dynamics
raw_lnAt1 <- filter(ecoevo_df_allWeeks_imputed_withWeek0, Aphids_sub == 1) %>% # df
  filter(week < 12) %>%
  select(week, Treatment, ID, log1pNt1 = lnAt1) %>%
  mutate(resp = "lnAt1")

# get raw parasitoid population dynamics
raw_lnPt1 <- filter(ecoevo_df_allWeeks_imputed_withWeek0) %>% 
  filter(week < 12) %>%
  select(week, Treatment, ID, log1pNt1 = lnPt1) %>% 
  mutate(resp = "lnPt1")

# organize equilibrium aphid abundance data for plotting
eco_equilib_df <- data.frame(Treatment = factor(c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid")),
                             aphid_equilibrium = c(eq_aphid[1], 0), # setting to 0 for aphid-parasitoid because the equilibrium is not feasible
                             week = c(12,12)) # set as last week for plot, but its a long-term prediction
round(expm1(eq_aphid[1]),0) # equilibrium aphid abundance

# make plot of ecological dynamics
eco_dynamics_plot <- bind_rows(epred_lnAt1 %>% filter(Treatment %in% c("RYG","RYGP")), 
                               epred_lnPt1 %>% filter(Treatment %in% c("RYG","RYGP"))) %>%
  #filter(drop_prediction == 0) %>%
  mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))) %>%
  ggplot(aes(x = week, y = estimate)) +
  geom_line(aes(color = resp), linewidth = 1) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = resp), alpha = 0.3) +
  facet_wrap(~Treatment, ncol = 1) +
  geom_line(data = bind_rows(raw_lnAt1, raw_lnPt1) %>% 
              filter(week != 4, Treatment %in% c("RYG","RYGP")) %>% 
              mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))), 
            aes(x = week, y = log1pNt1, color = resp, group = paste(resp,ID)),
            alpha = 0.2) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) + #, breaks = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,"\u221e")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(data = eco_equilib_df, aes(x = week, y = aphid_equilibrium), size = 2, color = cbbPalette[1])

# get expected value for red morph frequency dynamics in polycultures
epred_Rt1 <- predictions(ecoevo_dynamics_brm, 
                         by = c("week","Treatment"),
                         newdata = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1), 
                         resp = "Rt1", 
                         re_formula = NULL, 
                         type = "response") %>%
  mutate(resp = "Rt1")

# get expected value for yellow morph frequency dynamics in polycultures
epred_Yt1 <- predictions(ecoevo_dynamics_brm, 
                         by = c("week","Treatment"),
                         newdata = filter(ecoevo_df_allWeeks_MeasurementError, Freq_sub == 1), 
                         resp = "Yt1", 
                         re_formula = NULL, 
                         type = "response") %>%
  mutate(resp = "Yt1")

# expected mean frequences of red and yellow morph at week 11
filter(epred_Rt1, week == 11, Treatment == "RYGP")$estimate
filter(epred_Yt1, week == 11, Treatment == "RYGP")$estimate

# raw red morph frequency dynamics
raw_Rt1 <- filter(ecoevo_df_allWeeks_imputed_withWeek0, Freq_sub == 1) %>% 
  filter(week < 12) %>%
  select(week, Treatment, ID, Freq_t1 = Rt1) %>% 
  mutate(resp = "Rt1")#, drop_prediction = 0)

# raw yellow morph frequency dynamics
raw_Yt1 <- filter(ecoevo_df_allWeeks_imputed_withWeek0, Freq_sub == 1) %>% 
  filter(week < 12) %>%
  select(week, Treatment, ID, Freq_t1 = Yt1) %>% 
  mutate(resp = "Yt1")#, drop_prediction = 0)

# organize equilibrium frequency data for plot
evo_equilib_df <- data.frame(Treatment = factor(c("RYG","RYG","RYGP","RYGP"), 
                                                levels = c("RYG","RYGP"), 
                                                labels = c("Aphid","Aphid-Parasitoid")),
                             freq_equilibrium = c(eq_aphid[2], eq_aphid[3], NA, NA), # frequencies are not feasible in aphid-parasitoid system
                             resp = c("Rt1","Yt1","Rt1","Yt1"),
                             week = c(12, 12, 12, 12))
round(eq_aphid[2],2) # red morph equilibrium frequency
round(eq_aphid[3],2) # yellow morph equilibrium frequency

# make plot of evolutionary dynamics
evo_dynamics_plot <- bind_rows(epred_Rt1, epred_Yt1) %>%
  filter(Treatment %in% c("RYG","RYGP")) %>% # drop_prediction == 0, 
  mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))) %>%
  ggplot(aes(x = week, y = estimate)) +
  geom_line(aes(color = resp), linewidth = 1) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = resp), alpha = 0.3) +
  facet_wrap(~Treatment, ncol = 1) +
  geom_line(data = bind_rows(raw_Rt1, raw_Yt1) %>% 
              filter(week != 4, Treatment %in% c("RYG","RYGP")) %>% 
              mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))), 
            aes(x = week, y = Freq_t1, color = resp, group = paste(resp,ID)),
            alpha = 0.2) + 
  scale_color_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  scale_x_continuous(name = "Week", breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,"\u221e")) +
  ylab("Frequency (z)") +
  geom_hline(yintercept = c(0,1), linetype = "dotted") +
  geom_point(data = evo_equilib_df, aes(x = week, y = freq_equilibrium, color = resp), size = 2)

# generate figure 3
plot_ecoevo_dynamics <- plot_grid(eco_dynamics_plot, # + ggtitle("Ecological dynamics"), 
                                  evo_dynamics_plot, # + ggtitle("Evolutionary dynamics"), 
                                  ncol = 2,
                                  labels = "AUTO")
plot_ecoevo_dynamics
save_plot(filename = "figures/Fig3_ecoevo_dynamics.pdf",
         plot = plot_ecoevo_dynamics, base_height = 6)


# Plot of monoculture population dynamics ----
bind_rows(epred_lnAt1 %>% filter(Treatment %in% c("RP","YP","GP")), 
          epred_lnPt1 %>% filter(Treatment %in% c("RP","YP","GP"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("RP","YP","GP"), labels = c("Red","Yellow","Green"))) %>%
  ggplot(aes(x = week, y = estimate)) +
  geom_line(aes(color = Treatment), linewidth = 1) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = Treatment), alpha = 0.3) +
  facet_wrap(~resp) + 
  scale_color_manual(name = NULL, values = cbbPalette[c(7,5,4)], labels = c("Red morph","Yellow morph","Green morph")) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(7,5,4)], labels = c("Red morph","Yellow morph","Green morph")) +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) + #, breaks = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(0,2,4,6,8,10,12)) +
  geom_hline(yintercept = 0, linetype = "dotted")

