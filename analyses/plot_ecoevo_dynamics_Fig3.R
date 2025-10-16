## Setup ----

# load and view data
source('code/manage_data.R')
ecoevo_df 

# load other required libraries
library(brms)
library(cowplot)
library(marginaleffects)

# load posterior predictions (and bayesian MAR1 model of eco-evo dynamics)
load("analyses/posteriors_ecoevo_model.RData")

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Plot of eco-evolutionary dynamics in aphid and aphid-parasitoid systems (Fig. 3) ----

# get expected value of aphid population dynamics in polycultures
epred_lnAt1 <- predictions(ecoevo_dynamics_brm, 
                           by = c("week","Treatment"),
                           newdata = filter(ecoevo_df, Aphids_sub == 1), 
                           resp = "lnAt1", 
                           re_formula = NULL, 
                           type = "response") %>%
  mutate(resp = "lnAt1")

# get expected value of parasitoid population dynamics in polycultures
epred_lnPt1 <- predictions(ecoevo_dynamics_brm, 
                           by = c("week","Treatment"),
                           newdata = filter(ecoevo_df, Ptoids_sub == 1), 
                           resp = "lnPt1", 
                           re_formula = NULL, 
                           type = "response") %>%
  mutate(resp = "lnPt1")

# get raw aphid population dynamics
raw_lnAt1 <- filter(df, Treatment %in% c("RYG","RYGP")) %>% 
  select(week, Treatment, ID, Aphids_cage) %>%
  mutate(resp = "lnAt1", log1pNt1 = log1p(Aphids_cage)) %>%
  select(-Aphids_cage)

# get raw parasitoid population dynamics
raw_lnPt1 <- filter(df, Treatment %in% c("RYGP")) %>% 
  select(week, Treatment, ID, Ptoids_cage) %>% 
  mutate(resp = "lnPt1", log1pNt1 = log1p(Ptoids_cage)) %>%
  select(-Ptoids_cage)

# organize equilibrium aphid abundance data for plotting
eco_equilib_df <- data.frame(Treatment = factor(c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid")),
                             aphid_equilibrium = c(eq_aphid[1], 0), # setting to 0 for aphid-parasitoid because the equilibrium is not feasible
                             week = c(14,14)) # set as last week for plot, but its a long-term prediction
round(expm1(eq_aphid[1]),0) # equilibrium aphid abundance

# make plot of ecological dynamics
eco_dynamics_plot <- bind_rows(epred_lnAt1 %>% filter(Treatment %in% c("RYG","RYGP")), 
                               epred_lnPt1 %>% filter(Treatment %in% c("RYG","RYGP"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))) %>%
  ggplot(aes(x = week, y = estimate)) +
  geom_line(aes(color = resp), linewidth = 1) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = resp), alpha = 0.3) +
  facet_wrap(~Treatment, ncol = 1) +
  geom_line(data = bind_rows(raw_lnAt1, raw_lnPt1) %>% filter(week != 4) %>% mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))), 
            aes(x = week, y = log1pNt1, color = resp, group = paste(resp,ID)),
            alpha = 0.1) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) + #, breaks = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(0,2,4,6,8,10,12,14)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(data = eco_equilib_df, aes(x = week, y = aphid_equilibrium), size = 2, color = cbbPalette[1])

# get expected value for red morph frequency dynamics in polycultures
epred_Rt1 <- predictions(ecoevo_dynamics_brm, 
                         by = c("week","Treatment"),
                         newdata = filter(ecoevo_df, Freq_sub == 1), 
                         resp = "Rt1", 
                         re_formula = NULL, 
                         type = "response") %>%
  mutate(resp = "Rt1")

# get expected value for yellow morph frequency dynamics in polycultures
epred_Yt1 <- predictions(ecoevo_dynamics_brm, 
                         by = c("week","Treatment"),
                         newdata = filter(ecoevo_df, Freq_sub == 1), 
                         resp = "Yt1", 
                         re_formula = NULL, 
                         type = "response") %>%
  mutate(resp = "Yt1")

# raw red morph frequency dynamics
raw_Rt1 <- filter(df, Treatment %in% c("RYG","RYGP")) %>% 
  select(week, Treatment, ID, Freq_t1 = R_cage_freq) %>% 
  mutate(resp = "Rt1")

# raw yellow morph frequency dynamics
raw_Yt1 <- filter(df, Treatment %in% c("RYG","RYGP")) %>% 
  select(week, Treatment, ID, Freq_t1 = Y_cage_freq) %>% 
  mutate(resp = "Yt1")

# organize equilibrium frequency data for plot
evo_equilib_df <- data.frame(Treatment = factor(c("RYG","RYG","RYGP","RYGP"), 
                                                levels = c("RYG","RYGP"), 
                                                labels = c("Aphid","Aphid-Parasitoid")),
                             freq_equilibrium = c(eq_aphid[2], eq_aphid[3], NA, NA), # frequencies are not feasible in aphid-parasitoid system
                             resp = c("Rt1","Yt1","Rt1","Yt1"),
                             week = c(14,14, 14, 14))
round(eq_aphid[2],2) # red morph equilibrium frequency
round(eq_aphid[3],2) # yellow morph equilibrium frequency

# make plot of evolutionary dynamics
evo_dynamics_plot <- bind_rows(epred_Rt1, epred_Yt1) %>%
  mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))) %>%
  ggplot(aes(x = week, y = estimate)) +
  geom_line(aes(color = resp), linewidth = 1) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, fill = resp), alpha = 0.3) +
  facet_wrap(~Treatment, ncol = 1) +
  geom_line(data = bind_rows(raw_Rt1, raw_Yt1) %>% filter(week != 4) %>% mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP"), labels = c("Aphid","Aphid-Parasitoid"))), 
            aes(x = week, y = Freq_t1, color = resp, group = paste(resp,ID)),
            alpha = 0.1) +
  scale_color_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  scale_x_continuous(name = "Week", breaks = c(0,2,4,6,8,10,12,14)) +
  ylab("Frequency (z)") +
  geom_hline(yintercept = c(0,1), linetype = "dotted") +
  geom_point(data = evo_equilib_df, aes(x = week, y = freq_equilibrium, color = resp), size = 2)

# generate figure 3
plot_ecoevo_dynamics <- plot_grid(eco_dynamics_plot, # + ggtitle("Ecological dynamics"), 
                                  evo_dynamics_plot, # + ggtitle("Evolutionary dynamics"), 
                                  ncol = 2,
                                  labels = "AUTO")
plot_ecoevo_dynamics
save_plot(filename = "figures/ecoevo_dynamics_Fig3.pdf", 
          plot = plot_ecoevo_dynamics, base_height = 6) 