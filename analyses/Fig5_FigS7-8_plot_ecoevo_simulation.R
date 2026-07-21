
## Setup ----

# load posterior predictions (and bayesian MAR1 model of eco-evo dynamics)
load("analyses/posteriors_ecoevo_model.RData")

# source function for simulation
source('code/simulate_MAR1_dynamics.R')

# load and view data
source('code/manage_data.R')

# load other required libraries
library(brms)
library(cowplot)

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Setup conditions for simulation ----

# get initial conditions: polyculture with parasitoid at week 7
inits <- filter(ecoevo_df_allWeeks_imputed, Treatment == "RYGP", week == 7) %>%
  summarise(mean_lnAt1 = mean(lnAt1, na.rm = T),
            mean_lnPt1 = mean(lnPt1, na.rm = T),
            mean_Rt1 = mean(Rt1, na.rm = T),
            mean_Yt1 = mean(Yt1, na.rm = T)) %>%
  t()
inits

# duration in weeks (or time steps)
dur <- 150

# additions is a flexible matrix that would allow additions of species or traits
# after the first time point. Here we just use it for the first time point
additions <- matrix(data = c(inits, rep(0, (dur-1)*nrow(inits))), nrow = nrow(inits), ncol = dur)
additions

# if applicable, apply a press perturbation to the parasitoid (change in mortality rate)
# here, we increase the parasitoid's mortality rate to make the ecological-part of the system feasible
aphid_press <- 1.3

## Fig. 5: Simulate and plot eco-evolutionary dynamics for aphid-parasitoid system with a feasible ecological equilibrium ----

# simulation with all eco-evolutionary effects
ecoevo_sim <- EcoEvo_MAR1_Dynamics(Additions = additions,
                                   X.vector = x_aphid.ptoid_median + c(aphid_press,0,0,0),
                                   J.matrix = M_aphid.ptoid_median, 
                                   Abund.Positions = c(1,2),
                                   Freq.Positions = c(3,4),
                                   Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP","R","Y"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco", "evo"),
         sim_type = "Eco-Evo feedback") 

# check that evolutionary constraints (1 > freq > 0, and sum_freq < 1) is respected
ecoevo_sim %>%
  select(-ecoevo) %>%
  mutate(Week.adj = Week + 7) %>%
  pivot_wider(names_from = "response", values_from = "value") %>%
  mutate(sum_freq = R + Y) %>%
  filter(sum_freq > 1) 
# which it is

# simulation with only ecological effects.
# we do this by subsetting the relevant matrices and vectors to only the ecological parts
ecoonly_sim <- EcoEvo_MAR1_Dynamics(Additions = additions[c(1,2),],
                                    X.vector = as.matrix(x_aphid.ptoid_median[c("lnA","lnP"),]) + c(aphid_press,0), 
                                    J.matrix = M_aphid.ptoid_median[c("lnA","lnP"),c("lnA","lnP")], 
                                    Abund.Positions = c(1,2),
                                    Freq.Positions = NULL,
                                    Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco"),
         sim_type = "Ecological feedback") 

# Fig. 5
plot_nonequilib_ecoevo_effect <- bind_rows(ecoevo_sim, ecoonly_sim) %>% 
  filter(ecoevo == "eco") %>%
  mutate(Week.adj = Week + 7) %>%
  ggplot(aes(x = Week.adj, y = value, color = response, linetype = sim_type, group = paste(response, sim_type))) +
  geom_line() +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(0,25,50,75,100,125,150)) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_linetype_manual(name = NULL, values = c("dashed","solid"), labels = c("Eco-evo","Ecology only")) +
  geom_hline(yintercept = 0, linetype = "dotted")
plot_nonequilib_ecoevo_effect
save_plot("figures/Fig5_nonequilibrium_ecoevo_effect.pdf", plot_nonequilib_ecoevo_effect, base_height = 4)

# get ecological values for last week of simulation
ecoevo_sim %>% filter(ecoevo == "eco", Week == 149)
ecoonly_sim %>% filter(Week == 149)

# eco-evo effect on parasitoid abundances
expm1(4.87) - expm1(5.55) # reduced by 127 individuals
(expm1(4.87)/expm1(5.55)) - 1 # reduced by 50% on raw scale

# eco-evo effect on aphid abundances
expm1(7.21) - expm1(7.61) # reduced aphid abundance by ~665 individuals
(expm1(7.21)/expm1(7.61))-1 # reduced aphid abundance by 33%


## Fig. S7: Simulate and plot evolutionary dynamics for aphid-parasitoid system with a feasible ecological equilibrium ----
# using same simulation data as for Fig. 5, but only showing evolutionary dynamics
plot_nonequilib_evo_dyn <- bind_rows(ecoevo_sim, ecoonly_sim) %>% 
  filter(ecoevo == "evo") %>%
  mutate(Week.adj = Week + 7) %>% # based on initial conditions
  ggplot(aes(x = Week.adj, y = value, color = response, group = paste(response, sim_type))) +
  geom_line() +
  scale_color_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  scale_x_continuous("Week", breaks = c(0,25,50,75,100,125,150)) + 
  ylab("Frequency (z)") 
plot_nonequilib_evo_dyn
save_plot("figures/FigS7_nonequilibrium_evo_dyn.pdf", plot_nonequilib_evo_dyn, base_height = 5)

## Fig. S8: Simulate and plot eco-evolutionary dynamics for aphid-parasitoid system WITHOUT a feasible ecological equilibrium ----
# only real change from code for Fig. 5 is that we do not apply a press perturbation to parasitoid mortality rates
# i.e. this uses the estimates from the data.

# simulation with all eco-evolutionary effects
ecoevo_sim_nonFeasibileEquilib <- EcoEvo_MAR1_Dynamics(Additions = additions,
                                   X.vector = x_aphid.ptoid_median, # + c(aphid_press,0,0,0),
                                   J.matrix = M_aphid.ptoid_median, 
                                   Abund.Positions = c(1,2),
                                   Freq.Positions = c(3,4),
                                   Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP","R","Y"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco", "evo"),
         sim_type = "Eco-Evo feedback") 

# simulation with only ecological effects.
# we do this by subsetting the relevant matrices and vectors to only the ecological parts
ecoonly_sim_nonFeasibileEquilib <- EcoEvo_MAR1_Dynamics(Additions = additions[c(1,2),],
                                    X.vector = as.matrix(x_aphid.ptoid_median[c("lnA","lnP"),]), # + c(aphid_press,0),
                                    J.matrix = M_aphid.ptoid_median[c("lnA","lnP"),c("lnA","lnP")], 
                                    Abund.Positions = c(1,2),
                                    Freq.Positions = NULL,
                                    Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco"),
         sim_type = "Ecological feedback") 

# Fig. S8
plot_nonequilib_ecoevo_effect_nonFeasibleEquilib <- bind_rows(ecoevo_sim_nonFeasibileEquilib, ecoonly_sim_nonFeasibileEquilib) %>% 
  filter(ecoevo == "eco") %>%
  mutate(Week.adj = Week + 7) %>%
  ggplot(aes(x = Week.adj, y = value, color = response, linetype = sim_type, group = paste(response, sim_type))) +
  geom_line() +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(7:20)) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_linetype_manual(name = NULL, values = c("dashed","solid"), labels = c("Eco-evo","Ecology only")) +
  geom_hline(yintercept = 0, linetype = "dotted")
plot_nonequilib_ecoevo_effect_nonFeasibleEquilib
save_plot("figures/FigS8_nonequilibrium_ecoevo_effect.pdf", plot_nonequilib_ecoevo_effect_nonFeasibleEquilib, base_height = 4)


## Compare ecological dynamics of red and yellow morph monocultures ----
inits5 <- filter(ecoevo_df_allWeeks_imputed, Treatment %in% c("RP","YP"), week == 5) %>%
  summarise(mean_lnAt1 = mean(lnAt1, na.rm = T),
            mean_lnPt1 = mean(lnPt1, na.rm = T),
            mean_Rt1 = mean(Rt1, na.rm = T),
            mean_Yt1 = mean(Yt1, na.rm = T)) %>%
  t()
inits5

M_aphid.ptoid_median_RP.YP.GP <- M_aphid.ptoid_median
M_aphid.ptoid_median_RP.YP.GP[c("R","Y"),] <- 0
M_aphid.ptoid_median_RP.YP.GP["R","R"] <- 1
M_aphid.ptoid_median_RP.YP.GP["Y","Y"] <- 1
x_aphid.ptoid_median_RP.YP.GP <- x_aphid.ptoid_median
x_aphid.ptoid_median_RP.YP.GP[c("R","Y"),] <- 0

additions <- matrix(data = c(inits5, rep(0, (dur-1)*nrow(inits5))), nrow = nrow(inits5), ncol = dur)
additions_RP <- additions
additions_RP[3,1] <- 1 # R
additions_RP[4,1] <- 0 # Y
RP_sim <- EcoEvo_MAR1_Dynamics(Additions = additions_RP,
                     X.vector = x_aphid.ptoid_median_RP.YP.GP,
                     J.matrix = M_aphid.ptoid_median_RP.YP.GP, 
                     Abund.Positions = c(1,2),
                     Freq.Positions = c(3,4),
                     Duration = dur)  %>%
  pivot_longer(cols = c("lnA","lnP","R","Y"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco", "evo"),
         sim_type = "Eco-Evo feedback") 

additions_YP <- additions
additions_YP[3,1] <- 0 # R
additions_YP[4,1] <- 1 # Y
YP_sim <- EcoEvo_MAR1_Dynamics(Additions = additions_YP,
                               X.vector = x_aphid.ptoid_median_RP.YP.GP,
                               J.matrix = M_aphid.ptoid_median_RP.YP.GP, 
                               Abund.Positions = c(1,2),
                               Freq.Positions = c(3,4),
                               Duration = dur)  %>%
  pivot_longer(cols = c("lnA","lnP","R","Y"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco", "evo"),
         sim_type = "Eco-Evo feedback")

bind_rows(RP_sim %>% mutate(Treatment = "RP"), 
          YP_sim %>% mutate(Treatment = "YP")) %>% 
  filter(ecoevo == "eco") %>%
  mutate(Week.adj = Week + 7) %>%
  ggplot(aes(x = Week.adj, y = value, color = Treatment, group = paste(response, Treatment))) +
  geom_line() +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(7:20)) +
  scale_color_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~response)
# qualitatively matches the experimental analysis
