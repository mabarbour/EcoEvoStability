
## Setup ----

# source function for simulation
source('code/simulate_MAR1_dynamics.R')

# load and view data
source('code/manage_data.R')
ecoevo_df 

# load other required libraries
library(brms)
library(cowplot)

# load posterior predictions (and bayesian MAR1 model of eco-evo dynamics)
load("analyses/posteriors_ecoevo_model.RData")

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Setup conditions for simulation ----

# get initial conditions: polyculture with parasitoid at week 7
inits <- filter(ecoevo_df, Treatment == "RYGP", week == 7) %>%
  summarise(mean_lnAt1 = mean(lnAt1, na.rm = T),
            mean_lnPt1 = mean(lnPt1, na.rm = T),
            mean_Rt1 = mean(Rt1, na.rm = T),
            mean_Yt1 = mean(Yt1, na.rm = T)) %>%
  t()
inits

# duration in weeks (or time steps)
dur <- 100 

# additions is a flexible matrix that would allow additions of species or traits
# after the first time point. Here we just use it for the first time point
additions <- matrix(data = c(inits, rep(0, (dur-1)*nrow(inits))), nrow = nrow(inits), ncol = dur)
additions

# if applicable, apply a press perturbation to the parasitoid (change in mortality rate)
# here, we increase the parasitoid's mortality rate to make the ecological-part of the system feasible
ptoid_press <- -2

## Fig. 5: Simulate and plot eco-evolutionary dynamics for aphid-parasitoid system with a feasible ecological equilibrium ----

# simulation with all eco-evolutionary effects
I.matrix_ecoevo <- diag(1, nrow = nrow(M_aphid.ptoid_median)) # identity matrix for eco-evo effects
ecoevo_sim <- EcoEvo_MAR1_Dynamics(Additions = additions,
                                   X.vector = x_aphid.ptoid_median + c(0,ptoid_press,0,0),
                                   J.matrix = (I.matrix_ecoevo + M_aphid.ptoid_median), 
                                   Abund.Positions = c(1,2),
                                   Freq.Positions = c(3,4),
                                   Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP","R","Y"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco", "evo"),
         sim_type = "Eco-Evo feedback") 

# at week.adj 21 (start with initial conditions at week 7), 
# the predicted frequencies are no longer biologically plausible (sum of frequences > 1), 
# therefore, we'll restrict our comparisons to before this timepoint.
ecoevo_sim %>%
  select(-ecoevo) %>%
  mutate(Week.adj = Week + 7) %>%
  pivot_wider(names_from = "response", values_from = "value") %>%
  mutate(sum_freq = R + Y) %>%
  filter(sum_freq > 1) 

# simulation with only ecological effects.
# we do this by subsetting the relevant matrices and vectors to only the ecological parts
I.matrix_ecoonly <- diag(1, nrow = nrow(M_aphid.ptoid_median[c("lnA","lnP"),c("lnA","lnP")])) # identity matrix for ecological effects only
ecoonly_sim <- EcoEvo_MAR1_Dynamics(Additions = additions[c(1,2),],
                                    X.vector = as.matrix(x_aphid.ptoid_median[c("lnA","lnP"),]) + c(0,ptoid_press),
                                    J.matrix = (I.matrix_ecoonly + M_aphid.ptoid_median[c("lnA","lnP"),c("lnA","lnP")]), 
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
  filter(Week.adj < 21) %>% # impose evolutionary constraint
  ggplot(aes(x = Week.adj, y = value, color = response, linetype = sim_type, group = paste(response, sim_type))) +
  geom_line() +
  scale_y_continuous(name = "Abundance (N)", 
                     breaks = log1p(c(0,10,100,1000,10000)),
                     labels = c(0,10,100,1000,10000)) +
  scale_x_continuous(name = "Week", breaks = c(7:20)) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Parasitoid")) +
  scale_linetype_manual(name = NULL, values = c("dashed","solid"), labels = c("Eco-evo","Ecology only")) +
  geom_hline(yintercept = 0, linetype = "dotted")
plot_nonequilib_ecoevo_effect
save_plot("figures/nonequilibrium_ecoevo_effect_Fig5.pdf", plot_nonequilib_ecoevo_effect, base_height = 4)

# get values for last week of simulate before evolutionary constraint breaks
bind_rows(ecoevo_sim, ecoonly_sim) %>%
  mutate(Week.adj = Week + 7) %>%
  filter(Week.adj == 20) %>%
  select(-ecoevo) %>%
  pivot_wider(names_from = "response", values_from = "value")

# eco-evo effect on parasitoid abundances
expm1(2.57) - expm1(3.02) # reduced by 7-8 individuals
(expm1(2.57)/expm1(3.02)) - 1 # reduced by 38% on raw scale
(2.57/3.02)-1 # reduced by 15% on log scale

# eco-evo effect on aphid abundances
expm1(6.23) - expm1(6.70) # reduced aphid abundance by ~305 individuals
(expm1(6.23)/expm1(6.70))-1 # reduced aphid abundance by 38%
(6.23/6.70)-1 # reduced by 7% on log scale


## Fig. S1: Simulate and plot evolutionary dynamics for aphid-parasitoid system with a feasible ecological equilibrium ----
# using same simulation data as for Fig. 5, but only showing evolutionary dynamics
plot_nonequilib_evo_dyn <- bind_rows(ecoevo_sim, ecoonly_sim) %>% 
  filter(ecoevo == "evo") %>%
  mutate(Week.adj = Week + 7) %>% # based on initial conditions
  filter(Week.adj < 21) %>% # impose evolutionary constraint
  ggplot(aes(x = Week.adj, y = value, color = response, group = paste(response, sim_type))) +
  geom_line() +
  scale_color_manual(name = NULL, values = cbbPalette[c(7,5)], labels = c("Red morph","Yellow morph")) +
  scale_x_continuous("Week", breaks = 7:20) + 
  ylab("Frequency (z)") 
plot_nonequilib_evo_dyn
save_plot("figures/nonequilibrium_evo_dyn_FigS1.pdf", plot_nonequilib_evo_dyn, base_height = 5)

## Fig. S2: Simulate and plot eco-evolutionary dynamics for aphid-parasitoid system WITHOUT a feasible ecological equilibrium ----
# only real change from code for Fig. 5 is that we do not apply a press perturbation to parasitoid mortality rates
# i.e. this uses the estimates from the data.

# simulation with all eco-evolutionary effects
ecoevo_sim_nonFeasibileEquilib <- EcoEvo_MAR1_Dynamics(Additions = additions,
                                   X.vector = x_aphid.ptoid_median, # + c(0,ptoid_press,0,0),
                                   J.matrix = (I.matrix_ecoevo + M_aphid.ptoid_median), 
                                   Abund.Positions = c(1,2),
                                   Freq.Positions = c(3,4),
                                   Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP","R","Y"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco", "evo"),
         sim_type = "Eco-Evo feedback") 

# simulation with only ecological effects.
# we do this by subsetting the relevant matrices and vectors to only the ecological parts
ecoonly_sim_nonFeasibileEquilib <- EcoEvo_MAR1_Dynamics(Additions = additions[c(1,2),],
                                    X.vector = as.matrix(x_aphid.ptoid_median[c("lnA","lnP"),]), # + c(0,ptoid_press),
                                    J.matrix = (I.matrix_ecoonly + M_aphid.ptoid_median[c("lnA","lnP"),c("lnA","lnP")]), 
                                    Abund.Positions = c(1,2),
                                    Freq.Positions = NULL,
                                    Duration = dur) %>%
  pivot_longer(cols = c("lnA","lnP"), names_to = "response") %>%
  mutate(ecoevo = ifelse(response %in% c("lnA","lnP"), "eco"),
         sim_type = "Ecological feedback") 

# Fig. S2
plot_nonequilib_ecoevo_effect_nonFeasibleEquilib <- bind_rows(ecoevo_sim_nonFeasibileEquilib, ecoonly_sim_nonFeasibileEquilib) %>% 
  filter(ecoevo == "eco") %>%
  mutate(Week.adj = Week + 7) %>%
  filter(Week.adj < 21) %>% # impose evolutionary constraint
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
save_plot("figures/nonequilibrium_ecoevo_effect_FigS2.pdf", plot_nonequilib_ecoevo_effect_nonFeasibleEquilib, base_height = 4)
