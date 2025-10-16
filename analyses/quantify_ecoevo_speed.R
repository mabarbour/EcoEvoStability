## Setup ----

# load and view data
source('code/manage_data.R')
ecoevo_df 

# load other required libraries
library(lmerTest)
library(marginaleffects)
library(cowplot)

# set plot theme
theme_set(theme_cowplot())

# manage data for eco-evo speed analyses 
tG <- 7/10 # 7 day sampling period, 10 day generation time 
speed_df <- ecoevo_df %>%
  filter(Freq_sub == 1, week > 6, Aphids_cage_lag1 > 0, Aphids_cage > 0) %>%
  mutate(eco_speed = abs(log(Aphids_cage/Aphids_cage_lag1)),
         # method below is more appropriate for continuous time sampling, rather than discrete time
         # eco_speed_alt = abs((Aphids_cage - Aphids_cage_lag1)/Aphids_cage_lag1*tG),
         evo_speed_Y = abs((Yt1 - Yt)/Yt*tG),
         evo_speed_R = abs((Rt1 - Rt)/Rt*tG),
         max_evo_speed = pmax(evo_speed_Y, evo_speed_R),
         speed_ratio_Y = evo_speed_Y/eco_speed,
         speed_ratio_R = evo_speed_R/eco_speed,
         speed_ratio_max = max_evo_speed/eco_speed) 

## fit and analyse model ----

# fit model
log_speed_ratio_m1 <- lmer(log(speed_ratio_max) ~ week*ptoid_in + (1|ID), data = speed_df)
summary(log_speed_ratio_m1) # suggests random effect for cage ID could be dropped

# overall speed ratio
avg_predictions(log_speed_ratio_m1, re.form = NA, transform = exp) 

# change in speed ratio over time
avg_comparisons(log_speed_ratio_m1, variables = "week", re.form = NA, transform = exp) 

# effect of adding parasitoid
avg_comparisons(log_speed_ratio_m1, variables = "ptoid_in", re.form = NA, transform = exp)

# does the parasitoid alter the speed ratio over time (i.e. statistical interaction)?
avg_comparisons(log_speed_ratio_m1, variables = "week", by = "ptoid_in", re.form = NA, hypothesis = "b2 - b1 = 0", transform = exp)
