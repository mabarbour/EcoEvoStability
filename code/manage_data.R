## Setup ----

# load libraries
library(tidyverse)

# load and view data
aphid_data_final <- read_csv("data/aphid_data_final.csv")
aphid_data_final

# manipulate data for analysis
df <- aphid_data_final %>%
  mutate(# Estimate aphid abundance for the entire cage. First calculate aphids per plant then multiply by 16 (total number of plants per cage)
         R_cage = as.integer((R*aphids_per_click)/plants*16),
         Y_cage = as.integer((Y*aphids_per_click)/plants*16),
         G_cage = as.integer((G*aphids_per_click)/plants*16),
         M_cage = as.integer((M*1)/plants*16), # always 1 per click
         P_cage = P + ptoids_added, # adult ptoids always estimated per cage, includes 2 ptoids added on weeks 3 and 4
         W_cage = W, # mostly an estimate per cage, but a small fraction of these aphids come from plant counts...
         Aphids_cage = R_cage + Y_cage + G_cage, # non-winged
         Ptoids_cage = P_cage + M_cage, # add adult and mummy stages together to determine total parasitoid abundance
         R_cage_freq = R_cage/Aphids_cage,
         Y_cage_freq = Y_cage/Aphids_cage,
         G_cage_freq = G_cage/Aphids_cage,
         Treatment = ifelse(ID %in% c("3RP", "8RP", "13RP", "18RP", "23RP", "28RP"), "RP",
                            ifelse(ID %in% c("4YP", "9YP", "14YP", "19YP", "24YP", "29YP"), "YP",
                                   ifelse(ID %in% c("5GP", "10GP", "15GP", "20GP", "25GP", "30GP"), "GP",
                                          ifelse(ID %in% c("2RYGP", "7RYGP", "12RYGP", "17RYGP", "22RYGP", "27RYGP"), "RYGP",
                                                 ifelse(ID %in% c("1RYG", "6RYG", "11RYG", "16RYG", "21RYG", "26RYG"), "RYG", NA))))),
         ptoid_treatments = ifelse(Treatment %in% c("RP","YP","GP","RYGP"), "Ptoids", "No ptoids"),
         gen_div = ifelse(Treatment %in% c("RP","YP","GP"), 1, 3),
         f.week = factor(week, levels = sort(unique(week)))) %>%
  group_by(ID) %>%
  mutate(# determine lagged abundances for MAR(1) models
         R_cage_lag1 = lag(R_cage, n = 1, order_by = week),
         Y_cage_lag1 = lag(Y_cage, n = 1, order_by = week),
         G_cage_lag1 = lag(G_cage, n = 1, order_by = week),
         M_cage_lag1 = lag(M_cage, n = 1, order_by = week),
         P_cage_lag1 = lag(P_cage, n = 1, order_by = week),
         W_cage_lag1 = lag(W_cage, n = 1, order_by = week),
         Aphids_cage_lag1 = R_cage_lag1 + Y_cage_lag1 + G_cage_lag1,
         Ptoids_cage_lag1 = M_cage_lag1 + P_cage_lag1,
         # determine lag 2 abundances for MAR(P) models
         R_cage_lag2 = lag(R_cage, n = 2, order_by = week),
         Y_cage_lag2 = lag(Y_cage, n = 2, order_by = week),
         G_cage_lag2 = lag(G_cage, n = 2, order_by = week),
         M_cage_lag2 = lag(M_cage, n = 2, order_by = week),
         P_cage_lag2 = lag(P_cage, n = 2, order_by = week),
         W_cage_lag2 = lag(W_cage, n = 2, order_by = week),
         Aphids_cage_lag2 = R_cage_lag2 + Y_cage_lag2 + G_cage_lag2,
         Ptoids_cage_lag2 = M_cage_lag2 + P_cage_lag2,
         # lagged frequencies may be useful too
         R_cage_freq_lag1 = lag(R_cage_freq, n = 1, order_by = week),
         Y_cage_freq_lag1 = lag(Y_cage_freq, n = 1, order_by = week),
         G_cage_freq_lag1 = lag(G_cage_freq, n = 1, order_by = week)) %>%
  ungroup()
df # view data

# subset and rename columns for downstream analyses 
ecoevo_df <- df %>% 
  filter(week > 0, week < 4 | week > 6) %>% # exclude weeks of non-sampling and heat-wave
  # mutate and rename variables to facilitate model processing
  mutate(ptoid_in = ifelse(Treatment %in% c("RYGP","RP","YP","GP") & week > 3, 1, 0),
         Aphids_sub = ifelse(Aphids_cage_lag1 > 0 | aphids_extinct == 0, 1, 0), 
         Ptoids_sub = ifelse(Ptoids_cage_lag1 > 0 | ptoids_extinct == 0 & ptoid_in == 1, 1, 0), #  & ptoid_in == 1
         Freq_sub = ifelse(Treatment %in% c("RYG","RYGP"), 1, 0),
         # log1p abundance variables
         lnAt1 = log1p(Aphids_cage),
         lnAt = log1p(Aphids_cage_lag1),
         lnPt1 = log1p(Ptoids_cage),
         lnPt = log1p(Ptoids_cage_lag1), 
         # shorten names for frequency variables
         Rt1 = R_cage_freq,
         Rt = R_cage_freq_lag1,
         Gt1 = G_cage_freq,
         Gt = G_cage_freq_lag1,
         Yt1 = Y_cage_freq,
         Yt = Y_cage_freq_lag1) %>% 
  # focus on select data columns for subsequent analyses
  select(week, ID, Treatment, 
         lnAt1, lnAt, lnPt1, lnPt, 
         Rt1, Rt, Gt1, Gt, Yt1, Yt,
         Aphids_cage, Aphids_cage_lag1,
         ptoid_in, Aphids_sub, Ptoids_sub, Freq_sub) 

# set NA in subsets to zero for brms to work
ecoevo_df$Aphids_sub[is.na(ecoevo_df$Aphids_sub) == TRUE] <- 0
ecoevo_df$Ptoids_sub[is.na(ecoevo_df$Ptoids_sub) == TRUE] <- 0
ecoevo_df$Freq_sub[is.na(ecoevo_df$Freq_sub) == TRUE] <- 0

# view data
ecoevo_df
