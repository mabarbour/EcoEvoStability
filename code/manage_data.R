## Setup ----

# load libraries
library(tidyverse)

# load and view data
aphid_data_final <- read_csv("data/aphid_data_final.csv")
aphid_data_final

# ptoids_added should be zero for all RYG treatments (no parasitoids) in week 3 and 4, but it isn't (data entry error)
aphid_data_final <- aphid_data_final %>%
  mutate(ptoids_added = ifelse(ID %in% c("1RYG", "6RYG", "11RYG", "16RYG", "21RYG", "26RYG"), 0, ptoids_added))

filter(aphid_data_final, week %in% c(3,4), ID %in% c("1RYG", "6RYG", "11RYG", "16RYG", "21RYG", "26RYG")) %>%
  select(week, ID, ptoids_added) # double-checked that all is good

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
         ptoid_treatments = ifelse(Treatment %in% c("RP","YP","GP","RYGP"), "Ptoids", "No ptoids")) %>%
         # gen_div = ifelse(Treatment %in% c("RP","YP","GP"), 1, 3),
         # f.week = factor(week, levels = sort(unique(week)))) %>%
  group_by(ID) %>%
  mutate(# for Aphids and Ptoids _cage_lag1
         R_cage_lag1 = lag(R_cage, n = 1, order_by = week),
         Y_cage_lag1 = lag(Y_cage, n = 1, order_by = week),
         G_cage_lag1 = lag(G_cage, n = 1, order_by = week),
         M_cage_lag1 = lag(M_cage, n = 1, order_by = week),
         P_cage_lag1 = lag(P_cage, n = 1, order_by = week),
         W_cage_lag1 = lag(W_cage, n = 1, order_by = week),
         # for downstream filtering
         Aphids_cage_lag1 = R_cage_lag1 + Y_cage_lag1 + G_cage_lag1,
         Ptoids_cage_lag1 = M_cage_lag1 + P_cage_lag1) %>%
         # # determine lag 2 abundances for MAR(P) models
         # R_cage_lag2 = lag(R_cage, n = 2, order_by = week),
         # Y_cage_lag2 = lag(Y_cage, n = 2, order_by = week),
         # G_cage_lag2 = lag(G_cage, n = 2, order_by = week),
         # M_cage_lag2 = lag(M_cage, n = 2, order_by = week),
         # P_cage_lag2 = lag(P_cage, n = 2, order_by = week),
         # W_cage_lag2 = lag(W_cage, n = 2, order_by = week),
         # Aphids_cage_lag2 = R_cage_lag2 + Y_cage_lag2 + G_cage_lag2,
         # Ptoids_cage_lag2 = M_cage_lag2 + P_cage_lag2,
         # lagged frequencies may be useful too
         # R_cage_freq_lag1 = lag(R_cage_freq, n = 1, order_by = week),
         # Y_cage_freq_lag1 = lag(Y_cage_freq, n = 1, order_by = week),
         # G_cage_freq_lag1 = lag(G_cage_freq, n = 1, order_by = week),
         # R_cage_freq_lag2 = lag(R_cage_freq, n = 2, order_by = week),
         # Y_cage_freq_lag2 = lag(Y_cage_freq, n = 2, order_by = week),
         # G_cage_freq_lag2 = lag(G_cage_freq, n = 2, order_by = week)) %>%
  ungroup()
df # view data

# # subset and rename columns for downstream analyses 
# ecoevo_df <- df %>% 
#   filter(week > 0, week < 4 | week > 6) %>% # exclude weeks of non-sampling and heat-wave
#   # mutate and rename variables to facilitate model processing
#   mutate(ptoid_in = ifelse(Treatment %in% c("RYGP","RP","YP","GP") & week > 3, 1, 0),
#          Aphids_sub = ifelse(Aphids_cage_lag1 > 0 | aphids_extinct == 0, 1, 0), 
#          Ptoids_sub = ifelse(Ptoids_cage_lag1 > 0 | ptoids_extinct == 0 & ptoid_in == 1, 1, 0), #  & ptoid_in == 1
#          Freq_sub = ifelse(Treatment %in% c("RYG","RYGP"), 1, 0),
#          # log1p abundance variables
#          lnAt1 = log1p(Aphids_cage),
#          lnAt = log1p(Aphids_cage_lag1),
#          lnPt1 = log1p(Ptoids_cage),
#          lnPt = log1p(Ptoids_cage_lag1), 
#          # shorten names for frequency variables
#          Rt1 = R_cage_freq,
#          Rt = R_cage_freq_lag1,
#          Gt1 = G_cage_freq,
#          Gt = G_cage_freq_lag1,
#          Yt1 = Y_cage_freq,
#          Yt = Y_cage_freq_lag1,
#          # lag2 variables
#          lnAt0 = log1p(Aphids_cage_lag2),
#          lnPt0 = log1p(Ptoids_cage_lag2),
#          Rt0 = R_cage_freq_lag2,
#          Yt0 = Y_cage_freq_lag2) %>% 
#   # focus on select data columns for subsequent analyses
#   select(week, ID, Treatment, 
#          lnAt1, lnAt, lnPt1, lnPt, lnAt0, lnPt0, Rt0, Yt0,
#          Rt1, Rt, Gt1, Gt, Yt1, Yt,
#          Aphids_cage, Aphids_cage_lag1,
#          ptoid_in, Aphids_sub, Ptoids_sub, Freq_sub, plants_no_growth) 
# 
# # set NA in subsets to zero for brms to work
# ecoevo_df$Aphids_sub[is.na(ecoevo_df$Aphids_sub) == TRUE] <- 0
# ecoevo_df$Ptoids_sub[is.na(ecoevo_df$Ptoids_sub) == TRUE] <- 0
# ecoevo_df$Freq_sub[is.na(ecoevo_df$Freq_sub) == TRUE] <- 0
# 
# # view data
# ecoevo_df
# 
# 
# 
# ## ecoevo_df with missing data
# ecoevo_df_allWeeks <- df %>% 
#   filter(week > 0) %>% # include weeks of non-sampling and heat-wave
#   # mutate and rename variables to facilitate model processing
#   mutate(ptoid_in = ifelse(Treatment %in% c("RYGP","RP","YP","GP") & week > 3, 1, 0),
#          Aphids_sub = ifelse(Aphids_cage_lag1 > 0 | aphids_extinct == 0, 1, 0), 
#          Ptoids_sub = ifelse(Ptoids_cage_lag1 > 0 | ptoids_extinct == 0 & ptoid_in == 1, 1, 0), #  & ptoid_in == 1
#          Freq_sub = ifelse(Treatment %in% c("RYG","RYGP"), 1, 0),
#          # log1p abundance variables
#          lnAt1 = log1p(Aphids_cage),
#          lnAt = log1p(Aphids_cage_lag1),
#          lnPt1 = log1p(Ptoids_cage),
#          lnPt = log1p(Ptoids_cage_lag1), 
#          # shorten names for frequency variables
#          Rt1 = R_cage_freq,
#          Rt = R_cage_freq_lag1,
#          Gt1 = G_cage_freq,
#          Gt = G_cage_freq_lag1,
#          Yt1 = Y_cage_freq,
#          Yt = Y_cage_freq_lag1) %>%
#          # lag2 variables
#          #lnAt0 = log1p(Aphids_cage_lag2),
#          #lnPt0 = log1p(Ptoids_cage_lag2),
#          #Rt0 = R_cage_freq_lag2,
#          #Yt0 = Y_cage_freq_lag2) %>% 
#   # focus on select data columns for subsequent analyses
#   select(week, ID, Treatment, 
#          lnAt1, lnAt, lnPt1, lnPt, #lnAt0, lnPt0, Rt0, Yt0,
#          Rt1, Rt, Gt1, Gt, Yt1, Yt,
#          #Aphids_cage, Aphids_cage_lag1,
#          ptoid_in, Aphids_sub, Ptoids_sub, Freq_sub, plants_no_growth) 
# 
# # set NA in subsets to zero for brms to work
# ecoevo_df_allWeeks$Aphids_sub[is.na(ecoevo_df_allWeeks$Aphids_sub) == TRUE] <- 0
# ecoevo_df_allWeeks$Ptoids_sub[is.na(ecoevo_df_allWeeks$Ptoids_sub) == TRUE] <- 0
# ecoevo_df_allWeeks$Freq_sub[is.na(ecoevo_df_allWeeks$Freq_sub) == TRUE] <- 0
# 
# # view data
# ecoevo_df_allWeeks


### Imputed dataset ----

# first some prep steps
prep_df <- df %>% 
  # mutate and rename variables to facilitate model processing
  mutate(ptoid_in = ifelse(Treatment %in% c("RYGP","RP","YP","GP") & week > 3, 1, 0),
         Aphids_sub = ifelse(Aphids_cage_lag1 > 0 | aphids_extinct == 0, 1, 0), 
         Ptoids_sub = ifelse(Ptoids_cage_lag1 > 0 | ptoids_extinct == 0 & ptoid_in == 1, 1, 0), 
         # with abundance cutoff #Freq_sub = ifelse(Aphids_cage >= 26 & Treatment %in% c("RYG","RYGP"), 1, 0), # 26 chosen to maximize the number of data points, capture the drop in aphid abundances. Next lowest value is 6... 
         # without abundance cutoff, added Aphids_cage > 0 to avoid imputing missing value for when no aphids were detected. 
         Freq_sub = ifelse(aphids_extinct == 0 & Aphids_cage > 0 & Treatment %in% c("RYG","RYGP"), 1, 0), # is.na(R_cage_freq) == FALSE & is.na(Y_cage_freq) == FALSE
         # log1p abundance variables
         lnAt1 = log1p(Aphids_cage),
         lnPt1 = log1p(Ptoids_cage),
         # shorten names for raw abundance variables
         At1 = Aphids_cage,
         Pt1 = Ptoids_cage,
         R_At1 = R_cage,
         Y_At1 = Y_cage,
         G_At1 = G_cage,
         # shorten names for frequency variables
         Rt1 = R_cage_freq,
         Yt1 = Y_cage_freq,
         Gt1 = G_cage_freq) %>%
  # focus on select data columns for subsequent analyses
  select(week, ID, Treatment, 
         lnAt1, lnPt1,
         At1, Pt1,
         R_At1, Y_At1, G_At1,
         Rt1, Yt1, Gt1, 
         ptoid_in, Aphids_sub, Ptoids_sub, Freq_sub, plants_no_growth)

# for each cage, impute the mean between weeks 3 and 5
imputed_week4_data <- prep_df %>%
  filter(week %in% c(3,5)) %>%
  group_by(ID, Treatment) %>%
  summarise(lnAt1 = mean(lnAt1),
            lnPt1 = mean(lnPt1), 
            At1 = mean(At1), # not using downstream, but doing for consistency
            Pt1 = mean(Pt1), # not using downstream, but doing for consistency
            R_At1 = mean(R_At1), # not using downstream, but doing for consistency
            Y_At1 = mean(Y_At1), # not using downstream, but doing for consistency
            G_At1 = mean(G_At1), # not using downstream, but doing for consistency
            Rt1 = mean(Rt1),
            Yt1 = mean(Yt1),
            Gt1 = mean(Gt1)) %>%
  mutate(week = 4,
         plants_no_growth = 0, # assuming that none of the plants reached this level
         ptoid_in = ifelse(Treatment %in% c("RYGP","RP","YP","GP"), 1, 0),
         Aphids_sub = 1,
         Ptoids_sub = ifelse(ptoid_in == 1, 1, 0),
         Freq_sub = ifelse(Treatment %in% c("RYG","RYGP"), 1, 0)) %>%
  ungroup() 

# insert imputed week 4 into dataset, then calculate lags
ecoevo_df_allWeeks_imputed_withWeek0 <- bind_rows(
  filter(prep_df, week < 4), 
  imputed_week4_data, 
  filter(prep_df, week > 4)) %>%
  group_by(ID) %>%
  mutate(
    # determine lagged abundances for MAR(1) models
    # much easier to determine lags after imputing week 4
    lnAt = lag(lnAt1, n = 1, order_by = week),
    lnPt = lag(lnPt1, n = 1, order_by = week),
    At = lag(At1, n = 1, order_by = week), # not using downstream, but doing for consistency
    Pt = lag(Pt1, n = 1, order_by = week), # not using downstream, but doing for consistency
    R_At = lag(R_At1, n = 1, order_by = week), # not using downstream, but doing for consistency
    Y_At = lag(Y_At1, n = 1, order_by = week), # not using downstream, but doing for consistency
    G_At = lag(G_At1, n = 1, order_by = week), # not using downstream, but doing for consistency
    Rt = lag(Rt1, n = 1, order_by = week),
    Yt = lag(Yt1, n = 1, order_by = week),
    Gt = lag(Gt1, n = 1, order_by = week)) %>%
  ungroup() %>%
  mutate(f.rowid = as.character(row_number()))

# set NA in subsets to zero for brms to work
ecoevo_df_allWeeks_imputed_withWeek0$Aphids_sub[is.na(ecoevo_df_allWeeks_imputed_withWeek0$Aphids_sub) == TRUE] <- 0
ecoevo_df_allWeeks_imputed_withWeek0$Ptoids_sub[is.na(ecoevo_df_allWeeks_imputed_withWeek0$Ptoids_sub) == TRUE] <- 0
ecoevo_df_allWeeks_imputed_withWeek0$Freq_sub[is.na(ecoevo_df_allWeeks_imputed_withWeek0$Freq_sub) == TRUE] <- 0

# view data
ecoevo_df_allWeeks_imputed_withWeek0

# ecoevo_df_allWeeks_imputed_withWeek0 %>% filter(ID == "12RYGP", week %in% c(13,14))
# # Cage ID 12RYGP had a difficult to interpret pattern betweek weeks 13 and 14
# # Week 13: No aphids or parasitoids were detected (lnAt1 and lnPt1 = 0), therefore no frequencies were recorded.
# # Week 14: Two aphid individuals detected, 1 yellow and 1 green, and no parasitoids were found.
# # Our model is intended to focus on the deterministic dynamics.
# # Given the potential shift in parasitoid treatment (no longer any parasitoids), stochasticity, and the chance of an unintentional introduction
# # the week 14 datapoint will be excluded.
# ID_12RYGP_week_14 <- filter(ecoevo_df_allWeeks_imputed_withWeek0, ID == "12RYGP", week == 14)$rowid
# ecoevo_df_allWeeks_imputed_withWeek0$Aphids_sub[ID_12RYGP_week_14] <- 0
# ecoevo_df_allWeeks_imputed_withWeek0$Freq_sub[ID_12RYGP_week_14] <- 0

ecoevo_df_allWeeks_imputed <- ecoevo_df_allWeeks_imputed_withWeek0 %>%
  filter(week > 0) 
# keeping a version above with week 0 for plotting raw data in other scripts

# source and join Myzus degree data by week
source("code/Degday_Myzus.R")
ecoevo_df_allWeeks_imputed <- left_join(ecoevo_df_allWeeks_imputed, select(Myzus_dd_week, week = week_id, EGN_sine))

# check data points were certain criteria ----
sum(expm1(ecoevo_df_allWeeks_imputed$lnAt) >= 50, na.rm = T) / 
  sum(expm1(ecoevo_df_allWeeks_imputed$lnAt) > 0, na.rm = T) # 85% of data had greater than 50 aphids
sum(expm1(ecoevo_df_allWeeks_imputed$lnAt) >= 100, na.rm = T) / 
  sum(expm1(ecoevo_df_allWeeks_imputed$lnAt) > 0, na.rm = T) # 84% of data had greater than 100 aphids

ecoevo_df_allWeeks_imputed %>%
  filter(Freq_sub == 1) %>%
  summarise(`Rt_btw_0.1-0.9` = sum(Rt > 0.1 & Rt < 0.9, na.rm = T) / nrow(.),
            `Yt_btw_0.1-0.9` = sum(Yt > 0.1 & Yt < 0.9, na.rm = T) / nrow(.))
# 99% of data had Rt and Yt values between 0.1 and 0.9

ecoevo_df_allWeeks_imputed %>%
  filter(Freq_sub == 1) %>%
  summarise(`Rt_btw_0.2-0.8` = sum(Rt > 0.2 & Rt < 0.8, na.rm = T) / nrow(.),
            `Yt_btw_0.2-0.8` = sum(Yt > 0.2 & Yt < 0.8, na.rm = T) / nrow(.))
# 94% of data had Rt and Yt values between 0.2 and 0.8

ggplot(filter(ecoevo_df_allWeeks_imputed, Freq_sub == 1)) +
  geom_boxplot(aes(x = 1, y = abs(Rt1 - Rt)), color = "red") +
  geom_boxplot(aes(x = 2, y = abs(Yt1 - Yt)), color = "yellow") +
  ylab("abs(Freq_change)")
ecoevo_df_allWeeks_imputed %>%
  filter(Freq_sub == 1) %>%
  summarise(`Rt_abs_less_0.2` = sum(abs(Rt1 - Rt) < 0.2, na.rm = T) / nrow(.),
            `Yt_abs_less_0.2` = sum(abs(Yt1 - Yt) < 0.2, na.rm = T) / nrow(.)) 
# 95% and 92% of weeks had absolute change of < 0.2 for Rt and Yt, respectively