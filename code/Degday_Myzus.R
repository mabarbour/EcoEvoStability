# Degree days of Myzus persicae from degday function ----
source("code/Temp_data.R")

#Establish Lower and upper tresholds in F (Whalon, M. E., and Z. Smilowitz. 1979)
Thresh_low = 39.2
Thresh_up = 86

# Establish Mean generation time (Whalon, M. E., and Z. Smilowitz. 1979)
MGT = 129.78
date_ref <- as.Date("2024-06-03")

Myzus_dd <- temp_journaliere %>%
  mutate(dbl_tri= dd_dbl_tri(daily_min = temp_min_f, daily_max = temp_max_f, nextday_min = tmin_next,
                             thresh_low = Thresh_low, thresh_up = Thresh_up),
         dbl_sine = dd_dbl_sine(daily_min = temp_min_f, daily_max = temp_max_f, nextday_min = tmin_next,
                                thresh_low = Thresh_low, thresh_up = Thresh_up)) %>%
  mutate(
    week_id = ceiling(as.numeric(date - date_ref) / 7), # + 1, # calculate Monday - Sunday of the week before
    week_start = date_ref + week_id * 7 # (week_id - 1)
  ) %>%
  filter(week_id < 15, week_id > 0)

Myzus_dd %>%
  group_by(week_id, week_start) %>%
  summarise(count = n())

# explore temperature over time
Myzus_dd %>%
  mutate(week_6_highlight = ifelse(week_id == 6, 1, 0)) %>%
  ggplot(aes(x = date, y = temp_max_f)) +
  geom_line(aes(color = week_6_highlight)) +
  scale_color_viridis_c()

# Grouping by week and Including estimated generation number (EGN) ----
Myzus_dd_week <- Myzus_dd %>%
  filter(week_id < 15) %>%
  group_by(week_start, week_id) %>%
  summarise(
    DD_dbl_tri = sum(dbl_tri, na.rm = TRUE),
    DD_dbl_sine = sum(dbl_sine, na.rm = TRUE)
  )

Myzus_dd_week <- Myzus_dd_week %>%
  mutate(
    EGN_tri = DD_dbl_tri / MGT,
    EGN_sine = DD_dbl_sine / MGT
  )

# explore how degree days and estimated number of generations varies over time.
Myzus_dd_week %>%
  mutate(week_6_highlight = ifelse(week_id == 6, 1, 0)) %>%
  ggplot(aes(x = week_id, y = EGN_sine)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(color = factor(week_6_highlight))) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = 1:14)
