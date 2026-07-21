
## Setup ----

# load and view data
source('code/manage_data.R')

# data subsets
mono_noPtoid_df <- ecoevo_df_allWeeks_imputed %>% 
  filter(ptoid_in == 0, Treatment %in% c("RP","YP","GP")) %>%
  # converting to factor for gam fitting
  mutate(Treatment = factor(Treatment, levels = c("RP","YP","GP")),
         ID = factor(ID))

mono_yesPtoid_df <- ecoevo_df_allWeeks_imputed %>% 
  filter(ptoid_in == 1, week!= 4, Treatment %in% c("RP","YP","GP")) # same as week > 4

# note uneven sample sizes among treatments after week 11
mono_yesPtoid_df %>%
  select(Treatment, week, At1, Pt1) %>%
  drop_na() %>%
  group_by(Treatment, week) %>%
  summarise(sample_size = n()) %>%
  print(n = 28)

# data for fitting to monoculture treatments
mono_yesPtoid_df_toWeek11 <- mono_yesPtoid_df %>% 
  filter(week < 12) %>%
  # converting to factor for gam fitting
  mutate(Treatment = factor(Treatment, levels = c("RP","YP","GP")),
         ID = factor(ID))

poly_df <- ecoevo_df_allWeeks_imputed %>% 
  # restrict to weeks with the parasitoid
  filter(week > 4, Treatment %in% c("RYG","RYGP")) 

# during this period, there was only 1 extinction, otherwise it is completely balanced
# it also captures the period where aphid morph responses were clearest (in monoculture)
# see code snippet for sample_size below
poly_df %>%
  select(Treatment, week, Rt1) %>%
  drop_na() %>%
  group_by(Treatment, week) %>%
  summarise(sample_size = n())

# data for fitting to polyculture treatments
poly_df_toWeek11 <- poly_df %>%
  filter(week < 12) %>%
  # converting to factor for gam fitting
  mutate(Treatment = factor(Treatment, levels = c("RYG","RYGP")),
         ID = factor(ID))

# additional libraries
library(brms)
library(marginaleffects) # reminder: default is type = "response" and re_formula = NULL
library(cowplot)
library(lmerTest)
library(mgcv)
library(glmmTMB)

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
mono_pal <- c("#D55E00","#F0E442","#009E73")

# function to calculate percent effect size from a negative binomial model
calc_percent_effect_nb <- function(x) {(exp(x) - 1)*100}

## Monoculture: Test the hypothesis that aphids morphs vary in population size in the absence of parasitoids ----
# note data is log instead of log1p transformed here because it doesn't need to be

# GAM
gam_nb_aphid_noPtoid <- gam(At1 ~ Treatment + 
                      s(week, by = Treatment, k = 3) + 
                      s(week, ID, bs = "fs", k = 3), 
                    family = nb(theta = NULL, link = "log"),
                    data = mono_noPtoid_df,
                    method = "REML")
gam.check(gam_nb_aphid_noPtoid) # p-values not great
gratia::appraise(gam_nb_aphid_noPtoid) # residuals seem okay
gratia::draw(gam_nb_aphid_noPtoid)
summary(gam_nb_aphid_noPtoid)

# cage-level predictions
plot_predictions(gam_nb_aphid_noPtoid, by = c("week","Treatment","ID"),
                 type = "link") + # remove observation-level random effect
  geom_point(data = mono_noPtoid_df, aes(x = week, y = log(At1), color = Treatment), alpha = 0.3) +
  scale_color_manual(values = mono_pal) +
  scale_fill_manual(values = mono_pal)

# treatment-level effects
Fig_S1_mono_At1_noPtoid <- predictions(gam_nb_aphid_noPtoid,
            variables = "Treatment", 
            by = c("week","Treatment"),
            type = "link",
            exclude = "s(week,ID)") %>%
  ggplot(aes(x = week, y = exp(estimate))) +
  geom_line(aes(color = Treatment), position = position_dodge(width = 0.4), alpha = 0.5) +
  geom_pointrange(aes(ymax = exp(conf.high), ymin = exp(conf.low), color = Treatment, fill = Treatment),
                  position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymax = exp(estimate + std.error), ymin = exp(estimate - std.error), color = Treatment, fill = Treatment),
                  position = position_dodge(width = 0.4), linewidth = 1.5) +
  scale_color_manual(values = mono_pal) + 
  scale_fill_manual(values = mono_pal) +
  scale_x_continuous("Week of experiment", breaks = 1:3) +
  ggtitle("No Parasitoids") +
  scale_y_continuous(name = "Aphid abundance", transform = "log", breaks = c(1,10,100,1000,2000,4000,8000,10000))
Fig_S1_mono_At1_noPtoid
save_plot(file = "figures/Fig_S1_mono_At1_noPtoid.pdf", plot = Fig_S1_mono_At1_noPtoid)

plot_comparisons(gam_nb_aphid_noPtoid, 
                 variables = "Treatment", by = "week", 
                 exclude = "s(week,ID)", type = "link", transform = calc_percent_effect_nb) +
  geom_hline(yintercept = 0, linetype = "dotted")
avg_comparisons(gam_nb_aphid_noPtoid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link", transform = calc_percent_effect_nb)
avg_comparisons(gam_nb_aphid_noPtoid, variables = "Treatment", exclude = "s(week,ID)", type = "link", transform = calc_percent_effect_nb)
avg_comparisons(gam_nb_aphid_noPtoid, variables = "Treatment", 
                newdata = filter(mono_noPtoid_df, week < 3), 
                exclude = "s(week,ID)", type = "link", transform = calc_percent_effect_nb)
# trend for higher yellow morph abundance in absence of parasitoids across weeks 1 and 2,
# but differences disappear at week 3, likely due to density dependence.

# test robustness of gam (gam.check p-values) by fitting a linear model with glmmTMB
glmmTMB_nb_aphid_noPtoid <- glmmTMB(At1 ~ (week + I(week^2))*Treatment + (week + I(week^2)||ID),
                            family = nbinom1(link = "log"),
                            data = mono_noPtoid_df)
summary(glmmTMB_nb_aphid_noPtoid)
plot_predictions(glmmTMB_nb_aphid_noPtoid, 
                 by = c("week","Treatment"),
                 vcov = TRUE,
                 type = "link",
                 re.form = NA)
plot_comparisons(glmmTMB_nb_aphid_noPtoid, 
                 variables = "Treatment",
                 by = c("week"),
                 vcov = TRUE,
                 type = "link",
                 re.form = NA) +
  geom_hline(yintercept = 0, linetype = "dashed")
# similar inference, but perhaps an increase at week 2


## Monoculture: Test the hypothesis that aphids morphs vary in population size in the presence of parasitoids ----

# GAM
gam_nb_aphid <- gam(At1 ~ Treatment + 
                  s(week, by = Treatment, k = 6) + 
                  s(week, ID, bs = "fs", k = 6), 
                family = nb(theta = NULL, link = "log"),
                data = mono_yesPtoid_df_toWeek11, 
                method = "REML")
gam.check(gam_nb_aphid) # p-values okay
gratia::appraise(gam_nb_aphid) # residuals seem okay
summary(gam_nb_aphid)

# cage-level predictions
plot_predictions(gam_nb_aphid, by = c("week","Treatment","ID"),
                 type = "link") + # remove observation-level random effect
  geom_point(data = mono_yesPtoid_df_toWeek11, aes(x = week, y = log(At1), color = Treatment), alpha = 0.3) +
  scale_color_manual(values = mono_pal) +
  scale_fill_manual(values = mono_pal)

# treatment-level effects
Fig_S2_mono_At1_withPtoid <- predictions(gam_nb_aphid,
            variables = "Treatment", 
            by = c("week","Treatment"),
            type = "link",
            exclude = "s(week,ID)") %>%
  ggplot(aes(x = week, y = exp(estimate))) +
  geom_line(aes(color = Treatment), position = position_dodge(width = 0.4), alpha = 0.5) +
  geom_pointrange(aes(ymax = exp(conf.high), ymin = exp(conf.low), color = Treatment, fill = Treatment),
                  position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymax = exp(estimate + std.error), ymin = exp(estimate - std.error), color = Treatment, fill = Treatment),
                  position = position_dodge(width = 0.4), linewidth = 1.5) +
  scale_color_manual(values = mono_pal) + 
  scale_fill_manual(values = mono_pal) +
  scale_x_continuous("Week of experiment", breaks = 5:11) +
  ggtitle("With Parasitoids") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(name = "Aphid abundance", transform = "log", breaks = c(1,10,100,1000,10000))
Fig_S2_mono_At1_withPtoid
save_plot(file = "figures/Fig_S2_mono_At1_withPtoid.pdf", plot = Fig_S2_mono_At1_withPtoid)

plot_comparisons(gam_nb_aphid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link") +
  geom_hline(yintercept = 0, linetype = "dotted")
avg_comparisons(gam_nb_aphid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link")
avg_comparisons(gam_nb_aphid, variables = "Treatment", exclude = "s(week,ID)", type = "link", comparison = "difference")
avg_comparisons(gam_nb_aphid, variables = "Treatment", exclude = "s(week,ID)", type = "link", comparison = "difference", transform = exp)
# clear evidence that yellow morph is more abundant in the presence of parasitoids than the red morph
# this pattern is particularly pronounced in weeks 10 and 11 when relative parasitoid abundances were high. 
# yellow morph also had transiently higher abundance at week 6, when parasitoid abundances were low...

plot_predictions(gam_nb_aphid, variables = "Treatment", by = "Treatment", type = "link") # shows counterfactual predictions on link scale
avg_predictions(gam_nb_aphid, variables = "Treatment", by = "Treatment", type = "link")
exp(7.52)/exp(5.15) # reflects fold differences in abundances, which is the same as taking the exponent of the average treatment effect on the link (log) scale

## Monoculture: Test the hypothesis that parasitoid population size varies among aphid morphs ----

# GAM
gam_nb_ptoid <- gam(Pt1 ~ Treatment + 
                      # chose k based on fit to cage-level predictions
                      s(week, by = Treatment, k = 6) + 
                      s(week, ID, bs = "fs", k = 6), 
                    family = nb(theta = NULL, link = "log"),
                    data = mono_yesPtoid_df_toWeek11, # %>% filter(week > 5),
                    method = "REML")
gam.check(gam_nb_ptoid) # p-values off, but cage-level prediction look okay. p-values improves when filtering week > 5, but qualitative inference remains the same
gratia::appraise(gam_nb_ptoid) # residuals seem okay
summary(gam_nb_ptoid)

# cage-level predictions
plot_predictions(gam_nb_ptoid, by = c("week","Treatment","ID"),
                 type = "link") + # remove observation-level random effect
  geom_point(data = mono_yesPtoid_df_toWeek11, aes(x = week, y = log(Pt1), color = Treatment), alpha = 0.3) +
  scale_color_manual(values = mono_pal) +
  scale_fill_manual(values = mono_pal)

# treatment-level effects
Fig_S3_mono_Pt1 <- predictions(gam_nb_ptoid,
            variables = "Treatment", 
            by = c("week","Treatment"),
            type = "link",
            exclude = "s(week,ID)") %>%
  ggplot(aes(x = week, y = exp(estimate))) +
  geom_line(aes(color = Treatment), position = position_dodge(width = 0.4), alpha = 0.5) +
  geom_pointrange(aes(ymax = exp(conf.high), ymin = exp(conf.low), color = Treatment, fill = Treatment),
                  position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymax = exp(estimate + std.error), ymin = exp(estimate - std.error), color = Treatment, fill = Treatment),
                  position = position_dodge(width = 0.4), linewidth = 1.5) +
  scale_color_manual(values = mono_pal) + 
  scale_fill_manual(values = mono_pal) +
  scale_x_continuous("Week of experiment", breaks = 5:11) +
  #ggtitle("With Parasitoids") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(name = "Parasitoid abundance", transform = "log", breaks = c(1,10,100,1000,10000))
Fig_S3_mono_Pt1
save_plot(file = "figures/Fig_S3_mono_Pt1.pdf", plot = Fig_S3_mono_Pt1)

plot_comparisons(gam_nb_ptoid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link") +
  geom_hline(yintercept = 0, linetype = "dotted")
avg_comparisons(gam_nb_ptoid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link")
avg_comparisons(gam_nb_ptoid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link", 
                comparison = "difference", transform = exp)
avg_comparisons(gam_nb_ptoid, variables = "Treatment", exclude = "s(week,ID)", type = "link")
# no general difference across morphs, but the pattern was more dynamic
# parasitoid abundance on yellow morphs was lower than red morph at week 7 (despite greater aphid abundance at week 6),
# but parasitoid abundance was higher on yellow morphs at week 10 (trend) and 11 (clear)

# for week 7 we will reverse the comparison RP vs YP by taking the absolute value
avg_comparisons(gam_nb_ptoid, variables = "Treatment", by = "week", exclude = "s(week,ID)", type = "link") %>%
  filter(week == "7") %>%
  mutate(estimate = exp(abs(estimate)),
         conf.low = exp(abs(conf.low)),
         conf.high = exp(abs(conf.high)))

## Polyculture: Test the hypothesis that red morph frequency varies between parasitoid treatments ----

# Fit and check GAM
poly_Rt1_gam <- gam(Rt1 ~ Treatment + 
                      s(week, by = Treatment, k = 4) + 
                      s(week, ID, bs = "fs", k = 4),
                    family = betar(link = "logit"),
                    data = poly_df_toWeek11 %>% 
                      mutate(# nudged 2 data points to fit betabinomial assumptions
                             Rt1 = ifelse(Rt1 == 1, 0.99, Rt1)), 
                    method = "REML")
gam.check(poly_Rt1_gam) # p-values okay
gratia::appraise(poly_Rt1_gam) # residuals seem okay
gratia::draw(poly_Rt1_gam)
summary(poly_Rt1_gam)

# cage-level predictions
plot_predictions(poly_Rt1_gam, by = c("week","Treatment","ID"),
                 type = "response") + 
  geom_point(data = poly_df_toWeek11, aes(x = week, y = Rt1, color = Treatment), alpha = 0.3) +
  scale_color_manual(values = c("grey","black")) +
  scale_fill_manual(values = c("grey","black"))

# treatment effect for the average cage
Fig_S5_poly_Rt1_freq <- plot_predictions(poly_Rt1_gam, by = c("week","Treatment"), exclude = c("s(week,ID)","s(week_f,ID)")) +
  scale_color_manual(values = c("black","#0072B2"), labels = c("No parasitoid","With parasitoid")) + 
  scale_fill_manual(values = c("black","#0072B2"), labels = c("No parasitoid","With parasitoid")) +
  scale_x_continuous(name = "Week of experiment", breaks = 5:11) +
  scale_y_continuous(name = "Red morph frequency", lim = c(0,1)) +
  geom_hline(yintercept = 0.33, linetype = "dashed")
Fig_S5_poly_Rt1_freq
save_plot(file = "figures/Fig_S5_poly_Rt1_freq.pdf", plot = Fig_S5_poly_Rt1_freq)

plot_comparisons(poly_Rt1_gam, variables = "Treatment", by = "week", exclude = c("s(week,ID)")) +
  geom_hline(yintercept = 0, linetype = "dotted")
avg_comparisons(poly_Rt1_gam, variables = "Treatment", by = "week", exclude = c("s(week,ID)"))
avg_comparisons(poly_Rt1_gam, variables = "Treatment", exclude = c("s(week,ID)"))
# clear evidence of positive selection on red morph frequency by end of experiment at week 11


## Polyculture: Test the hypothesis that yellow morph frequency varies between parasitoid treatments ----

# Fit and check GAM
poly_Yt1_gam <- gam(Yt1 ~ Treatment + 
                  s(week, by = Treatment, k = 7) + 
                  s(week, ID, bs = "fs", k = 7),
                  family = betar(link = "logit"),
                data = poly_df_toWeek11 %>% 
                  mutate(# nudged 2 data points to fit betabinomial assumptions
                         Yt1 = ifelse(Yt1 == 0, 0.01, Yt1)), 
                method = "REML")
gam.check(poly_Yt1_gam) # p-values okay with k=7, worse with k=6
gratia::appraise(poly_Yt1_gam) # residuals look okay
gratia::draw(poly_Yt1_gam)
summary(poly_Yt1_gam)

# cage-level predictions
plot_predictions(poly_Yt1_gam, by = c("week","Treatment","ID"),
                 type = "response") +
  geom_point(data = poly_df_toWeek11, aes(x = week, y = Yt1, color = Treatment), alpha = 0.3) +
  scale_color_manual(values = c("black","#0072B2")) +
  scale_fill_manual(values = c("black","#0072B2"))

# treatment effect for the average cage
Fig_S4_poly_Yt1_freq <- plot_predictions(poly_Yt1_gam, by = c("week","Treatment"), type = "response", exclude = c("s(week,ID)","s(week_f,ID)")) +
  scale_color_manual(values = c("black","#0072B2"), labels = c("No parasitoid","With parasitoid")) + 
  scale_fill_manual(values = c("black","#0072B2"), labels = c("No parasitoid","With parasitoid")) +
  scale_x_continuous(name = "Week of experiment", breaks = 5:11) +
  scale_y_continuous(name = "Yellow morph frequency", lim = c(0,1)) +
  geom_hline(yintercept = 0.33, linetype = "dashed")
Fig_S4_poly_Yt1_freq
save_plot(file = "figures/Fig_S4_poly_Yt1_freq.pdf", plot = Fig_S4_poly_Yt1_freq)

plot_comparisons(poly_Yt1_gam, variables = "Treatment", by = "week", exclude = c("s(week,ID)")) +
  geom_hline(yintercept = 0, linetype = "dotted")
avg_comparisons(poly_Yt1_gam, variables = "Treatment", by = "week", exclude = c("s(week,ID)"))
avg_comparisons(poly_Yt1_gam, variables = "Treatment", exclude = c("s(week,ID)","s(week_f,ID)"))
# clear evidence of selection to reduce yellow morph frequency by week 8
avg_comparisons(poly_Yt1_gam, variables = "Treatment", newdata = filter(poly_df_toWeek11, week > 7), exclude = c("s(week,ID)"))


# test robustness of gam (gam.check p-values) by fitting a linear model with glmmTMB
poly_Yt1_glmmTMB <- glmmTMB(cbind(Y_At1, At1-Y_At1) ~ week*Treatment + (week||ID),
                            family = betabinomial(link = "logit"),
                            data = poly_df_toWeek11 %>% 
                              mutate(Yt1 = ifelse(Yt1 == 0, 0.01, Yt1)))
summary(poly_Yt1_glmmTMB)
plot_predictions(poly_Yt1_glmmTMB, 
                 by = c("week","Treatment"),
                 vcov = TRUE,
                 type = "response",
                 re.form = NA)
# same qualitative answer, so we report the gam

## Polyculture: Test the hypothesis that green morph frequency varies between parasitoid treatments ----

# Fit and check GAM
poly_Gt1_gam <- gam(Gt1 ~ Treatment + 
                      s(week, by = Treatment, k = 6) + 
                      s(week, ID, bs = "fs", k = 6),
                    family = betar(link = "logit"),
                    data = poly_df_toWeek11 %>% 
                      mutate(# nudged 5 data points to fit betabinomial assumptions
                             Gt1 = ifelse(Gt1 == 0, 0.01, Gt1)), 
                    method = "REML")
gam.check(poly_Gt1_gam) # p-values okay
gratia::appraise(poly_Gt1_gam) # residuals seem okay
gratia::draw(poly_Gt1_gam)
summary(poly_Gt1_gam)

# cage-level predictions
plot_predictions(poly_Gt1_gam, by = c("week","Treatment","ID"),
                 type = "response") +
  geom_point(data = poly_df_toWeek11, aes(x = week, y = Gt1, color = Treatment), alpha = 0.3) +
  scale_color_manual(values = c("grey","black")) +
  scale_fill_manual(values = c("grey","black"))

# treatment effect for the average cage
Fig_S6_poly_Gt1_freq <- plot_predictions(poly_Gt1_gam, by = c("week","Treatment"), exclude = c("s(week,ID)","s(week_f,ID)")) +
  scale_color_manual(values = c("black","#0072B2"), labels = c("No parasitoid","With parasitoid")) + 
  scale_fill_manual(values = c("black","#0072B2"), labels = c("No parasitoid","With parasitoid")) +
  scale_x_continuous(name = "Week of experiment", breaks = 5:11) +
  scale_y_continuous(name = "Green morph frequency", lim = c(0,1)) +
  geom_hline(yintercept = 0.33, linetype = "dashed")
Fig_S6_poly_Gt1_freq
save_plot(file = "figures/Fig_S6_poly_Gt1_freq.pdf", plot = Fig_S6_poly_Gt1_freq)


plot_comparisons(poly_Gt1_gam, variables = "Treatment", by = "week", exclude = c("s(week,ID)")) +
  geom_hline(yintercept = 0, linetype = "dotted")
avg_comparisons(poly_Gt1_gam, variables = "Treatment", by = "week", exclude = c("s(week,ID)"))
avg_comparisons(poly_Gt1_gam, variables = "Treatment", exclude = c("s(week,ID)"))
# clear evidence of selection to decrease green morph frequency by the end of the experiment,
# but there was also a transient trend for an increase at week 9
