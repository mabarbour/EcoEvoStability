
## Setup ----

# load libraries
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)

# load posterior predictions
load("analyses/posteriors_ecoevo_model.RData")

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Plot eco-evolutionary stability for aphid and aphid-parasitoid systems (Fig. 2) ----

plot_ecoevo_stability <- bind_rows(select(draws_ecoevo_aphid_df, stability = M_stability_aphid) %>% mutate(type = "Eco-evolutionary\neffects (J)", interaction = "Aphid"),
                                   select(draws_ecoevo_aphid.ptoid_df, stability = M_stability_aphid.ptoid) %>% mutate(type = "Eco-evolutionary\neffects (J)", interaction = "Aphid-Parasitoid"),
                                   select(draws_ecoevo_aphid_df, stability = A_stability_aphid) %>% mutate(type = "Ecological\nfeedback (A)", interaction = "Aphid"),
                                   select(draws_ecoevo_aphid.ptoid_df, stability = A_stability_aphid.ptoid) %>% mutate(type = "Ecological\nfeedback (A)", interaction = "Aphid-Parasitoid"),
                                   select(draws_ecoevo_aphid_df, stability = D_evoecoevo_stability_aphid) %>% mutate(type = "Evolutionary +\nEvo-Eco-Evo feedback\n(D+CA^-1(-B))", interaction = "Aphid"),
                                   select(draws_ecoevo_aphid.ptoid_df, stability = D_evoecoevo_stability_aphid.ptoid) %>% mutate(type = "Evolutionary +\nEvo-Eco-Evo feedback\n(D+CA^-1(-B))", interaction = "Aphid-Parasitoid")) %>%
  ggplot(aes(x = type, y = stability, color = interaction, fill = interaction)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_pointinterval(position = position_dodge(width = 0.25)) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Aphid-Parasitoid")) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)], labels = c("Aphid","Aphid-Parasitoid")) +
  ylab(expression("Stability (max(" * lambda * "))")) +
  xlab("")
plot_ecoevo_stability
save_plot("figures/ecoevo_stability_Fig2.pdf", plot = plot_ecoevo_stability, base_height = 5)

## Calculate differences in stability between aphid and aphid-parasitoid systems ----

# total eco-evolutionary stability
diff_total_stability <- cbind(draws_ecoevo_aphid_df, draws_ecoevo_aphid.ptoid_df) %>%
  mutate(diff_M_stability = M_stability_aphid.ptoid - M_stability_aphid)
round(median(diff_total_stability$diff_M_stability),2) # median difference
1-sum(diff_total_stability$diff_M_stability > 0)/nrow(diff_total_stability) # one-sided probability 

# ecological stability (A matrix)
diff_eco_stability <- cbind(draws_ecoevo_aphid_df, draws_ecoevo_aphid.ptoid_df) %>%
  mutate(diff_eco_stability = A_stability_aphid.ptoid - A_stability_aphid)
round(median(diff_eco_stability$diff_eco_stability),2) # median difference
1-sum(diff_eco_stability$diff_eco_stability > 0)/nrow(diff_eco_stability) # one-sided probability

# eco-evo stability (D + evo-eco-evo feedback)
diff_D_evoecoevo_stability <- cbind(draws_ecoevo_aphid_df, draws_ecoevo_aphid.ptoid_df) %>%
  mutate(diff_D_evoecoevo_stability = D_evoecoevo_stability_aphid.ptoid - D_evoecoevo_stability_aphid)
round(median(diff_D_evoecoevo_stability$diff_D_evoecoevo_stability),2) # median difference
1-sum(diff_D_evoecoevo_stability$diff_D_evoecoevo_stability > 0)/nrow(diff_D_evoecoevo_stability) # one-sided probability

