
## Setup ----

# load libraries
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)
library(patchwork)

# load posterior predictions
load("analyses/posteriors_ecoevo_model.RData")

# set plot theme
theme_set(theme_cowplot())

# color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## Plot eco-evolutionary effects for aphid and aphid-parasitoid systems (Fig. 4) ----

# A-matrix (eco-to-eco effects)
A_mat_data <- bind_rows(
  plyr::ldply(M_aphid, .fun = function(x) x["lnA","lnA"]) %>% mutate(interaction = "Aphid", resp = "ln(N[A])", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnA","lnA"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[A])", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnA","lnP"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[A])", effect = "ln(N[P])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnP","lnA"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[P])", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnP","lnP"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[P])", effect = "ln(N[P])", draw = 1:4000))
A_mat_plot <- A_mat_data %>%
  ggplot(aes(x = 0, y = V1, color = interaction, fill = interaction)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_pointinterval(position = position_dodge(width = 0.25)) + # , show.legend = F
  coord_flip() +
  facet_grid(resp ~ effect, switch = "y", labeller = label_parsed) +
  ylab(expression(partialdiff*ln(N[i])/partialdiff*ln(N[j]))) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)]) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)]) + 
  ggtitle("Ecological feedback (A)") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# B-matrix (evo-to-eco effects)
B_mat_data <- bind_rows(
  plyr::ldply(M_aphid, .fun = function(x) x["lnA","R"]) %>% mutate(interaction = "Aphid", resp = "ln(N[A])", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid, .fun = function(x) x["lnA","Y"]) %>% mutate(interaction = "Aphid", resp = "ln(N[A])", effect = "z[Y]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnA","R"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[A])", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnA","Y"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[A])", effect = "z[Y]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnP","R"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[P])", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["lnP","Y"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "ln(N[P])", effect = "z[Y]", draw = 1:4000))
B_mat_plot <- B_mat_data %>%
  ggplot(aes(x = 0, y = V1, color = interaction, fill = interaction)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_pointinterval(position = position_dodge(width = 0.25)) + # , show.legend = F
  coord_flip() +
  facet_grid(resp ~ effect, switch = "y", labeller = label_parsed) +
  ylab(expression(partialdiff*ln(N[j])/partialdiff*z[i])) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)]) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)]) + 
  ggtitle("Evo-to-eco effect (B)") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# C-matrix (eco-to-evo effects)
C_mat_data <- bind_rows(
  plyr::ldply(M_aphid, .fun = function(x) x["R","lnA"]) %>% mutate(interaction = "Aphid", resp = "z[R]", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid, .fun = function(x) x["Y","lnA"]) %>% mutate(interaction = "Aphid", resp = "z[Y]", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["R","lnA"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[R]", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["R","lnP"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[R]", effect = "ln(N[P])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["Y","lnA"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[Y]", effect = "ln(N[A])", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["Y","lnP"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[Y]", effect = "ln(N[P])", draw = 1:4000))
C_mat_plot <- C_mat_data %>%
  ggplot(aes(x = 0, y = V1, color = interaction, fill = interaction)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_pointinterval(position = position_dodge(width = 0.25)) + # , show.legend = F
  coord_flip() +
  facet_grid(resp ~ effect, switch = "y", labeller = label_parsed) +
  ylab(expression(partialdiff*z[i]/partialdiff*ln(N[j]))) +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)]) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)]) + 
  scale_y_continuous(breaks = c(-0.05,0,0.05)) +
  ggtitle("Eco-to-evo effect (C)") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# D-matrix (evo-to-evo effects)
D_mat_data <- bind_rows(
  plyr::ldply(M_aphid, .fun = function(x) x["R","R"]) %>% mutate(interaction = "Aphid", resp = "z[R]", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid, .fun = function(x) x["R","Y"]) %>% mutate(interaction = "Aphid", resp = "z[R]", effect = "z[Y]", draw = 1:4000),
  plyr::ldply(M_aphid, .fun = function(x) x["Y","R"]) %>% mutate(interaction = "Aphid", resp = "z[Y]", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid, .fun = function(x) x["Y","Y"]) %>% mutate(interaction = "Aphid", resp = "z[Y]", effect = "z[Y]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["R","R"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[R]", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["R","Y"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[R]", effect = "z[Y]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["Y","R"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[Y]", effect = "z[R]", draw = 1:4000),
  plyr::ldply(M_aphid.ptoid, .fun = function(x) x["Y","Y"]) %>% mutate(interaction = "Aphid-Parasitoid", resp = "z[Y]", effect = "z[Y]", draw = 1:4000))
D_mat_plot <- D_mat_data %>%
  ggplot(aes(x = 0, y = V1, color = interaction, fill = interaction)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_pointinterval(position = position_dodge(width = 0.25)) + 
  #stat_slabinterval(slab_alpha = 0.2, position = position_dodge(width = 0.25)) + # , show.legend = F
  coord_flip() +
  facet_grid(resp ~ effect, switch = "y", labeller = label_parsed) +
  ylab(expression(partialdiff*z[i]/partialdiff*z[j])) +
  xlab("") +
  scale_color_manual(name = NULL, values = cbbPalette[c(1,6)]) +
  scale_fill_manual(name = NULL, values = cbbPalette[c(1,6)]) +
  ggtitle("Evolutionary feedback (D)") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# put together
plot_ecoevo_effects <- (A_mat_plot + B_mat_plot + C_mat_plot + D_mat_plot) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
plot_ecoevo_effects
save_plot("figures/ecoevo_effects_Fig4.pdf", plot = plot_ecoevo_effects, base_height = 6)

## Calculate differences in effects (parameters) between aphid and aphid-parasitoid systems ----

# direct ecological effects
A_mat_data %>%
  select(-draw) %>%
  unite(col = "term", resp:effect) %>%
  group_by(term, interaction) %>%
  median_qi() # same as reported in figure

A_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid', term == "ln(N[A])_ln(N[A])") %>%
  # negative frequency dependence in aphid only system
  summarise(median = median(V1),
            p_onetail = 1 - sum(V1 < 0)/length(V1))

# evo-to-eco effects 
B_mat_data %>%
  select(-draw) %>%
  unite(col = "term", resp:effect) %>%
  group_by(term, interaction) %>%
  median_qi()  # same as reported in figure

B_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid') %>%
  group_by(term) %>%
  # no evidence of evo-to-eco effects in aphid-only system
  summarise(median = median(V1),
            p_onetail = 1 - sum(V1 > 0)/length(V1))

B_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid-Parasitoid', term == "ln(N[A])_z[Y]") %>%
  summarise(median = median(V1),
            p_onetail = 1 - sum(V1 < 0)/length(V1))

B_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid-Parasitoid', term %in% c("ln(N[P])_z[R]", "ln(N[P])_z[Y]")) %>%
  pivot_wider(names_from = "term", values_from = "V1") %>%
  mutate(diff = `ln(N[P])_z[Y]` - `ln(N[P])_z[R]`) %>%
  select(diff) %>%
  # difference in resistance to parasitoid between yellow and red morph 
  summarise(diff_median = median(diff),
            p_onetail = 1 - sum(diff < 0)/length(diff))

# eco-to-evo effects
C_mat_data %>%
  select(-draw) %>%
  unite(col = "term", resp:effect) %>%
  group_by(term, interaction) %>%
  median_qi() # same as reported in figure

C_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid') %>%
  group_by(term) %>%
  # no evidence of eco-to-evo effects in aphid-only system
  summarise(median = median(V1),
            p_onetail = 1 - sum(V1 < 0)/length(V1))

# direct evolutionary effects
D_mat_data %>%
  select(-draw) %>%
  unite(col = "term", resp:effect) %>%
  group_by(term, interaction) %>%
  median_qi() # same as reported in figure

D_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid', term == "z[R]_z[R]") %>%
  # evidence for negative frequency dependent selection on red morph in aphid-only system
  summarise(median = median(V1),
            p_onetail = 1 - sum(V1 < 0)/length(V1))

D_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(interaction == 'Aphid', term == "z[Y]_z[Y]") %>%
  # evidence for negative frequency-dependent selection in aphid-only system
  summarise(median = median(V1),
            p_onetail = 1 - sum(V1 < 0)/length(V1))

D_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(term == "z[R]_z[R]") %>%
  pivot_wider(names_from = "interaction", values_from = "V1") %>%
  mutate(diff = `Aphid-Parasitoid` - `Aphid`) %>%
  select(diff) %>%
  # change in frequency dependent selection on red morph between aphid and aphid-parasitoid systems
  summarise(diff_median = median(diff),
            p_onetail = 1 - sum(diff > 0)/length(diff))

D_mat_data %>%
  unite(col = "term", resp:effect) %>%
  filter(term == "z[Y]_z[Y]") %>%
  pivot_wider(names_from = "interaction", values_from = "V1") %>%
  mutate(diff = `Aphid-Parasitoid` - `Aphid`) %>%
  select(diff) %>%
  # change in frequency dependent selection on yellow morph between aphid and aphid-parasitoid systems
  summarise(diff_median = median(diff),
            p_onetail = 1 - sum(diff < 0)/length(diff))

