# Folder contents:

- `00_demonstrate_ecoevo_FigS1-6.R` - R script to reproduce analyses that quantify eco-to-evo and evo-to-eco effects using traditional statistical analyses.

- `01_quantify_ecoevo_speed.R` - R script to reproduce analyses that quantify the relative speed of evolutionary vs. ecological dynamics.

- `02_fit_ecoevo_model.R` - R script to fit eco-evolutionary dynamics using a Bayesian first-order multivariate autoregressive model (MAR(1)) with the R package `brms`.

- `ecoevo_dynamics_brms.rds` - Output of Bayesian MAR(1) model, saved to speed-up downstream analyses.

- `03_measure_ecoevo_stability.R` - R script to generate posterior predictions (from `ecoevo_dynamics_brms.rds`) of eco-evolutionary effects and higher-order properties such as stability.

- `posteriors_ecoevo_model.RData` - Output of `measure_ecoevo_stability.R` that is reused in scripts to generate figures.

- `Fig2_plot_ecoevo_stability.R` - R script to generate Figure 2 and associated statistical analyses (e.g., median effect sizes and 95% credible intervals).

- `Fig3_plot_ecoevo_dynamics.R` - R script to generate Figure 3 and associated statistical analyses.

- `Fig4_plot_ecoevo_effects.R` - R script to generate Figure 4 and associated statistical analyses.

- `Fig5_FigS7-8_plot_ecoevo_simulation.R` - R script to generate Figure 5 and associated statistical analyses as well as Figures S1 and S2.
