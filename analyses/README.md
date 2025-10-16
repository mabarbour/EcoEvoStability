# Folder contents:

-   `quantify_ecoevo_speed.R` - R script to reproduce analyses that quantify the relative speed of evolutionary vs. ecological dynamics.

-   `fit_ecoevo_model.R` - R script to fit eco-evolutionary dynamics using a Bayesian first-order multivariate autoregressive model (MAR(1)) with the R package `brms`.

-   `ecoevo_dynamics_brms.rds` - Output of Bayesian MAR(1) model, saved to speed-up downstream analyses.

-   `measure_ecoevo_stability.R` - R script to generate posterior predictions (from `ecoevo_dynamics_brms.rds`) of eco-evolutionary effects and higher-order properties such as stability.

-   `posteriors_ecoevo_model.RData` - Output of `measure_ecoevo_stability.R` that is reused in scripts to generate figures.

-   `plot_ecoevo_stability_Fig2.R` - R script to generate Figure 2 and associated statistical analyses (e.g. median effect sizes, 95% credible intervals, and one-sided probability of effect direction).

-   `plot_ecoevo_dynamics_Fig3.R` - R script to generate Figure 3 and associated statistical analyses.

-   `plot_ecoevo_effects_Fig4.R` - R script to generate Figure 4 and associated statistical analyses.

-   `plot_ecoevo_simulation_Fig5_FigS1-2.R` - R script to generate Figure 5 and associated statistical analyses as well as Figures S1 and S2.
