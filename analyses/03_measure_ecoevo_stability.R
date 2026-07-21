## Setup ----

# start a new session, see https://stackoverflow.com/a/66127391/2554330 for library(matlib)
options(rgl.useNULL = TRUE)
library(rgl)
library(matlib) # for inverse

# load and view data
source('code/manage_data.R')

# load other required libraries
library(brms)

# load bayesian MAR1 model of eco-evo dynamics
ecoevo_dynamics_brm <- readRDS("analyses/ecoevo_dynamics_brm.rds")

# get all posterior draws for model parameters
draws <- as_draws_df(ecoevo_dynamics_brm, variable = "^b", regex = T) %>%
  rename(# simplify names of some aphid terms
          b_lnAt1_lnAt = bsp_lnAt1_milnAtidxEQind,
         `b_lnAt1_lnAt:ptoid_in` = `bsp_lnAt1_milnAtidxEQind:ptoid_in`,
         `b_lnAt1_ptoid_in:lnPt` = `bsp_lnAt1_ptoid_in:milnPtidxEQind`,
         b_lnAt1_Rt = bsp_lnAt1_miRtidxEQind,
         `b_lnAt1_Rt:ptoid_in` = `bsp_lnAt1_miRtidxEQind:ptoid_in`,
         b_lnAt1_Yt = bsp_lnAt1_miYtidxEQind,
         `b_lnAt1_Yt:ptoid_in` = `bsp_lnAt1_miYtidxEQind:ptoid_in`,
         # simplify names of some parasitoid terms
         b_lnPt1_lnAt = bsp_lnPt1_milnAtidxEQind,
         b_lnPt1_lnPt = bsp_lnPt1_milnPtidxEQind,
         b_lnPt1_Rt = bsp_lnPt1_miRtidxEQind,
         b_lnPt1_Yt = bsp_lnPt1_miYtidxEQind,
         # simplify names of some red morph terms
         b_Rt1_lnAt = bsp_Rt1_milnAtidxEQind,
         `b_Rt1_lnAt:ptoid_in` = `bsp_Rt1_milnAtidxEQind:ptoid_in`, 
         `b_Rt1_ptoid_in:lnPt` = `bsp_Rt1_ptoid_in:milnPtidxEQind`,
         b_Rt1_Rt = bsp_Rt1_miRtidxEQind,
         `b_Rt1_Rt:ptoid_in` = `bsp_Rt1_miRtidxEQind:ptoid_in`,
         b_Rt1_Yt = bsp_Rt1_miYtidxEQind,
         `b_Rt1_Yt:ptoid_in` = `bsp_Rt1_miYtidxEQind:ptoid_in`,
         # simplify names of some yellow morph terms
         b_Yt1_lnAt = bsp_Yt1_milnAtidxEQind,
         `b_Yt1_lnAt:ptoid_in` = `bsp_Yt1_milnAtidxEQind:ptoid_in`,
         `b_Yt1_ptoid_in:lnPt` = `bsp_Yt1_ptoid_in:milnPtidxEQind`,
         b_Yt1_Rt = bsp_Yt1_miRtidxEQind,
         `b_Yt1_Rt:ptoid_in` = `bsp_Yt1_miRtidxEQind:ptoid_in`,
         b_Yt1_Yt = bsp_Yt1_miYtidxEQind,
         `b_Yt1_Yt:ptoid_in` = `bsp_Yt1_miYtidxEQind:ptoid_in`)
  
## Get aphid eco-evo parameters (median values) ----
# these median values are useful for getting an overall picture of the dynamics
# and used for subsequent simulations

# vector of aphid intrinsic growth rates and baseline selection on morph frequencies (intercept terms)
# note that this combines vectors 'r' and 'x' from the manuscript.
x_aphid_median <- matrix(
  data = c(median(draws$b_lnAt1_Intercept), 
           median(draws$b_Rt1_Intercept),
           median(draws$b_Yt1_Intercept)),
  nrow = 3, ncol = 1, byrow = T,
  dimnames = list(c("lnA","R","Y"),c("x")))

# matrix of all eco-evolutionary effects
M_aphid_median <- matrix(
  # no offset along diagonals
  data = c(median(draws$b_lnAt1_lnAt), median(draws$b_lnAt1_Rt), median(draws$b_lnAt1_Yt),
           median(draws$b_Rt1_lnAt), median(draws$b_Rt1_Rt), median(draws$b_Rt1_Yt),
           median(draws$b_Yt1_lnAt), median(draws$b_Yt1_Rt), median(draws$b_Yt1_Yt)), 
  nrow = 3, ncol = 3, byrow = T,
  dimnames = list(c("lnA","R","Y"),c("lnA","R","Y")))

# identity matrix for aphid-only system
I_aphid <- diag(x = 1, nrow = nrow(M_aphid_median), ncol = ncol(M_aphid_median))

# subset matrix by ecological, eco-to-evo, evo-to-eco, and evolutionary effects
A_aphid_median <- M_aphid_median["lnA","lnA"]
B_aphid_median <- M_aphid_median["lnA",c("R","Y")]
C_aphid_median <- M_aphid_median[c("R","Y"),"lnA"]
D_aphid_median <- M_aphid_median[c("R","Y"),c("R","Y")]

# calculate evo-eco-evo feedback loop. 
# Note different formulas for discrete-time model!
inv_A_aphid_median <- 1/(1 - A_aphid_median) # note just simple multiplication because of the one value
evoecoevo_aphid_median <- C_aphid_median * inv_A_aphid_median * B_aphid_median # evo-eco-evo feedback loop
evoecoevo_aphid_median

# equilibrium values for aphid abundance and morph frequencies
# note that since I used an offset in the regression, this is the correct formulation
# since I added 1 back to the diagonals after modelling the offset, 
# I need to subtract the identity matrix first (I-B)^-1, see Ives et al. 2003, Ecol. Monog.
eq_aphid <- inv(I_aphid - M_aphid_median) %*% x_aphid_median
eq_aphid
1 - eq_aphid[2,] - eq_aphid[3, ] # green morph

# ecological equilibrium in absence of evolution
inv_A_aphid_median * x_aphid_median[c("lnA"),]
# lnA = 8.83
# exp(8.83) = 6,835
# this suggests the evolutionary dynamics could have an appreciable effect on aphid abundances

# evolutionary equilibrium in absence of ecology
inv(I_aphid[c(2,3),c(2,3)] - D_aphid_median) %*% x_aphid_median[c("R","Y"),]
# R = 0.40
# Y = 0.37
# this suggests the ecological dynamics could have an appreciable effect on morph frequencies...

# calculate stability of relevant matrices
max(Re(eigen(M_aphid_median)$values)) # eco-evo stability
max(Re(eigen(A_aphid_median)$values)) # eco stability
max(Re(eigen(D_aphid_median + evoecoevo_aphid_median)$values)) # evo + evoecoevo feedback stability

# exploratory: may be possible to further partition contribution of evo-eco-evo feedback and direct evolutionary effects
max(Re(eigen(D_aphid_median)$values)) # evo stability
max(Re(eigen(D_aphid_median + evoecoevo_aphid_median)$values)) - max(Re(eigen(D_aphid_median)$values))
# difference between D and D + evoecoevo feedback should gives insight to the unique role of the feedback on stability


## Get aphid-parasitoid eco-evo parameters (median values) ----
# these median values are useful for getting an overall picture of the dynamics
# and used for subsequent simulations

# vector of aphid and parasitoid intrinsic growth rates and baseline selection on morph frequencies (intercept terms)
# note that this combines vectors 'r' and 'x' from the manuscript.
x_aphid.ptoid_median <- matrix(
  # adding effect of 'ptoid_in' accounts for the general effect of adding the parasitoid on parameter estimates
  data = c(median(draws$b_lnAt1_Intercept + draws$b_lnAt1_ptoid_in),
           median(draws$b_lnPt1_Intercept),
           median(draws$b_Rt1_Intercept + draws$b_Rt1_ptoid_in),
           median(draws$b_Yt1_Intercept + draws$b_Yt1_ptoid_in)),
  nrow = 4, ncol = 1, byrow = T,
  dimnames = list(c("lnA","lnP","R","Y"),c("x")))
x_aphid.ptoid_median

# matrix of all eco-evolutionary effects
M_aphid.ptoid_median <- matrix(
  data = c(
    # note that we add relevant coefficients that are modified by the parasitoid
    # no offset along diagonal
    # row 1
    median(draws$b_lnAt1_lnAt + draws$`b_lnAt1_lnAt:ptoid_in`), 
    median(draws$`b_lnAt1_ptoid_in:lnPt`), 
    median(draws$b_lnAt1_Rt + draws$`b_lnAt1_Rt:ptoid_in`), 
    median(draws$b_lnAt1_Yt + draws$`b_lnAt1_Yt:ptoid_in`),
    # row 2, no need to add coeficients since these were only estimated with 'ptoid_in'
    median(draws$b_lnPt1_lnAt), median(draws$b_lnPt1_lnPt), median(draws$b_lnPt1_Rt), median(draws$b_lnPt1_Yt),
    # row 3
    median(draws$b_Rt1_lnAt + draws$`b_Rt1_lnAt:ptoid_in`), 
    median(draws$`b_Rt1_ptoid_in:lnPt`), 
    median(draws$b_Rt1_Rt + draws$`b_Rt1_Rt:ptoid_in`), 
    median(draws$b_Rt1_Yt + draws$`b_Rt1_Yt:ptoid_in`),
    # row 4
    median(draws$b_Yt1_lnAt + draws$`b_Yt1_lnAt:ptoid_in`), 
    median(draws$`b_Yt1_ptoid_in:lnPt`), 
    median(draws$b_Yt1_Rt + draws$`b_Yt1_Rt:ptoid_in`), 
    median(draws$b_Yt1_Yt + draws$`b_Yt1_Yt:ptoid_in`)),  
  nrow = 4, ncol = 4, byrow = T,
  dimnames = list(c("lnA","lnP","R","Y"),c("lnA","lnP","R","Y")))
M_aphid.ptoid_median

# identity matrix for aphid-ptoid system
I_M_aphid.ptoid <- diag(x = 1, nrow = nrow(M_aphid.ptoid_median), ncol = ncol(M_aphid.ptoid_median))

# get relevant sub matrices
A_aphid.ptoid_median <- M_aphid.ptoid_median[c("lnA","lnP"),c("lnA","lnP")]
B_aphid.ptoid_median <- M_aphid.ptoid_median[c("lnA","lnP"),c("R","Y")]
C_aphid.ptoid_median <- M_aphid.ptoid_median[c("R","Y"),c("lnA","lnP")]
D_aphid.ptoid_median <- M_aphid.ptoid_median[c("R","Y"),c("R","Y")]

# identity matrix for ecological and evolutionary subsystems
I_A_aphid.ptoid <- diag(x = 1, nrow = nrow(A_aphid.ptoid_median), ncol = ncol(A_aphid.ptoid_median))
I_D_aphid.ptoid <- diag(x = 1, nrow = nrow(D_aphid.ptoid_median), ncol = ncol(D_aphid.ptoid_median))

# calculate evo-eco-evo feedback loop
# Note different formulas for discrete-time model!
evoecoevo_aphid.ptoid_median <- C_aphid.ptoid_median %*% inv(I_A_aphid.ptoid - A_aphid.ptoid_median) %*% B_aphid.ptoid_median 
evoecoevo_aphid.ptoid_median

# if you only want to calculate the evo-to-eco effect or eco-to-evo effect (assuming stable equilibrium that takes into account ecological or feedbacks, respectively)
inv(I_A_aphid.ptoid - A_aphid.ptoid_median) %*% (B_aphid.ptoid_median) # evo-to-eco-effect
# interpretation: increasing Y from 0 to 1 increases aphid and parasitoid log densities by 1.5 and 1.3, respectively,
# in contrast, increasing R from 0 to 1 decreases aphid and parasitoid log densities by 0.4 and 0.7, respectively.
# whether these are small or large depend on the log densities...
inv(I_D_aphid.ptoid - D_aphid.ptoid_median) %*% (C_aphid.ptoid_median) # eco-to-evo effect
# interpretation: increasing parasitoid log density by 1 unit will decrease frequency of red morph (first row) and increase frequency of yellow morph (second row),
# which doesn't make sense, suggesting this equilibrium interpretation is not valid (it is unstable in fact, which makes sense)

# if evolutionary dynamics are faster than ecological dynamics, then the relevant feedback loop
# for stability is the eco-evo-eco feedback loop. We don't focus on this feedback loop in the manuscript
ecoevoeco_aphid.ptoid_median <- B_aphid.ptoid_median %*% inv(I_D_aphid.ptoid - D_aphid.ptoid_median) %*% C_aphid.ptoid_median # evo-eco-evo feedback loop
ecoevoeco_aphid.ptoid_median # note that this is unlikely to be accurate given the evolutionary feedback itself was unstable


# equilibrium values for discrete time model, see Ives et al. 2003, Ecol. Monog.
eq_aphid.ptoid <- inv(I_M_aphid.ptoid - M_aphid.ptoid_median) %*% x_aphid.ptoid_median
eq_aphid.ptoid # note that the equilibrium is not feasible (i.e. values outside plausible range)
# for example, red and yellow morph frequencies are beyond 0-1 interval
inv(I_M_aphid.ptoid - M_aphid.ptoid_median) %*% (x_aphid.ptoid_median + c(1,0,0,0))

# ecological equilibrium in absence of evolution
inv(I_A_aphid.ptoid - A_aphid.ptoid_median) %*% x_aphid.ptoid_median[c("lnA","lnP"),]
# lower aphid and lower ptoid
# this equilibrium is already unfeasible
inv(I_A_aphid.ptoid - A_aphid.ptoid_median) %*% (x_aphid.ptoid_median[c("lnA","lnP"),] + c(1,0)) # with aphid press perturbation


# evolutionary equilibrium in absence of ecology
inv(I_D_aphid.ptoid - D_aphid.ptoid_median) %*% x_aphid.ptoid_median[c("R","Y"),]
# greater differences between R and Y, with R because much smaller (but why?)
# but this evolutionary equilibrium is not feasible

# eco and evo equilibriums after dropping eco-to-evo and evo-to-eco effects
M_aphid.ptoid_median_dropEcoEvo <- M_aphid.ptoid_median
M_aphid.ptoid_median_dropEcoEvo[c("lnA","lnP"),c("R","Y")] <- 0
M_aphid.ptoid_median_dropEcoEvo[c("R","Y"),c("lnA","lnP")] <- 0
M_aphid.ptoid_median_dropEcoEvo
inv(I_M_aphid.ptoid - M_aphid.ptoid_median_dropEcoEvo) %*% x_aphid.ptoid_median
# matches separate analyses
max(Re(eigen(M_aphid.ptoid_median_dropEcoEvo)$values))

# stability analysis
max(Re(eigen(M_aphid.ptoid_median)$values)) # eco-evolutionary stability
max(Re(eigen(A_aphid.ptoid_median)$values)) # eco stability
max(Re(eigen(D_aphid.ptoid_median + evoecoevo_aphid.ptoid_median)$values)) # evo + evoecoevo feedback stability
max(Re(eigen(D_aphid.ptoid_median)$values)) # same as dropping eco-evo and evo-eco submatrices
# extra stability
max(Re(eigen(A_aphid.ptoid_median + ecoevoeco_aphid.ptoid_median)$values))

## Get aphid eco-evo parameters (all posteriors) ----
# the code below uses the posterior distribution for all parameters
# to propogate uncertainty to emergent properties, such as the equilibrium values and stability

# create empty matrices to be filled by the for-loop
M_aphid <- list()
M_stability_aphid <- list()
A_stability_aphid <- list()
D_stability_aphid <- list()
D_evoecoevo_stability_aphid <- list()

# create empty vectors to be fille by the for-loop
x_aphid <- list()
x_equilib_aphid <- list()
x_equilib_aphid_output <- list()

# for-loop
for(i in 1:nrow(draws)){
  # get vector of intrinsic growth rates for populations and selection for morph frequencies at low frequencies
  x_aphid[[i]] <- matrix(
    data = c(draws$b_lnAt1_Intercept[i], 
             draws$b_Rt1_Intercept[i],
             draws$b_Yt1_Intercept[i]),
    nrow = 3, ncol = 1, byrow = T,
    dimnames = list(c("lnA","R","Y"),c("x")))
  
  # get complete eco-evo effects matrix for aphids
  M_aphid[[i]] <- matrix(
    # no offset along diagonal
    data = c(draws$b_lnAt1_lnAt[i], draws$b_lnAt1_Rt[i], draws$b_lnAt1_Yt[i],
             draws$b_Rt1_lnAt[i], draws$b_Rt1_Rt[i], draws$b_Rt1_Yt[i],
             draws$b_Yt1_lnAt[i], draws$b_Yt1_Rt[i], draws$b_Yt1_Yt[i]),
    nrow = 3, ncol = 3, byrow = T,
    dimnames = list(c("lnA","R","Y"),c("lnA","R","Y")))
  
  # calculate equilibrium population sizes (log scale) and morph frequencies
  # equilibrium values for discrete time model, see Ives et al. 2003, Ecol. Monog.
  x_equilib_aphid[[i]] <- inv(I_aphid - M_aphid[[i]]) %*% x_aphid[[i]]
  x_equilib_aphid[[i]]
  rownames(x_equilib_aphid[[i]]) <- c("lnA","R","Y")
  x_equilib_aphid_output[[i]] <- data.frame(
    lnA_aphid = x_equilib_aphid[[i]][1],
    R_aphid = x_equilib_aphid[[i]][2],
    Y_aphid = x_equilib_aphid[[i]][3]) %>%
    mutate(G_aphid = 1 - R_aphid - Y_aphid)
  
  # extract relevant feedback loop
  A <- M_aphid[[i]][c("lnA"),c("lnA")] # eco
  C <- M_aphid[[i]][c("R","Y"),c("lnA")] # eco-to-evo
  B <- M_aphid[[i]][c("lnA"),c("R","Y")] # evo-to-eco
  D <- M_aphid[[i]][c("R","Y"),c("R","Y")] # evo
  evoecoevo <- C * (1/(1-A)) * B # evo-eco-evo feedback loop
  
  # calculate relevant stability metrics
  M_stability_aphid[[i]] <- max(Re(eigen(M_aphid[[i]])$values))
  A_stability_aphid[[i]] <- max(Re(eigen(A)$values))
  D_stability_aphid[[i]] <- max(Re(eigen(D)$values))
  D_evoecoevo_stability_aphid[[i]] <- max(Re(eigen(D + evoecoevo)$values))
}

# focus on key output
draws_ecoevo_aphid_df <- cbind(
  plyr::ldply(M_stability_aphid) %>% rename(M_stability_aphid = V1),
  plyr::ldply(A_stability_aphid) %>% rename(A_stability_aphid = V1),
  plyr::ldply(D_stability_aphid) %>% rename(D_stability_aphid = V1),
  plyr::ldply(D_evoecoevo_stability_aphid) %>% rename(D_evoecoevo_stability_aphid = V1),
  plyr::ldply(x_equilib_aphid_output))

## Get aphid-parasitoid eco-evo parameters (all posteriors) ----

# create matrices to be filled by for-loop
M_aphid.ptoid <- list()
M_stability_aphid.ptoid <- list()
A_stability_aphid.ptoid <- list()
D_stability_aphid.ptoid <- list()
D_evoecoevo_stability_aphid.ptoid <- list()

# create empty vectors to be filled by for-loop
x_aphid.ptoid <- list()
x_equilib_aphid.ptoid <- list()
x_equilib_aphid.ptoid_output <- list()

# for-loop
for(i in 1:nrow(draws)){
  # get vector of intrinsic growth rates for populations and selection for morph frequencies at low frequencies
  x_aphid.ptoid[[i]] <- matrix(
    data = c(draws$b_lnAt1_Intercept[i] + draws$b_lnAt1_ptoid_in[i], 
             draws$b_lnPt1_Intercept[i],
             draws$b_Rt1_Intercept[i] + draws$b_Rt1_ptoid_in[i],
             draws$b_Yt1_Intercept[i] + draws$b_Yt1_ptoid_in[i]),
    nrow = 4, ncol = 1, byrow = T,
    dimnames = list(c("lnA","lnP","R","Y"),c("x")))
  
  # get complete eco-evo matrix for aphid.ptoids
  M_aphid.ptoid[[i]] <- matrix(
    data = c(
      # note that we add relevant coefficients that are modified by the parasitoid
      # no offset along diagonal
      # row 1
      draws$b_lnAt1_lnAt[i] + draws$`b_lnAt1_lnAt:ptoid_in`[i], 
      draws$`b_lnAt1_ptoid_in:lnPt`[i],
      draws$b_lnAt1_Rt[i] + draws$`b_lnAt1_Rt:ptoid_in`[i],
      draws$b_lnAt1_Yt[i] + draws$`b_lnAt1_Yt:ptoid_in`[i],
      # row 2, no need to add coeficients since these were only estimated with 'ptoid_in'
      draws$b_lnPt1_lnAt[i], draws$b_lnPt1_lnPt[i], draws$b_lnPt1_Rt[i], draws$b_lnPt1_Yt[i],
      # row 3
      draws$b_Rt1_lnAt[i] + draws$`b_Rt1_lnAt:ptoid_in`[i],
      draws$`b_Rt1_ptoid_in:lnPt`[i],
      draws$b_Rt1_Rt[i] + draws$`b_Rt1_Rt:ptoid_in`[i], 
      draws$b_Rt1_Yt[i] + draws$`b_Rt1_Yt:ptoid_in`[i],
      # row 4
      draws$b_Yt1_lnAt[i] + draws$`b_Yt1_lnAt:ptoid_in`[i],
      draws$`b_Yt1_ptoid_in:lnPt`[i],
      draws$b_Yt1_Rt[i] + draws$`b_Yt1_Rt:ptoid_in`[i],
      draws$b_Yt1_Yt[i] + draws$`b_Yt1_Yt:ptoid_in`[i]),
    nrow = 4, ncol = 4, byrow = T,
    dimnames = list(c("lnA","lnP","R","Y"),c("lnA","lnP","R","Y")))
  
  # calculate equilibrium population sizes (log scale) and morph frequencies
  # equilibrium values for discrete time model, see Ives et al. 2003, Ecol. Monog.
  x_equilib_aphid.ptoid[[i]] <- inv(I_M_aphid.ptoid - M_aphid.ptoid[[i]]) %*% x_aphid.ptoid[[i]]
  x_equilib_aphid.ptoid[[i]]
  rownames(x_equilib_aphid.ptoid[[i]]) <- c("lnA","lnP","R","Y")
  x_equilib_aphid.ptoid_output[[i]] <- data.frame(
    lnA_aphid.ptoid = x_equilib_aphid.ptoid[[i]][1],
    lnP_aphid.ptoid = x_equilib_aphid.ptoid[[i]][2],
    R_aphid.ptoid = x_equilib_aphid.ptoid[[i]][3],
    Y_aphid.ptoid = x_equilib_aphid.ptoid[[i]][4]) %>%
    mutate(G_aphid.ptoid = 1 - R_aphid.ptoid - Y_aphid.ptoid)
  
  # extract relevant feedback loop
  A <- M_aphid.ptoid[[i]][c("lnA","lnP"),c("lnA","lnP")] # eco
  C <- M_aphid.ptoid[[i]][c("R","Y"),c("lnA","lnP")] # eco-to-evo
  B <- M_aphid.ptoid[[i]][c("lnA","lnP"),c("R","Y")] # evo-to-eco
  D <- M_aphid.ptoid[[i]][c("R","Y"),c("R","Y")] # evo
  evoecoevo <- C %*% inv(I_A_aphid.ptoid - A) %*% B # evo-eco-evo feedback loop

  # calculate relevant stability metrics
  M_stability_aphid.ptoid[[i]] <- max(Re(eigen(M_aphid.ptoid[[i]])$values))
  A_stability_aphid.ptoid[[i]] <- max(Re(eigen(A)$values))
  D_stability_aphid.ptoid[[i]] <- max(Re(eigen(D)$values))
  D_evoecoevo_stability_aphid.ptoid[[i]] <- max(Re(eigen(D + evoecoevo)$values))
}

# focus on key output
draws_ecoevo_aphid.ptoid_df <- cbind(
  plyr::ldply(M_stability_aphid.ptoid) %>% rename(M_stability_aphid.ptoid = V1),
  plyr::ldply(A_stability_aphid.ptoid) %>% rename(A_stability_aphid.ptoid = V1),
  plyr::ldply(D_stability_aphid.ptoid) %>% rename(D_stability_aphid.ptoid = V1),
  plyr::ldply(D_evoecoevo_stability_aphid.ptoid) %>% rename(D_evoecoevo_stability_aphid.ptoid = V1),
  plyr::ldply(x_equilib_aphid.ptoid_output))

## Save output for plotting in separate R scripts ----
save.image("analyses/posteriors_ecoevo_model.RData")
