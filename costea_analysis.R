
################# Load libraries ##################
library(tidyverse)
library(magrittr)
library(devtools)
load_all()

################# Load data ##################
load("costea2017_metaphlan2_profiles.rda")
load("costea2017_mock_composition.rda")
load("costea2017_sample_data.rda")

################# Filter data (species level, mock taxa only) ##################
costea2017_metaphlan2_profiles_species <- costea2017_metaphlan2_profiles %>%
  filter(sapply(Clade,function(x) grepl("s__",x,fixed = TRUE))) %>%
  filter(sapply(Clade, function(x) !grepl("t__",x,fixed = TRUE))) %>%
  filter(sapply(Clade,function(x)
    grepl("saccharolyticum",x,fixed = TRUE) |
      grepl("perfringens",x,fixed = TRUE) |
      grepl("melaninogenica",x,fixed = TRUE)|
      grepl("difficile",x,fixed = TRUE)|
      grepl("enterica",x,fixed = TRUE)|
      grepl("plantarum",x,fixed = TRUE)|
      grepl("cholerae",x,fixed = TRUE)|
      grepl("pseudotuberculosis",x,fixed = TRUE)|
      grepl("Blautia_hansenii",x,fixed = TRUE)|
      grepl("nucleatum",x,fixed = TRUE)
    ))

mock_cols <- costea2017_mock_composition$Taxon

mock_mat <- t(as.matrix(costea2017_mock_composition[,-1]))
colnames(mock_mat) <- mock_cols

################# Construct observation matrix W ##################

cell_concentration <- costea2017_mock_composition

W <- costea2017_mock_composition %>%
  as.data.frame() %>%
  (function(x){ colnames(x)[2] <- "cells_per_ml"
   return(x)}) %>%
  group_by(Taxon) %>%
  summarize(cells_per_ml = mean(cells_per_ml)) %>%
  (function(x){ m <- matrix(x$cells_per_ml,nrow = 1)
  colnames(m) <- x$Taxon
  return(m)})

W_reorder <- sapply(colnames(W),
       function(x) which(sapply(costea2017_metaphlan2_profiles_species$Clade,
                                function(d) grepl(x,d,fixed = TRUE)))) %>%
  as.numeric()

W <- rbind(W,
           costea2017_metaphlan2_profiles_species[W_reorder,-1] %>%
             as.matrix() %>%
             t)

################# Specify Designs, etc. to fit full and null models ##################

protocols <- sapply(rownames(W)[-1],
       function(x) costea2017_sample_data$Protocol[
         costea2017_sample_data$Run_accession == x
       ])

individuals <- sapply(rownames(W)[-1],
                    function(x) costea2017_sample_data$Individual[
                      costea2017_sample_data$Run_accession == x
                    ])

X <- lapply(protocols,
            function(x) as.numeric(
              c(x == "H",x == "Q", x == "W")
            ))
X <- do.call(rbind,X)
X[,1] <- 1

X <- rbind(0,X)

Z <- matrix(1,nrow = 30, ncol = 1)

Z_tilde <- matrix(0,nrow= 30, ncol = 1)
Z_tilde_gamma_cols <- 1

gammas <- apply(W,1,function(x) log(sum(x)))
gammas_fixed_indices <- rep(FALSE, length(gammas))

P <- matrix(W[1,]/sum(W[1,]),nrow =1, ncol = 10)

P_fixed_indices <- matrix(FALSE,nrow = 1, ncol = 10)

X_tilde <- matrix(0,ncol = 3, nrow= 1)

B <- matrix(0,
            ncol = 10,
            nrow = 3)

B_fixed_indices <- matrix(FALSE,ncol = 10,nrow = 3)
B_fixed_indices[,10] <- TRUE


P_tilde <- P*0
P_tilde_fixed_indices <- !P_fixed_indices
gamma_tilde <- matrix(0,ncol = 1, nrow = 1)
gamma_tilde_fixed_indices <- TRUE

################# Fit full and null models ##################

full_model  <-
  estimate_parameters(W = W,
                    X = X,
                    Z = Z,
                    Z_tilde = Z_tilde,
                    Z_tilde_gamma_cols = 1,
                    gammas = gammas,
                    gammas_fixed_indices = gammas_fixed_indices,
                    P = P,
                    P_fixed_indices = P_fixed_indices,
                    B = B,
                    B_fixed_indices = B_fixed_indices,
                    X_tilde = X_tilde,
                    P_tilde = P_tilde,
                    P_tilde_fixed_indices = P_tilde_fixed_indices,
                    gamma_tilde = gamma_tilde,
                    gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                    alpha_tilde = NULL,
                    Z_tilde_list = NULL,
                    barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                    barrier_scale = 10, #increments for value of barrier penalty
                    max_barrier = 1e12, #maximum value of barrier_t
                    initial_conv_tol = 1000,
                    final_conv_tol = 0.1,
                    final_f = 1e-6,
                    constraint_tolerance = 1e-10,
                    hessian_regularization = 0.01,
                    criterion = "Poisson",
                    subproblem_method = "Newton",
                    profile_P = FALSE,
                    profiling_maxit = 25,
                    wts = NULL,
                    verbose = FALSE)

### Create null model specification for test
null_model <- full_model

### set second and third rows of B equal to zero
null_model$B[2:3,] <- 0
null_model$B_fixed_indices[2:3,] <- TRUE

bootstrap_test <-
  bootstrap_lrt(W = W,
                fitted_model= full_model,
                null_param = null_model,
                n_boot = 1000,
                parallelize = TRUE)

### fit reparametrization with X_i = [1_H 1_Q 1_W]
### to (more easily) get CIs for protocol-specific effects
### (also can bootstrap from original model, but currently a bit more
### involved to pull out cis for quantities of form A%*%beta than for
### beta itself)
X_repar <- X
X_repar[,1] <- X_repar[,1] - X_repar[,2] - X_repar[,3]
full_reparam  <-
  estimate_parameters(W = W,
                      X = X_repar,
                      Z = Z,
                      Z_tilde = Z_tilde,
                      Z_tilde_gamma_cols = 1,
                      gammas = gammas,
                      gammas_fixed_indices = gammas_fixed_indices,
                      P = P,
                      P_fixed_indices = P_fixed_indices,
                      B = B,
                      B_fixed_indices = B_fixed_indices,
                      X_tilde = X_tilde,
                      P_tilde = P_tilde,
                      P_tilde_fixed_indices = P_tilde_fixed_indices,
                      gamma_tilde = gamma_tilde,
                      gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                      alpha_tilde = NULL,
                      Z_tilde_list = NULL,
                      barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                      barrier_scale = 10, #increments for value of barrier penalty
                      max_barrier = 1e12, #maximum value of barrier_t
                      initial_conv_tol = 1000,
                      final_conv_tol = 0.1,
                      final_f = 1e-6,
                      constraint_tolerance = 1e-10,
                      hessian_regularization = 0.01,
                      criterion = "Poisson",
                      subproblem_method = "Newton",
                      profile_P = FALSE,
                      profiling_maxit = 25,
                      wts = NULL,
                      verbose = TRUE)

full_cis <- bootstrap_ci(W,
                         fitted_model = full_reparam,
                         n_boot=1000,
                         m = NULL,
                         alpha = 0.05,
                         parallelize = TRUE,
                         ncores = 5,
                         seed = 3423,
                         return_models = FALSE,
                         verbose = FALSE

)

taxa <- colnames(W)
taxa <- lapply(taxa,
               function(x) strsplit(x,"_") %>%
                 (function(y) paste(substr(y[[1]][1],1,1),". ",
                                    y[[1]][2],sep = "",collapse = ""))) %>%
  do.call(c,.)

full_cis$ci %>%
  filter(param == "B") %>%
  mutate(Protocol = c("H","Q","W")[k],
         Estimate = round(value,2)) %>%
  mutate(Upper = round(upper_ci,2),
         Lower = round(lower_ci,2),
         Taxon = rep(taxa[1:9],3)) %>%
  dplyr::select(Protocol, Taxon, Estimate, Lower, Upper) %>%
  mutate(Estimate = apply(cbind(Estimate,Lower,Upper),
                          1, function(x) paste(x[1]," (", x[2], " - ", x[3], ")",
                                               sep = "",
                                               collapse = ""))) %>%
  select(c(Protocol, Taxon, Estimate)) %>%
  pivot_wider(id_cols = Taxon, names_from = Protocol,
              values_from = Estimate) %>%
  knitr::kable(format = "latex")


################# 10-fold cross-validation ##################

# construct folds
folds <- vector(10, mode = "list")
available <- 2:30
unique_individuals <- c(1:8,"M")

for(i in 1:9){
  folds[[i]] <- which(individuals == unique_individuals[i]) +1
}

folds[[10]] <- which(individuals %in% c("A","B")) +1

full_cv <- vector(10, mode = "list")
null_cv <- vector(10, mode = "list")
for(whichfoldout in 1:10){
  print(whichfoldout)
  heldout <- folds[[whichfoldout]]
  nheldout <- length(heldout)
  Z_cv <- Z
  Z_cv <- cbind(Z_cv - Z_cv*(1:30 %in% heldout))
  P_cv <- P
  P_fixed_indices_cv <- P_fixed_indices
  for(k in 1:nheldout){
    Z_cv <- cbind(Z_cv,as.numeric(1:30 == heldout[k]))
    P_cv <- rbind(P_cv,P)
    P_fixed_indices_cv <- rbind(P_fixed_indices_cv,
                                P_fixed_indices)
  }

  full_cv[[whichfoldout]]  <-
    estimate_parameters(W = W,
                        X = X,
                        Z = Z_cv,
                        Z_tilde = Z_tilde,
                        Z_tilde_gamma_cols = 1,
                        gammas = gammas,
                        gammas_fixed_indices = gammas_fixed_indices,
                        P = P_cv,
                        P_fixed_indices = P_fixed_indices_cv,
                        B = B,
                        B_fixed_indices = B_fixed_indices,
                        X_tilde = X_tilde,
                        P_tilde = P_tilde,
                        P_tilde_fixed_indices = P_tilde_fixed_indices,
                        gamma_tilde = gamma_tilde,
                        gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                        alpha_tilde = NULL,
                        Z_tilde_list = NULL,
                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                        barrier_scale = 10, #increments for value of barrier penalty
                        max_barrier = 1e12, #maximum value of barrier_t
                        initial_conv_tol = 1000,
                        final_conv_tol = 0.1,
                        final_f = 1e-6,
                        constraint_tolerance = 1e-10,
                        hessian_regularization = 0.01,
                        criterion = "Poisson",
                        subproblem_method = "Newton",
                        profile_P = FALSE,
                        profiling_maxit = 25,
                        wts = NULL,
                        verbose = FALSE)

  null_cv[[whichfoldout]]  <-
    estimate_parameters(W = W,
                        X = X[,1,drop = FALSE],
                        Z = Z_cv,
                        Z_tilde = Z_tilde,
                        Z_tilde_gamma_cols = 1,
                        gammas = gammas,
                        gammas_fixed_indices = gammas_fixed_indices,
                        P = P_cv,
                        P_fixed_indices = P_fixed_indices_cv,
                        B = B[1,,drop = FALSE],
                        B_fixed_indices = B_fixed_indices[1,,drop = FALSE],
                        X_tilde = X_tilde[,1,drop = FALSE],
                        P_tilde = P_tilde,
                        P_tilde_fixed_indices = P_tilde_fixed_indices,
                        gamma_tilde = gamma_tilde,
                        gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                        alpha_tilde = NULL,
                        Z_tilde_list = NULL,
                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                        barrier_scale = 10, #increments for value of barrier penalty
                        max_barrier = 1e12, #maximum value of barrier_t
                        initial_conv_tol = 1000,
                        final_conv_tol = 0.1,
                        final_f = 1e-6,
                        constraint_tolerance = 1e-10,
                        hessian_regularization = 0.01,
                        criterion = "Poisson",
                        subproblem_method = "Newton",
                        profile_P = FALSE,
                        profiling_maxit = 25,
                        wts = NULL,
                        verbose = FALSE)

}

full_cv_predictions <- lapply(1:10,
                              function(x)
                                full_cv[[x]]$varying[
                                  full_cv[[x]]$varying$param == "P"&
                                    full_cv[[x]]$varying$k>1,])

for(i in 1:10){
  full_cv_predictions[[i]]$k <- sapply(full_cv_predictions[[i]]$k,
                                       function(x) folds[[i]][x - 1])

}

full_cv_predictions <- do.call(rbind,full_cv_predictions)
fc_values <- W[1,]/sum(W[1,])
full_cv_predictions$fc_value <- sapply(full_cv_predictions$j,
                                       function(d) fc_values[d])
full_cv_predictions$protocol <-
  sapply(full_cv_predictions$k,
         function(d) protocols[d-1]) #  d - 1 bc k starts at 2
                                     #  (k = 1 is fc data)


full_cv_predictions$specimen <-
  sapply(full_cv_predictions$k,
         function(d) individuals[d - 1]) #  d - 1 bc k starts at 2
                                         #  (k = 1 is fc data)


null_cv_predictions <- lapply(1:10,
                              function(x)
                                null_cv[[x]]$varying[
                                  null_cv[[x]]$varying$param == "P"&
                                    null_cv[[x]]$varying$k>1,])

for(i in 1:10){
  null_cv_predictions[[i]]$k <- sapply(null_cv_predictions[[i]]$k,
                                       function(x) folds[[i]][x - 1])

}

null_cv_predictions <- do.call(rbind,null_cv_predictions)
fc_values <- W[1,]/sum(W[1,])
null_cv_predictions$fc_value <- sapply(null_cv_predictions$j,
                                       function(d) fc_values[d])
null_cv_predictions$protocol <-
  sapply(null_cv_predictions$k,
         function(d) protocols[d -1 ]) #  d - 1 bc k starts at 2 (k = 1 is fc data)

null_cv_predictions$specimen <-
  sapply(null_cv_predictions$k,
         function(d) individuals[d - 1]) #  d - 1 bc k starts at 2 (k = 1 is fc data)
null_cv_predictions$model <- "Null Model"
full_cv_predictions$model <- "Full Model"

W_prop <- W[-1,]
for(i in 1:nrow(W_prop)){
  W_prop[i,] <- W_prop[i,]/sum(W_prop[i,])
}

naive_predictions <- null_cv_predictions[numeric(0),]

for(i in 1:nrow(W_prop)){
  protocol <- protocols[i]
  specimen <- individuals[i]
  for(j in 1:ncol(W)){
    naive_predictions <- rbind(naive_predictions,
                               data.frame(value = W_prop[i,j],
                                          param = "P",
                                          k = i,
                                          j  = j,
                                          fc_value = fc_values[j],
                                          protocol = protocol,
                                          specimen = specimen,
                                          model = "Plug-in"))
  }
}




rbind(null_cv_predictions,full_cv_predictions,
      naive_predictions) %>%
  mutate(protocol = sapply(protocol, function(x) paste("Protocol ",x,
                                                       sep = "",
                                                       collapse = ""))) %>%
  # filter(!is.na(protocol)) %>% #why would protocol be NA? Check!
  ggplot() +
  geom_point(aes(x = fc_value, y = value, color= specimen),
             size = .5) +
  geom_line(aes(x = fc_value, y = value, color = specimen,
                group = as.factor(k)), size = .5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  facet_grid(model~protocol,scales = "free_y") +
  scale_color_grey() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  guides(color=guide_legend(title="Specimen")) +
  xlab("Relative Abundance Measured by Flow Cytometry") +
  ylab("Cross-Validated Estimated Relative Abundance")






