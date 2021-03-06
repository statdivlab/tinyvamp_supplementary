###### Load Libraries ######

library(tidyverse)
library(devtools)
library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(Matrix)
library(phyloseq)
library(gridExtra)
library(logsum)
library(fastnnls)
library(here)

load_all()

###### Load preprocessed Karstens et al. (2019) data ######
load(here("data/karstens_phyloseq.Rdata"))

### Inspect sample data
sample_data(ps_pooled_species)

### Pull out column indices corresponding to taxa in the mock
Pseudomonas <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Pseudomonas",x,
                                                                                          fixed = TRUE)) %>%
  which()

Escherichia <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Escherichia",x,
                                                                                          fixed = TRUE)) %>%
  which()

Salmonella <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Salmonella",x,
                                                                                         fixed = TRUE)) %>%
  which()

Limosilactobacillus <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Limosilactobacillus",x,
                                                                                                  fixed = TRUE)) %>%
  which()

Enterococcus <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Enterococcus",x,
                                                                                           fixed = TRUE)) %>%
  which()

Staphylococcus <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Staphylococcus",x,
                                                                                             fixed = TRUE)) %>%
  which()

Listeria <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Listeria",x,
                                                                                       fixed = TRUE)) %>%
  which()

Bacillus <- tax_table(ps_pooled_species)@.Data[,"Genus"]  %>% sapply(function(x) grepl("Bacillus",x,
                                                                                       fixed = TRUE)) %>%
  which()

mock_taxa <- c(Pseudomonas,
               Escherichia,
               Salmonella[1],
               Enterococcus,
               Staphylococcus,
               Listeria,
               Bacillus,
               Limosilactobacillus)

###### Construct count data matrix W ######
W <- ps_pooled_species@otu_table
dim(W)

W <- as.matrix(W)

# reorder columns so taxa in mock come last

W <- cbind(W[,!(1:ncol(W) %in% mock_taxa)],
           W[,mock_taxa])

#reorder rows so dilutions in increasing order

W <- W[order(rownames(W)),]


###### Construct helper function to specify models ######
specify_dilution_model <- function(
  samples_treated_as_known,
  samples_treated_as_unknown,
  intercept_in_Z_tilde = FALSE,
  separate_dilution_series = FALSE,
  use_alpha_tilde = FALSE,
  beta_for_dilution = FALSE,
  W){
  known_indicator <- as.numeric(1:9 %in% samples_treated_as_known)


  dilutions <- 0:8

  ###
  Z <- cbind(known_indicator,
             1 - known_indicator)

  colnames(Z) <- NULL

  ### Z_tilde
  # set up dilutions as covariate
  Z_tilde <- matrix(3^(dilutions)/exp(mean(log(3^(0:8)))),ncol = 1)

  # separate dilution covariate into two columns if needed
  if(separate_dilution_series){
    Z_tilde <- cbind(Z_tilde[,1]*known_indicator,
                     Z_tilde[,1]*(1 - known_indicator))
  }

  # add intercept column if needed
  if(intercept_in_Z_tilde){
    Z_tilde <- cbind(1, Z_tilde)
  }

  Z_tilde_gamma_cols <- 1:ncol(Z_tilde)

  if(use_alpha_tilde){

    Z_tilde_list <- list(
      known_indicator*Z_tilde[,ncol(Z_tilde),drop = FALSE],
      (1 - known_indicator)*Z_tilde[,ncol(Z_tilde),drop = FALSE]
    )
    if(intercept_in_Z_tilde){
      Z_tilde_list[[1]] <- cbind(1, Z_tilde_list[[1]])
      Z_tilde_list[[2]] <- cbind(0, Z_tilde_list[[2]])
    }
    Z_tilde <- NULL
    alpha_tilde <- rnorm(1)
  } else{
    Z_tilde_list <- NULL
    alpha_tilde <- NULL
  }


  ### gamma_tilde
  K_tilde <- max(Z_tilde_gamma_cols)
  gamma_tilde <- matrix(-5, nrow = K_tilde,ncol = 1)
  gamma_tilde[] <- rnorm(K_tilde,-2)

  gamma_tilde_fixed_indices <- matrix(FALSE, ncol = ncol(gamma_tilde),
                                      nrow = nrow(gamma_tilde))

  ### P_tilde
  P_tilde <- matrix(1/ncol(W),ncol = ncol(W), nrow = K_tilde)
  P_tilde_fixed_indices <- matrix(FALSE,
                                  ncol = ncol(P_tilde),
                                  nrow = nrow(P_tilde))

  # fix contaminants in mock to be zero
  P_tilde[,ncol(P_tilde_fixed_indices)] <- 0

  P_tilde_fixed_indices[,ncol(P_tilde_fixed_indices)] <- TRUE

  for(k_tilde in 1:nrow(P_tilde)){
    P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]] <-
      rexp(sum(!P_tilde_fixed_indices[k_tilde,]))
    P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]] <-
      P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]]/sum(
        P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]]
      )
  }

  ### B
  B <- matrix(0,ncol = ncol(W),
              nrow = 1)

  if(beta_for_dilution){
    B <- rbind(B,B)
  }


  B_fixed_indices <- matrix(TRUE,ncol = ncol(B),nrow = nrow(B))
  B_fixed_indices[,ncol(B) - 7:1] <- FALSE

  B[!B_fixed_indices] <- rnorm(sum(!B_fixed_indices))

  ### X
  X <- matrix(1,nrow = nrow(W),ncol = 1)
  if(beta_for_dilution){
    X <- cbind(X,0:8)
  }

  ### X_tilde
  X_tilde <- matrix(0, nrow = K_tilde,
                    ncol = nrow(B))

  ### gammas

  gammas <- log(apply(W,1,sum))
  gammas <- gammas + rnorm(length(gammas))

  gammas_fixed_indices <-rep(F,length(gammas))

  ### P
  P <- rbind(c(rep(0,ncol(W) - 8), rep(1/8,8)),
             rep(1/ncol(W),ncol(W)))

  P_fixed_indices <- matrix(FALSE,ncol = ncol(W),nrow = 2)
  P_fixed_indices[1,] <- TRUE

  P[2,!P_fixed_indices[2,]] <- rexp(sum(!P_fixed_indices[2,]))
  P[2,!P_fixed_indices[2,]] <-   P[2,!P_fixed_indices[2,]]/sum(
    P[2,!P_fixed_indices[2,]]
  )


  return(list(
    "Z" = Z,
    "Z_tilde" = Z_tilde,
    "Z_tilde_gamma_cols" = Z_tilde_gamma_cols,
    "P" = P,
    "P_fixed_indices" = P_fixed_indices,
    "gammas" = gammas,
    "gammas_fixed_indices" = gammas_fixed_indices,
    "X_tilde" = X_tilde,
    "X" = X,
    "B" = B,
    "B_fixed_indices" = B_fixed_indices,
    "P_tilde" = P_tilde,
    "P_tilde_fixed_indices" = P_tilde_fixed_indices,
    "gamma_tilde" = gamma_tilde,
    "gamma_tilde_fixed_indices" = gamma_tilde_fixed_indices,
    "alpha_tilde" = alpha_tilde,
    "Z_tilde_list" = Z_tilde_list
  ))

}
######## Fit Model (No folds) without betas ########

Z <- matrix(1,ncol = 1, nrow = 9)

colnames(Z) <- NULL

### Z_tilde
# set up dilutions as covariate
Z_tilde <- matrix(3^(0:8)/exp(mean(log(3^(0:8)))),ncol = 1)

Z_tilde_gamma_cols <- 1
Z_tilde_list <- NULL
alpha_tilde <- NULL



### gamma_tilde
K_tilde <- max(Z_tilde_gamma_cols)
gamma_tilde <- matrix(0, nrow = K_tilde,ncol = 1)


gamma_tilde_fixed_indices <- matrix(FALSE, ncol = ncol(gamma_tilde),
                                    nrow = nrow(gamma_tilde))

### P_tilde
P_tilde <- matrix(1/ncol(W),ncol = ncol(W), nrow = K_tilde)
P_tilde_fixed_indices <- matrix(FALSE,
                                ncol = ncol(P_tilde),
                                nrow = nrow(P_tilde))

# fix contaminants in mock to be zero
P_tilde[,ncol(P_tilde_fixed_indices)] <- 0

P_tilde_fixed_indices[,ncol(P_tilde_fixed_indices)] <- TRUE

for(k_tilde in 1:nrow(P_tilde)){
  P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]] <-
    rexp(sum(!P_tilde_fixed_indices[k_tilde,]))
  P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]] <-
    P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]]/sum(
      P_tilde[k_tilde,!P_tilde_fixed_indices[k_tilde,]]
    )
}

### B
B <- matrix(0,ncol = ncol(W),
            nrow = 1)




B_fixed_indices <- matrix(TRUE,ncol = ncol(B),nrow = nrow(B))

### X
X <- matrix(0,nrow = nrow(W),ncol = 1)


### X_tilde
X_tilde <- matrix(0, nrow = K_tilde,
                  ncol = nrow(B))

### gammas

gammas <- log(apply(W,1,sum))
gammas <- gammas + rnorm(length(gammas))

gammas_fixed_indices <-rep(F,length(gammas))

### P
P <- matrix(c(rep(0,240),rep(1/8,8)),nrow = 1, ncol = ncol(W))

P_fixed_indices <- matrix(TRUE,ncol = ncol(W),nrow = 1)

B_fixed_indices <- matrix(TRUE,ncol = 248, nrow = 1)
B_fixed_indices[1,241:247] <- FALSE

full_karstens_model <-
  estimate_parameters(W = W,
                    X = matrix(1,ncol = 1, nrow = 9),
                    Z = Z,
                    Z_tilde = Z_tilde,
                    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
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
                    alpha_tilde = alpha_tilde,
                    Z_tilde_list = Z_tilde_list,
                    barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                    barrier_scale = 10, #increments for value of barrier penalty
                    max_barrier = 1e12, #maximum value of barrier_t
                    initial_conv_tol = 1000,
                    final_conv_tol = 0.1,
                    final_f = 1e-6,
                    constraint_tolerance = 1e-10,
                    hessian_regularization = 0.01,
                    criterion = "reweighted_Poisson",
                    subproblem_method = "Newton",
                    profile_P = TRUE,
                    profiling_maxit = 25,
                    verbose = TRUE,
                    wts = NULL)

full_cis <- bootstrap_ci(W,
                         fitted_model = full_karstens_model,
                         n_boot = 100,
                         parallelize = TRUE,
                         seed  = 4323)

full_cis$ci %>%
  filter(param == "B")

saveRDS(full_karstens_model,"full_karstens_model")

P <- matrix(1/ncol(W),nrow = 1, ncol = ncol(W))

P_fixed_indices <- matrix(FALSE,ncol = ncol(W),nrow = 1)

######## Three-fold CV ########



folds <- cbind(c(1,4,7),c(2,5,8),c(3,6,9))

fits <- vector(ncol(folds),
                      mode = "list")

tries <- fits

for(try in 1:length(tries)){
  tries[[try]] <- vector(1,mode = "list")
}
set.seed(829343)
for(fold in 1:ncol(folds)){
  message(paste("Fitting fold", fold))
  unknown_samples <- folds[,fold]
  known_samples <- setdiff(1:9,unknown_samples)

  for(try in 1:1){
    print(paste("Try ", try, sep = "", collapse = ""))
inputs <- specify_dilution_model(
  samples_treated_as_known = known_samples,
  samples_treated_as_unknown = unknown_samples,
  intercept_in_Z_tilde = FALSE,
  separate_dilution_series = FALSE,
  beta_for_dilution = FALSE,
  use_alpha_tilde = TRUE,
  W = W)

tries[[fold]][[try]] <-
  estimate_parameters(W = W,
                    X = inputs$X,
                    Z = inputs$Z,
                    Z_tilde = inputs$Z_tilde,
                    Z_tilde_gamma_cols = inputs$Z_tilde_gamma_cols,
                    gammas = inputs$gammas,
                    gammas_fixed_indices = inputs$gammas_fixed_indices,
                    P = inputs$P,
                    P_fixed_indices = inputs$P_fixed_indices,
                    B = inputs$B,
                    B_fixed_indices = inputs$B_fixed_indices,
                    X_tilde = inputs$X_tilde,
                    P_tilde = inputs$P_tilde,
                    P_tilde_fixed_indices = inputs$P_tilde_fixed_indices,
                    gamma_tilde = inputs$gamma_tilde,
                    gamma_tilde_fixed_indices = inputs$gamma_tilde_fixed_indices,
                    alpha_tilde = inputs$alpha_tilde,
                    Z_tilde_list = inputs$Z_tilde_list,
                    barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                    barrier_scale = 10, #increments for value of barrier penalty
                    max_barrier = 1e12, #maximum value of barrier_t
                    initial_conv_tol = 1000,
                    final_conv_tol = 0.1,
                    final_f = 1e-6,
                    constraint_tolerance = 1e-10,
                    hessian_regularization = 0.01,
                    criterion = "reweighted_Poisson",
                    subproblem_method = "Newton",
                    verbose = TRUE,
                    profile_P = TRUE,
                    profiling_maxit = 25,
                    wts = NULL)

  }
  lls <- sapply(1:1,
                function(x) tries[[fold]][[x]]$objective)
  best_fit <- which.min(lls)

  fits[[fold]] <- tries[[fold]][[best_fit]]
}


rbind(
do.call(rbind,lapply(1:3,
       function(d)
         fits[[d]]$varying %>%
         # filter(param == "P",j>240) %>%
         mutate(Fold = d,
                model = "Beta in Dilution")))#,
# do.call(rbind,
# lapply(1:3,
#        function(d)
#          gmm_fits[[d]]$varying %>%
#          # filter(param == "P",j>240) %>%
#          mutate(Fold = d,
#                 model = "Beta Varying in Dilution"))
# )
# ) %>%
#
#   do.call(rbind,lapply(4:6,
#                        function(d)
#                          simple_fits[[d]]$varying %>%
#                          # filter(param == "P",j>240) %>%
#                          mutate(Fold = d,
#                                 model = "Beta in Dilution"))) %>%
#   mutate(Fold_type = sapply(Fold, function(x) ifelse(x<4,"Unbalanced",
#                                                      "Balanced"))) %>%
) %>%
  mutate(Fold = as.factor(Fold)) %>%
  # dplyr::filter(model == "Beta in Dilution") %>%
  # dplyr::filter(Fold_type == "Balanced") %>%
  saveRDS(file = "karstens_3fold")




read_proportions <- vector(9, mode = "list")
names(read_proportions) <- paste("D",c(0:8), sep ="")
for(dilution in 0:8){
  read_counts <- W[paste("D",dilution, sep = "", collapse = ""),]
  read_prop <- read_counts/sum(read_counts)
  read_proportions[[dilution + 1]] <- read_prop
}

mock_taxa_names <- sapply(mock_taxa,
                          function(x) paste(as.character(
                            tax_table(ps_pooled_species)[x,6:7]),
                            sep = " ",
                            collapse = " "))

mock_taxa_names["ASV_1"] <- "Limosilactobacillus (species unclassified)"
values <- c(do.call(c,read_proportions), c(rep(0,240),rep(1/8,8)))
taxa_names <- c("Salmonella (species unclassified)",
                rep("Other",239),
                mock_taxa_names)
taxa_names <- factor(taxa_names,
                     levels = c(mock_taxa_names,
                                "Salmonella (species unclassified)",
                                "Other"))
dilutions <- rep(c(0:8,"Theoretical"),each = 248)

observed_read_proportions <-
  data.frame("Read_Proportions" = values,
             "Taxon" = rep(taxa_names,10),
             "Dilution" = as.character(dilutions))

observed_read_proportions$Type <- sapply(observed_read_proportions$Dilution,
                                         function(x) ifelse(x == "Theoretical",
                                                            "Theoretical",
                                                            "Observed"))

observed_props <- observed_read_proportions %>%
  group_by(Dilution, Taxon, Type) %>%
  summarize(Proportion = sum(Read_Proportions))

karstens_3fold <- readRDS("karstens_3fold")


fitted_proportions_3fold <- karstens_3fold %>%
  filter(param == "P") %>%
  select(c("value","Fold","j")) %>%
  mutate("Taxon" = sapply(j,function(x)
    ifelse(x<241 & x>1,"Other",
           ifelse(x == 241,
                  "P. aeruginosa",
                  ifelse(x == 242,
                         "E. coli",
                         ifelse(x == 243,
                                "S. enterica",
                                ifelse(x == 244,
                                       "E. faecium",
                                       ifelse(x == 245,
                                              "S. aureus",
                                              ifelse(x == 246,
                                                     "L. monocytogenes",
                                                     ifelse(x == 247,
                                                            "B. halotolerans",
                                                            ifelse(x == 248,
                                                                   "L. fermentum",
                                                                   "S. [unclassified]")))))))))))


fitted_proportions_3fold$Taxon %<>%
  factor(levels = c("P. aeruginosa",
                    "E. coli",
                    "S. enterica",
                    "L. fermentum",
                    "E. faecium",
                    "S. aureus",
                    "L. monocytogenes",
                    "B. halotolerans",
                    "S. [unclassified]",
                    "Other"))

taxon_conversion <- cbind(c("P. aeruginosa",
                            "E. coli",
                            "S. enterica",
                            "E. faecium",
                            "S. aureus",
                            "L. monocytogenes",
                            "B. halotolerans",
                            "L. fermentum",
                            "S. [unclassified]",
                            "Other"),
                          c(
                            "Pseudomonas aeruginosa",
                            "Escherichia-Shigella coli",
                            "Salmonella enterica",
                            "Enterococcus faecium",
                            "Staphylococcus aureus",
                            "Listeria monocytogenes",
                            "Bacillus halotolerans",
                            "Limosilactobacillus (species unclassified)",
                            "Salmonella (species unclassified)",
                            "Other"))



observed_props$Taxon <- sapply(observed_props$Taxon,
                               function(d) taxon_conversion[taxon_conversion[,2]==d,1])

observed_props <- cbind(observed_props[,1:3],
                        data.frame("Fold" =
                                     sapply(observed_props$Type,
                                            function(x)
                                              ifelse(x == "Observed",
                                                     x,
                                                     "Theoretical \nComposition"))),

                        observed_props[,4])

fitted_props_3fold <- data.frame("Dilution" = "Not Applicable",
                                 "Taxon" = fitted_proportions_3fold$Taxon,
                                 "Type" = "Fitted",
                                 "Fold" = fitted_proportions_3fold$Fold,
                                 "Proportion" = fitted_proportions_3fold$value)

results_taxa_lumped_3fold <- rbind(observed_props,
                                   fitted_props_3fold)

results_taxa_lumped_3fold$Taxon %<>%
  factor(levels = c("P. aeruginosa",
                    "E. coli",
                    "S. enterica",
                    "L. fermentum",
                    "E. faecium",
                    "S. aureus",
                    "L. monocytogenes",
                    "B. halotolerans",
                    "S. [unclassified]",
                    "Other"))

lumped_taxa <- unique(results_taxa_lumped_3fold$Taxon)
names(lumped_taxa) <- c("P. aeruginosa",
                        "E. coli",
                        "S. enterica",
                        "E. faecium",
                        "S. aureus",
                        "L. monocytogenes",
                        "B. halotolerans",
                        "L. fermentum",
                        "S. [unclassified]",
                        "Other")
results_taxa_lumped_3fold$nice_taxon <- factor(sapply(results_taxa_lumped_3fold$Taxon,
                                                function(x) names(lumped_taxa)[lumped_taxa ==x] ),
                                         levels = c("P. aeruginosa",
                                                    "E. coli",
                                                    "S. enterica",
                                                    "L. fermentum",
                                                    "E. faecium",
                                                    "S. aureus",
                                                    "L. monocytogenes",
                                                    "B. halotolerans",
                                                    "S. [unclassified]",
                                                    "Other"))

# results_taxa_lumped_3fold %>% View()
# observed_proportions <-


fitted_prop_plot <-
  results_taxa_lumped_3fold %>%
  mutate(Fold = sapply(Fold, function(x) ifelse(
    x%in% 1:3,
    paste("Fold ",x,
                                               sep = "",
                                               collapse = ""),
    x))) %>%
  # filter(Fitted=="Observed") %>%
  filter(Fold != "Observed") %>%
  # filter(Type %in% c("Fitted","Theoretical")) %>%
  ggplot() +
  geom_point(aes(x= Taxon, y = Proportion),
             size = .5) +
  geom_line(aes(x= Taxon, y = Proportion,
                group = interaction(Dilution,Fold),
                linetype = Type,
                color = Fold))+
  theme_bw() +
  guides(linetype = guide_legend(ncol = 1,
                                 keyheight=0.15,
                                 default.unit="inch"),
         color = guide_legend(ncol = 1,
                              keyheight=0.1,
                              default.unit="inch")) +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),

        legend.spacing.y = unit(0.25, 'cm')) +
  scale_color_manual(values = sequential_hcl(4, palette = "LightGrays")[4:1]) +
  guides(linetype=guide_legend(title="Type")) +
  xlab("Taxon") +
  ylab("Proportion") +
  ylim(c(0,.8))


observed_prop_plot <-
  results_taxa_lumped_3fold %>%
  # filter(Fitted=="Observed") %>%
  filter(!(Fold %in% 1:3)) %>%
  # filter(Type %in% c("Fitted","Theoretical")) %>%
  ggplot() +
  geom_point(aes(x= Taxon, y = Proportion),
             size = .5) +
  geom_line(aes(x= Taxon, y = Proportion,
                group = interaction(Dilution,Fold),
                linetype = Type,
                color = Dilution))+
  theme_bw() +
  guides(linetype = guide_legend(ncol = 1,
                                 keyheight=0.15,
                                 default.unit="inch"),
         color = guide_legend(ncol = 1,
                              keyheight=0.1,
                              default.unit="inch")) +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),

        legend.spacing.y = unit(0.25, 'cm')) +
  scale_color_manual(values = sequential_hcl(12, palette = "LightGrays")[c(1:9,1)]) +
  guides(linetype=guide_legend(title="")) +
  xlab("Taxon") +
  ylab("Proportion") +
  ylim(c(0,.8))

###### Log Ratio Target Taxon reads to Contaminant reads ######

for_lr_plot <- observed_read_proportions %>%
  group_by(Dilution) %>%
  summarize(specimen_total = sum(Read_Proportions[!(Taxon %in% c("Salmonella (species unclassified)",
                                                                 "Other"))]),
            contam_total = sum(Read_Proportions[(Taxon %in% c(
              "Other"))]),
            unclass_salm_total = sum(Read_Proportions[(Taxon %in% c("Salmonella (species unclassified)"))])) %>%
  filter(Dilution != "Theoretical")

lr_plot_intercept <- optim(0,
                           function(x) with(for_lr_plot,
                                            sum((log(contam_total/specimen_total) -
                                                   log(3)*as.numeric(Dilution) - x)^2)),
                           method = "Brent",
                           lower = -100,
                           upper = 100)
lr_dil_plot <- for_lr_plot %>%
  mutate(linear_fit = as.numeric(Dilution)*log(3) + lr_plot_intercept$par,
         arbitrary_group = "a",
         legend_text = "Theoretical Slope \n(Intercept Fit by \nLeast Squares)") %>%
  ggplot() +
  geom_point(aes(x = Dilution, y = log(contam_total/specimen_total))) +
  # geom_abline(aes(intercept = -7.507174, slope = log(3)), linetype = 2) +
  geom_line(aes(x = Dilution, y = linear_fit, group = arbitrary_group, linetype = legend_text)) +
  scale_linetype_manual(values = 2,name = "") +
  ylab("Log Ratio Contaminant to \nSpecimen Reads") +
  xlab("Number Three-fold Dilutions") +
  # coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom")

###### Log Read Ratios (L. fermentum denominator) Vs. Dilution ######
lumped_taxa <- unique(observed_read_proportions$Taxon)[-(1:2)]
names(lumped_taxa) <- c("P. aeruginosa",
                        "E. coli",
                        "S. enterica",
                        "E. faecium",
                        "S. aureus",
                        "L. monocytogenes",
                        "B. halotolerans",
                        "L. fermentum")

lrs <- observed_read_proportions %>%
  filter(!(Taxon %in% c("Other","Salmonella (species unclassified)"))) %>%
  filter(Dilution != "Theoretical") %>%
  group_by(Dilution) %>%
  mutate(lr_read_proportions = log(Read_Proportions) -
           log(Read_Proportions[Taxon == "Limosilactobacillus (species unclassified)"])) %>%
  ungroup %>%
  group_by(Taxon) %>%
  mutate(mean_lr = mean(lr_read_proportions)) %>%
  ungroup
lr_tax_plot <- lrs %>%
  filter(Taxon != "Limosilactobacillus (species unclassified)") %>%
  group_by(Taxon, Dilution) %>%
  mutate(Taxon = names(lumped_taxa[lumped_taxa == Taxon])) %>%
  ggplot(aes(x = Dilution, y = lr_read_proportions, group = Taxon,linetype = Taxon)) +
  geom_point() +
  geom_line() +
  # geom_smooth(alpha = .3) +
  # geom_line(aes(x = Dilution, y = mean_lr, group = Taxon, color = Taxon),
  #           linetype = 2) +
  scale_color_brewer(palette = "Dark2") +
  ylab("Log Read Ratio \n(Denominator Taxon L. fermentum)") +
  xlab("Number Three-fold Dilutions") +
  # facet_grid(Taxon~.)+
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0))


pdf("karstens_figure.pdf", width = 7,
    height = 4)
grid.arrange(observed_prop_plot, lr_dil_plot,
             fitted_prop_plot, lr_tax_plot,nrow = 2 )
dev.off()


###### summaries for paper #######

### RMSE
cv_rmse <- fitted_props_3fold %>%
  mutate(Theoretical = sapply(Taxon,
                              function(x)
                                if(x %in% c("S. [unclassified]","Other")){
                                  0
                                }else{0.125})) %>%
  group_by(Fold) %>%
  summarize(RMSE = sqrt(mean((Proportion - Theoretical)^2)))

mean(cv_rmse$RMSE)

plugins <- lapply(1:9, function(x) W[x,]/sum(W[x,]))

pooled <- apply(W,2,sum)/sum(W)

plugin_rmse <- sapply(plugins,
                      function(x) sqrt(mean((x - c(rep(0,240),rep(0.125,8)))^2)))

median(plugin_rmse)
range(plugin_rmse)

pooled_mse <- sqrt(mean((pooled - c(rep(0,240),rep(0.125,8)))^2))




cv_nonzero <- fitted_props_3fold %>%
  mutate(Theoretical = sapply(Taxon,
                              function(x)
                                if(x %in% c("S. [unclassified]","Other")){
                                  0
                                }else{0.125})) %>%
  group_by(Fold) %>%
  summarize(nonzero = sum(Proportion[Theoretical == 0] !=0))

mean(cv_nonzero$nonzero)
plugin_nonzero <- sapply(plugins,
                         function(x)
                           sum(x[1:240]!=0)
                         )

mean(plugin_nonzero)
range(plugin_nonzero)

################################## Full Karstens Without Betas ######################################

B[] <- 0
B_fixed_indices[] <- TRUE
P <- matrix(1/248,nrow = 1,ncol = 248)
P_fixed_indices <- matrix(FALSE,ncol = 248, nrow = 1)
no_eff_karstens_model <-
  estimate_parameters(W = W,
                      X = X,
                      Z = Z,
                      Z_tilde = Z_tilde,
                      Z_tilde_gamma_cols = Z_tilde_gamma_cols,
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
                      alpha_tilde = alpha_tilde,
                      Z_tilde_list = Z_tilde_list,
                      barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                      barrier_scale = 10, #increments for value of barrier penalty
                      max_barrier = 1e12, #maximum value of barrier_t
                      initial_conv_tol = 1000,
                      final_conv_tol = 0.1,
                      final_f = 1e-6,
                      constraint_tolerance = 1e-10,
                      hessian_regularization = 0.01,
                      criterion = "reweighted_Poisson",
                      subproblem_method = "Newton",
                      profile_P = TRUE,
                      profiling_maxit = 25,
                      verbose = TRUE,
                      wts = NULL)

saveRDS(no_eff_karstens_model,"karstens_model_no_betas")

set.seed(432343)
no_eff_cis <- bootstrap_ci(W = W,
                           fitted_model = no_eff_karstens_model,
                           n_boot = 1000,
                           parallelize = TRUE,
                           ncores = 5,
                           verbose = TRUE)

# saveRDS(no_eff_cis,"karstens_cis_no_betas")
no_eff_cis <- readRDS("karstens_cis_no_betas")

no_eff_cis$ci %>%
  filter(param == "P") %>%
  mutate(covered = as.numeric(lower_ci ==0)) %>%
  mutate(covered = sapply(covered, function(x)
    ifelse(x == 1, "Yes","No"))) %>%
  ggplot() +
  geom_errorbar(aes(x = j, ymin = lower_ci,ymax = upper_ci,color =
                      as.factor(covered)),
                width = .1) +
  scale_y_log10() +
  scale_color_manual(values = c("black","gray")) +
  labs(color = "Zero Included in
Marginal 95% CI")+
  xlab("Taxon Index") +
  ylab("Estimated Relative Abundance \n(Up to Detection Effects)")+
  theme_bw()

no_eff_cis$ci %>%
  filter(param == "P") %>%
  mutate(covered = as.numeric(lower_ci ==0)) %>%
  filter(j<241) %>%
  with(c(mean(covered),
       sum(covered),
       which(!covered)))


#### For simulation #####

log_read_depths <- log(apply(W[,241:248],1,sum))
dilutions <- 0:8
truncated_dilutions <- dilutions[5:9]
truncated_log_read_depths <- log_read_depths[5:9]
plot(dilutions, log_read_depths,ylim = c(9,15))
lines(lowess(log_read_depths~dilutions), col = "red")
trunc_lm <- lm(truncated_log_read_depths~truncated_dilutions)
lines(truncated_dilutions, predict(trunc_lm))
lines(dilutions[0:5],rep(mean(log_read_depths[0:4]),5))
### Some attempt at dispersion:
var(log_read_depths[0:4])


means <- meaninate(gammas = full_karstens$gammas,
                   B = full_karstens$B,
                   P = full_karstens$P,
                   X = full_karstens$X,
                   Z = full_karstens$Z,
                   X_tilde = full_karstens$X_tilde,
                   Z_tilde = full_karstens$Z_tilde,
                   Z_tilde_gamma_cols = 1,
                   P_tilde = full_karstens$P_tilde,
                   gamma_tilde = full_karstens$gamma_tilde,
                   alpha_tilde = NULL,
                   Z_tilde_list = NULL,
                   return_separate = FALSE)


sq_errors <- as.numeric((W - means)^2)
means <- as.numeric(means)
sq_errors <- sq_errors[]
means_sq <- means^2
sq_errors_minus_means <- sq_errors - means
moderately_silly_model <- lm(sq_errors_minus_means~means_sq - 1 )
floor(1/moderately_silly_model$coef)
#some mild sanity checks -- i.e. that this is not entirely silly
par(mfrow = c(1,2))
plot(log(means), log(sq_errors), pch = ".")
lines(log(means)[order(means)], log(predict(moderately_silly_model) + means)[order(means)],col = "red")
plot(sqrt(means), sqrt(sq_errors),pch = ".")
lines(sqrt(means[order(means)]), sqrt((predict(moderately_silly_model) + means)[order(means)]),col = "red")

