################ Data load-in ################
library(tidyverse)
library(devtools)
library(data.table)
library(magrittr)
library(Matrix)
# library(sparseinv)
# install_github("https://github.com/ailurophilia/logsum")
# install_github("https://github.com/ailurophilia/fastnnls")
library(logsum)
library(fastnnls)
library(scales)
load_all()

brooks_counts <- read.csv("brooks_counts.csv")
brooks_truth <- read.csv("brooks_truth.csv")


brooks_counts$Plate <- sapply(brooks_counts$Sample,
                              function(x) substr(x, 6,6))




brooks_counts$Barcode <- sapply(as.character(brooks_counts$Sample),
                                function(x) strsplit(x, "_", fixed = TRUE)[[1]][2])

brooks_truth$convenience_id <- apply(brooks_truth, 1,
                                     function(x) paste("sample",
                                                       as.numeric(x["Plate"]),
                                                       as.numeric(x["Barcode"]),
                                                       sep  = "_",
                                                       collapse = "_"))

### rename brooks_truth cols so they match brooks_counts
colnames(brooks_truth)[colnames(brooks_truth) == "Lcrispatus"] <-
  "Lactobacillus.crispatus_cluster"

colnames(brooks_truth)[colnames(brooks_truth) == "Liners"] <-
  "Lactobacillus.iners"

colnames(brooks_truth)[colnames(brooks_truth) == "Gvaginalis"] <-
  "Gardnerella.vaginalis"

colnames(brooks_truth)[colnames(brooks_truth) == "Avaginae"] <-
  "Atopobium.vaginae"

colnames(brooks_truth)[colnames(brooks_truth) == "Pbivia"] <-
  "Prevotella.bivia"

colnames(brooks_truth)[colnames(brooks_truth) == "Samnii"] <-
  "Sneathia.amnii"

colnames(brooks_truth)[colnames(brooks_truth) == "GroupBStrep"] <-
  "Streptococcus.agalactiae"


brooks_counts$convenience_id <- apply(brooks_counts, 1,
                                      function(x) paste("sample",
                                                        x["Plate"],
                                                        x["Barcode"],
                                                        sep = "_",
                                                        collapse = "_"))

######## Brooks Plates 1 and 2 ########
# get counts
brooks_counts_experiment1 <- brooks_counts %>% filter(Plate %in% c(1,2))
# get true (or at least intended) composition
brooks_truth_experiment1 <- brooks_truth %>% filter(Plate %in% c(1,2))


taxon_counts <- brooks_counts_experiment1
taxon_counts <- taxon_counts[,colnames(taxon_counts) %in% colnames(brooks_truth_experiment1)]

W <- taxon_counts %>% dplyr::select(Atopobium.vaginae,
                                    Gardnerella.vaginalis,
                                    Lactobacillus.crispatus_cluster,
                                    Lactobacillus.iners,
                                    Prevotella.bivia,
                                    Sneathia.amnii,
                                    Streptococcus.agalactiae)

saveRDS(W,"brooks_W")




truth_compressed <- round(brooks_truth_experiment1[,5:11], 2) %>%
  apply(1, function(x) paste(x, sep = "",collapse = "_"))
specimens <- character(0)
for(i in 1:80){
  if(truth_compressed[i] %in% specimens){

  } else{
    specimens <- c(specimens,truth_compressed[i])
  }
}


samples_to_specimens <- lapply(1:80, function(i) as.numeric(truth_compressed[i] == specimens)) %>%
  (function(x) do.call(rbind,x)) %>% as.matrix

specimen_indices <- apply(samples_to_specimens,2,function(x) min((1:80)[x==1]))

P <- brooks_truth_experiment1[specimen_indices,] %>%
  dplyr::select(Atopobium.vaginae,
                Gardnerella.vaginalis,
                Lactobacillus.crispatus_cluster,
                Lactobacillus.iners,
                Prevotella.bivia,
                Sneathia.amnii,
                Streptococcus.agalactiae) %>%
  as.matrix
saveRDS(P,"brooks_P")


#check that samples_to_specimens is correct
max(abs(samples_to_specimens%*%P - brooks_truth_experiment1[,c("Atopobium.vaginae",
                                                               "Gardnerella.vaginalis",
                                                               "Lactobacillus.crispatus_cluster",
                                                               "Lactobacillus.iners",
                                                               "Prevotella.bivia",
                                                               "Sneathia.amnii",
                                                               "Streptococcus.agalactiae")]))

counts_to_samples <- lapply(1:80, function(i)
  as.numeric(taxon_counts$convenience_id[i] == brooks_truth_experiment1$convenience_id)) %>%
  (function(x) do.call(rbind,x))

rbind(W[1,1:7],
      (counts_to_samples%*%samples_to_specimens%*%P)[1,1:7])
Z <- counts_to_samples%*%samples_to_specimens
saveRDS(Z,"brooks_Z")


get_overlap_matrix <- function(x){
  overlap_matrix <- matrix(0,nrow = length(x),ncol = length(x))
  for(i in 1:length(x)){
    for(j in i:length(x)){
      if(x[i]>0 & x[j]>0){
        overlap_matrix[i,j] <- overlap_matrix[j,i] <- 1
      }
    }
  }
  return(overlap_matrix)
}

get_laplacian <- function(P_subset){
  overlap_list <-
    lapply(1:nrow(P_subset),function(i) get_overlap_matrix(P_subset[i,]))
  overlap_matrix <- overlap_list[[1]]*0
  for(i in 1:length(overlap_list)){
    overlap_matrix <- overlap_matrix + overlap_list[[i]]
  }
  laplacian <- 0*overlap_matrix
  for(i in 1:nrow(laplacian)){
    for(j in i:nrow(laplacian)){
      if(i == j){
        laplacian[i,j] <- sum(overlap_matrix[i,-i]>0)
      } else{
        laplacian[i,j] <- laplacian[j,i] <- ifelse(overlap_matrix[i,j]>0,-1,0)
      }
    }
  }
  return(laplacian)
}

is_connected <- function(P_subset,cutoff = 1e-5){
  laplacian <- get_laplacian(P_subset)
  fiedler <- eigen(laplacian)$value
  fiedler <- fiedler[order(fiedler)][2]
  return(fiedler>1e-5)
}

get_acceptable_sample <- function(n_per_plate,
                      plate_info,
                      P,
                      Z){
  unique_plates <- unique(plate_info)
  sample_indices <- numeric(0)
  for(plate_no in 1:length(unique_plates)){
    connected <- FALSE
    plate_indices <- which(plate_info==unique_plates[plate_no])
    while(!connected){
      potential_Z_indices <- sample(plate_indices,n_per_plate,replace = FALSE)
      potential_P_indices <- sapply(potential_Z_indices,
                                  function(x) which(Z[x,]>0))
      connected <- is_connected(P[potential_P_indices,])
    }
  sample_indices <- c(sample_indices,potential_Z_indices)
  }
  return(sample_indices)
}

do_one_simulation_brooks <- function(n_per_plate,
                              plate_info,
                              P,
                              Z,
                              W,
                              criterion){
known_indices <- get_acceptable_sample(n_per_plate,
                                       plate_info,
                                       P,
                                       Z)
unique_plates <- unique(plate_info)
Z_tilde <- do.call(rbind,lapply(1:length(plate_info),
                  function(i) as.numeric(unique_plates == plate_info[i])))

W_reordered <- rbind(W[known_indices,],
                     W[-known_indices,])

Z_tilde_reordered <- rbind(Z_tilde[known_indices,],
                           Z_tilde[-known_indices,])

Z_tilde_gamma_cols <- 1:ncol(Z_tilde_reordered)

P_known <- lapply(known_indices,
                  function(x) P[which(Z[x,]==1),])
P_known <- do.call(rbind,P_known)

P_unknown <- matrix(1/ncol(P),nrow = nrow(Z) - length(known_indices),
                    ncol = ncol(P))

P_for_fit <- rbind(P_known,
                   P_unknown)
P_fixed_indices <- matrix(FALSE,nrow = nrow(P_for_fit),
                          ncol = ncol(P_for_fit))

P_fixed_indices[1:length(known_indices),] <-TRUE

Z_for_fit <- diag(nrow(Z))

temp_model <- estimate_parameters(W = W_reordered,
                    X = matrix(1,nrow = nrow(W),ncol = 1),
                    Z = Z_for_fit,
                    Z_tilde = Z_tilde_reordered,
                    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                    gammas = apply(W_reordered,1,function(x) log(sum(x))),
                    gammas_fixed_indices = rep(F,nrow(W_reordered)),
                    P = P_for_fit,
                    P_fixed_indices = P_fixed_indices,
                    B = matrix(0,ncol = ncol(W),nrow = 1),
                    B_fixed_indices = matrix(c(rep(FALSE,ncol(W) - 1),
                                               TRUE),
                                             nrow = 1),
                    X_tilde = matrix(0,nrow = 2,ncol = 1),
                    P_tilde = matrix(1/ncol(W),nrow = 2,ncol= ncol(W)),
                    P_tilde_fixed_indices = matrix(FALSE, ncol = ncol(W),
                                                   nrow = 2),
                    gamma_tilde = matrix(rep(0,2),ncol = 1),
                    gamma_tilde_fixed_indices = matrix(rep(FALSE,2), ncol = 1),
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
                    criterion = criterion,
                    profile_P = TRUE,
                    verbose = FALSE)

return(list(Z_for_original_P = Z[-known_indices,],
            estimated_P = temp_model$P[-(1:length(known_indices)),],
            known = known_indices))

}

nsims <- 100

sims_3_poisson <- vector(nsims,mode = "list")
sims_5_poisson <- vector(nsims,mode = "list")
sims_10_poisson <- vector(nsims,mode = "list")
sims_20_poisson <- vector(nsims,mode = "list")

sims_3_reweighted <- vector(nsims,mode = "list")
sims_5_reweighted <- vector(nsims,mode = "list")
sims_10_reweighted <- vector(nsims,mode = "list")
sims_20_reweighted <- vector(nsims,mode = "list")

for(i in 1:nsims){
  print(i)

  print("Poisson 3 known")
  set.seed(i)
  sims_3_poisson[[i]] <-
    do_one_simulation_brooks(n_per_plate = 3,
                      plate_info = brooks_counts$Plate[1:80],
                      P = P,
                      Z = Z,
                      W = W,
                      criterion = "Poisson")
  print("Reweighted Poisson 3 known")
  set.seed(i)
  sims_3_reweighted[[i]] <-
    do_one_simulation_brooks(n_per_plate = 3,
                      plate_info = brooks_counts$Plate[1:80],
                      P = P,
                      Z = Z,
                      W = W,
                      criterion = "reweighted_Poisson")

  print("Poisson 5 known")
  set.seed(i)
 sims_5_poisson[[i]] <-
   do_one_simulation_brooks(n_per_plate = 5,
                     plate_info = brooks_counts$Plate[1:80],
                     P = P,
                     Z = Z,
                     W = W,
                     criterion = "Poisson")
 print("Reweighted Poisson 5 known")
 set.seed(i)
 sims_5_reweighted[[i]] <-
   do_one_simulation_brooks(n_per_plate = 5,
                     plate_info = brooks_counts$Plate[1:80],
                     P = P,
                     Z = Z,
                     W = W,
                     criterion = "reweighted_Poisson")
 print("Poisson 10 known")
 set.seed(i)
 sims_10_poisson[[i]] <-
   do_one_simulation_brooks(n_per_plate = 10,
                     plate_info = brooks_counts$Plate[1:80],
                     P = P,
                     Z = Z,
                     W = W,
                     criterion = "Poisson")
 print("Reweighted Poisson 10 known")
 set.seed(i)
 sims_10_reweighted[[i]] <-
   do_one_simulation_brooks(n_per_plate = 10,
                     plate_info = brooks_counts$Plate[1:80],
                     P = P,
                     Z = Z,
                     W = W,
                     criterion = "reweighted_Poisson")
 print("Poisson 20 known")
 set.seed(i)
 sims_20_poisson[[i]] <-
   do_one_simulation_brooks(n_per_plate = 20,
                     plate_info = brooks_counts$Plate[1:80],
                     P = P,
                     Z = Z,
                     W = W,
                     criterion = "Poisson")
 print("Reweighted Poisson 20 known")
 set.seed(i)
 sims_20_reweighted[[i]] <-
   do_one_simulation_brooks(n_per_plate = 20,
                     plate_info = brooks_counts$Plate[1:80],
                     P = P,
                     Z = Z,
                     W = W,
                     criterion = "reweighted_Poisson")


}

sim_summary_poisson_3 <- lapply(1:length(sims_3_poisson),
       function(x) cbind(rbind(sims_3_poisson[[x]]$Z_for_original_P%*%P,
                         sims_3_poisson[[x]]$estimated_P),
                         "sample" = rep((1:80)[(!(1:80 %in% sims_3_poisson[[x]]$known))],2)) %>%
         data.frame() %>%
         mutate(type = c(rep("Theoretical",nrow(sims_3_poisson[[x]]$Z_for_original_P)),
                rep("Estimated",nrow(sims_3_poisson[[x]]$Z_for_original_P)))) %>%
         mutate(simulation_index = x,
                method = "Poisson",
                n_known = 3))


sim_summary_reweighted_3 <- lapply(1:length(sims_3_reweighted),
                                   function(x) cbind(rbind(sims_3_reweighted[[x]]$Z_for_original_P%*%P,
                                                           sims_3_reweighted[[x]]$estimated_P),
                                                     "sample" = rep((1:80)[(!(1:80 %in% sims_3_reweighted[[x]]$known))],2)) %>%
                                     data.frame() %>%
                                     mutate(type = c(rep("Theoretical",nrow(sims_3_reweighted[[x]]$Z_for_original_P)),
                                                     rep("Estimated",nrow(sims_3_reweighted[[x]]$Z_for_original_P)))) %>%
                                     mutate(simulation_index = x,
                                            method = "reweighted",
                                            n_known = 3))



sim_summary_poisson_5 <- lapply(1:length(sims_5_poisson),
                                function(x) cbind(rbind(sims_5_poisson[[x]]$Z_for_original_P%*%P,
                                                        sims_5_poisson[[x]]$estimated_P),
                                                  "sample" = rep((1:80)[(!(1:80 %in% sims_5_poisson[[x]]$known))],2)) %>%
                                  data.frame() %>%
                                  mutate(type = c(rep("Theoretical",nrow(sims_5_poisson[[x]]$Z_for_original_P)),
                                                  rep("Estimated",nrow(sims_5_poisson[[x]]$Z_for_original_P)))) %>%
                                  mutate(simulation_index = x,
                                         method = "Poisson",
                                         n_known = 5))


sim_summary_reweighted_5 <- lapply(1:length(sims_5_reweighted),
                                   function(x) cbind(rbind(sims_5_reweighted[[x]]$Z_for_original_P%*%P,
                                                           sims_5_reweighted[[x]]$estimated_P),
                                                     "sample" = rep((1:80)[(!(1:80 %in% sims_5_reweighted[[x]]$known))],2)) %>%
                                     data.frame() %>%
                                     mutate(type = c(rep("Theoretical",nrow(sims_5_reweighted[[x]]$Z_for_original_P)),
                                                     rep("Estimated",nrow(sims_5_reweighted[[x]]$Z_for_original_P)))) %>%
                                     mutate(simulation_index = x,
                                            method = "reweighted",
                                            n_known = 5))




sim_summary_poisson_10 <- lapply(1:length(sims_10_poisson),
                                function(x) cbind(rbind(sims_10_poisson[[x]]$Z_for_original_P%*%P,
                                                        sims_10_poisson[[x]]$estimated_P),
                                                  "sample" = rep((1:80)[(!(1:80 %in% sims_10_poisson[[x]]$known))],2)) %>%
                                  data.frame() %>%
                                  mutate(type = c(rep("Theoretical",nrow(sims_10_poisson[[x]]$Z_for_original_P)),
                                                  rep("Estimated",nrow(sims_10_poisson[[x]]$Z_for_original_P)))) %>%
                                  mutate(simulation_index = x,
                                         method = "Poisson",
                                         n_known = 10))


sim_summary_reweighted_10 <- lapply(1:length(sims_10_reweighted),
                                   function(x) cbind(rbind(sims_10_reweighted[[x]]$Z_for_original_P%*%P,
                                                           sims_10_reweighted[[x]]$estimated_P),
                                                     "sample" = rep((1:80)[(!(1:80 %in% sims_10_reweighted[[x]]$known))],2)) %>%
                                     data.frame() %>%
                                     mutate(type = c(rep("Theoretical",nrow(sims_10_reweighted[[x]]$Z_for_original_P)),
                                                     rep("Estimated",nrow(sims_10_reweighted[[x]]$Z_for_original_P)))) %>%
                                     mutate(simulation_index = x,
                                            method = "reweighted",
                                            n_known = 10))


sim_summary_poisson_20 <- lapply(1:length(sims_20_poisson),
                                function(x) cbind(rbind(sims_20_poisson[[x]]$Z_for_original_P%*%P,
                                                        sims_20_poisson[[x]]$estimated_P),
                                                  "sample" = rep((1:80)[(!(1:80 %in% sims_20_poisson[[x]]$known))],2)) %>%
                                  data.frame() %>%
                                  mutate(type = c(rep("Theoretical",nrow(sims_20_poisson[[x]]$Z_for_original_P)),
                                                  rep("Estimated",nrow(sims_20_poisson[[x]]$Z_for_original_P)))) %>%
                                  mutate(simulation_index = x,
                                         method = "Poisson",
                                         n_known = 20))


sim_summary_reweighted_20 <- lapply(1:length(sims_20_reweighted),
                                   function(x) cbind(rbind(sims_20_reweighted[[x]]$Z_for_original_P%*%P,
                                                           sims_20_reweighted[[x]]$estimated_P),
                                                     "sample" = rep((1:80)[(!(1:80 %in% sims_20_reweighted[[x]]$known))],2)) %>%
                                     data.frame() %>%
                                     mutate(type = c(rep("Theoretical",nrow(sims_20_reweighted[[x]]$Z_for_original_P)),
                                                     rep("Estimated",nrow(sims_20_reweighted[[x]]$Z_for_original_P)))) %>%
                                     mutate(simulation_index = x,
                                            method = "reweighted",
                                            n_known = 20))


sim_summary <- do.call(rbind,c(sim_summary_poisson_3,
                               sim_summary_reweighted_3,
                               sim_summary_poisson_5,
                               sim_summary_reweighted_5,
                               sim_summary_poisson_10,
                               sim_summary_reweighted_10,
                               sim_summary_poisson_20,
                               sim_summary_reweighted_20))

# saveRDS(sim_summary, "Brooks_simulation_summary")
sim_summary <- readRDS("Brooks_simulation_summary")
# sim_summary <- readRDS("Brooks_simulation_summary")
# sim_summary %>%
#   pivot_longer(-c(type,sample_index,simulation_index,method,n_known)) %>%
#   pivot_wider(names_from = type,
#               values_from = value,
#               id_cols = c(simulation_index,sample_index,name,
#                           method,n_known)) %>%
#   mutate(name = sapply(name,
#                        function(x)
#                          c("A. vaginae",
#                            "G. vaginalis",
#                            "L. crispatus",
#                            "L. iners",
#                            "P. bivia",
#                            "S. amnii",
#                            "S. agalactiae")[
#                              which(c("Atopobium.vaginae",
#                                      "Gardnerella.vaginalis",
#                                      "Lactobacillus.crispatus_cluster",
#                                      "Lactobacillus.iners",
#                                      "Prevotella.bivia",
#                                      "Sneathia.amnii",
#                                      "Streptococcus.agalactiae") ==x)
#                            ]
#   )) %>%
#   ggplot() +
#   geom_point(aes(x = Theoretical, y= Estimated,color = name),
#              position = position_dodge(.05),
#              size = .5)+
#   geom_abline(aes(intercept = 0, slope = 1),linetype = 2) +
#   theme_bw() +
#   scale_color_brewer(palette = "Dark2") +
#   facet_grid(method~n_known) +
#   scale_y_sqrt() +
#   scale_x_sqrt()

# sim_summary %>%
#   pivot_longer(-c(type,sample_index,simulation_index,method,n_known)) %>%
#   pivot_wider(names_from = method,
#               values_from = value,
#               id_cols = c(simulation_index,sample_index,name,n_known,type)) %>%
#   filter(type == "Estimated") %>%
#   ggplot() +
#   geom_point(aes(x = Poisson, y = Reweighted,color= name),
#              size = .5) +
#   geom_abline(aes(intercept = 0, slope = 1),
#               linetype = 2,
#               alpha = .75) +
#   facet_wrap(~n_known) +
#   scale_y_sqrt() +
#   guides(color = guide_legend(title = "Taxon")) +
#   scale_x_sqrt() +
#   theme_bw() +
#   ylab("Reweighted Estimator") +
#   xlab("Poisson Estimator")

# sim_summary %>%
#   pivot_longer(-c(type,sample_index,simulation_index,method,n_known)) %>%
#   pivot_wider(names_from = c(method,type),
#               values_from = value,
#               id_cols = c(simulation_index,sample_index,name,n_known)) %>%
#   ggplot() +
#   geom_point(aes(x = (Poisson_Estimated - Poisson_Theoretical)^2,
#                  y = (Reweighted_Estimated- Reweighted_Theoretical)^2),
#              size = .5) +
#   geom_smooth(aes(x = (Poisson_Estimated - Poisson_Theoretical)^2,
#                   y = (Reweighted_Estimated- Reweighted_Theoretical)^2),
#              size = .5,
#              se = FALSE) +
#   geom_abline(aes(intercept = 0, slope = 1),
#               linetype = 2,
#               alpha = .75) +
#   facet_wrap(~n_known) +
#   scale_y_sqrt() +
#   scale_x_sqrt()

# sim_summary %>%
#   pivot_longer(-c(type,sample_index,simulation_index,method,n_known)) %>%
#   pivot_wider(names_from = c(method,type),
#               values_from = value,
#               id_cols = c(simulation_index,sample_index,name,n_known)) %>%
#   group_by(n_known,name) %>%
#   summarize(Poisson_mse = mean((Poisson_Estimated - Poisson_Theoretical)^2),
#             Poisson_se = sd((Poisson_Estimated - Poisson_Theoretical)^2),
#             Reweighted_mse = mean((Reweighted_Estimated - Reweighted_Theoretical)^2),
#             Reweighted_se = sd((Reweighted_Estimated - Reweighted_Theoretical)^2)
#             ) %>%
#   ggplot() +
#   geom_point(aes(x = Poisson_mse, y = Reweighted_mse, color = as.factor(n_known))) +
#   geom_abline(aes(intercept = 0, slope = 1),linetype = 2) +
#   theme_bw()

### Calculate plugin estimator RMSE
plugin <- W
for(i in 1:80){
  plugin[i,] <- plugin[i,]/sum(plugin[i,])
}

plugin_rmse <- sqrt(mean((as.matrix(plugin) - Z%*%P)^2))

### Calculate RMSE by n_known and estimator

sim_summary %>%
  pivot_longer(-c(type,sample,simulation_index,method,n_known)) %>%
  pivot_wider(names_from = type,
              values_from = value,
              id_cols = c(simulation_index,sample,name,
                          method,n_known)) %>%
  group_by(method,n_known,sample) %>%
  summarize(mse = mean((Estimated - Theoretical)^2)) %>%
  group_by(method, n_known) %>%
  summarize(rmse = sqrt(mean(mse)))

### Calculate plugin prop_zero
plugin_nonzero <- mean(W[Z%*%P == 0] ==0)

### Calculated estimated prop_zero
zero_overall_summary <-
  sim_summary %>%
  pivot_longer(-c(type,sample,simulation_index,method,n_known)) %>%
  pivot_wider(names_from = type,
              values_from = value,
              id_cols = c(simulation_index,sample,name,
                          method,n_known)) %>%
  mutate(name = sapply(name,
                       function(x)
                         c("A. vaginae",
                           "G. vaginalis",
                           "L. crispatus",
                           "L. iners",
                           "P. bivia",
                           "S. amnii",
                           "S. agalactiae")[
                             which(c("Atopobium.vaginae",
                                     "Gardnerella.vaginalis",
                                     "Lactobacillus.crispatus_cluster",
                                     "Lactobacillus.iners",
                                     "Prevotella.bivia",
                                     "Sneathia.amnii",
                                     "Streptococcus.agalactiae") ==x)
                           ]
  )) %>%
  group_by(Theoretical,method,n_known,sample) %>%
  filter(Theoretical == 0) %>%
  summarize(prop_zero = mean(Estimated == 0)) %>%
  group_by(method, n_known) %>%
  summarize(prop_zero = mean(prop_zero))


#### Plots for Supplement
rmse_plot <- sim_summary %>%
  pivot_longer(-c(type,sample,simulation_index,method,n_known)) %>%
  pivot_wider(names_from = type,
              values_from = value,
              id_cols = c(simulation_index,sample,name,
                          method,n_known)) %>%
  mutate(name = sapply(name,
                       function(x)
                         c("A. vaginae",
                           "G. vaginalis",
                           "L. crispatus",
                           "L. iners",
                           "P. bivia",
                           "S. amnii",
                           "S. agalactiae")[
                             which(c("Atopobium.vaginae",
                                     "Gardnerella.vaginalis",
                                     "Lactobacillus.crispatus_cluster",
                                     "Lactobacillus.iners",
                                     "Prevotella.bivia",
                                     "Sneathia.amnii",
                                     "Streptococcus.agalactiae") ==x)
                           ]
  )) %>%
  group_by(Theoretical,method,n_known,name,sample) %>%
  summarize(mse = mean((Theoretical - Estimated)^2)) %>%
  group_by(Theoretical,method, n_known, name) %>%
  summarize(mse = mean(mse)) %>%
  mutate(method = sapply(method,
                         function(x)
                           ifelse(x == "Poisson",
                                  "Unweighted",
                                  "Reweighted"))) %>%
  ggplot() +
  geom_point(aes(x = 0.5-abs(Theoretical - .5),
                 y = sqrt(mse),
                 color = as.factor(n_known)),
             size = .5,
             position = position_dodge(.03)) +
  geom_line(aes(x = 0.5-abs(Theoretical - .5),
                 y = sqrt(mse),
                group = interaction(as.factor(n_known),
                                    method),
                 color = as.factor(n_known),
                linetype = method),
             position = position_dodge(.03)) +
  facet_grid(.~name)+
  scale_x_continuous(breaks = c(0,.2,.4))+
  guides(color= guide_legend(title = "Number Known \nper Plate"),
         linetype = FALSE) +
  scale_color_manual(values = colorspace::sequential_hcl(6, palette = "LightGrays")[4:1])  +
  scale_linetype_manual(values = c(1,4)) +
  theme_bw() +
  ylab("\n \nRoot Mean Square Error") +
  xlab("Theoretical Abundance")

plugin_nonzero <- sapply(1:7, function(j) mean(W[Z%*%P[,j] == 0,j] ==0))

names(plugin_nonzero) <- colnames(W)

plugin_nonzero <- data.frame(prop_zero = plugin_nonzero,
           name =  c("A. vaginae",
                     "G. vaginalis",
                     "L. crispatus",
                     "L. iners",
                     "P. bivia",
                     "S. amnii",
                     "S. agalactiae"))


  zero_summary <-
  sim_summary %>%
  pivot_longer(-c(type,sample,simulation_index,method,n_known)) %>%
  pivot_wider(names_from = type,
              values_from = value,
              id_cols = c(simulation_index,sample,name,
                          method,n_known)) %>%
  mutate(name = sapply(name,
                       function(x)
                         c("A. vaginae",
                           "G. vaginalis",
                           "L. crispatus",
                           "L. iners",
                           "P. bivia",
                           "S. amnii",
                           "S. agalactiae")[
                             which(c("Atopobium.vaginae",
                                     "Gardnerella.vaginalis",
                                     "Lactobacillus.crispatus_cluster",
                                     "Lactobacillus.iners",
                                     "Prevotella.bivia",
                                     "Sneathia.amnii",
                                     "Streptococcus.agalactiae") ==x)
                           ]
  )) %>%
  group_by(Theoretical,method,n_known,name,sample) %>%
  filter(Theoretical == 0) %>%
  summarize(prop_zero = mean(Estimated == 0)) %>%
    group_by(method, n_known, name) %>%
    summarize(prop_zero = mean(prop_zero)) %>%
    mutate(method = sapply(method,
                           function(x)
                             ifelse(x == "Poisson",
                                    "Unweighted",
                                    "Reweighted")))

zero_plot <- ggplot() +
  geom_point(aes(x = as.factor(n_known), y = prop_zero
                ),
             size = 0.5,
             data = zero_summary) +
  geom_line(aes(x = as.factor(n_known), y = prop_zero,

                group = interaction(name,method),
                linetype = method),
            data = zero_summary) +
  geom_abline(aes(intercept = prop_zero,slope = 0),
              data = plugin_nonzero,
              linetype = 3) +
  scale_y_continuous(
    limits= c(.25,.85))+
  scale_linetype_manual(values = c(1,4)) +
  facet_grid(.~name) +
  guides(linetype = guide_legend(title ="Estimator")) +
  xlab("Number Known per Plate") +
  ylab("Proportion of Samples with \nTaxon Theoretically Absent \nEstimated not to Contain Taxon") +
  theme_bw()

rmse_plot <- egg::set_panel_size(rmse_plot,
                            width  = unit(2.75, "cm"),
                            height = unit(5.5, "cm"))

zero_plot <- egg::set_panel_size(zero_plot,
                                 width  = unit(2.75, "cm"),
                                 height = unit(5.5, "cm"))

gridExtra::grid.arrange(rmse_plot,zero_plot)

sim_summary %>%
  pivot_longer(-c(type,sample,simulation_index,method,n_known)) %>%
  pivot_wider(names_from = type,
              values_from = value,
              id_cols = c(simulation_index,sample,name,
                          method,n_known)) %>%
  mutate(name = sapply(name,
                       function(x)
                         c("A. vaginae",
                           "G. vaginalis",
                           "L. crispatus",
                           "L. iners",
                           "P. bivia",
                           "S. amnii",
                           "S. agalactiae")[
                             which(c("Atopobium.vaginae",
                                     "Gardnerella.vaginalis",
                                     "Lactobacillus.crispatus_cluster",
                                     "Lactobacillus.iners",
                                     "Prevotella.bivia",
                                     "Sneathia.amnii",
                                     "Streptococcus.agalactiae") ==x)
                           ]
  )) %>%
  group_by(method,n_known,sample) %>%
  summarize(mse  = (mean((Theoretical - Estimated)^2))) %>%
  group_by(method, n_known) %>%
  summarize(rmse = sqrt(mean(mse))) %>%
  ggplot() +
  geom_point(aes(x = n_known, y = rmse)) +
  geom_line(aes(x = n_known, y = rmse,
                group = method,
                linetype = method)) +
  # ylim(c(0,.045))+
  xlab("Number Known per Plate") +
  ylab("Root Mean Square Error") +
  theme_bw()




