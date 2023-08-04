library(tidyverse)
library(mvgimme)
library(cowplot)
source("scripts/aux_funs.R")
set.seed(35037)

data <- mvgimme::simData

# multicore
mv_res <- multiverse.gimme(data = data,
                 groupcutoffs = c(.5, .6, .7),
                 subcutoffs = c(.5, .6, .7),
                 n.cores = 10,
                 subgroup = TRUE,
                 save_output = TRUE,
                 save_dir = "output/test",
                 prune_output = TRUE)

ref_fit <- mvgimme::gimme(data = data, 
                 groupcutoff = .75,
                 subcutoff   = .75,
                 subgroup = TRUE)

# saveRDS(mv_res, "mv_res.RDS")
# saveRDS(ref_fit, "ref_fit.RDS")
mv_res <- readRDS("mv_res.RDS")
ref_fit <- readRDS("ref_fit.RDS")



test_combs <- expand.grid(
  groupcutoffs = c(.5, .6, .7),
  subcutoffs = .51,
  rmsea.cuts = .05,
  srmr.cuts = .05,
  nnfi.cuts = .95,
  cfi.cuts = .95,
  n.excellent = 2
)


test_all <- multiverse.compare(l_res = mv_res, 
                               ref_model = ref_fit)

test_group <- multiverse.compare.group(l_res = mv_res,
                                 ref_model = ref_fit)


test_subgroup <- multiverse.compare.subgroup(l_res = mv_res,
                                             ref_model = ref_fit)

test_individual <- multiverse.compare.individual(l_res = mv_res,
                                                 ref_model = ref_fit)





# New try in data generation

# contemp.
a_mat <- matrix(data = runif(n = 36, min = .15, max = .5), nrow = 6, ncol = 6)
a_prob <- matrix(rbinom(36, size = 1, 0.2), nrow = 6, ncol = 6)
a_mat <- a_mat * a_prob
diag(a_mat) <- rep(0, 6)

# lagged
phi_mat <- matrix(data = runif(n = 36), nrow = 6, ncol = 6)
phi_prob <- matrix(rbinom(36, size = 1, 0.2), nrow = 6, ncol = 6)
phi_mat <- phi_mat * phi_prob
diag(phi_mat) <- rep(.3, 6)

psi_mat <- matrix(data = rep(0, 36), nrow = 6, ncol = 6)
diag(psi_mat) <- rep(1, 6)
# psi_mat <- as.matrix(Matrix::forceSymmetric(psi_mat))
group_ind <- c(rep(1, 10), rep(2, 10))
n_ind <- length(group_ind)
n_tp <- 150

data <- simulateVARtest(A = a_mat, 
                   Phi = phi_mat,
                   # Psi = psi_mat,
                   # subAssign = group_ind,
                   N = n_ind,
                   Obs = n_tp)



# With subgroups of equal size


# change position of one effect, add small noise to other effects
# for contemporaneous matrix
a_1 <- a_mat
a_2 <- a_1
a_2[5,3] <- a_2[6,3]
a_2[6,3] <- 0
nonzero <- which(a_2 != 0)
a_2[nonzero] <- a_2[nonzero] + rnorm(length(nonzero), mean = 0, sd = .05)

a_list <- list(a_1, a_2)


# for temporal matrix
phi_1 <- phi_mat
phi_2 <- phi_1
phi_2[5,1] <- phi_2[6,1]
phi_2[6,1] <- 0
nonzero <- which(phi_2 != 0)
phi_2[nonzero] <- phi_2[nonzero] + rnorm(length(nonzero), mean = 0, sd = .05)

phi_list <- list(phi_1, phi_2)

group_ind <- c(rep(1, 10), rep(2, 10))
n_ind <- length(group_ind)
n_tp <- 150


psi_list <- list(psi_mat, psi_mat)
data_sub <- simulateVARtest(A = a_list, 
                                    Phi = phi_list,
                                    Psi = psi_list,
                                    subAssign = group_ind,
                                    N = n_ind,
                                    Obs = n_tp)


# Test fitting a model on it
test_res <- mvgimme::gimme(data_sub$dataList, 
               subgroup = TRUE,
               ar = TRUE)



# -------------------------------------------------------------------------
# Visualizations ----------------------------------------------------------
# -------------------------------------------------------------------------



# Group-Level Viz ---------------------------------------------------------
#--- Number of Group-Level Edges


#--- Which edges are present/absent


#--- Heterogeneity
# Specification Curve Style
# order results by their heterogeneity
# create specification plot which colors based on low, mid, high
# inspired by https://github.com/masurp/specr/blob/33b12bd005c44076711ffd0fb71774703f28ab3c/R/plot_choices.r


# create bogus dataframe
group_cuts <- c(.50, .60, .75, .80)
sub_cuts <- c(.50, .60, .75, .80)
rmsea_cuts <- c(.03, .05, .08)
srmr_cuts <- c(.03, .05, .08)
nnfi_cuts <- c(.90, .95, .97)
cfi_cuts <- c(.90, .95, .97)
n_excels <- c(1, 2, 3)

df_mv <- expand.grid(
  groupcutoffs = group_cuts,
  subcutoffs = sub_cuts,
  rmsea.cuts = rmsea_cuts,
  srmr.cuts = srmr_cuts,
  nnfi.cuts = nnfi_cuts,
  cfi.cuts = cfi_cuts,
  n.excellent = n_excels
)

df_mv$heterogeneity_g <- rnorm(nrow(df_mv), mean = .3, sd = .005)



# Color palette
library(RColorBrewer)
palette_3 <- colorRampPalette(brewer.pal(9, "RdBu"))(3)
palette_4 <- colorRampPalette(brewer.pal(9, "RdBu"))(4)

# TODO change blue color values

# Set values explicitly 
palette_full <- c(palette_4[1:2], palette_3[3], palette_4[3:4])
names(palette_full) <- c("liberal", "medium-liberal", "medium", "medium-strict", "strict")

# Define the desired levels for each factor
group_levels <- c("liberal", "medium-liberal", "medium-strict", "strict")
sub_levels <- c("liberal", "medium-liberal", "medium-strict", "strict")
rmsea_levels <- c("strict", "medium", "liberal")
srmr_levels <- c("strict", "medium", "liberal")
nnfi_levels <- c("liberal", "medium", "strict")
cfi_levels <- c("liberal", "medium", "strict")
n_excels_levels <- c("strict", "medium", "liberal")



# For lineplot
plot_outcome <- function(mv_res, 
                         var,
                         specs = NULL){   # hard-coded for nor
  
  var <- enquo(var)
  
  mv_res %>% 
    dplyr::mutate(variable = as.numeric(!!var)) %>% 
    dplyr::arrange(variable) %>% 
    dplyr::mutate(iteration = dplyr::row_number()) %>%
    ggplot(aes(x = .data$iteration,
               y = variable)) + 
    geom_point(size = 0.8)+
    theme_multiverse()
}

# For Specification Plot
plot_specification <- function(mv_res,
                               var,
                               specs = NULL){    # hard-coded for now
    
  var <- enquo(var)
  
  mv_res %>% 
    dplyr::mutate(variable = as.numeric(!!var)) %>% 
    dplyr::arrange(variable) %>% 
    dplyr::mutate(iteration = dplyr::row_number()) %>%
    dplyr::mutate(across(c(groupcutoffs, subcutoffs,
                           rmsea.cuts, srmr.cuts,
                           cfi.cuts, nnfi.cuts,
                           n.excellent), ~ as.factor(.))) %>% 
    dplyr::mutate(groupcutoffs = fct_recode(groupcutoffs, !!!setNames(as.character(group_cuts), group_levels)),
                  subcutoffs = fct_recode(subcutoffs, !!!setNames(as.character(sub_cuts), sub_levels)),
                  rmsea.cuts = fct_recode(rmsea.cuts, !!!setNames(as.character(rmsea_cuts), rmsea_levels)),
                  srmr.cuts = fct_recode(srmr.cuts, !!!setNames(as.character(srmr_cuts), srmr_levels)),
                  cfi.cuts = fct_recode(cfi.cuts, !!!setNames(as.character(cfi_cuts), cfi_levels)),
                  nnfi.cuts = fct_recode(nnfi.cuts, !!!setNames(as.character(nnfi_cuts), nnfi_levels)),
                  n.excellent = fct_recode(n.excellent, !!!setNames(as.character(n_excels), n_excels_levels))) %>% 
    tidyr::pivot_longer(cols = c(groupcutoffs,
                                 subcutoffs,
                                 rmsea.cuts,
                                 srmr.cuts,
                                 cfi.cuts,
                                 nnfi.cuts,
                                 n.excellent),
                        values_to = "value", names_to = "specification") %>% 
  ggplot(aes(x = .data$iteration,
             y = 1,
             color = .data$value)) + 
    geom_point(shape = 124, size = 15
               #pch='.'   #for faster plotting
                )+
    scale_y_continuous(limits = c(0.99, 1.01), expand = c(0,0))+
    theme_multiverse()+
    scale_color_manual(values = palette_full)+
    facet_wrap(specification~., 
               ncol = 1, 
               strip.position = "left")+
    labs(y = "",
         color = "Specification")+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.text = element_text(size = rel(1.2)))
}

# for testing
df_mv_test <- dplyr::slice_sample(df_mv, n = 400)

test_out <- plot_outcome(mv_res = df_mv, var = heterogeneity_g)
# test_out

test_spec <- plot_specification(mv_res = df_mv, var = heterogeneity_g)
# test_spec
# ggsave("test.svg", test_spec, device = "svg", height = 11, width = 11)


# Combine Plots
comb_plot <- cowplot::plot_grid(test_out, test_spec, 
                   nrow = 2, align = "h",
                   rel_heights = c(0.5,1))
ggsave("comb_plot.svg", comb_plot, device = "svg", height = 11, width = 11)


# Subgroup-Level ----------------------------------------------------------
# Number of subgroups



# Visualize ARI/VI/Modularity across datasets



# Individual-Level --------------------------------------------------------
# Difference in edge weights
# we could also do this with specification curve analysis?


# difference in network summaries



# difference in fit measures




