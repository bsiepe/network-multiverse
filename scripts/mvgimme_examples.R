library(tidyverse)
library(mvgimme)
library(cowplot)
library(tictoc)
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


# Create longer mv_res for testing
# Create an empty list to store the duplicated elements
mv_long <- list()

# Loop through each element in the initial list
mv_long <- c(rep(mv_res, 100))

test_combs <- expand.grid(
  groupcutoffs = c(.5, .6, .7),
  subcutoffs = .51,
  rmsea.cuts = .05,
  srmr.cuts = .05,
  nnfi.cuts = .95,
  cfi.cuts = .95,
  n.excellent = 2
)


tic()
test_all <- multiverse.compare(l_res = mv_long, 
                               ref_model = ref_fit)
toc()


test_group <- multiverse.compare.group(l_res = mv_res,
                                 ref_model = ref_fit)


test_subgroup <- multiverse.compare.subgroup(l_res = mv_res,
                                             ref_model = ref_fit)

test_individual <- multiverse.compare.individual(l_res = mv_res,
                                                 ref_model = ref_fit)





# New try in data generation

# contemp.
# a_mat <- matrix(data = runif(n = 36, min = .15, max = .5), nrow = 6, ncol = 6)
# a_prob <- matrix(rbinom(36, size = 1, 0.2), nrow = 6, ncol = 6)
# a_mat <- a_mat * a_prob
# diag(a_mat) <- rep(0, 6)

# manually
a_mat <- matrix(c(
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  0.00, 0.00, 0.25, 0.40, 0.00, 0.00,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.15,
  0.00, 0.00, 0.35, 0.00, 0.00, 0.00
), nrow = 6, byrow = TRUE)


# Write to supplemental material (fn from https://www.r-bloggers.com/2020/08/matrix-to-latex/)
array_to_LaTeX <- function(arr){
  rows <- apply(arr, MARGIN=1, paste, collapse = " & ")
  matrix_string <- paste(rows, collapse = " \\\\ ")
  return(paste("\\begin{bmatrix}", matrix_string, "\\end{bmatrix}"))
  }

# lagged
# phi_mat <- matrix(data = runif(n = 36, min = .05, max = .4), nrow = 6, ncol = 6)
# 
# phi_prob <- matrix(rbinom(36, size = 1, 0.2), nrow = 6, ncol = 6)
# phi_mat <- phi_mat * phi_prob
# diag(phi_mat) <- rep(.3, 6)

ar_e <- 0.3
# manually
phi_mat <- matrix(c(
  ar_e, 0.00, 0.00, 0.00, 0.00, 0.00,
  0.00, ar_e, 0.00, 0.00, 0.10, 0.00,
  0.00, 0.00, ar_e, 0.00, 0.00, 0.20,
  0.00, 0.00, 0.00, ar_e, 0.00, 0.00,
  0.00, 0.00, 0.00, 0.00, ar_e, 0.00,
  0.00, 0.00, 0.00, 0.00, 0.00, ar_e
), nrow = 6, byrow = TRUE)


psi_mat <- matrix(data = rep(0, 36), nrow = 6, ncol = 6)
diag(psi_mat) <- rep(1, 6)

# With subgroups of equal size
# change position of one effect, add small noise to other effects
# for contemporaneous matrix
a_1 <- a_mat
a_2 <- a_1
a_2[5,3] <- a_2[6,3]
a_2[6,3] <- 0
nonzero <- which(a_2 != 0)
a_2[nonzero] <- a_2[nonzero] + rnorm(length(nonzero), mean = 0, sd = .1)
a_2 

# Create a fixed change matrix
change_a <- matrix(c(
  0, 0,  0.00000000,  0.00000000, 0,  0.00000000,
  0, 0, -0.01010176, -0.06727818, 0,  0.00000000,
  0, 0,  0.00000000,  0.00000000, 0,  0.00000000,
  0, 0,  0.00000000,  0.00000000, 0,  0.00000000,
  0, 0,  0.20420106,  0.00000000, 0, -0.04994893,
  0, 0, -0.35000000,  0.00000000, 0,  0.00000000
), nrow = 6, byrow = TRUE)

a_list <- list(a_1, a_2)


# for temporal matrix
phi_1 <- phi_mat
phi_2 <- phi_1
phi_2[5,1] <- phi_2[6,1]
phi_2[6,1] <- 0
nonzero <- which(phi_2 != 0)
phi_2[nonzero] <- phi_2[nonzero] + rnorm(length(nonzero), mean = 0, sd = .1)

# Create a fixed change matrix
phi_2 - phi_1
change_phi <- matrix(c(
  0.0299806, 0.0000000, 0.0000000,  0.00000000,  0.0000000,  0.00000000,
  0.0000000, 0.1384933, 0.0000000,  0.00000000,  0.2270527,  0.00000000,
  0.0000000, 0.0000000, 0.2266066,  0.00000000,  0.0000000, -0.02793758,
  0.0000000, 0.0000000, 0.0000000, -0.06865572,  0.0000000,  0.00000000,
  0.0000000, 0.0000000, 0.0000000,  0.00000000, -0.1298429,  0.00000000,
  0.0000000, 0.0000000, 0.0000000,  0.00000000,  0.0000000,  0.05377093
), nrow = 6, byrow = TRUE)



phi_list <- list(phi_1, phi_2)

group_ind <- c(rep(1, 20), rep(2, 20))
n_ind <- length(group_ind)
n_tp <- 100


psi_list <- list(psi_mat, psi_mat)
data_sub <- simulateVARtest(A = a_list, 
                                    Phi = phi_list,
                                    Psi = psi_list,
                                    subAssign = group_ind,
                                    N = n_ind,
                                    Obs = n_tp,
                                    indA = 0.3,
                                    indPhi = 0.3)


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





# for testing
df_mv_test <- dplyr::slice_sample(df_mv, n = 400)
test_new <- dplyr::slice_sample(test, n = 500)

test_out <- plot_outcome(mv_res = df_mv, var = heterogeneity_g)
# test_out
test_out <- plot_outcome(mv_res = test_new, var = mean_diff_dens_temp_i) 


test_spec <- plot_specification(mv_res = df_mv, var = heterogeneity_g)
# test_spec
test_spec <- plot_specification(mv_res = test_new, var = mean_diff_dens_temp_i)

# ggsave("test.svg", test_spec, device = "svg", height = 11, width = 11)


# Combine Plots
comb_plot <- cowplot::plot_grid(test_out, test_spec, 
                   nrow = 2, align = "h",
                   rel_heights = c(0.5,1))
ggsave("comb_plot.svg", comb_plot, device = "svg", height = 16, width = 16)


# Subgroup-Level ----------------------------------------------------------
# Number of subgroups



# Visualize ARI/VI/Modularity across datasets



# Individual-Level --------------------------------------------------------
# Difference in edge weights
# we could also do this with specification curve analysis?


# difference in network summaries



# difference in fit measures


# Difference in most central edge
# Could do lollipop plot, two facets
# edge in reference model colored 






