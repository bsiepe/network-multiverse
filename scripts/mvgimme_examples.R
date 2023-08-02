library(tidyverse)
library(mvgimme)
library(cowplot)
source("scripts/aux_funs.R")

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
    geom_point()+
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
               #pch='.'
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

test_out <- plot_outcome(mv_res = df_mv_test, var = heterogeneity_g)
test_out

test_spec <- plot_specification(mv_res = df_mv_test, var = heterogeneity_g)
test_spec
ggsave("test.svg", test_spec, device = "svg", height = 11, width = 11)


# Combine Plots
cowplot::plot_grid(test_out, test_spec, 
                   nrow = 2, align = "h",
                   rel_heights = c(0.5,1))



# Subgroup-Level ----------------------------------------------------------




# Individual-Level --------------------------------------------------------









