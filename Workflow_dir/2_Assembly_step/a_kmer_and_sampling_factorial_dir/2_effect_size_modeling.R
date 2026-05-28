# =============================================================================
# effect_size_modeling.R
# Add-on module for transrate_EDA.R
#
# PURPOSE
#   Quantify and visualise the effect of k-mer size vs. sample size (sampling_set)
#   on assembly Accuracy, using the factorial output of `benchmark.sh -m both`.
#
# HOW TO USE
#   Source this AFTER `metricsdf` is created in transrate_EDA.R (~line 235):
#       source("effect_size_modeling.R")
#   Requires `metricsdf` with columns:
#       vfold_set, kmer, sampling_set, TP, FP, FN, Accuracy
#   and `opt$out` (output directory).
#
# DEPENDENCIES (install once)
#   install.packages(c("glmmTMB","car","performance","mgcv","gratia",
#                       "ggeffects","ranger","pdp","patchwork","broom.mixed"))
# =============================================================================


rm(list = ls())

if(!is.null(dev.list())) dev.off()


# ── Theme ────────────────────────────────────────────────────────────────────
my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  # base_family <- if ("GillSans" %in% systemfonts::system_fonts()$family) "GillSans" else ""
  base_family <- "GillSans"
  theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      legend.position   = legend_pos,
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "gray90", color = "white"),
      strip.text        = element_text(angle = 0, size = base_size, hjust = 0),
      axis.text         = element_text(size = rel(0.7), color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      ...
    )
}

source("~/Documents/GitHub/ConotoxinBenchmark/3_kmer_dir/kmer_and_sampling_dir/transrate_EDA.R")

opt$out <-  "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)      # beta-binomial GLMM
  library(car)          # Type-II/III Wald ANOVA
  library(performance)  # r2() for (G)LMMs
  library(mgcv)         # GAM smooth surface
  library(ggeffects)    # tidy marginal effects
  library(ranger)       # random forest
  library(pdp)          # partial dependence
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 0. Prepare the modelling frame
# -----------------------------------------------------------------------------
# kmer arrives as character ("19", "21", ...) -> numeric.
# denom is the binomial "number of trials"; failures = FN + FP.
# We keep BOTH a numeric and a factor version of each predictor: numeric for the
# GLMM/GAM/RF, factor for descriptive interaction plots and for any anova that
# wants discrete levels.

model_df <- metricsdf %>%
  mutate(
    kmer_n  = as.numeric(as.character(kmer)),
    samp_n  = as.numeric(sampling_set),
    vfold   = factor(vfold_set),
    denom   = TP + FN + FP,
    failure = FN + FP
  ) %>%
  filter(denom > 0, !is.na(kmer_n), !is.na(samp_n)) %>%
  # centre + scale numeric predictors so coefficients are standardised
  # (a 1-unit change = 1 SD), which makes them directly comparable as
  # effect sizes and stabilises the interaction term.
  mutate(
    kmer_z = as.numeric(scale(kmer_n)),
    samp_z = as.numeric(scale(samp_n))
  )

message("Modelling frame: ", nrow(model_df), " rows, ",
        n_distinct(model_df$kmer_n), " kmer levels x ",
        n_distinct(model_df$samp_n), " sample levels x ",
        n_distinct(model_df$vfold), " folds")

# =============================================================================
# 1. Beta-binomial GLMM  (primary inferential model)
# =============================================================================
# Accuracy is a proportion of a known count -> binomial-type response.
# betabinomial family absorbs overdispersion. Random intercept (1|vfold)
# accounts for the repeated-measures resampling design.

glmm_full <- glmmTMB(
  cbind(TP, failure) ~ kmer_z * samp_z + (1 | vfold),
  family = betabinomial, data = model_df
)

# Reduced models for semipartial effect sizes (Type-II style)
glmm_noint  <- update(glmm_full, . ~ kmer_z + samp_z + (1 | vfold))
glmm_nokmer <- update(glmm_full, . ~ samp_z          + (1 | vfold))
glmm_nosamp <- update(glmm_full, . ~ kmer_z          + (1 | vfold))

cat("\n=== Beta-binomial GLMM: fixed-effect summary ===\n")
print(summary(glmm_full)$coefficients$cond)

cat("\n=== Type-III Wald ANOVA (which terms matter) ===\n")
print(car::Anova(glmm_full, type = 3))

# =============================================================================
# 2. Effect size: semipartial Delta-R^2 (marginal R^2)
# =============================================================================
# marginal R^2 = variance explained by FIXED effects only (Nakagawa).
# Delta-R^2 for a term = R2(model with term) - R2(model without it).
# This reads directly as "k-mer explains X% of the variance in Accuracy".

r2m <- function(m) {
  v <- tryCatch(performance::r2_nakagawa(m)$R2_marginal,
                error = function(e) NA_real_)
  as.numeric(v)
}

R2_full  <- r2m(glmm_full)
R2_noint <- r2m(glmm_noint)
R2_nokmer<- r2m(glmm_nokmer)
R2_nosamp<- r2m(glmm_nosamp)

effect_size <- tibble(
  term  = c("k-mer size", "sample size", "kmer x sample"),
  dR2   = c(R2_noint - R2_nokmer,      # kmer | sample
            R2_noint - R2_nosamp,      # sample | kmer
            R2_full  - R2_noint),      # interaction
) %>%
  mutate(
    dR2     = pmax(dR2, 0),                       # guard tiny negatives
    pct_var = 100 * dR2 / sum(dR2, na.rm = TRUE)  # share of explained variance
  )

cat("\n=== Effect size (semipartial Delta-R^2 on Accuracy) ===\n")
print(effect_size)

write_tsv(effect_size, file.path(opt$out, "effect_size_glmm.tsv"))

# =============================================================================
# 3. GAM smooth surface  (shape of each effect + formal interaction test)
# =============================================================================
# quasibinomial absorbs overdispersion; weights = denom.
# ti(kmer, sampling) is the pure interaction term -> its p-value tests
# whether the kmer effect depends on sample size.

gam_fit <- mgcv::gam(
  Accuracy ~ s(kmer_n, k = 5) + s(samp_n, k = 6) +
    ti(kmer_n, samp_n, k = c(5, 5)) +
    s(vfold, bs = "re"),
  family  = quasibinomial,
  weights = denom,
  data    = model_df, method = "REML"
)

cat("\n=== GAM summary (deviance explained, per-smooth EDF & p) ===\n")
print(summary(gam_fit))
# Heads-up: a non-significant ti() term => effects are essentially additive.

# =============================================================================
# 4. Random-forest cross-check  (model-agnostic importance)
# =============================================================================
# Permutation importance answers "kmer vs sample" without assuming a
# functional form; PDP gives the marginal shape.

set.seed(123)

rf_fit <- ranger(
  Accuracy ~ kmer_n + samp_n + vfold,
  data = model_df, num.trees = 1000,
  importance = "permutation", respect.unordered.factors = TRUE
)

rf_imp <- enframe(ranger::importance(rf_fit),
                  name = "term", value = "importance") %>%
  dplyr::filter(term %in% c("kmer_n", "samp_n")) %>%
  dplyr::mutate(term = dplyr::case_match(term, "kmer_n" ~ "k-mer size", "samp_n" ~ "sample size", .default = term))

cat("\n=== Random-forest permutation importance ===\n")
print(rf_imp)

# 1-D partial dependence for each predictor
pdp_kmer <- pdp::partial(rf_fit, pred.var = "kmer_n", train = model_df)
pdp_samp <- pdp::partial(rf_fit, pred.var = "samp_n", train = model_df)

# =============================================================================
# 5. Visualisation
# =============================================================================
es_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

## 5a. Interaction plot --------------------------------------------------------
# Parallel lines => additive; fanning => interaction.
int_summary <- model_df %>%
  group_by(kmer_n, samp_n) %>%
  summarise(Accuracy = mean(Accuracy), .groups = "drop")

p_interaction <- ggplot(int_summary,
                        aes(kmer_n, Accuracy,
                            group = samp_n, colour = samp_n)) +
  geom_line() + geom_point(size = 1) +
  scale_colour_viridis_c(name = "sample size") +
  labs(title = "A) Interaction plot",
       x = "k-mer size", y = "Mean Accuracy") +
  es_theme

## 5b. Model-based marginal effects -------------------------------------------
# Predicted Accuracy from the GAM, with CI ribbons.
me_kmer <- ggeffects::ggpredict(gam_fit, terms = "kmer_n [all]")
me_samp <- ggeffects::ggpredict(gam_fit, terms = "samp_n [all]")

p_me_kmer <- plot(me_kmer) +
  labs(title = "B) Marginal effect of k-mer", x = "k-mer size",
       y = "Predicted Accuracy") + es_theme
p_me_samp <- plot(me_samp) +
  labs(title = "C) Marginal effect of sample size", x = "sample size",
       y = "Predicted Accuracy") + es_theme

## 5c. Effect-size bar ---------------------------------------------------------
p_effect <- ggplot(effect_size,
                   aes(reorder(term, pct_var), pct_var, fill = term)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct_var)), hjust = -0.15, size = 3) +
  coord_flip(clip = "off") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "D) Effect size (share of explained variance)",
       x = NULL, y = "% of explained variance in Accuracy") +
  es_theme

## 5d. RF importance + partial-dependence -------------------------------------
p_rf_imp <- ggplot(rf_imp, aes(reorder(term, importance), importance,
                               fill = term)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_flip() + scale_fill_brewer(palette = "Set2") +
  labs(title = "E) RF permutation importance", x = NULL,
       y = "Importance") + es_theme

p_pdp <- bind_rows(
  pdp_kmer %>% rename(x = kmer_n) %>% mutate(predictor = "k-mer size"),
  pdp_samp %>% rename(x = samp_n) %>% mutate(predictor = "sample size")
) %>%
  ggplot(aes(x, yhat)) +
  geom_line(colour = "steelblue") +
  facet_wrap(~predictor, scales = "free_x") +
  labs(title = "F) Partial dependence (random forest)",
       x = NULL, y = "Predicted Accuracy") + es_theme

## 5e. Fitted-surface heatmap (model version of your left panel) --------------
grid <- expand_grid(
  kmer_n = seq(min(model_df$kmer_n), max(model_df$kmer_n), length.out = 60),
  samp_n = seq(min(model_df$samp_n), max(model_df$samp_n), length.out = 60),
  vfold  = factor(levels(model_df$vfold)[1], levels = levels(model_df$vfold))
)
grid$Accuracy_hat <- as.numeric(
  predict(gam_fit, newdata = grid, type = "response",
          exclude = "s(vfold)")            # average over the random effect
)

p_surface <- ggplot(grid, aes(kmer_n, samp_n, fill = Accuracy_hat)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = Accuracy_hat), colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Accuracy") +
  labs(title = "G) GAM-fitted Accuracy surface",
       x = "k-mer size", y = "sample size") + es_theme

# -----------------------------------------------------------------------------
# Save individual panels and one combined figure
# -----------------------------------------------------------------------------
combined <- (p_interaction | p_me_kmer | p_me_samp) /
  (p_effect      | p_rf_imp  | p_pdp)      /
  (p_surface     + plot_spacer() + plot_spacer()) &
  my_custom_theme()

ggsave(file.path(opt$out, "kmer_vs_sample_effectsize.png"),
       combined, width = 14, height = 11, dpi = 300)

message("\nWrote effect-size figure -> ",
        file.path(opt$out, "kmer_vs_sample_effectsize.png"))
message("Wrote effect-size table  -> ",
        file.path(opt$out, "effect_size_glmm.tsv"))

# -----------------------------------------------------------------------------
# Console interpretation hint
# -----------------------------------------------------------------------------
cat("\n--- INTERPRETATION GUIDE ---\n",
    "* effect_size$pct_var  : relative effect of kmer vs sample vs interaction.\n",
    "* car::Anova p-values  : whether each term is statistically significant.\n",
    "* GAM ti() term        : significant => effects interact; n.s. => additive.\n",
    "* RF importance + PDP  : model-agnostic confirmation of the above.\n",
    "If kmer's pct_var >> sample's, the vertical banding in your heatmap is\n",
    "confirmed: k-mer size is the dominant driver of assembly Accuracy.\n",
    sep = "")



quit()

# =============================================================================
# 4. Quick descriptive plots
# =============================================================================

my_theme <- function(base_size = 10, legend_pos = "top", ...) {
  theme_bw(base_family = "GillSans", base_size = base_size) +
    theme(legend.position = legend_pos,
          strip.placement = "outside", 
          strip.background = element_rect(fill = 'gray90', color = 'white'),
          strip.text = element_text(angle = 0, size = base_size, hjust = 0), 
          axis.text = element_text(size = rel(0.7), color = "black"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          ...
    )
}
cols_to <- c("Precision", "Sensitivity", "Accuracy", "Fscore")
recode_to <- c("A) Precision", "B) Sensitivity", "C) Accuracy", "D) F-score")
names(recode_to) <- cols_to

long_metrics <- metricsdf %>%
  pivot_longer(any_of(cols_to), names_to = "metric", values_to = "y") %>%
  mutate(metric = recode_factor(metric, !!!recode_to))

# Choose the x-axis: kmer if we swept it, otherwise sampling_set
xvar <- if (has_kmer) "kmer" else "sampling_set"

color_var <- if (length(unique(metricsdf$Assembly)) > 1) "Assembly" else NULL

# xvar <- "sampling_set"

var_ <- "Accuracy" #Accuracy Precision Sensitivity Fscore

breaks <- pretty(metricsdf[[var_]])
k      <- length(breaks)

labels_auto <- c(
  as.character(breaks[1]),                                    # "0"
  if (k > 2) paste0(breaks[1:(k-2)], "-", breaks[2:(k-1)]), # middle intervals
  paste0(">", breaks[k-1])                                    # open upper bin
)

# Dynamic palette: one colour per factor level
n_levels   <- length(labels_auto)
fill_vals  <- rev(scales::pal_brewer(direction = -1, palette = "YlGnBu")(n_levels))
fill_named <- setNames(fill_vals, labels_auto)


p_kmer_sampling <- metricsdf |>
  group_by(kmer, sampling_set) |> rstatix::get_summary_stats(all_of(var_), type = "mean_se") |>
  mutate(countfactor   = cut(mean,
                             breaks = c(-1, breaks),
                             labels = labels_auto)) |>
  ggplot(aes(kmer, sampling_set, fill = countfactor)) +
  geom_tile(color = "white", linewidth = 0.4, na.rm = FALSE) +
  scale_fill_manual(
    values   = fill_named,
    na.value = "grey90",
    breaks   = labels_auto          # legend ordered low → high
  ) +
  my_theme(base_size = 10, legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 1,
                             title          = var_,
                             title.position = "top",
                             title.hjust    = 0,
                             keywidth       = unit(0.35, "cm"),
                             keyheight      = unit(0.35, "cm"),
                             override.aes   = list(size = 5)
  ))


ggsave(file.path(opt$out, "kmer_vs_sampling_metric.png"),
       p_kmer_sampling, width = 4, height = 5, dpi = 300)

p_main <- long_metrics %>%
  filter(sampling_set == 1) |>
  ggplot(aes(x = .data[[xvar]], y = y,
             group = if (!is.null(color_var)) .data[[color_var]] else 1,
             color = if (!is.null(color_var)) .data[[color_var]] else NULL)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "pointrange", shape = 1) +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = xvar, y = NULL,
       caption = sprintf("transrate_EDA.R | rcov >= %s | min length %d",
                         opt$rcov, opt$minlen)) +
  my_theme()

ggsave(file.path(opt$out, "benchmark_metrics.png"),
       p_main, width = 8, height = 6, dpi = 300)

message("Wrote ", file.path(opt$out, "benchmark_metrics.png"))

# Precision-vs-Sensitivity scatter
p_scatter <- metricsdf %>%
  ggplot(aes(x = Sensitivity, y = Precision,
             color = if (has_kmer) factor(kmer) else factor(sampling_set))) +
  geom_point(aes(alpha = sampling_set)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4,
              aes(group = 1), color = "black") +
  labs(color = if (has_kmer) "kmer" else "sampling_set",
       caption = "transrate_EDA.R") +
  my_theme() +
  guides(
    alpha = guide_legend(nrow = 2,
                         title          = "Sample",
                         title.position = "top",
                         title.hjust    = 0,
                         keywidth       = unit(0.35, "cm"),
                         keyheight      = unit(0.35, "cm"),
                         override.aes   = list(size = 2)),
    
    color = guide_legend(nrow = 2,
                         title          = "Kmer",
                         title.position = "top",
                         title.hjust    = 0,
                         keywidth       = unit(0.35, "cm"),
                         keyheight      = unit(0.35, "cm"),
                         override.aes   = list(size = 2))
    
  )

ggsave(file.path(opt$out, "precision_vs_sensitivity.png"),
       p_scatter, width = 5, height = 5, dpi = 300)

# Stratified contig-quality counts (mimicking the "% alignment" plot)
align_strata <- transratedf %>%
  mutate(strat = case_when(
    reference_coverage == 1   ~ "100% alignment",
    reference_coverage >= 0.95 ~ ">= 95% alignment",
    reference_coverage >= 0.90 ~ ">= 90% alignment",
    reference_coverage >= 0.80 ~ ">= 80% alignment",
    TRUE                      ~ "< 80% alignment"
  )) %>%
  count(!!!syms(group_vars), strat) %>%
  mutate(strat = factor(strat, levels = c(
    "< 80% alignment", ">= 80% alignment", ">= 90% alignment",
    ">= 95% alignment", "100% alignment")))


n_pallet <- length(unique(align_strata$strat))

scale_col <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_col <- structure(scale_col, names = sort(levels(unique(align_strata$strat))))

scale_fill <- ggsci::pal_uchicago(alpha = 0.8)(n_pallet) 

scale_fill <- structure(scale_fill, names = sort(levels(unique(align_strata$strat))))


p_strata <- align_strata %>%
  ggplot(aes(x = .data[[xvar]], y = n, color = strat, fill = strat, group = strat)) +
  # Add connecting lines
  stat_summary(fun.data = mean_se, geom = "line") +
  # Add error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, size = 0.5) +
  # Add points for means
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, color = "white", stroke = 0.5) +
  # stat_summary(fun.data = mean_se, geom = "pointrange", shape = 1, size = 0.3) +
  { if (length(unique(metricsdf$Assembly)) > 1) facet_wrap(~Assembly) else NULL } +
  labs(y = "Number of assembled conotoxins", x = "Kmer size",color = NULL,
       caption = "transrate_EDA.R") +
  my_theme(legend.text = element_text(size = 3), legend_pos = "top") +
  scale_color_manual("", values = scale_col) +
  scale_fill_manual("", values = scale_col) 

opt$out <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

ggsave(file.path(opt$out, "alignment_strata.png"),
       p_strata, width = 4, height = 4, dpi = 300)

message("\nDone. Outputs in ", normalizePath(opt$out))

library(patchwork)

PSAVE <- wrap_plots(
  C = p_strata,
  B = p_scatter,
  A = p_kmer_sampling,
  design = "ABC"
) +
  plot_layout(
    heights = 1,
    widths  = 0.7
  ) &
  theme(legend.position = "top", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.spacing = unit(0, "cm"))


ggsave(
  PSAVE,
  filename = "kmer_to_sample.png",
  path     = opt$out,
  width    = 12,
  height   = 5,
  device   = png,
  dpi      = 800
)