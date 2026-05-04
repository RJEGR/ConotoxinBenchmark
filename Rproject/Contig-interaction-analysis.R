# =============================================================================
#  ANÁLISIS DE INTERACCIONES LINEALES Y NO LINEALES
#  Contexto: Clasificación de proteínas/DNA por verosimilitud de ensamblaje
#  Dataset:  contigbaseddatalonger.csv
#  Autor:    [Editar con tu nombre]
#  Fecha:    2026-04
# =============================================================================
#
#  ESTRUCTURA DEL SCRIPT:
#  BLOQUE 0 — Configuración y carga de librerías
#  BLOQUE 1 — Carga y preprocesamiento
#  BLOQUE 2 — Análisis LINEAL (correlaciones, ANOVA, modelos lineales)
#  BLOQUE 3 — Análisis NO LINEAL (Random Forest, GBM, SHAP, MI)
#  BLOQUE 4 — Interacciones Numéricas × Categóricas (boxplots, PCA, UMAP)
#  BLOQUE 5 — Modelado predictivo (Multinomial + RF + XGBoost)
#  BLOQUE 6 — Exportar resultados
# =============================================================================


# =============================================================================
# BLOQUE 0 — CONFIGURACIÓN
# =============================================================================


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# ── 0.1  Instalar paquetes faltantes ─────────────────────────────────────────
pkgs_needed <- c(
  # Datos y utilidades
  "tidyverse", "data.table", "janitor",
  # Estadística lineal
  "corrplot", "ggcorrplot", "psych", "car", "effectsize",
  # Machine Learning
  "randomForest", "ranger", "xgboost", "caret",
  # SHAP / Importancia
  "shapviz", "SHAPforxgboost",
  # Información mutua / no lineal
  "infotheo", "minerva",
  # Visualización avanzada
  "ggplot2", "ggpubr", "patchwork", "viridis", "RColorBrewer",
  "GGally", "ggridges", "ggforce",
  # Reducción dimensional
  "FactoMineR", "factoextra",
  # Modelos multinomiales
  "nnet", "broom", "broom.mixed"
)

# Instalar solo los que faltan (descomenta si es la primera vez)
# install.packages(setdiff(pkgs_needed, installed.packages()[,"Package"]))

# ── 0.2  Cargar librerías ────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(corrplot)
  library(ggcorrplot)
  library(psych)
  library(car)
  library(effectsize)
  library(randomForest)
  library(ranger)
  library(xgboost)
  library(caret)
  library(infotheo)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(viridis)
  library(GGally)
  library(ggridges)
  library(FactoMineR)
  library(factoextra)
  library(nnet)
  library(broom)
})

# ── 0.3  Parámetros globales (EDITAR AQUÍ) ───────────────────────────────────
setwd("~/Documents/GitHub/ConotoxinBenchmark/INPUTS/")
RUTA_CSV        <- "contig-based-data-longer.csv"   # Ruta al archivo CSV
DIR_SALIDA      <- "resultados_interacciones"     # Carpeta de resultados
SEMILLA         <- 42                             # Semilla aleatoria
N_ARBOLES_RF    <- 500                            # Árboles en Random Forest
TOP_N_FEATURES  <- 15                             # Top features a visualizar
PALETA_CLASES   <- c(                             # Colores por clase final
  "full"     = "#10b981",
  "multi"    = "#00d4ff",
  "chimera"  = "#ef4444",
  "error"    = "#f59e0b",
  "fragment" = "#a78bfa"
)

# Crear directorio de salida
if (!dir.exists(DIR_SALIDA)) dir.create(DIR_SALIDA, recursive = TRUE)

set.seed(SEMILLA)
cat("\n✓ Configuración completada\n")


# =============================================================================
# BLOQUE 1 — CARGA Y PREPROCESAMIENTO
# =============================================================================

cat("── Cargando datos...\n")
df_raw <- fread(RUTA_CSV, na.strings = c("", "NA", "N/A", "NaN"))

# ── 1.1  Definir columnas por tipo ───────────────────────────────────────────
COLS_NUMERIC <- c(
  "Cys_number", "Hydrophobicity", "seq_freq", "Protein_width", "Score_sf",
  "at_skew", "coverage", "cpg_count", "cpg_ratio", "eff_count", "eff_length",
  "gc_skew", "in_bridges", "length", "linguistic_complexity_6", "orf_length",
  "p_bases_covered", "p_good", "p_not_segmented", "p_seq_true",
  "prop_gc", "score", "tpm", "reference_coverage"
)

COLS_CATEG_BINARIA <- c(
  "sf_evidence",   # Concordance / Ambiguous
  "tab"            # Regex / pHMM
)

COLS_CATEG_CUALIT <- c(
  "final_annotation",  # TARGET: full/multi/chimera/error/fragment
  "prelim_cat",        # m->n, 1->n, m->1, 1->1
  "Region",            # (Mature)_(Pro)_(Signal) etc.
  "summarise",         # % alignment quality
  "Method",            # Assembler
  "gs_conoSorter",     # Superfamilia conotoxina
  "gs_conoServer"      # Servidor de predicción
)

# ── 1.2  Preprocesamiento ────────────────────────────────────────────────────
df <- df_raw %>%
  as_tibble() %>%
  # Convertir numéricas
  mutate(across(all_of(COLS_NUMERIC), as.numeric)) %>%
  # Convertir categóricas
  mutate(across(all_of(c(COLS_CATEG_BINARIA, COLS_CATEG_CUALIT)), as.factor)) %>%
  # Transformaciones logarítmicas para variables muy sesgadas
  mutate(
    log1p_tpm       = log1p(tpm),
    log1p_coverage  = log1p(coverage),
    log1p_eff_count = log1p(eff_count),
    log1p_in_bridges= log1p(in_bridges),
    log1p_cpg_count = log1p(cpg_count),
    log1p_ref_coverage = log1p(reference_coverage),
    
  ) %>%
  # Filtrar filas sin target
  filter(!is.na(final_annotation))

# Dataset limpio solo con numericas (sin NA en cols base)
COLS_NUM_BASE <- c("at_skew", "coverage", "cpg_count", "cpg_ratio", "eff_count",
                   "eff_length", "gc_skew", "in_bridges", "length",
                   "linguistic_complexity_6", "orf_length", "p_bases_covered",
                   "p_good", "p_not_segmented", "p_seq_true", "prop_gc",
                   "score", "tpm", "reference_coverage","log1p_tpm", "log1p_coverage",
                   "log1p_eff_count", "log1p_in_bridges", "log1p_cpg_count", "log1p_ref_coverage")

df_num_completo <- df %>% select(all_of(COLS_NUM_BASE), final_annotation)

cat(sprintf("✓ Dataset: %d filas × %d columnas\n", nrow(df), ncol(df)))
cat(sprintf("  Clases en final_annotation: %s\n",
            paste(levels(df$final_annotation), collapse = ", ")))
cat(sprintf("  Filas con todas las numéricas base: %d\n",
            sum(complete.cases(df[, COLS_NUM_BASE]))))


# =============================================================================
# BLOQUE 2 — ANÁLISIS LINEAL
# =============================================================================

cat("\n══ BLOQUE 2: Análisis Lineal ══\n")

# ── 2.1  Matriz de Correlación de Pearson (numéricas) ────────────────────────
cat("  [2.1] Correlación de Pearson...\n")

mat_cor <- df %>%
  select(all_of(COLS_NUM_BASE)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson")

p1_corrplot <- ggcorrplot(
  mat_cor,
  method    = "square",
  type      = "upper",
  lab       = TRUE,
  lab_size  = 2.2,
  colors    = c("#ef4444", "#0f1117", "#10b981"),
  outline.color = "grey20",
  title     = "Correlación de Pearson — Métricas numéricas de calidad de contigs",
  ggtheme   = theme_minimal(base_size = 11)
) +
  theme(
    plot.title       = element_text(size = 12, face = "bold", color = "grey10"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    legend.position  = "bottom"
  )

ggsave(file.path(DIR_SALIDA, "2_1_correlacion_pearson.png"),
       p1_corrplot, width = 14, height = 12, dpi = 150)

# Pares con alta correlación (|r| > 0.7)
cor_df <- mat_cor %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
  filter(var1 < var2, abs(r) > 0.70) %>%
  arrange(desc(abs(r)))

cat("  Pares con |r| > 0.70:\n")
print(cor_df)


# ── 2.2  Correlación de Spearman (robusta a no normalidad) ───────────────────
cat("  [2.2] Correlación de Spearman...\n")

mat_spear <- df %>%
  select(all_of(COLS_NUM_BASE)) %>%
  cor(use = "pairwise.complete.obs", method = "spearman")

p2_spearman <- ggcorrplot(
  mat_spear,
  method    = "circle",
  type      = "upper",
  lab       = FALSE,
  colors    = c("#7c3aed", "#f8fafc", "#f59e0b"),
  title     = "Correlación de Spearman — Métricas numéricas (robusta a outliers)",
  ggtheme   = theme_minimal(base_size = 10)
)

ggsave(file.path(DIR_SALIDA, "2_2_correlacion_spearman.png"),
       p2_spearman, width = 12, height = 10, dpi = 150)


# ── 2.3  ANOVA: Numéricas ~ final_annotation ────────────────────────────────
cat("  [2.3] ANOVA uno a uno (numérica ~ final_annotation)...\n")

resultados_anova <- map_dfr(COLS_NUM_BASE, function(col_nm) {
  datos_col <- df %>%
    select(y = !!sym(col_nm), grupo = final_annotation) %>%
    drop_na()
  
  if (nrow(datos_col) < 50) return(NULL)
  
  # Modelo ANOVA
  mod <- aov(y ~ grupo, data = datos_col)
  sm  <- summary(mod)[[1]]
  
  # Tamaño de efecto (eta cuadrado)
  eta2 <- eta_squared(mod, partial = FALSE)$Eta2[1]
  
  tibble(
    variable  = col_nm,
    F_stat    = sm[["F value"]][1],
    p_valor   = sm[["Pr(>F)"]][1],
    eta2      = eta2,
    sig       = case_when(
      sm[["Pr(>F)"]][1] < 0.001 ~ "***",
      sm[["Pr(>F)"]][1] < 0.01  ~ "**",
      sm[["Pr(>F)"]][1] < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )
})

resultados_anova <- resultados_anova %>%
  arrange(desc(eta2)) %>%
  mutate(p_ajustado = p.adjust(p_valor, method = "BH"))  # Corrección FDR

cat("  Top 10 variables por eta² (efecto sobre final_annotation):\n")
print(head(resultados_anova, 10))

# Plot de eta² ordenado
p3_anova <- resultados_anova %>%
  mutate(variable = fct_reorder(variable, eta2)) %>%
  ggplot(aes(x = eta2, y = variable, fill = eta2)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sig), hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "ANOVA: Efecto de cada variable numérica sobre final_annotation",
    subtitle = "Tamaño de efecto η² (Eta cuadrado) — valores más altos = mayor diferencia entre clases",
    x        = "η² (Eta cuadrado)",
    y        = NULL,
    caption  = "Significancia: *** p<0.001 · ** p<0.01 · * p<0.05 · n.s. no significativo (corrección FDR Benjamini-Hochberg)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(DIR_SALIDA, "2_3_anova_eta2.png"),
       p3_anova, width = 12, height = 8, dpi = 150)


# ── 2.4  Modelo Lineal Generalizado (Multinomial Logístico) ──────────────────
cat("  [2.4] Regresión Multinomial Logística...\n")

# Usar solo columnas sin NA para este modelo
df_multinom <- df %>%
  select(final_annotation, all_of(COLS_NUM_BASE), sf_evidence, tab, prelim_cat) %>%
  drop_na(score, p_good, coverage, log1p_tpm, orf_length, prop_gc,
          linguistic_complexity_6, p_not_segmented, p_seq_true, final_annotation) %>%
  mutate(
    # Referenciar contra "full" (clase más "buena")
    final_annotation = relevel(final_annotation, ref = "full"),
    sf_evidence      = relevel(factor(sf_evidence), ref = "Concordance")
  )

# Modelo reducido con features clave
mod_multinom <- multinom(
  final_annotation ~ reference_coverage + score + p_good + log1p_coverage + log1p_tpm +
    orf_length + prop_gc + linguistic_complexity_6 +
    p_not_segmented + p_seq_true + p_bases_covered +
    at_skew + gc_skew + cpg_ratio,
  data    = df_multinom,
  MaxNWts = 5000,
  maxit   = 200,
  trace   = FALSE
)

cat("  Coeficientes del modelo multinomial (top):\n")

coef_tabla <- tidy(mod_multinom, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  arrange(p.value) %>%
  mutate(
    OR  = exp(estimate),
    sig = ifelse(p.value < 0.05, "✓", "")
  )

print(head(coef_tabla, 20))

# Guardar tabla completa
write.csv(coef_tabla,
          file.path(DIR_SALIDA, "2_4_coeficientes_multinomial.csv"),
          row.names = FALSE)

# Plot de coeficientes (Odds Ratios)
p4_OR <- coef_tabla %>%
  filter(abs(estimate) > 0.1) %>%
  mutate(term = fct_reorder(term, estimate)) %>%
  ggplot(aes(x = estimate, y = term, color = y.level)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3) +
  facet_wrap(~ y.level, scales = "free_x", nrow = 2) +
  scale_color_manual(values = PALETA_CLASES) +
  labs(
    title    = "Regresión Multinomial: Coeficientes log-odds vs. clase 'full'",
    subtitle = "Barras = IC 95% · Referencia: full (transcripto completo de alta calidad)",
    x = "Coeficiente (log-odds)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "2_4_multinomial_coeficientes.png"),
       p4_OR, width = 15, height = 8, dpi = 150)


# ── 2.5  Kruskal-Wallis (no paramétrico por clase) ───────────────────────────
cat("  [2.5] Kruskal-Wallis (alternativa no paramétrica a ANOVA)...\n")

resultados_kw <- map_dfr(COLS_NUM_BASE, function(col_nm) {
  datos_col <- df %>% select(y = !!sym(col_nm), grupo = final_annotation) %>% drop_na()
  if (nrow(datos_col) < 30) return(NULL)
  kw <- kruskal.test(y ~ grupo, data = datos_col)
  tibble(
    variable  = col_nm,
    statistic = kw$statistic,
    p_valor   = kw$p.value,
    p_ajust   = p.adjust(kw$p.value, method = "BH")
  )
}) %>% arrange(p_valor)

write.csv(resultados_kw,
          file.path(DIR_SALIDA, "2_5_kruskal_wallis.csv"),
          row.names = FALSE)

cat("  Variables con p-ajustado < 0.05:", sum(resultados_kw$p_ajust < 0.05, na.rm = TRUE), "\n")


# =============================================================================
# BLOQUE 3 — ANÁLISIS NO LINEAL
# =============================================================================

cat("\n══ BLOQUE 3: Análisis No Lineal ══\n")

# ── 3.1  Información Mutua (MI) ──────────────────────────────────────────────
cat("  [3.1] Información Mutua Numérica × final_annotation...\n")

# Discretizar numéricas para MI
df_disc <- df %>%
  select(all_of(COLS_NUM_BASE)) %>%
  drop_na() %>%
  mutate(across(everything(), ~ discretize(., disc = "equalfreq", nbins = 10)[[1]]))

target_disc <- df %>%
  filter(!is.na(final_annotation)) %>%
  slice(which(complete.cases(df[, COLS_NUM_BASE]))) %>%
  pull(final_annotation) %>%
  as.integer()

mi_scores <- map_dbl(COLS_NUM_BASE, function(col_nm) {
  x_vec <- df_disc[[col_nm]]
  if (length(x_vec) != length(target_disc)) return(NA)
  mutinformation(x_vec, target_disc)
})

mi_df <- tibble(
  variable = COLS_NUM_BASE,
  MI       = mi_scores
) %>% arrange(desc(MI))

cat("  Top 10 por Información Mutua con final_annotation:\n")

print(head(mi_df, 10))

p5_MI <- mi_df %>%
  mutate(variable = fct_reorder(variable, MI)) %>%
  ggplot(aes(x = MI, y = variable, fill = MI)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(option = "magma", direction = -1, begin = 0.2) +
  labs(
    title    = "Información Mutua: Variables numéricas → final_annotation",
    subtitle = "MI captura relaciones lineales Y no lineales. Mayor MI = mayor dependencia estadística",
    x = "Información Mutua (bits)", y = NULL,
    caption  = "Discretización: 10 bins de frecuencia igual (equalfreq)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "3_1_informacion_mutua.png"),
       p5_MI, width = 11, height = 8, dpi = 150)


# ── 3.2  Random Forest: Importancia de Features (Gini + Permutación) ─────────
cat("  [3.2] Random Forest — Importancia de features...\n")

df_rf <- df %>%
  select(final_annotation, all_of(COLS_NUM_BASE)) %>%
  drop_na(score, p_good, log1p_coverage, log1p_tpm, orf_length, prop_gc,
          linguistic_complexity_6, p_not_segmented, p_seq_true,
          p_bases_covered, at_skew, gc_skew, cpg_ratio, final_annotation)

# Usar ranger (más rápido que randomForest para 65K rows)
modelo_rf <- ranger(
  formula          = final_annotation ~ .,
  data             = df_rf,
  num.trees        = N_ARBOLES_RF,
  importance       = "permutation",   # más robusto que "impurity"
  classification   = TRUE,
  class.weights    = table(df_rf$final_annotation) / nrow(df_rf),  # pesos por imbalance
  seed             = SEMILLA,
  num.threads      = parallel::detectCores() - 1
)

cat(sprintf("  RF OOB accuracy: %.3f\n", 1 - modelo_rf$prediction.error))

# Importancia
imp_df <- tibble(
  variable   = names(modelo_rf$variable.importance),
  importancia = modelo_rf$variable.importance
) %>% arrange(desc(importancia))

cat("  Top 10 features por importancia de permutación RF:\n")

print(head(imp_df, 10))

p6_RF_imp <- imp_df %>%
  head(TOP_N_FEATURES) %>%
  mutate(variable = fct_reorder(variable, importancia)) %>%
  ggplot(aes(x = importancia, y = variable, fill = importancia)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(option = "cividis", direction = -1) +
  labs(
    title    = "Random Forest: Importancia por permutación (Top features)",
    subtitle = sprintf("ranger · %d árboles · class.weights compensando imbalance 11× · OOB acc = %.3f",
                       N_ARBOLES_RF, 1 - modelo_rf$prediction.error),
    x = "Importancia (caída en accuracy por permutación)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "3_2_rf_importancia_permutacion.png"),
       p6_RF_imp, width = 11, height = 7, dpi = 150)


# ── 3.3  XGBoost + SHAP Values ───────────────────────────────────────────────
cat("  [3.3] XGBoost + SHAP interaction values...\n")

# Preparar matriz para XGBoost
features_xgb <- c("score", "p_good", "log1p_coverage", "log1p_tpm",
                  "orf_length", "prop_gc", "linguistic_complexity_6",
                  "p_not_segmented", "p_seq_true", "p_bases_covered",
                  "at_skew", "gc_skew", "cpg_ratio", "log1p_in_bridges",
                  "log1p_eff_count", "log1p_cpg_count", "eff_length", "reference_coverage")

df_xgb <- df %>%
  select(final_annotation, all_of(features_xgb)) %>%
  drop_na()

X_mat <- df_xgb %>% select(-final_annotation) %>% as.matrix()
y_vec <- as.integer(df_xgb$final_annotation) - 1  # 0-indexed

dtrain <- xgb.DMatrix(data = X_mat, label = y_vec)

params_xgb <- list(
  objective        = "multi:softprob",
  num_class        = length(levels(df_xgb$final_annotation)),
  eta              = 0.1,
  max_depth        = 6,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  eval_metric      = "mlogloss",
  seed             = SEMILLA
)

modelo_xgb <- xgb.train(
  params   = params_xgb,
  data     = dtrain,
  nrounds  = 300,
  verbose  = 0
)

# Importancia XGBoost (Gain, Cover, Frequency)
xgb_imp <- xgb.importance(model = modelo_xgb) %>% as_tibble()

p7_xgb_imp <- xgb_imp %>%
  head(TOP_N_FEATURES) %>%
  mutate(Feature = fct_reorder(Feature, Gain)) %>%
  pivot_longer(cols = c(Gain, Cover, Frequency), names_to = "tipo", values_to = "valor") %>%
  ggplot(aes(x = valor, y = Feature, fill = tipo)) +
  geom_col(position = "dodge", width = 0.6) +
  facet_wrap(~tipo, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "XGBoost: Importancia de features — Gain, Cover, Frequency",
    subtitle = "Gain = reducción de impureza total (más relevante biológicamente)",
    x = "Valor", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "3_3_xgboost_importancia.png"),
       p7_xgb_imp, width = 14, height = 7, dpi = 150)

# ── SHAP Values (requiere paquete shapviz o SHAPforxgboost) ──────────────────
# Descomenta si tienes SHAPforxgboost instalado:
# library(SHAPforxgboost)
# shap_vals <- shap.values(xgb_model = modelo_xgb, X_train = X_mat)
# shap_long  <- shap.prep(shap_contrib = shap_vals$shap_score, X_train = X_mat)
# p_shap     <- shap.plot.summary(shap_long)
# ggsave(file.path(DIR_SALIDA, "3_3_shap_summary.png"), p_shap, width=12, height=8, dpi=150)


# ── 3.4  Correlación No Lineal — Maximal Information Coefficient (MIC) ───────
cat("  [3.4] MIC entre pares de numéricas clave...\n")

cols_mic <- c("score", "p_good", "log1p_coverage", "log1p_tpm",
              "orf_length", "prop_gc", "linguistic_complexity_6",
              "p_not_segmented", "at_skew", "gc_skew","reference_coverage")

df_mic_input <- df %>%
  select(all_of(cols_mic)) %>%
  drop_na() %>%
  slice_sample(n = min(5000, nrow(.)))  # MIC es lento en >10K filas

# Calcular MIC (mine() de minerva)
# Descomenta si tienes minerva instalado:
library(minerva)
mic_result <- mine(df_mic_input)
mic_mat    <- mic_result$MIC
p_mic <- ggcorrplot(mic_mat, method="square", lab=TRUE, lab_size=2.5,
                    colors=c("#0f1117","#7c3aed","#f59e0b"),
                    title="MIC — Maximal Information Coefficient (no lineal)")
# ggsave(file.path(DIR_SALIDA, "3_4_MIC_nonlinear.png"), p_mic, width=10, height=9, dpi=150)

cat("  [MIC: descomenta el bloque si tienes minerva instalado]\n")


# =============================================================================
# BLOQUE 4 — INTERACCIONES NUMÉRICAS × CATEGÓRICAS
# =============================================================================

cat("\n══ BLOQUE 4: Interacciones Numéricas × Categóricas ══\n")

# ── 4.1  Ridge plots: distribución por clase ─────────────────────────────────
cat("  [4.1] Ridge plots — distribución de features por final_annotation...\n")

# Top 6 features por eta² (más discriminativas)
top6_vars <- resultados_anova$variable[1:6]

df_ridge <- df %>%
  select(final_annotation, all_of(top6_vars)) %>%
  drop_na() %>%
  pivot_longer(-final_annotation, names_to = "variable", values_to = "valor") %>%
  filter(is.finite(valor))

p8_ridge <- ggplot(df_ridge, aes(x = valor, y = final_annotation,
                                 fill = final_annotation)) +
  geom_density_ridges(alpha = 0.75, scale = 1.5, rel_min_height = 0.01,
                      bandwidth = NULL, show.legend = FALSE) +
  facet_wrap(~ variable, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = PALETA_CLASES) +
  labs(
    title    = "Distribución de features discriminativas por clase de anotación",
    subtitle = "Top 6 variables por ANOVA η² — separación = poder discriminativo",
    x = "Valor", y = "Clase final_annotation"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "4_1_ridge_plots_por_clase.png"),
       p8_ridge, width = 14, height = 9, dpi = 150)


# ── 4.2  Boxplots: Numéricas × sf_evidence (binaria) ─────────────────────────
cat("  [4.2] Boxplots Numéricas × sf_evidence (binaria)...\n")

vars_box <- c("score", "p_good", "log1p_coverage", "linguistic_complexity_6",
              "prop_gc", "orf_length", "reference_coverage")

intercept <- "sf_evidence" # 

df_box_sf <- df %>%
  select(all_of(c(intercept, vars_box))) %>%
  drop_na(all_of(intercept)) %>%
  pivot_longer(-all_of(intercept), names_to = "variable", values_to = "valor") 

names(df_box_sf)[names(df_box_sf) %in% intercept] <- "Intercept"


p9_box_sf <- ggplot(df_box_sf, aes(x = Intercept, y = valor, fill = Intercept)) +
  geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4,
                     vjust = -0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  # scale_fill_manual(values = c("Concordance" = "#10b981", "Ambiguous" = "#ef4444")) +
  labs(
    title    = "Violin + Boxplot: Métricas numéricas × Intercept (Concordance vs Ambiguous)",
    subtitle = "Comparación de distribuciones — prueba de Wilcoxon bilateral",
    x = NULL, y = "Valor", fill = intercept
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(DIR_SALIDA, "4_2_violin_sf_evidence.png"),
       p9_box_sf, width = 14, height = 9, dpi = 150)


# ── 4.3  Boxplots: Numéricas × tab (Regex vs pHMM) ───────────────────────────
cat("  [4.3] Boxplots Numéricas × tab (Regex vs pHMM)...\n")

df_box_tab <- df %>%
  select(tab, all_of(vars_box)) %>%
  drop_na(tab) %>%
  pivot_longer(-tab, names_to = "variable", values_to = "valor")

p10_box_tab <- ggplot(df_box_tab, aes(x = tab, y = valor, fill = tab)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Regex" = "#00d4ff", "pHMM" = "#7c3aed")) +
  labs(
    title    = "Violin + Boxplot: Métricas numéricas × tab (Regex vs pHMM)",
    subtitle = "Método de detección de señal peptídica",
    x = NULL, y = "Valor", fill = "tab"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(DIR_SALIDA, "4_3_violin_tab_regex_pHMM.png"),
       p10_box_tab, width = 14, height = 9, dpi = 150)


# ── 4.4  Heatmap: Mediana de numéricas por clase (normalizado) ───────────────
cat("  [4.4] Heatmap de perfiles de clase...\n")

heatmap_data <- df %>%
  select(final_annotation, all_of(COLS_NUM_BASE)) %>%
  drop_na(final_annotation) %>%
  group_by(final_annotation) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-final_annotation, names_to = "variable", values_to = "mediana") %>%
  group_by(variable) %>%
  mutate(z_score = scale(mediana)[, 1]) %>%  # normalizar por variable
  ungroup()

p11_heatmap <- ggplot(heatmap_data, aes(x = final_annotation, y = variable,
                                        fill = z_score)) +
  geom_tile(color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = round(mediana, 2)), size = 2.5, color = "white") +
  scale_fill_gradient2(low = "#1e40af", mid = "#0f1117", high = "#dc2626",
                       midpoint = 0, name = "Z-score") +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  labs(
    title    = "Perfiles de clase: Mediana de métricas numéricas (Z-score por variable)",
    subtitle = "Rojo = alto relativo · Azul = bajo relativo — números = valor crudo de mediana",
    x = "Clase final_annotation", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

ggsave(file.path(DIR_SALIDA, "4_4_heatmap_perfiles_clase.png"),
       p11_heatmap, width = 12, height = 10, dpi = 150)


# ── 4.5  PCA: Reducción dimensional coloreada por clase ──────────────────────
cat("  [4.5] PCA sobre numéricas completas...\n")

# -----------------------------------------------------------------------------
# FIX: Error "infinite or missing values in 'x'" en prcomp/svd
#
# CAUSA RAÍZ: drop_na() parcial + COLS_NUM_BASE contiene columnas con 54% NA
#   (Cys_number, Hydrophobicity, seq_freq, Protein_width, Score_sf,
#    reference_coverage). Esas columnas sobreviven el drop_na() parcial,
#    luego scale() produce NaN/Inf, y svd() falla.
#
# ESTRATEGIA A — PCA sobre columnas SIEMPRE disponibles (sin NA, sin Inf)
#   Usa las 17 features de calidad de ensamblaje que tienen 0% de missingness.
#   Es la más reproducible y estable.
#
# ESTRATEGIA B — PCA sobre TODAS las columnas (incluyendo conotoxina)
#   Requiere imputación con la mediana antes de escalar.
#   Descomenta el bloque B si quieres incluir Cys_number, Hydrophobicity, etc.
# -----------------------------------------------------------------------------

# ── Columnas sin missingness estructural (0% NA en dataset) ──────────────────
COLS_PCA_COMPLETAS <- c(
  "score", "p_good", "p_not_segmented", "p_seq_true", "p_bases_covered",
  "log1p_tpm", "log1p_coverage", "log1p_eff_count", "log1p_in_bridges",
  "log1p_cpg_count", "orf_length", "prop_gc", "linguistic_complexity_6",
  "at_skew", "gc_skew", "cpg_ratio", "eff_length", "reference_coverage"
)

# ─────────────────────────────────────────────────────────────────────────────
# ESTRATEGIA A (DEFAULT): solo columnas con 0% NA
# ─────────────────────────────────────────────────────────────────────────────

df_pca_input <- df %>%
  select(final_annotation, all_of(COLS_PCA_COMPLETAS)) %>%
  filter(!is.na(final_annotation)) %>%
  # Paso 1: eliminar filas con NA o Inf en CUALQUIER columna seleccionada
  filter(if_all(all_of(COLS_PCA_COMPLETAS), ~ is.finite(.x)))

# Paso 2: construir matriz numérica — eliminar columnas de varianza cero
X_pca_raw <- df_pca_input %>% select(-final_annotation)

cols_var_ok <- X_pca_raw %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE) > 1e-10)) %>%
  unlist()

X_pca_raw <- X_pca_raw[, cols_var_ok, drop = FALSE]

# Paso 3: estandarizar (z-score) — usar scale() con center=TRUE, scale=TRUE
#   de forma explícita para evitar NaN cuando sd=0
X_pca <- as.data.frame(
  scale(X_pca_raw, center = TRUE, scale = TRUE)
)

# Paso 4: verificación defensiva — eliminar filas/cols con NA o Inf residuales
cols_ok  <- colSums(!is.finite(as.matrix(X_pca))) == 0
X_pca    <- X_pca[, cols_ok, drop = FALSE]
rows_ok  <- complete.cases(X_pca) & apply(X_pca, 1, function(r) all(is.finite(r)))
X_pca    <- X_pca[rows_ok, , drop = FALSE]
labels_pca <- df_pca_input$final_annotation[rows_ok]

cat(sprintf("  PCA input: %d filas × %d columnas (estrategia A: 0%% NA)\n",
            nrow(X_pca), ncol(X_pca)))
cat(sprintf("  Columnas incluidas: %s\n", paste(names(X_pca), collapse = ", ")))

res_pca <- prcomp(X_pca, center = FALSE, scale. = FALSE)

# Varianza explicada
var_exp <- round(100 * res_pca$sdev^2 / sum(res_pca$sdev^2), 1)
cat(sprintf("  PC1: %.1f%% · PC2: %.1f%% · PC3: %.1f%% de varianza\n",
            var_exp[1], var_exp[2], var_exp[3]))

pca_coords <- as.data.frame(res_pca$x[, 1:3]) %>%
  mutate(clase = df_pca_input$final_annotation)

p12_pca <- ggplot(pca_coords %>% slice_sample(n = min(10000, nrow(.))),
                  aes(x = PC1, y = PC2, color = clase)) +
  geom_point(alpha = 0.3, size = 0.8) +
  stat_ellipse(aes(group = clase), type = "t", level = 0.90,
               linewidth = 1.2, linetype = "solid") +
  scale_color_manual(values = PALETA_CLASES) +
  labs(
    title    = "PCA — Espacio latente de métricas de calidad de contigs",
    subtitle = sprintf("PC1 (%.1f%%) × PC2 (%.1f%%) — Elipses = 90%% IC por clase",
                       var_exp[1], var_exp[2]),
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color = "Clase"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

ggsave(file.path(DIR_SALIDA, "4_5_PCA_por_clase.png"),
       p12_pca, width = 11, height = 8, dpi = 150)

# Biplot de loadings
p13_loadings <- fviz_pca_var(res_pca,
                             col.var = "contrib",
                             gradient.cols = c("#00d4ff", "#7c3aed", "#ef4444"),
                             repel = TRUE,
                             title = "PCA Biplot — Contribución de variables a PC1 y PC2") +
  theme_minimal(base_size = 11)

ggsave(file.path(DIR_SALIDA, "4_5_PCA_biplot_loadings.png"),
       p13_loadings, width = 11, height = 9, dpi = 150)


# ── 4.6  Scatter: Interacciones clave (score × p_good por clase) ─────────────
cat("  [4.6] Scatterplots de interacciones clave...\n")

p14_scatter <- df %>%
  select(score, p_good, log1p_coverage, log1p_tpm, final_annotation) %>%
  drop_na() %>%
  slice_sample(n = min(8000, nrow(.))) %>%
  ggplot(aes(x = score, y = p_good, color = final_annotation)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.5) +
  scale_color_manual(values = PALETA_CLASES) +
  facet_wrap(~ final_annotation, nrow = 1) +
  labs(
    title    = "score × p_good — Espacio de verosimilitud por clase",
    subtitle = "score = calidad del ensamblaje · p_good = proporción de bases buenas · Curva = LOESS",
    x = "score (calidad ensamblaje)", y = "p_good (proporción bases buenas)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", color = "grey10"))

ggsave(file.path(DIR_SALIDA, "4_6_scatter_score_pgood.png"),
       p14_scatter, width = 16, height = 5, dpi = 150)


# ── 4.7  GGally Pairs: Matrix de dispersión clave ────────────────────────────
cat("  [4.7] GGally pairs matrix (top 5 features)...\n")

top5 <- resultados_anova$variable[1:5]

p15_pairs <- df %>%
  select(all_of(top5), final_annotation) %>%
  drop_na() %>%
  slice_sample(n = min(3000, nrow(.))) %>%
  ggpairs(
    columns = 1:5,
    mapping = aes(color = final_annotation, alpha = 0.4),
    upper   = list(continuous = "cor"),
    lower   = list(continuous = "points"),
    diag    = list(continuous = "densityDiag")
  ) +
  scale_color_manual(values = PALETA_CLASES) +
  scale_fill_manual(values = PALETA_CLASES) +
  theme_minimal(base_size = 9) +
  labs(title = "GGally — Pares de variables más discriminativas")

ggsave(file.path(DIR_SALIDA, "4_7_ggpairs_top5.png"),
       p15_pairs, width = 14, height = 12, dpi = 150)


# ── 4.8  Interacción prelim_cat × final_annotation (alluvial-like) ───────────
cat("  [4.8] Tabla cruzada prelim_cat × final_annotation...\n")

tabla_cruzada <- df %>%
  count(prelim_cat, final_annotation) %>%
  drop_na() %>%
  group_by(prelim_cat) %>%
  mutate(prop = n / sum(n))

p16_alluvial <- ggplot(tabla_cruzada,
                       aes(x = prelim_cat, y = prop, fill = final_annotation)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(prop > 0.05, scales::percent(prop, 1), "")),
            position = position_stack(vjust = 0.5), size = 3, color = "white",
            fontface = "bold") +
  scale_fill_manual(values = PALETA_CLASES) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title    = "Composición de final_annotation por prelim_cat",
    subtitle = "prelim_cat codifica la cardinalidad de la relación contig→transcripto",
    x = "prelim_cat (cardinalidad)", y = "Proporción", fill = "final_annotation"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(DIR_SALIDA, "4_8_prelim_cat_vs_annotation.png"),
       p16_alluvial, width = 11, height = 7, dpi = 150)


# =============================================================================
# BLOQUE 5 — COMPARACIÓN FINAL DE MODELOS
# =============================================================================

cat("\n══ BLOQUE 5: Comparación de modelos (CV = Fold) ══\n")

# ── 5.1  Matriz de confusión RF ───────────────────────────────────────────────
cat("  [5.1] Matriz de confusión OOB (RF)...\n")

# Predicciones OOB de ranger
pred_oob  <- modelo_rf$predictions  # factor con niveles
true_lab  <- df_rf$final_annotation

if (!is.null(pred_oob) && length(pred_oob) == length(true_lab)) {
  conf_mat <- table(Predicho = pred_oob, Real = true_lab)
  cat("\n  Matriz de confusión OOB:\n")
  print(conf_mat)
  
  # Métricas por clase
  for (cls in levels(true_lab)) {
    tp  <- conf_mat[cls, cls]
    fp  <- sum(conf_mat[cls, ]) - tp
    fn  <- sum(conf_mat[, cls]) - tp
    prec <- tp / (tp + fp)
    rec  <- tp / (tp + fn)
    f1   <- 2 * prec * rec / (prec + rec)
    cat(sprintf("  %8s: Precision=%.3f | Recall=%.3f | F1=%.3f\n",
                cls, prec, rec, f1))
  }
}


# ── 5.2  Tabla resumen comparativa ───────────────────────────────────────────
cat("\n  [5.2] Tabla resumen de resultados...\n")

resumen_final <- tibble(
  Análisis  = c("Pearson > 0.70", "Kruskal-Wallis sig.", "ANOVA η² máx",
                "MI máx (bits)", "RF OOB accuracy", "Top predictor (RF)",
                "Top predictor (MI)", "Top predictor (XGBoost)"),
  Resultado = c(
    paste(nrow(cor_df), "pares"),
    paste(sum(resultados_kw$p_ajust < 0.05, na.rm=TRUE), "variables"),
    sprintf("%.3f (%s)", resultados_anova$eta2[1], resultados_anova$variable[1]),
    sprintf("%.4f (%s)", mi_df$MI[1], mi_df$variable[1]),
    sprintf("%.3f", 1 - modelo_rf$prediction.error),
    imp_df$variable[1],
    mi_df$variable[1],
    xgb_imp$Feature[1]
  )
)

print(resumen_final)
write.csv(resumen_final, file.path(DIR_SALIDA, "5_resumen_resultados.csv"), row.names = FALSE)


# =============================================================================
# BLOQUE 6 — EXPORTAR RESULTADOS
# =============================================================================

cat("\n══ BLOQUE 6: Exportando resultados ══\n")

write.csv(resultados_anova,  file.path(DIR_SALIDA, "6_anova_eta2_ranking.csv"),  row.names = FALSE)
write.csv(mi_df,             file.path(DIR_SALIDA, "6_informacion_mutua.csv"),   row.names = FALSE)
write.csv(imp_df,            file.path(DIR_SALIDA, "6_rf_importancia.csv"),      row.names = FALSE)
write.csv(xgb_imp,           file.path(DIR_SALIDA, "6_xgboost_importancia.csv"), row.names = FALSE)
write.csv(cor_df,            file.path(DIR_SALIDA, "6_pares_alta_correlacion.csv"), row.names = FALSE)

# Guardar modelos
saveRDS(modelo_rf,  file.path(DIR_SALIDA, "modelo_rf.rds"))
xgb.save(modelo_xgb, file.path(DIR_SALIDA, "modelo_xgb.model"))

cat(sprintf("\n✅ Análisis completo. Archivos guardados en: %s/\n", DIR_SALIDA))
cat("   Plots generados:\n")
for (f in list.files(DIR_SALIDA, pattern = "\\.png$")) cat(sprintf("   · %s\n", f))
cat("   CSVs generados:\n")
for (f in list.files(DIR_SALIDA, pattern = "\\.csv$")) cat(sprintf("   · %s\n", f))


# =============================================================================
# NOTA FINAL: Extensiones opcionales
# =============================================================================
#
#  ► UMAP no lineal (reemplaza PCA para visualización): Hard-memory
#     library(umap)
#     um <- umap(X_pca)
#     umap_df <- data.frame(um$layout, clase = df_pca_input$final_annotation)
#     ggplot(umap_df, aes(X1, X2, color=clase)) + geom_point(alpha=0.3)
#
#  ► GAM para relaciones no lineales suavizadas:
#     library(mgcv)
#     mod_gam <- gam(as.numeric(final_annotation=="full") ~ s(score) + s(p_good) +
#                    s(log1p_coverage) + s(prop_gc), data=df, family=binomial)
#     plot(mod_gam, pages=1)
#
#  ► Regresión Lasso/Ridge para selección de features lineales:
#     library(glmnet)
#     cv_lasso <- cv.glmnet(X_mat, y_vec, family="multinomial", alpha=1, nfolds=10)
#     plot(cv_lasso)
#
#  ► Causal Discovery (PC algorithm): errors running
#     library(pcalg)
#     suf_stat <- list(C = cor(X_pca), n = nrow(X_pca))
#     pc_fit   <- pc(suf_stat, indepTest=gaussCItest, alpha=0.01, labels = names(X_pca))
#     plot(pc_fit, main = "")
#
# =============================================================================