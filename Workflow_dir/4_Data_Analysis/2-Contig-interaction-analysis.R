# =============================================================================
#  ANÁLISIS DE INTERACCIONES LINEALES Y NO LINEALES — v2.1 (CONFIGURABLE)
#  Contexto: Clasificación de proteínas/DNA por verosimilitud de ensamblaje
#  Dataset:  contigbaseddatalonger.csv
#  Autor:    [Editar con tu nombre]
#  Fecha:    2026-04
# =============================================================================
#
#  v2.1 NOVEDADES (vs v2.0 — auditoría end-to-end):
#  • DIR_SALIDA se sufija automáticamente con INTERCEPT_VAR
#    -> evita sobrescribir resultados al cambiar de variable de contraste.
#  • Defensa contra colisión: si ya existe una col 'Intercept' se renombra.
#  • top_n adaptativo: ridge plots y ggpairs no fallan si ANOVA tiene <6 filas.
#  • Skip automático de ggpairs cuando Intercept tiene >MAX_NIVELES_GGPAIRS.
#  • Heatmap robusto: z-score = 0 si sd<1e-10, filtra medianas NA.
#  • Matriz de confusión: niveles forzados, métricas con guardas (div/0).
#
#  v2.0 NOVEDADES:
#  • Variable de contraste ("Intercept") configurable globalmente:
#    Permite probar diferentes categóricas (final_annotation, prelim_cat,
#    sf_evidence, tab, Region, summarise, Method...) sin tocar el código.
#  • PALETA_CLASES dinámica: se selecciona automáticamente según INTERCEPT_VAR
#    (predefinida o auto-generada con viridis).
#  • PCA con Estrategia B por defecto (imputación por mediana):
#    Incluye Cys_number, Hydrophobicity, seq_freq, Protein_width, Score_sf,
#    reference_coverage en el análisis.
#  • Nueva variable derivada: log1p_ref_coverage.
#  • reference_coverage incluida en COLS_NUM_BASE.
#  • Validación automática de niveles, varianza y NA/Inf.
#  • Detección automática binario/multi-clase (XGBoost, multinomial).
#  • Bloque 4.8 dinámico: cross-tabulación contra TODAS las secundarias.
#  • Skip inteligente de bloques 4.2/4.3 si Intercept es sf_evidence/tab.
#
#  ESTRUCTURA:
#  BLOQUE 0 — Configuración global y librerías
#  BLOQUE 1 — Carga, transformaciones, renombrado de Intercept
#  BLOQUE 2 — Análisis LINEAL
#  BLOQUE 3 — Análisis NO LINEAL
#  BLOQUE 4 — Interacciones Numéricas × Categóricas
#  BLOQUE 5 — Comparación de modelos
#  BLOQUE 6 — Exportación
# =============================================================================



rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


# =============================================================================
# BLOQUE 0 — CONFIGURACIÓN
# =============================================================================

# ── 0.1  Paquetes necesarios ─────────────────────────────────────────────────
pkgs_needed <- c(
  "tidyverse", "data.table", "janitor",
  "corrplot", "ggcorrplot", "psych", "car", "effectsize",
  "randomForest", "ranger", "xgboost", "caret",
  "infotheo",
  "ggplot2", "ggpubr", "patchwork", "viridis", "RColorBrewer",
  "GGally", "ggridges", "ggforce",
  "FactoMineR", "factoextra",
  "nnet", "broom"
)
# install.packages(setdiff(pkgs_needed, installed.packages()[, "Package"]))

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(janitor)
  library(corrplot); library(ggcorrplot); library(psych)
  library(car); library(effectsize)
  library(randomForest); library(ranger); library(xgboost); library(caret)
  library(infotheo)
  library(ggplot2); library(ggpubr); library(patchwork)
  library(viridis); library(RColorBrewer)
  library(GGally); library(ggridges)
  library(FactoMineR); library(factoextra)
  library(nnet); library(broom)
})


# ── 0.2  PARÁMETROS GLOBALES (EDITAR AQUÍ) ──────────────────────────────────
setwd("~/Documents/GitHub/ConotoxinBenchmark/INPUTS/")

RUTA_CSV        <- "contig-based-data-longer.csv"   # Ruta al archivo CSV
DIR_SALIDA_BASE<- "resultados_interacciones"   # se sufija con INTERCEPT_VAR
SEMILLA        <- 42
N_ARBOLES_RF   <- 500
TOP_N_FEATURES <- 15
N_SAMPLE_PLOT  <- 10000   # nº máx. de puntos para scatter/PCA (rendimiento)
MAX_NIVELES_GGPAIRS <- 8  # ggpairs ilegible si Intercept tiene más niveles
MAX_NIVELES_FACETS  <- 12 # facets ilegibles si Intercept tiene más niveles

# ── INTERCEPT (variable de contraste) ────────────────────────────────────────
# Cambia INTERCEPT_VAR y INTERCEPT_REF para reorientar TODO el análisis.
# Niveles válidos para INTERCEPT_REF según la variable elegida:
#   final_annotation : full · multi · chimera · error · fragment
#   prelim_cat       : 1->1 · 1->n · m->1 · m->n
#   sf_evidence      : Concordance · Ambiguous
#   tab              : Regex · pHMM
#   Region           : (Mature) · (Pro-region) · (Signal) · combinaciones
#   summarise        : 100% / >=95% / >=90% / >=80% / <80% alignment
#   Method           : 11 ensambladores (TRINITY, SOAPDENOVO, etc.)

INTERCEPT_VAR <- "final_annotation"   # <-- variable categórica de contraste
INTERCEPT_REF <- "full"               # <-- nivel de referencia para relevel()

# ── PCA: estrategia ──────────────────────────────────────────────────────────
# "A" = solo columnas con 0% NA  (más estable, menos features)
# "B" = imputación por mediana    (incluye Cys_number, Hydrophobicity, etc.)  ★ DEFAULT
PCA_STRATEGY <- "B"

# ── Paletas predefinidas por variable (auto-generadas si falta una) ──────────
PALETAS_PREDEFINIDAS <- list(
  final_annotation = c(
    "Full"     = "#10b981",
    "Multi"    = "#00d4ff",
    "Chimera"  = "#ef4444",
    "Error"    = "#f59e0b",
    "Fragment" = "#a78bfa",
    "No-hit" = "#f59e0b"
  ),
  prelim_cat = c(
    "1->1" = "#10b981",
    "1->n" = "#00d4ff",
    "m->1" = "#f59e0b",
    "m->n" = "#a78bfa"
  ),
  sf_evidence = c(
    "Concordance" = "#10b981",
    "Ambiguous"   = "#ef4444"
  ),
  tab = c(
    "Regex" = "#00d4ff",
    "pHMM"  = "#7c3aed"
  ),
  summarise = c(
    "100% alignment"   = "#10b981",
    ">= 95% alignment" = "#84cc16",
    ">= 90% alignment" = "#f59e0b",
    ">= 80% alignment" = "#fb923c",
    "< 80 % alignment" = "#ef4444"
  )
  # Region, Method, gs_conoSorter, gs_conoServer → auto-generadas con viridis
)

recode_col <- c("full","multi", "fragment","chimera", "NaN")

recode_col <- structure(c("Full","Multi", "Fragment", "Chimera", "No-hit"), names = recode_col)


# MEJORA v2.1: subdirectorio por INTERCEPT_VAR + INTERCEPT_REF para no
# sobrescribir runs. Cada experimento tendrá su propia carpeta:
#   resultados_interacciones/final_annotation_full/
#   resultados_interacciones/final_annotation_multi/
#   resultados_interacciones/prelim_cat_1to1/
ref_safe   <- gsub("[^A-Za-z0-9]+", "", INTERCEPT_REF)  # sanitiza "1->1" -> "11"
DIR_SALIDA <- file.path(DIR_SALIDA_BASE,
                        sprintf("%s_%s", INTERCEPT_VAR, ref_safe))
if (!dir.exists(DIR_SALIDA)) dir.create(DIR_SALIDA, recursive = TRUE)
set.seed(SEMILLA)

cat(strrep("=", 70), "\n", sep = "")
cat("  ANÁLISIS DE INTERACCIONES — v2.1\n")
cat(sprintf("  INTERCEPT      : %s  (ref = %s)\n", INTERCEPT_VAR, INTERCEPT_REF))
cat(sprintf("  PCA strategy   : %s  (A=sin imputar · B=imputar mediana)\n", PCA_STRATEGY))
cat(sprintf("  Salida         : %s/\n", DIR_SALIDA))
cat(strrep("=", 70), "\n\n", sep = "")


# =============================================================================
# BLOQUE 1 — CARGA Y PREPROCESAMIENTO
# =============================================================================

cat("-- BLOQUE 1: Carga y preprocesamiento --\n")
df_raw <- fread(RUTA_CSV, na.strings = c("", "NA", "N/A", "NaN"))


df_raw |>
  drop_na(p_good) |>
  mutate(final_annotation = ifelse(is.na(final_annotation), "No-hit", final_annotation)) |>
  dplyr::mutate(final_annotation = dplyr::recode_factor(final_annotation, !!!recode_col)) |>
  ggplot(aes(x = p_good, y = final_annotation, fill=final_annotation)) + # model: p_good+score * reference_coverage
  # facet_wrap(~ Method) +ggplot2::stat_ecdf()
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.5, point_alpha = 1, alpha = 0.7) +
  scale_fill_manual(values = PALETAS_PREDEFINIDAS$final_annotation)

df_raw |>
  drop_na(p_good) |>
  mutate(final_annotation = ifelse(is.na(final_annotation), "No-hit", final_annotation)) |>
  dplyr::mutate(final_annotation = dplyr::recode_factor(final_annotation, !!!recode_col)) |>
  ggplot(aes(x = log10(in_bridges+1), y = final_annotation, fill=final_annotation)) + # model: p_good+score * reference_coverage
  # facet_wrap(~ Method) +ggplot2::stat_ecdf()
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.5, point_alpha = 1, alpha = 0.7) +
  scale_fill_manual(values = PALETAS_PREDEFINIDAS$final_annotation)

# ── 1.1  Definir columnas por tipo ───────────────────────────────────────────
COLS_NUMERIC <- c(
  # Quality metrics (0% NA)
  "n_hits" ,"n_unique_subjects","max_contigs_per_subject","best_bitscore","best_pident","best_max_coverage",
  "at_skew", "coverage", "cpg_count", "cpg_ratio", "eff_count", "eff_length",
  "gc_skew", "in_bridges", "length", "linguistic_complexity_6", "orf_length",
  "p_bases_covered", "p_good", "p_not_segmented", "p_seq_true",
  "prop_gc", "score", "tpm",
  # ~29% NA
  "reference_coverage",
  # ~55% NA (conotoxin-specific)
  "Cys_number", "Hydrophobicity", "seq_freq", "Protein_width", "Score_sf"
)

COLS_CATEG_TODAS <- c(
  "final_annotation", "prelim_cat", "Region", "summarise",
  "organismlatin","organismdiet","organismregion",
  "Method", "gs_conoSorter", "gs_conoServer", "sf_evidence", "tab"
)

 
# ── 1.2  VALIDACIÓN del INTERCEPT (early-fail) ───────────────────────────────
if (!INTERCEPT_VAR %in% COLS_CATEG_TODAS) {
  stop(sprintf(
    "INTERCEPT_VAR = '%s' no es una columna categórica válida.\n  Opciones: %s",
    INTERCEPT_VAR, paste(COLS_CATEG_TODAS, collapse = ", ")
  ))
}
if (!INTERCEPT_VAR %in% names(df_raw)) {
  stop(sprintf("Columna '%s' no existe en %s", INTERCEPT_VAR, RUTA_CSV))
}

niveles_disponibles <- as.character(df_raw[[INTERCEPT_VAR]])
niveles_disponibles <- unique(niveles_disponibles[!is.na(niveles_disponibles)])

if (!INTERCEPT_REF %in% niveles_disponibles) {
  stop(sprintf(
    "INTERCEPT_REF = '%s' no existe en columna '%s'.\n  Niveles válidos: %s",
    INTERCEPT_REF, INTERCEPT_VAR, paste(niveles_disponibles, collapse = " · ")
  ))
}

cat(sprintf("  Niveles de %s: %s\n",
            INTERCEPT_VAR, paste(niveles_disponibles, collapse = " · ")))


# ── 1.3  Conversión de tipos + transformaciones logarítmicas ─────────────────
df <- df_raw |>
  as_tibble() |>
  mutate(reference_coverage = ifelse(is.na(reference_coverage), 0, reference_coverage)) |>
  mutate(across(all_of(intersect(COLS_NUMERIC, names(df_raw))), as.numeric)) |>
  mutate(across(all_of(intersect(COLS_CATEG_TODAS, names(df_raw))), as.factor)) |>
  # log1p para variables muy sesgadas (todos los conteos/expresión)
  mutate(
    log1p_tpm          = log1p(tpm),
    log1p_coverage     = log1p(coverage),
    log1p_eff_count    = log1p(eff_count),
    log1p_in_bridges   = log1p(in_bridges),
    log1p_cpg_count    = log1p(cpg_count),
    log1p_ref_coverage = log1p(reference_coverage)   # NUEVO
  )


# ── 1.4  RENOMBRAR INTERCEPT_VAR -> "Intercept" ─────────────────────────────
# Esto permite que TODOS los bloques usen `Intercept` sin importar
# qué variable categórica eligió el usuario.

# Defensa: evitar colisión si ya existiera una columna llamada "Intercept"
if ("Intercept" %in% names(df) && INTERCEPT_VAR != "Intercept") {
  names(df)[names(df) == "Intercept"] <- "Intercept_orig"
  warning("Columna preexistente 'Intercept' renombrada a 'Intercept_orig'")
}

names(df)[names(df) == INTERCEPT_VAR] <- "Intercept"
df$Intercept <- relevel(factor(df$Intercept), ref = INTERCEPT_REF)
df <- df |> filter(!is.na(Intercept))

cat(sprintf("  Renombrado: %s -> 'Intercept' (ref = %s)\n", INTERCEPT_VAR, INTERCEPT_REF))
cat(sprintf("  Dataset post-filtro: %d filas x %d columnas\n", nrow(df), ncol(df)))


# ── 1.5  Construir PALETA_CLASES dinámica ────────────────────────────────────
niveles_intercept <- levels(df$Intercept)
n_niveles         <- length(niveles_intercept)
ES_BINARIO        <- (n_niveles == 2)

if (INTERCEPT_VAR %in% names(PALETAS_PREDEFINIDAS)) {
  PALETA_CLASES <- PALETAS_PREDEFINIDAS[[INTERCEPT_VAR]]
  faltantes <- setdiff(niveles_intercept, names(PALETA_CLASES))
  if (length(faltantes) > 0) {
    extras <- viridis::viridis(length(faltantes), option = "plasma", end = 0.85)
    names(extras) <- faltantes
    PALETA_CLASES <- c(PALETA_CLASES, extras)
  }
  PALETA_CLASES <- PALETA_CLASES[niveles_intercept]
} else {
  PALETA_CLASES <- viridis::viridis(n_niveles, option = "plasma", end = 0.85)
  names(PALETA_CLASES) <- niveles_intercept
}

cat(sprintf("  PALETA: %d colores (%s)\n",
            length(PALETA_CLASES),
            ifelse(INTERCEPT_VAR %in% names(PALETAS_PREDEFINIDAS),
                   "predefinida", "viridis-plasma auto")))

if (n_niveles > 12) {
  warning(sprintf(
    "Intercept tiene %d niveles. Algunas visualizaciones (ridge, scatter facetado) podrían ser ilegibles.\n",
    n_niveles))
}


# ── 1.6  COLS_NUM_BASE (todas las numéricas + log1p) ─────────────────────────
COLS_NUM_BASE <- c(
  # "at_skew", "coverage", "cpg_count", "cpg_ratio", "eff_count",
  # "eff_length", "gc_skew", "in_bridges", "length",
  # "linguistic_complexity_6", "orf_length", "p_bases_covered",
  # "p_good", "p_not_segmented", "p_seq_true", "prop_gc",
  # "score", "tpm", "reference_coverage",
  COLS_NUMERIC,
  "log1p_tpm", "log1p_coverage", "log1p_eff_count",
  "log1p_in_bridges", "log1p_cpg_count", "log1p_ref_coverage"
)

# Columnas SIN NA en el dataset (verificar empíricamente)
COLS_SIN_NA <- COLS_NUM_BASE[
  vapply(COLS_NUM_BASE,
         function(c) !any(is.na(df[[c]])) && !any(is.infinite(df[[c]])),
         logical(1))
]
cat(sprintf("  Numéricas: %d totales · %d sin NA\n",
            length(COLS_NUM_BASE), length(COLS_SIN_NA)))


# =============================================================================
# BLOQUE 2 — ANÁLISIS LINEAL
# =============================================================================

cat("\n== BLOQUE 2: Análisis Lineal ==\n")

# ── 2.1  Correlación de Pearson ──────────────────────────────────────────────
cat("  [2.1] Correlación de Pearson...\n")

mat_cor <- df |> select(all_of(COLS_NUM_BASE)) |>
  cor(use = "pairwise.complete.obs", method = "pearson")

p1_corrplot <- ggcorrplot(
  mat_cor, method = "square", type = "upper",
  lab = TRUE, lab_size = 2.2,
  colors = c("#ef4444", "#0f1117", "#10b981"),
  outline.color = "grey20",
  title = "Correlación de Pearson - Métricas numéricas",
  ggtheme = theme_minimal(base_size = 11)
) +
  theme(
    plot.title = element_text(size = 12, face = "bold", color = "grey10"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8), legend.position = "bottom"
  )

ggsave(file.path(DIR_SALIDA, "2_1_correlacion_pearson.png"),
       p1_corrplot, width = 14, height = 12, dpi = 150)

cor_df <- mat_cor |> as.data.frame() |> rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  filter(var1 < var2, abs(r) > 0.70) |> arrange(desc(abs(r)))

cat(sprintf("    %d pares con |r| > 0.70\n", nrow(cor_df)))


# ── 2.2  Correlación de Spearman ─────────────────────────────────────────────
cat("  [2.2] Correlación de Spearman...\n")

mat_spear <- df |> select(all_of(COLS_NUM_BASE)) |>
  cor(use = "pairwise.complete.obs", method = "spearman")

p2_spearman <- ggcorrplot(
  mat_spear, method = "circle", type = "upper", lab = FALSE,
  colors = c("#7c3aed", "#f8fafc", "#f59e0b"),
  title = "Correlación de Spearman (robusta a outliers)",
  ggtheme = theme_minimal(base_size = 10)
)

ggsave(file.path(DIR_SALIDA, "2_2_correlacion_spearman.png"),
       p2_spearman, width = 12, height = 10, dpi = 150)


# ── 2.3  ANOVA: Numéricas ~ Intercept ────────────────────────────────────────
cat("  [2.3] ANOVA Numéricas ~ Intercept...\n")

resultados_anova <- map_dfr(COLS_NUM_BASE, function(col_nm) {
  datos_col <- df |>
    select(y = !!sym(col_nm), grupo = Intercept) |> drop_na()
  
  if (nrow(datos_col) < 50 || length(unique(datos_col$grupo)) < 2) return(NULL)
  
  mod <- aov(y ~ grupo, data = datos_col)
  sm  <- summary(mod)[[1]]
  eta2 <- tryCatch(eta_squared(mod, partial = FALSE)$Eta2[1],
                   error = function(e) NA_real_)
  
  tibble(
    variable = col_nm,
    F_stat   = sm[["F value"]][1],
    p_valor  = sm[["Pr(>F)"]][1],
    eta2     = eta2,
    sig      = case_when(
      sm[["Pr(>F)"]][1] < 0.001 ~ "***",
      sm[["Pr(>F)"]][1] < 0.01  ~ "**",
      sm[["Pr(>F)"]][1] < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )
}) |>
  arrange(desc(eta2)) |>
  mutate(p_ajustado = p.adjust(p_valor, method = "BH"))

cat(sprintf("    Top 1: %s (eta2 = %.3f)\n",
            resultados_anova$variable[1], resultados_anova$eta2[1]))

p3_anova <- resultados_anova |>
  mutate(variable = fct_reorder(variable, eta2)) |>
  ggplot(aes(x = eta2, y = variable, fill = eta2)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sig), hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = sprintf("ANOVA: Numéricas vs. %s", INTERCEPT_VAR),
    subtitle = "Tamaño de efecto eta cuadrado - separación entre clases",
    x = "Eta cuadrado", y = NULL,
    caption = "*** p<0.001 · ** p<0.01 · * p<0.05 · n.s. = no significativo (FDR Benjamini-Hochberg)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        panel.grid.major.y = element_blank())

ggsave(file.path(DIR_SALIDA, "2_3_anova_eta2.png"),
       p3_anova, width = 12, height = 8, dpi = 150)


# ── 2.4  Regresión Multinomial Logística ─────────────────────────────────────
# (multinom() también funciona con 2 niveles = regresión logística clásica)
cat("  [2.4] Regresión Multinomial Logística...\n")

df_multinom <- df |>
  select(Intercept, all_of(COLS_NUM_BASE)) |>
  drop_na(Intercept, score, p_good, log1p_coverage, log1p_tpm, orf_length,
          prop_gc, linguistic_complexity_6, p_not_segmented, p_seq_true,
          p_bases_covered, at_skew, gc_skew, cpg_ratio, log1p_ref_coverage)

mod_multinom <- multinom(
  Intercept ~ score + p_good + log1p_coverage + log1p_tpm +
    orf_length + prop_gc + linguistic_complexity_6 +
    p_not_segmented + p_seq_true + p_bases_covered +
    at_skew + gc_skew + cpg_ratio + log1p_ref_coverage,
  data    = df_multinom,
  MaxNWts = 5000, maxit = 200, trace = FALSE
)

coef_tabla <- tidy(mod_multinom, conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  arrange(p.value) |>
  mutate(OR = exp(estimate), sig = ifelse(p.value < 0.05, "*", ""))

# Para regresión binaria, broom no incluye y.level → agregar manualmente
if (!"y.level" %in% names(coef_tabla)) {
  no_ref_level <- setdiff(levels(df_multinom$Intercept), INTERCEPT_REF)[1]
  coef_tabla$y.level <- no_ref_level
}

write.csv(coef_tabla,
          file.path(DIR_SALIDA, "2_4_coeficientes_multinomial.csv"),
          row.names = FALSE)
cat(sprintf("    %d coeficientes (%d significativos)\n",
            nrow(coef_tabla), sum(coef_tabla$p.value < 0.05, na.rm = TRUE)))

p4_OR <- coef_tabla |>
  filter(abs(estimate) > 0.1, !is.na(conf.low), !is.na(conf.high)) |>
  mutate(term = fct_reorder(term, estimate)) |>
  ggplot(aes(x = estimate, y = term, color = y.level)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3) +
  facet_wrap(~ y.level, scales = "free_x",
             nrow = ifelse(ES_BINARIO, 1, 2)) +
  scale_color_manual(values = PALETA_CLASES, na.value = "grey50") +
  labs(
    title    = sprintf("Multinomial: log-odds vs. clase de referencia '%s'", INTERCEPT_REF),
    subtitle = sprintf("IC 95%% · Variable de contraste: %s", INTERCEPT_VAR),
    x = "Coeficiente (log-odds)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "2_4_multinomial_coeficientes.png"),
       p4_OR, width = ifelse(ES_BINARIO, 9, 15), height = 8, dpi = 150)


# ── 2.5  Kruskal-Wallis ──────────────────────────────────────────────────────
cat("  [2.5] Kruskal-Wallis...\n")

resultados_kw <- map_dfr(COLS_NUM_BASE, function(col_nm) {
  datos_col <- df |> select(y = !!sym(col_nm), grupo = Intercept) |> drop_na()
  if (nrow(datos_col) < 30 || length(unique(datos_col$grupo)) < 2) return(NULL)
  kw <- kruskal.test(y ~ grupo, data = datos_col)
  tibble(variable = col_nm,
         statistic = kw$statistic,
         p_valor   = kw$p.value)
}) |>
  arrange(p_valor) |>
  mutate(p_ajust = p.adjust(p_valor, method = "BH"))

write.csv(resultados_kw,
          file.path(DIR_SALIDA, "2_5_kruskal_wallis.csv"),
          row.names = FALSE)
cat(sprintf("    %d variables con p-ajustado < 0.05\n",
            sum(resultados_kw$p_ajust < 0.05, na.rm = TRUE)))


# =============================================================================
# BLOQUE 3 — ANÁLISIS NO LINEAL
# =============================================================================

cat("\n== BLOQUE 3: Análisis No Lineal ==\n")

# ── 3.1  Información Mutua ────────────────────────────────────────────────────
cat("  [3.1] Información Mutua Numéricas x Intercept...\n")

# Solo columnas sin NA para discretización limpia
df_mi_input <- df |>
  select(all_of(COLS_SIN_NA), Intercept) |>
  drop_na()

df_disc <- df_mi_input |>
  select(all_of(COLS_SIN_NA)) |>
  mutate(across(everything(), ~ {
    res <- discretize(., disc = "equalfreq", nbins = 10)
    as.integer(res[[1]])
  }))

target_disc <- as.integer(df_mi_input$Intercept)

mi_scores <- map_dbl(COLS_SIN_NA, function(col_nm) {
  mutinformation(df_disc[[col_nm]], target_disc)
})

mi_df <- tibble(variable = COLS_SIN_NA, MI = mi_scores) |> arrange(desc(MI))

p5_MI <- mi_df |>
  mutate(variable = fct_reorder(variable, MI)) |>
  ggplot(aes(x = MI, y = variable, fill = MI)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(option = "magma", direction = -1, begin = 0.2) +
  labs(
    title    = sprintf("Información Mutua: Numéricas -> %s", INTERCEPT_VAR),
    subtitle = "MI captura relaciones lineales Y no lineales",
    x = "Información Mutua (bits)", y = NULL,
    caption = "Discretización: 10 bins de frecuencia igual"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "3_1_informacion_mutua.png"),
       p5_MI, width = 11, height = 8, dpi = 150)
cat(sprintf("    Top 1 MI: %s (%.4f bits)\n", mi_df$variable[1], mi_df$MI[1]))


# ── 3.2  Random Forest (ranger) ──────────────────────────────────────────────
cat("  [3.2] Random Forest - importancia por permutación...\n")

df_rf <- df |>
  select(Intercept, all_of(COLS_NUM_BASE)) |>
  drop_na(Intercept, score, p_good, log1p_coverage, log1p_tpm, orf_length,
          prop_gc, linguistic_complexity_6, p_not_segmented, p_seq_true,
          p_bases_covered, at_skew, gc_skew, cpg_ratio, log1p_ref_coverage)

# Pesos por imbalance (inversamente proporcionales a frecuencia)
class_w <- as.numeric(table(df_rf$Intercept))
class_w <- (1 / class_w) / sum(1 / class_w)
names(class_w) <- levels(df_rf$Intercept)

modelo_rf <- ranger(
  formula        = Intercept ~ .,
  data           = df_rf,
  num.trees      = N_ARBOLES_RF,
  importance     = "permutation",
  classification = TRUE,
  class.weights  = class_w,
  seed           = SEMILLA,
  num.threads    = max(1, parallel::detectCores() - 1)
)

cat(sprintf("    OOB accuracy = %.3f\n", 1 - modelo_rf$prediction.error))

imp_df <- tibble(
  variable    = names(modelo_rf$variable.importance),
  importancia = modelo_rf$variable.importance
) |> arrange(desc(importancia))

p6_RF_imp <- imp_df |>
  head(TOP_N_FEATURES) |>
  mutate(variable = fct_reorder(variable, importancia)) |>
  ggplot(aes(x = importancia, y = variable, fill = importancia)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(option = "cividis", direction = -1) +
  labs(
    title    = "Random Forest: importancia por permutación",
    subtitle = sprintf("ranger · %d árboles · OOB acc = %.3f · target = %s",
                       N_ARBOLES_RF, 1 - modelo_rf$prediction.error, INTERCEPT_VAR),
    x = "Caída en accuracy por permutación", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "3_2_rf_importancia_permutacion.png"),
       p6_RF_imp, width = 11, height = 7, dpi = 150)


# ── 3.3  XGBoost (auto-detecta binario vs multi-clase) ───────────────────────
cat("  [3.3] XGBoost...\n")

features_xgb <- c("score", "p_good", "log1p_coverage", "log1p_tpm",
                  "log1p_ref_coverage", "reference_coverage",
                  "orf_length", "prop_gc", "linguistic_complexity_6",
                  "p_not_segmented", "p_seq_true", "p_bases_covered",
                  "at_skew", "gc_skew", "cpg_ratio", "log1p_in_bridges",
                  "log1p_eff_count", "log1p_cpg_count", "eff_length")

df_xgb <- df |> select(Intercept, all_of(features_xgb)) |> drop_na()
X_mat  <- df_xgb |> select(-Intercept) |> as.matrix()
y_vec  <- as.integer(df_xgb$Intercept) - 1   # 0-indexed
dtrain <- xgb.DMatrix(data = X_mat, label = y_vec)

n_clases <- length(levels(df_xgb$Intercept))

params_xgb <- list(
  objective        = ifelse(ES_BINARIO, "binary:logistic", "multi:softprob"),
  eta              = 0.1, max_depth = 6, subsample = 0.8,
  colsample_bytree = 0.8,
  eval_metric      = ifelse(ES_BINARIO, "logloss", "mlogloss"),
  seed             = SEMILLA
)
if (!ES_BINARIO) params_xgb$num_class <- n_clases

modelo_xgb <- xgb.train(params = params_xgb, data = dtrain,
                        nrounds = 300, verbose = 0)

xgb_imp <- xgb.importance(model = modelo_xgb) |> as_tibble()

p7_xgb_imp <- xgb_imp |>
  head(TOP_N_FEATURES) |>
  mutate(Feature = fct_reorder(Feature, Gain)) |>
  pivot_longer(c(Gain, Cover, Frequency), names_to = "tipo", values_to = "valor") |>
  ggplot(aes(x = valor, y = Feature, fill = tipo)) +
  geom_col(position = "dodge", width = 0.6) +
  facet_wrap(~tipo, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = sprintf("XGBoost - predicción de %s", INTERCEPT_VAR),
    subtitle = sprintf("Modo: %s · Gain = reducción de impureza total (más relevante)",
                       ifelse(ES_BINARIO, "binary:logistic", "multi:softprob")),
    x = "Valor", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "3_3_xgboost_importancia.png"),
       p7_xgb_imp, width = 14, height = 7, dpi = 150)

# ── SHAP (opcional, ver NOTA FINAL) ──────────────────────────────────────────


# =============================================================================
# BLOQUE 4 — INTERACCIONES NUMÉRICAS × CATEGÓRICAS
# =============================================================================

cat("\n== BLOQUE 4: Interacciones Numéricas x Categóricas ==\n")

# ── 4.1  Ridge plots por Intercept ───────────────────────────────────────────
cat("  [4.1] Ridge plots - top features por eta2...\n")

# Defensa: tomar como máximo las filas disponibles
n_top_ridge <- min(6, nrow(resultados_anova))
top6_vars   <- resultados_anova$variable[seq_len(n_top_ridge)]

df_ridge <- df |> select(Intercept, all_of(top6_vars)) |> drop_na() |>
  pivot_longer(-Intercept, names_to = "variable", values_to = "valor") |>
  filter(is.finite(valor))

p8_ridge <- ggplot(df_ridge, aes(x = valor, y = Intercept, fill = Intercept)) +
  geom_density_ridges(alpha = 0.75, scale = 1.5, rel_min_height = 0.01,
                      show.legend = FALSE) +
  facet_wrap(~ variable, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = PALETA_CLASES) +
  labs(
    title    = sprintf("Distribución de features discriminativas por %s", INTERCEPT_VAR),
    subtitle = "Top 6 variables por ANOVA eta2 - separación = poder discriminativo",
    x = "Valor", y = INTERCEPT_VAR
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "4_1_ridge_plots_por_clase.png"),
       p8_ridge, width = 14, height = 9, dpi = 150)


# ── 4.2  Boxplots × sf_evidence (omitir si sf_evidence es Intercept) ─────────
cat("  [4.2] Violins x sf_evidence...\n")
vars_box <- c("score", "p_good", "log1p_coverage", "linguistic_complexity_6",
              "prop_gc", "orf_length")

if (INTERCEPT_VAR == "sf_evidence") {
  cat("    OMITIDO - Intercept ya es sf_evidence (ver Block 4.1)\n")
} else if (!"sf_evidence" %in% names(df)) {
  cat("    OMITIDO - sf_evidence no disponible\n")
} else {
  df_box_sf <- df |> select(sf_evidence, all_of(vars_box)) |>
    drop_na(sf_evidence) |>
    pivot_longer(-sf_evidence, names_to = "variable", values_to = "valor")
  
  p9_box_sf <- ggplot(df_box_sf, aes(x = sf_evidence, y = valor, fill = sf_evidence)) +
    geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, vjust = -0.5) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = PALETAS_PREDEFINIDAS$sf_evidence) +
    labs(title = "Métricas numéricas x sf_evidence",
         subtitle = "Wilcoxon · Concordance vs Ambiguous",
         x = NULL, y = "Valor", fill = "sf_evidence") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave(file.path(DIR_SALIDA, "4_2_violin_sf_evidence.png"),
         p9_box_sf, width = 14, height = 9, dpi = 150)
}


# ── 4.3  Boxplots × tab (omitir si tab es Intercept) ─────────────────────────
cat("  [4.3] Violins x tab...\n")

if (INTERCEPT_VAR == "tab") {
  cat("    OMITIDO - Intercept ya es tab (ver Block 4.1)\n")
} else if (!"tab" %in% names(df)) {
  cat("    OMITIDO - tab no disponible\n")
} else {
  df_box_tab <- df |> select(tab, all_of(vars_box)) |>
    drop_na(tab) |>
    pivot_longer(-tab, names_to = "variable", values_to = "valor")
  
  p10_box_tab <- ggplot(df_box_tab, aes(x = tab, y = valor, fill = tab)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = PALETAS_PREDEFINIDAS$tab) +
    labs(title = "Métricas numéricas x tab (Regex vs pHMM)",
         subtitle = "Método de detección de señal peptídica",
         x = NULL, y = "Valor", fill = "tab") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave(file.path(DIR_SALIDA, "4_3_violin_tab_regex_pHMM.png"),
         p10_box_tab, width = 14, height = 9, dpi = 150)
}


# ── 4.4  Heatmap perfiles de clase (Z-score) ─────────────────────────────────
cat("  [4.4] Heatmap de perfiles de clase...\n")

COLS_NUM_BASE <- COLS_NUM_BASE %in% 

heatmap_data <- df |>
  select(Intercept, all_of(COLS_NUM_BASE)) |>
  group_by(Intercept) |>
  summarise(across(everything(), ~ median(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(-Intercept, names_to = "variable", values_to = "mediana") |>
  group_by(variable) |>
  mutate(
    # Defensa: scale() devuelve NaN si sd=0 o todo NA
    z_score = {
      v   <- mediana
      mu  <- mean(v, na.rm = TRUE)
      sd_ <- stats::sd(v, na.rm = TRUE)
      if (is.na(sd_) || sd_ < 1e-10) rep(0, length(v)) else (v - mu) / sd_
    }
  ) |>
  ungroup() |>
  filter(!is.na(mediana))   # quitar combinaciones (clase × variable) sin datos

p11_heatmap <- 
  heatmap_data |>
  dplyr::mutate(Intercept = dplyr::recode_factor(Intercept, !!!recode_col)) |>
  ggplot(aes(x = Intercept, y = variable, fill = z_score)) +
  geom_tile(color = "grey20", linewidth = 0.3) +
  geom_text(aes(label = round(mediana, 2)), size = 2.5, color = "white") +
  scale_fill_gradient2(low = "#1e40af", mid = "#0f1117", high = "#dc2626",
                       midpoint = 0, name = "Z-score") +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  labs(
    title    = sprintf("Perfiles por %s - mediana de métricas (Z-score)", INTERCEPT_VAR),
    subtitle = "Rojo = alto relativo · Azul = bajo · números = mediana cruda",
    x = INTERCEPT_VAR, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 9))


ggsave(file.path(DIR_SALIDA, "4_4_heatmap_perfiles_clase.png"),
       p11_heatmap, width = 5, height = 7, dpi = 150)

heatmap_data |>
  dplyr::mutate(Intercept = dplyr::recode_factor(Intercept, !!!recode_col)) |>
  write.csv(file.path(DIR_SALIDA, "4_4_heatmap_perfiles_clase.csv"),
            row.names = FALSE)


# ── 4.5  PCA — Estrategia B (DEFAULT, con imputación por mediana) ────────────
cat(sprintf("  [4.5] PCA - Estrategia %s...\n", PCA_STRATEGY))

if (PCA_STRATEGY == "B") {
  # ESTRATEGIA B: incluye Cys_number, Hydrophobicity, etc. con imputación
  COLS_PCA <- c(
    "score", "p_good", "p_not_segmented", "p_seq_true", "p_bases_covered",
    "log1p_tpm", "log1p_coverage", "log1p_eff_count", "log1p_in_bridges",
    "log1p_cpg_count", "orf_length", "prop_gc", "linguistic_complexity_6",
    "at_skew", "gc_skew", "cpg_ratio", "eff_length",
    # Conotoxina (~55% NA - imputar por mediana)
    "Cys_number", "Hydrophobicity", "seq_freq", "Protein_width", "Score_sf",
    # reference_coverage (~29% NA - imputar)
    "log1p_ref_coverage"
  )
  
  df_pca_input <- df |>
    select(Intercept, all_of(COLS_PCA)) |>
    filter(!is.na(Intercept))
  
  # Imputación por mediana columna a columna (con manejo de Inf)
  X_pca_raw <- df_pca_input |> select(-Intercept) |>
    mutate(across(everything(), ~ {
      .x[!is.finite(.x)] <- NA
      med <- median(.x, na.rm = TRUE)
      ifelse(is.na(.x), med, .x)
    }))
  
} else {
  # ESTRATEGIA A: solo columnas con 0% NA
  COLS_PCA <- intersect(c(
    "score", "p_good", "p_not_segmented", "p_seq_true", "p_bases_covered",
    "log1p_tpm", "log1p_coverage", "log1p_eff_count", "log1p_in_bridges",
    "log1p_cpg_count", "orf_length", "prop_gc", "linguistic_complexity_6",
    "at_skew", "gc_skew", "cpg_ratio", "eff_length"
  ), COLS_SIN_NA)
  
  df_pca_input <- df |>
    select(Intercept, all_of(COLS_PCA)) |>
    filter(!is.na(Intercept)) |>
    filter(if_all(all_of(COLS_PCA), ~ is.finite(.x)))
  
  X_pca_raw <- df_pca_input |> select(-Intercept)
}

# Eliminar columnas de varianza ~0
cols_var_ok <- vapply(X_pca_raw, function(x) sd(x, na.rm = TRUE) > 1e-10, logical(1))
X_pca_raw <- X_pca_raw[, cols_var_ok, drop = FALSE]

# Escalar (z-score)
X_pca <- as.data.frame(scale(X_pca_raw, center = TRUE, scale = TRUE))

# Verificación defensiva
cols_ok <- colSums(!is.finite(as.matrix(X_pca))) == 0
X_pca   <- X_pca[, cols_ok, drop = FALSE]
rows_ok <- complete.cases(X_pca) & apply(X_pca, 1, function(r) all(is.finite(r)))
X_pca   <- X_pca[rows_ok, , drop = FALSE]
labels_pca <- df_pca_input$Intercept[rows_ok]

cat(sprintf("    Input: %d filas x %d cols\n", nrow(X_pca), ncol(X_pca)))
cat(sprintf("    Cols incluidas: %s\n",
            paste(names(X_pca), collapse = ", ")))

# PCA
res_pca <- prcomp(X_pca, center = TRUE, scale. = FALSE)
var_exp <- round(100 * res_pca$sdev^2 / sum(res_pca$sdev^2), 1)
cat(sprintf("    PC1: %.1f%% · PC2: %.1f%% · PC3: %.1f%%\n",
            var_exp[1], var_exp[2], var_exp[3]))

pca_coords <- as.data.frame(res_pca$x[, 1:3]) |> mutate(clase = labels_pca)

set.seed(SEMILLA)

p12_pca <- ggplot(pca_coords |> slice_sample(n = min(N_SAMPLE_PLOT, nrow(pca_coords))),
                  aes(x = PC1, y = PC2, color = clase)) +
  geom_point(alpha = 0.3, size = 0.8) +
  stat_ellipse(aes(group = clase), type = "t", level = 0.90, linewidth = 1.2) +
  scale_color_manual(values = PALETA_CLASES) +
  labs(
    title    = sprintf("PCA - Espacio latente · color = %s", INTERCEPT_VAR),
    subtitle = sprintf("Estrategia %s · PC1 (%.1f%%) x PC2 (%.1f%%) · Elipses 90%%",
                       PCA_STRATEGY, var_exp[1], var_exp[2]),
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]), color = INTERCEPT_VAR
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(DIR_SALIDA, "4_5_PCA_por_clase.png"),
       p12_pca, width = 11, height = 8, dpi = 150)

# Biplot
p13_loadings <- fviz_pca_var(
  res_pca, col.var = "contrib",
  gradient.cols = c("#00d4ff", "#7c3aed", "#ef4444"),
  repel = TRUE,
  title = sprintf("PCA Biplot - Estrategia %s", PCA_STRATEGY)
) + theme_minimal(base_size = 11)

ggsave(file.path(DIR_SALIDA, "4_5_PCA_biplot_loadings.png"),
       p13_loadings, width = 11, height = 9, dpi = 150)


# ── 4.6  Scatter score × p_good por Intercept ────────────────────────────────
cat("  [4.6] Scatter score x p_good...\n")

# nrow de facets dinámico
n_facet_rows <- ifelse(n_niveles <= 5, 1, ceiling(n_niveles / 5))

set.seed(SEMILLA)
p14_scatter <- df |>
  select(score, p_good, log1p_coverage, log1p_tpm, Intercept) |> drop_na() |>
  slice_sample(n = min(N_SAMPLE_PLOT, nrow(df))) |>
  ggplot(aes(x = score, y = p_good, color = Intercept)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.5) +
  scale_color_manual(values = PALETA_CLASES) +
  facet_wrap(~ Intercept, nrow = n_facet_rows) +
  labs(
    title    = sprintf("score x p_good por %s", INTERCEPT_VAR),
    subtitle = "score = calidad ensamblaje · p_good = bases buenas · curva = LOESS",
    x = "score", y = "p_good"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_SALIDA, "4_6_scatter_score_pgood.png"),
       p14_scatter, width = 16, height = 5 * n_facet_rows, dpi = 150)


# ── 4.7  GGally pairs (top 5 features) ───────────────────────────────────────
cat("  [4.7] GGally pairs...\n")

if (n_niveles > MAX_NIVELES_GGPAIRS) {
  cat(sprintf("    OMITIDO - Intercept tiene %d niveles (>%d), ggpairs sería ilegible\n",
              n_niveles, MAX_NIVELES_GGPAIRS))
} else {
  n_top_pairs <- min(5, nrow(resultados_anova))
  top5 <- resultados_anova$variable[seq_len(n_top_pairs)]
  
  set.seed(SEMILLA)
  # ggpairs es O(n²) en pares — limitar para rendimiento con dataset grande
  n_sample_pairs <- min(round(N_SAMPLE_PLOT / 3), 3000)
  p15_pairs <- df |>
    select(all_of(top5), Intercept) |> drop_na() |>
    slice_sample(n = min(n_sample_pairs, nrow(df))) |>
    ggpairs(
      columns = seq_len(n_top_pairs),
      mapping = aes(color = Intercept, alpha = 0.4),
      upper   = list(continuous = "cor"),
      lower   = list(continuous = "points"),
      diag    = list(continuous = "densityDiag")
    ) +
    scale_color_manual(values = PALETA_CLASES) +
    scale_fill_manual(values = PALETA_CLASES) +
    theme_minimal(base_size = 9) +
    labs(title = sprintf("Pares de variables más discriminativas (color = %s)", INTERCEPT_VAR))
  
  ggsave(file.path(DIR_SALIDA, "4_7_ggpairs_top5.png"),
         p15_pairs, width = 14, height = 12, dpi = 150)
}


# ── 4.8  Cross-tabulación Intercept × otras categóricas (LOOP) ───────────────
cat("  [4.8] Cross-tabulación con otras categóricas...\n")

# Categóricas secundarias: todas excepto Intercept y las de alta cardinalidad
secundarias_candidatas <- setdiff(
  c("prelim_cat", "sf_evidence", "tab", "summarise", "Region"),
  INTERCEPT_VAR
)
secundarias_disponibles <- intersect(secundarias_candidatas, names(df))

for (sec in secundarias_disponibles) {
  tabla_cruzada <- df |>
    count(.data[[sec]], Intercept) |>
    drop_na() |>
    group_by(.data[[sec]]) |>
    mutate(prop = n / sum(n)) |>
    ungroup()
  
  if (nrow(tabla_cruzada) == 0) next
  
  p_alluv <- ggplot(tabla_cruzada,
                    aes(x = .data[[sec]], y = prop, fill = Intercept)) +
    geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(prop > 0.05, scales::percent(prop, 1), "")),
              position = position_stack(vjust = 0.5), size = 3,
              color = "white", fontface = "bold") +
    scale_fill_manual(values = PALETA_CLASES) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(guide = guide_axis(angle = 25)) +
    labs(
      title    = sprintf("%s x %s - composición por categoría", sec, INTERCEPT_VAR),
      subtitle = sprintf("Distribución de %s dentro de cada nivel de %s",
                         INTERCEPT_VAR, sec),
      x = sec, y = "Proporción", fill = INTERCEPT_VAR
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(DIR_SALIDA, sprintf("4_8_crosstab_%s_vs_Intercept.png", sec)),
         p_alluv, width = 11, height = 7, dpi = 150)
  cat(sprintf("    OK %s x Intercept generado\n", sec))
}


# =============================================================================
# BLOQUE 5 — COMPARACIÓN FINAL DE MODELOS
# =============================================================================

cat("\n== BLOQUE 5: Comparación de modelos ==\n")

# ── 5.1  Matriz de confusión OOB del RF ──────────────────────────────────────
pred_oob <- modelo_rf$predictions
true_lab <- df_rf$Intercept

if (!is.null(pred_oob) && length(pred_oob) == length(true_lab)) {
  # Forzar mismos niveles para que table() devuelva matriz cuadrada
  pred_oob <- factor(pred_oob, levels = levels(true_lab))
  conf_mat <- table(Predicho = pred_oob, Real = true_lab)
  cat("  Matriz de confusión OOB:\n")
  print(conf_mat)
  
  metricas_clase <- map_dfr(levels(true_lab), function(cls) {
    tp   <- conf_mat[cls, cls]
    fp   <- sum(conf_mat[cls, ]) - tp
    fn   <- sum(conf_mat[, cls]) - tp
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
      2 * prec * rec / (prec + rec) else NA_real_
    tibble(clase = cls, n_real = sum(conf_mat[, cls]),
           precision = prec, recall = rec, f1 = f1)
  })
  print(metricas_clase)
  write.csv(metricas_clase,
            file.path(DIR_SALIDA, "5_1_metricas_por_clase.csv"),
            row.names = FALSE)
}


# ── 5.2  Tabla resumen comparativa ───────────────────────────────────────────
resumen_final <- tibble(
  Análisis = c(
    "INTERCEPT_VAR", "INTERCEPT_REF", "PCA strategy", "Filas usadas",
    "Pearson |r| > 0.70", "Kruskal-Wallis sig.",
    "ANOVA eta2 max", "MI max (bits)",
    "RF OOB accuracy", "Top predictor (RF)",
    "Top predictor (MI)", "Top predictor (XGBoost)"
  ),
  Resultado = c(
    INTERCEPT_VAR, INTERCEPT_REF, PCA_STRATEGY, as.character(nrow(df)),
    paste(nrow(cor_df), "pares"),
    paste(sum(resultados_kw$p_ajust < 0.05, na.rm = TRUE), "variables"),
    sprintf("%.3f (%s)", resultados_anova$eta2[1], resultados_anova$variable[1]),
    sprintf("%.4f (%s)", mi_df$MI[1], mi_df$variable[1]),
    sprintf("%.3f", 1 - modelo_rf$prediction.error),
    imp_df$variable[1], mi_df$variable[1], xgb_imp$Feature[1]
  )
)

print(resumen_final)

write.csv(resumen_final, file.path(DIR_SALIDA, "5_2_resumen_resultados.csv"),
          row.names = FALSE)


# =============================================================================
# BLOQUE 6 — EXPORTAR RESULTADOS
# =============================================================================

cat("\n== BLOQUE 6: Exportando ==\n")

# MEJORA v2.1: manifest del run — resumen mínimo para comparar experimentos
manifest <- tibble(
  parameter = c("INTERCEPT_VAR", "INTERCEPT_REF", "PCA_STRATEGY",
                "SEMILLA", "N_ARBOLES_RF", "N_SAMPLE_PLOT",
                "n_filas", "n_niveles", "ES_BINARIO",
                "RF_OOB_accuracy", "top_anova", "top_RF", "top_XGB",
                "fecha", "DIR_SALIDA"),
  value     = c(INTERCEPT_VAR, INTERCEPT_REF, PCA_STRATEGY,
                as.character(SEMILLA), as.character(N_ARBOLES_RF),
                as.character(N_SAMPLE_PLOT),
                as.character(nrow(df)), as.character(n_niveles),
                as.character(ES_BINARIO),
                sprintf("%.4f", 1 - modelo_rf$prediction.error),
                resultados_anova$variable[1],
                imp_df$variable[1],
                xgb_imp$Feature[1],
                as.character(Sys.time()),
                DIR_SALIDA)
)
write.csv(manifest, file.path(DIR_SALIDA, "0_manifest_run.csv"), row.names = FALSE)

write.csv(resultados_anova, file.path(DIR_SALIDA, "6_anova_eta2_ranking.csv"),  row.names = FALSE)
write.csv(mi_df,            file.path(DIR_SALIDA, "6_informacion_mutua.csv"),   row.names = FALSE)
write.csv(imp_df,           file.path(DIR_SALIDA, "6_rf_importancia.csv"),      row.names = FALSE)
write.csv(xgb_imp,          file.path(DIR_SALIDA, "6_xgboost_importancia.csv"), row.names = FALSE)
write.csv(cor_df,           file.path(DIR_SALIDA, "6_pares_alta_correlacion.csv"), row.names = FALSE)

saveRDS(modelo_rf, file.path(DIR_SALIDA, "modelo_rf.rds"))
xgb.save(modelo_xgb, file.path(DIR_SALIDA, "modelo_xgb.model"))

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("ANÁLISIS COMPLETO - INTERCEPT = %s (ref = %s)\n",
            INTERCEPT_VAR, INTERCEPT_REF))
cat(sprintf("   Carpeta de salida : %s/\n", DIR_SALIDA))
cat(sprintf("   PNG generados     : %d\n",
            length(list.files(DIR_SALIDA, pattern = "\\.png$"))))
cat(sprintf("   CSV generados     : %d\n",
            length(list.files(DIR_SALIDA, pattern = "\\.csv$"))))
cat(strrep("=", 70), "\n", sep = "")


# =============================================================================
# NOTA FINAL: Extensiones opcionales (revisiones)
# =============================================================================
#
#  Todas las extensiones usan ahora 'Intercept' (la variable renombrada) y
#  PALETA_CLASES (la paleta dinámica). Adaptan automáticamente el análisis
#  a la categórica elegida en INTERCEPT_VAR.
#
#  > UMAP no lineal (reemplaza PCA para visualización):
#     library(umap)
#     um <- umap(X_pca)
#     umap_df <- data.frame(um$layout, clase = labels_pca)
#     ggplot(umap_df, aes(X1, X2, color = clase)) +
#       geom_point(alpha = 0.3) +
#       scale_color_manual(values = PALETA_CLASES) +
#       labs(title = sprintf("UMAP por %s", INTERCEPT_VAR))
#
#  > GAM para relaciones no lineales suavizadas (modela P(clase = referencia)):
#     library(mgcv)
#     mod_gam <- gam(as.numeric(Intercept == INTERCEPT_REF) ~
#                    s(score) + s(p_good) + s(log1p_coverage) +
#                    s(log1p_ref_coverage) + s(prop_gc),
#                    data = df, family = binomial)
#     plot(mod_gam, pages = 1)
#
#  > Regresión Lasso/Ridge para selección de features lineales:
# Da coeficientes muy diferentes debido a los nfolds que varian entre corrida 
# library(glmnet)
# cv_lasso <- cv.glmnet(X_mat, y_vec,
#                       family = ifelse(ES_BINARIO, "binomial", "multinomial"),
#                       alpha = 1, nfolds = 10)
# plot(cv_lasso)
# 
# coefdf <- coef(cv_lasso, s = "lambda.min")
# 
# names(coefdf) <- levels(df_xgb$Intercept)
# 
# coefdf

#
#  > Causal Discovery (PC algorithm):
#     library(pcalg)
#     suf_stat <- list(C = cor(X_pca), n = nrow(X_pca))
#     pc_fit   <- pc(suf_stat, indepTest = gaussCItest,
#                    alpha = 0.01, p = ncol(X_pca))
#     plot(pc_fit, main = sprintf("Grafo causal - PC algorithm (target: %s)",
#                                 INTERCEPT_VAR))
#
# > SHAP values para XGBoost (interpretabilidad por instancia):
# library(SHAPforxgboost)
# shap_vals <- shap.values(xgb_model = modelo_xgb, X_train = X_mat)
# shap_long <- shap.prep(shap_contrib = shap_vals$shap_score, X_train = X_mat)
# shap.plot.summary(shap_long)
# shap.plot.dependence(shap_long, x = "score", color_feature = "p_good")
#
#  > MIC (Maximal Information Coefficient - relaciones no lineales):
#     library(minerva)
#     mic_input <- df_pca_input |> select(-Intercept) |>
#       slice_sample(n = min(5000, nrow(.)))
#     mic_result <- mine(mic_input)
#     ggcorrplot(mic_result$MIC, method = "square", lab = TRUE, lab_size = 2.5,
#                colors = c("#0f1117", "#7c3aed", "#f59e0b"),
#                title = "MIC - Maximal Information Coefficient (no lineal)")
#
#  > Conditional Random Forest (evita sesgo de Gini en features con muchos niveles):
#     library(party)
#     cf <- cforest(Intercept ~ ., data = df_rf, controls = cforest_unbiased())
#     varimp(cf, conditional = TRUE)
#
#  > Test post-hoc de Dunn (después de Kruskal-Wallis, comparaciones múltiples):
#     library(FSA)
#     dunnTest(score ~ Intercept, data = df, method = "bh")
#
# =============================================================================
