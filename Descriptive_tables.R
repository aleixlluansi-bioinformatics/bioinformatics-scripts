library(dplyr)
library(tibble)
library(survival)
library(gt)

# ============================================================
# 1. Preparar datos de supervivencia
# ============================================================

# CAMBIAR a 6 si PFS_time está expresado en meses.
tiempo_6m <- 182

datos_surv <- df_patients %>%
  mutate(
    Evento_PFS = as.numeric(as.character(Evento_PFS)),
    PFS_time = as.numeric(as.character(PFS_time))
  )

# Comprobar codificación
eventos_validos <- na.omit(unique(datos_surv$Evento_PFS))

if (!all(eventos_validos %in% c(0, 1))) {
  stop("Evento_PFS debe estar codificado como 0 = censurado y 1 = evento.")
}

# Eliminar registros sin tiempo, evento o grupo
datos_surv <- datos_surv %>%
  filter(
    !is.na(PFS_time),
    PFS_time >= 0,
    !is.na(Evento_PFS),
    !is.na(grupo_tratamiento)
  )

if (nrow(datos_surv) == 0) {
  stop("No hay pacientes con datos válidos de supervivencia.")
}

# ============================================================
# 2. Función para resumir supervivencia
# ============================================================

resumen_surv <- function(data, tiempo_6m = 182) {
  
  # Excluir datos no válidos
  data <- data %>%
    filter(
      !is.na(PFS_time),
      PFS_time >= 0,
      !is.na(Evento_PFS),
      Evento_PFS %in% c(0, 1)
    )
  
  n <- nrow(data)
  
  if (n == 0) {
    return(tibble(
      `Median follow-up` = NA_character_,
      `PFS events, n (%)` = NA_character_,
      `6-month PFS, % (95% CI)` = NA_character_
    ))
  }
  
  # Número y porcentaje de eventos
  ev_pfs <- sum(data$Evento_PFS == 1)
  
  pct <- 100 * ev_pfs / n
  
  # ----------------------------------------------------------
  # Mediana de seguimiento mediante reverse Kaplan-Meier
  # ----------------------------------------------------------
  
  fu_fit <- survfit(
    Surv(PFS_time, 1 - Evento_PFS) ~ 1,
    data = data
  )
  
  mediana_fu <- unname(fu_fit$table["median"])
  
  if (length(mediana_fu) == 0 ||
      !is.finite(mediana_fu)) {
    mediana_fu <- NA_real_
  }
  
  # ----------------------------------------------------------
  # PFS a 6 meses
  # ----------------------------------------------------------
  
  pfs_fit <- survfit(
    Surv(PFS_time, Evento_PFS) ~ 1,
    data = data
  )
  
  # Evitar estimar la tasa más allá del seguimiento observado
  if (max(data$PFS_time) < tiempo_6m) {
    
    tasa_txt <- "NE"
    
  } else {
    
    s6 <- summary(
      pfs_fit,
      times = tiempo_6m
    )
    
    if (length(s6$surv) == 0 || is.na(s6$surv[1])) {
      
      tasa_txt <- "NE"
      
    } else {
      
      tasa <- 100 * s6$surv[1]
      lci <- 100 * s6$lower[1]
      uci <- 100 * s6$upper[1]
      
      tasa_txt <- sprintf(
        "%.1f (%.1f–%.1f)",
        tasa, lci, uci
      )
    }
  }
  
  tibble(
    `Median follow-up` = ifelse(
      is.na(mediana_fu),
      "NE",
      sprintf("%.1f", mediana_fu)
    ),
    `PFS events, n (%)` = sprintf(
      "%d (%.1f%%)",
      ev_pfs, pct
    ),
    `6-month PFS, % (95% CI)` = tasa_txt
  )
}

# ============================================================
# 3. Resúmenes por grupo
# ============================================================

tab_overall <- resumen_surv(
  datos_surv,
  tiempo_6m = tiempo_6m
)

tab_IT <- resumen_surv(
  filter(datos_surv, grupo_tratamiento == "IT"),
  tiempo_6m = tiempo_6m
)

tab_ITCT <- resumen_surv(
  filter(datos_surv, grupo_tratamiento == "IT+CT"),
  tiempo_6m = tiempo_6m
)

# ============================================================
# 4. Comparaciones estadísticas
# ============================================================

# Fisher: proporción de eventos de PFS
datos_eventos <- datos_surv %>%
  filter(!is.na(Evento_PFS))

tabla_eventos <- table(
  datos_eventos$grupo_tratamiento,
  datos_eventos$Evento_PFS
)

p_events <- if (
  nrow(tabla_eventos) >= 2 &&
  ncol(tabla_eventos) >= 2
) {
  fisher.test(tabla_eventos)$p.value
} else {
  NA_real_
}

# Log-rank: comparación de curvas de PFS
datos_logrank <- datos_surv %>%
  filter(!is.na(grupo_tratamiento)) %>%
  mutate(grupo_tratamiento = droplevels(grupo_tratamiento))

if (nlevels(datos_logrank$grupo_tratamiento) >= 2) {
  
  fit_logrank <- survdiff(
    Surv(PFS_time, Evento_PFS) ~ grupo_tratamiento,
    data = datos_logrank
  )
  
  p_logrank <- pchisq(
    fit_logrank$chisq,
    df = length(fit_logrank$n) - 1,
    lower.tail = FALSE
  )
  
} else {
  p_logrank <- NA_real_
}

# ============================================================
# 5. Construir tabla
# ============================================================

tabla2 <- tibble(
  Variable = c(
    "Median follow-up",
    "PFS events, n (%)",
    "6-month PFS, % (95% CI)"
  ),
  Overall = unlist(tab_overall, use.names = FALSE),
  IT = unlist(tab_IT, use.names = FALSE),
  `IT+CT` = unlist(tab_ITCT, use.names = FALSE),
  p.value = c(
    NA_real_,
    p_events,
    p_logrank
  )
)

# ============================================================
# 6. Presentar tabla
# ============================================================

tabla2 %>%
  gt() %>%
  tab_header(
    title = "Survival outcomes"
  ) %>%
  cols_label(
    Variable = "Outcome",
    Overall = "Overall",
    IT = "IT",
    `IT+CT` = "IT+CT",
    p.value = "p-value"
  ) %>%
  fmt_number(
    columns = p.value,
    decimals = 3
  ) %>%
  sub_missing(
    columns = everything(),
    missing_text = ""
  )