
# ============================================================
# Kaplan-Meier: Progression-Free Survival (PFS)
# ============================================================

library(survival)
library(survminer)
library(ggplot2)

# 1. Cargar datos
NEWDATA <- read.csv(
  "~/Documentos/Neuroblastoma/Segunda secuenciacion/swisstransfer_28ecf44c-3133-4677-b5ae-28de6b147272/Resultados NGS_UCA-GTB_08092025/Results_rmd/NEWDATA.csv",
  stringsAsFactors = FALSE
)

# 2. Convertir fechas
NEWDATA$Fecha_inicio_IT <- as.Date(
  NEWDATA$Fecha_inicio_IT,
  format = "%d/%m/%y"
)

NEWDATA$Fecha_evento <- as.Date(
  NEWDATA$Fecha_evento,
  format = "%d/%m/%y"
)

# 3. Calcular tiempo de seguimiento en días
# IMPORTANTE:
# Fecha_evento debe representar la fecha de progresión,
# recaída o muerte, según la definición de PFS.
NEWDATA$PFS_time <- as.numeric(
  NEWDATA$Fecha_evento - NEWDATA$Fecha_inicio_IT
)

# 4. Revisar los datos antes del análisis
summary(NEWDATA$PFS_time)
table(NEWDATA$Evento_PFS, useNA = "ifany")
table(NEWDATA$NE_bin, useNA = "ifany")

# 5. Filtrar la población de análisis
NEWDATA_eval <- subset(
  NEWDATA,
  NE_bin == 0 &
    !is.na(PFS_time) &
    PFS_time >= 0 &
    !is.na(Evento_PFS) &
    Evento_PFS %in% c(0, 1)
)

# Comprobar que hay pacientes y eventos
if (nrow(NEWDATA_eval) == 0) {
  stop("No hay pacientes válidos para el análisis.")
}

if (sum(NEWDATA_eval$Evento_PFS == 1) == 0) {
  warning("No hay eventos de PFS en la población seleccionada.")
}

# 6. Kaplan-Meier global
fit_pfs <- survfit(
  Surv(PFS_time, Evento_PFS) ~ 1,
  data = NEWDATA_eval
)

x_max <- max(NEWDATA_eval$PFS_time, na.rm = TRUE)

KaplanMeier_PFS <- ggsurvplot(
  fit_pfs,
  data = NEWDATA_eval,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = FALSE,
  xlab = "Time (days)",
  ylab = "Progression-free survival probability",
  surv.median.line = "none",
  ggtheme = theme_classic(base_size = 14),
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  xlim = c(0, x_max),
  ylim = c(0, 1),
  font.x = 14,
  font.y = 14,
  font.tickslab = 14
)

# 7. Kaplan-Meier estratificado por cohorte
NEWDATA_eval$Cohorte <- factor(NEWDATA_eval$Cohorte)

fit_strat <- survfit(
  Surv(PFS_time, Evento_PFS) ~ Cohorte,
  data = NEWDATA_eval
)

KaplanMeier_PFS_stratifiedCohorte <- ggsurvplot(
  fit_strat,
  data = NEWDATA_eval,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time (days)",
  ylab = "Progression-free survival probability",
  ggtheme = theme_classic(base_size = 14),
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  xlim = c(0, x_max),
  ylim = c(0, 1),
  font.x = 14,
  font.y = 14,
  font.tickslab = 14,
  font.legend = 12
)

# 8. Guardar las figuras en PDF
pdf(
  "KaplanMeier_PFS_combined.pdf",
  width = 14,
  height = 8
)

combined <- arrange_ggsurvplots(
  list(
    KaplanMeier_PFS,
    KaplanMeier_PFS_stratifiedCohorte
  ),
  nrow = 1,
  ncol = 2
)

print(combined)
dev.off()

