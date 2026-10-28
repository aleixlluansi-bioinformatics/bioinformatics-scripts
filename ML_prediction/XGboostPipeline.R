library(caret)
library(randomForest)


set.seed(12345)
train_index <- createDataPartition(y=RF_input_generos_2$navm_confirm, p=0.75, list=FALSE)

## Subset the data
training_set <- RF_input_generos_2[train_index, ]
testing_set <- RF_input_generos_2[-train_index, ]

table(training_set$navm_confirm)
table(testing_set$navm_confirm)

#########
# SMOTE #
#########

# 🔹 Paso 1: forzar a 0/1 en train y test
training_set$navm_confirm <- ifelse(training_set$navm_confirm == 2, 1, 0)
testing_set$navm_confirm  <- ifelse(testing_set$navm_confirm == 2, 1, 0)

table(training_set$navm_confirm)
table(testing_set$navm_confirm)

# 🔹 Paso 2: aplicar SMOTE en train
feature_cols <- setdiff(names(training_set), "navm_confirm")

smote_res <- SMOTE(
  X = training_set[, feature_cols],
  target = training_set$navm_confirm,
  K = 5
)

smote_data <- as.data.frame(smote_res$data)
colnames(smote_data)[ncol(smote_data)] <- "navm_confirm"

# 🔹 Paso 3: asegurarse que ya está en 0/1
smote_data$navm_confirm <- as.numeric(as.character(smote_data$navm_confirm))
table(smote_data$navm_confirm)



# ---------------------------
# 1. Preparar datos
# ---------------------------
feature_cols <- core_common_all
label_col <- "navm_confirm"

# aseguramos que sea factor con niveles consistentes
smote_data[[label_col]] <- factor(smote_data[[label_col]], levels= c(0,1), labels = c("No","Yes"))
# ---------------------------
# 2. Definir control de CV
# ---------------------------
set.seed(123)
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"   # guardar para calcular métricas después
)

# ---------------------------
# 3. Grid de hiperparámetros
# ---------------------------
grid <- expand.grid(
  nrounds = c(50, 100,150, 200),         # número de boosting rounds
  max_depth = c(3, 5, 7),        # profundidad de árboles
  eta = c(0.01, 0.1, 0.3),       # learning rate
  gamma = c(0, 1),               # regularización mínima
  colsample_bytree = c(0.8, 1),  # columnas usadas por árbol
  min_child_weight = c(1, 5),    # mínimo peso de hijos
  subsample = c(0.8, 1)          # fracción de datos usada por cada árbol
)

# ---------------------------
# 4. Entrenar modelo
# ---------------------------
xgb_tuned <- train(
  x = smote_data[, feature_cols],
  y = smote_data[[label_col]],
  method = "xgbTree",
  metric = "ROC",
  trControl = control,
  tuneGrid = grid
)

# ---------------------------
# 5. Resultados
# ---------------------------
print(xgb_tuned)
print(xgb_tuned$bestTune)



feature_cols <- core_common_all
label_col <- "navm_confirm"



set.seed(123)

# CV 5x10
repeat_cv <- createMultiFolds(smote_data$navm_confirm, k = 5, times = 10)

# Inicializar métricas
auc_train_results  <- c()
auc_val_results    <- c()
acc_val_results    <- c()
balacc_val_results <- c()
sens_val_results   <- c()
spec_val_results   <- c()
feature_importance_list <- list()

for (i in seq_along(repeat_cv)) {
  
  train_fold <- smote_data[repeat_cv[[i]], ]
  val_fold   <- smote_data[-repeat_cv[[i]], ]
  
  # Saltar folds con solo una clase
  if(length(unique(train_fold[[label_col]])) < 2 || length(unique(val_fold[[label_col]])) < 2){
    auc_train_results  <- c(auc_train_results, NA)
    auc_val_results    <- c(auc_val_results, NA)
    acc_val_results    <- c(acc_val_results, NA)
    balacc_val_results <- c(balacc_val_results, NA)
    sens_val_results   <- c(sens_val_results, NA)
    spec_val_results   <- c(spec_val_results, NA)
    feature_importance_list[[i]] <- NA
    next
  }
  
  y_train <- train_fold[[label_col]]
  y_val   <- val_fold[[label_col]]
  
  dtrain <- xgb.DMatrix(data = as.matrix(train_fold[, feature_cols]), label = y_train)
  dval   <- xgb.DMatrix(data = as.matrix(val_fold[, feature_cols]), label = y_val)
  
  model <- xgboost(
    data = dtrain,
    objective = "binary:logistic",
    nrounds = xgb_tuned$bestTune$nrounds,
    max_depth = xgb_tuned$bestTune$max_depth,
    eta = xgb_tuned$bestTune$eta,
    gamma = xgb_tuned$bestTune$gamma,
    colsample_bytree = xgb_tuned$bestTune$colsample_bytree,
    min_child_weight = xgb_tuned$bestTune$min_child_weight,
    subsample = xgb_tuned$bestTune$subsample,
    verbose = 0
  )
  
  # Predicciones
  train_pred_prob <- predict(model, dtrain)
  val_pred_prob   <- predict(model, dval)
  
  # AUC
  roc_train <- pROC::roc(y_train, train_pred_prob, quiet = TRUE)
  roc_val   <- pROC::roc(y_val, val_pred_prob, quiet = TRUE)
  
  auc_train_results <- c(auc_train_results, as.numeric(pROC::auc(roc_train)))
  auc_val_results   <- c(auc_val_results, as.numeric(pROC::auc(roc_val)))
  
  # --- Threshold óptimo ---
  best_thr_df <- pROC::coords(roc_val, "best", ret = c("threshold"), best.method = "youden")
  best_thr <- as.numeric(best_thr_df$threshold)
  
  # --- Clasificación ---
  val_pred_class <- factor(ifelse(val_pred_prob > best_thr, 1, 0), levels = c(0,1))
  val_true       <- factor(val_fold[[label_col]], levels = c(0,1))
  
  # --- Confusion Matrix ---
  if(length(val_pred_class) == length(val_true)){
    cm <- caret::confusionMatrix(val_pred_class, val_true, positive = "1")
    acc_val_results    <- c(acc_val_results, cm$overall[["Accuracy"]])
    balacc_val_results <- c(balacc_val_results, cm$byClass[["Balanced Accuracy"]])
    sens_val_results   <- c(sens_val_results, cm$byClass[["Sensitivity"]])
    spec_val_results   <- c(spec_val_results, cm$byClass[["Specificity"]])
  } else {
    acc_val_results    <- c(acc_val_results, NA)
    balacc_val_results <- c(balacc_val_results, NA)
    sens_val_results   <- c(sens_val_results, NA)
    spec_val_results   <- c(spec_val_results, NA)
  }
  
}

metrics_summary <- data.frame(
  AUC_train   = mean(auc_train_results, na.rm = TRUE),
  AUC_val     = mean(auc_val_results, na.rm = TRUE),
  Accuracy    = mean(acc_val_results, na.rm = TRUE),
  Bal_Acc     = mean(balacc_val_results, na.rm = TRUE),
  Sensitivity = mean(sens_val_results, na.rm = TRUE),
  Specificity = mean(spec_val_results, na.rm = TRUE)
)

print(metrics_summary)







###############################
# Modelo final 100 reps y test #
###############################

library(xgboost)
library(shapviz)
library(dplyr)
library(caret)
library(pROC)

set.seed(123)
n_reps <- 100
seeds <- sample(1:10000, n_reps)

# Inicializar métricas
auc_test_results    <- c()
acc_test_results    <- c()
balacc_test_results <- c()
sens_test_results   <- c()
spec_test_results   <- c()
shap_train_list <- list()
shap_test_list  <- list()

for (s in seeds) {
  set.seed(s)
  
  # Entrenar modelo final sobre training set SMOTE
  xgb_model <- xgboost(
    data = as.matrix(smote_data[, core_common_all]),
    label = smote_data[[label_col]],
    nrounds = xgb_tuned$bestTune$nrounds,
    max_depth = xgb_tuned$bestTune$max_depth,
    eta = xgb_tuned$bestTune$eta,
    gamma = xgb_tuned$bestTune$gamma,
    colsample_bytree = xgb_tuned$bestTune$colsample_bytree,
    min_child_weight = xgb_tuned$bestTune$min_child_weight,
    subsample = xgb_tuned$bestTune$subsample,
    objective = "binary:logistic",
    verbose = 0
  )
  
  # Predicción sobre test set
  test_pred_prob <- predict(xgb_model, as.matrix(testing_set[, core_common_all]))
  test_true      <- testing_set[[label_col]]
  
  if(length(unique(test_true)) == 2){
    roc_test <- pROC::roc(test_true, test_pred_prob, levels = c("0","1"), quiet = TRUE)
    best_thresh <- 0.14
    
    test_true <- factor(test_true, levels = c(0,1))
    test_pred_class <- factor(ifelse(test_pred_prob > best_thresh, 1, 0), levels = c(0,1))
    cm <- caret::confusionMatrix(test_pred_class, test_true, positive = "1")
    
    auc_test_results    <- c(auc_test_results, as.numeric(pROC::auc(roc_test)))
    acc_test_results    <- c(acc_test_results, cm$overall["Accuracy"])
    balacc_test_results <- c(balacc_test_results, cm$byClass["Balanced Accuracy"])
    sens_test_results   <- c(sens_test_results, cm$byClass["Sensitivity"])
    spec_test_results   <- c(spec_test_results, cm$byClass["Specificity"])
  } else {
    auc_test_results    <- c(auc_test_results, NA)
    acc_test_results    <- c(acc_test_results, NA)
    balacc_test_results <- c(balacc_test_results, NA)
    sens_test_results   <- c(sens_test_results, NA)
    spec_test_results   <- c(spec_test_results, NA)
  }
  
  # ====================
  # SHAP sobre training set
  # ====================
  sv_train <- shapviz(xgb_model, X_pred = as.matrix(smote_data[, core_common_all]))
  shap_train_list[[s]] <- sv_train
  
  # ====================
  # SHAP sobre test set
  # ====================
  sv_test <- shapviz(xgb_model, X_pred = as.matrix(testing_set[, core_common_all]))
  shap_test_list[[s]] <- sv_test
}

# ---- Métricas promedio ----
final_metrics <- data.frame(
  AUC               = mean(auc_test_results, na.rm = TRUE),
  Accuracy          = mean(acc_test_results, na.rm = TRUE),
  Balanced_Accuracy = mean(balacc_test_results, na.rm = TRUE),
  Sensitivity       = mean(sens_test_results, na.rm = TRUE),
  Specificity       = mean(spec_test_results, na.rm = TRUE)
)
print(final_metrics)

# ---- Feature importance promedio training ----
feature_importance_train <- map_dfr(shap_train_list, function(sv){
  importance <- sv_importance(sv, kind = "bar", plot = FALSE)
  as.data.frame(importance)
}) %>%
  group_by(Feature) %>%
  summarise(mean_abs_shap = mean(Mean_abs_SHAP)) %>%
  arrange(desc(mean_abs_shap))

print(feature_importance_train)



#################################
# SHAP values sobre train y test #
#################################

# ---- Beeswarm plots (última repetición) ----
sv_train_last <- shap_train_list[[length(shap_train_list)]]
sv_test_last  <- shap_test_list[[length(shap_test_list)]]

# Training set
sv_importance(sv_train_last, kind = "beeswarm", max_display = 15, show_numbers = TRUE)
sv_importance(sv_train_last, show_numbers = TRUE)
sv_dependence(sv_train_last, v = feature_cols, share_y = TRUE)
sv_waterfall(sv_train_last, row_id = 2) +
  ggtitle("Waterfall plot")
sv_force(sv_train_last, row_id = 2) +
  ggtitle("Force plot")

# Test set
sv_importance(sv_test_last, kind = "beeswarm", max_display = 15, show_numbers = TRUE)
sv_importance(sv_test_last, show_numbers = TRUE)
sv_dependence(sv_test_last, v = feature_cols, share_y = TRUE)
sv_waterfall(sv_test_last, row_id = 4) +
  ggtitle("Waterfall plot")
sv_force(sv_test_last, row_id = 4) +
  ggtitle("Force plot")
sv_waterfall(sv_test_last, row_id = 5) +
  ggtitle("Waterfall plot")
sv_force(sv_test_last, row_id = 8) +
  ggtitle("Force plot")




#########################################################################################################
# Entrenar modelo 100 reps y distintas semillas y testar sobre test_set (25%) con remplazo sistematico #
##########################################################################################################


feature_cols<-core_xgb_features
set.seed(123)

# ----------------------------
# Parámetros y listas vacías
# ----------------------------
n_reps_final <- 100
seeds_final  <- sample(1:10000, n_reps_final)

final_auc    <- c()
final_aucpr  <- c()
final_f1     <- c()
final_acc    <- c()
final_balacc <- c()
final_sens   <- c()
final_spec   <- c()

shap_train_list <- list()
shap_test_list  <- list()


# ----------------------------
# Loop principal
# ----------------------------
for (i in seq_along(seeds_final)) {
  
  s <- seeds_final[i]
  set.seed(s)
  
  # Split estratificado
  train_idx <- createDataPartition(RF_input_generos_2[[label_col]], p = 0.75, list = FALSE)
  train_set <- RF_input_generos_2[train_idx, ]
  test_set  <- RF_input_generos_2[-train_idx, ]
  
  # SMOTE si clase minoritaria >= 2
  class_counts <- table(train_set[[label_col]])
  if(min(class_counts) >= 2){
    K_adj <- min(5, min(class_counts)-1)
    smote_res <- SMOTE(train_set[, feature_cols], train_set[[label_col]], K = K_adj)
    smote_data <- as.data.frame(smote_res$data)
    colnames(smote_data)[ncol(smote_data)] <- label_col
    smote_data[[label_col]] <- as.numeric(as.character(smote_data[[label_col]]))
  } else {
    smote_data <- train_set
  }
  
  # Entrenar XGBoost
  xgb_model <- xgboost(
    data = as.matrix(smote_data[, feature_cols]),
    label = smote_data[[label_col]],
    nrounds = 50,
    max_depth = 4,
    eta = 0.2,
    objective = "binary:logistic",
    verbose = 0
  )
  
  # Predicciones
  test_prob <- predict(xgb_model, as.matrix(test_set[, feature_cols]))
  test_true <- factor(test_set[[label_col]], levels=c(0,1))
  
  test_class <- factor(ifelse(test_prob > global_f1_threshold, 1, 0), levels=c(0,1))
  cm <- caret::confusionMatrix(test_class, test_true, positive="1")
  
  f1_val <- ifelse((cm$byClass["Sensitivity"] + cm$byClass["Pos Pred Value"]) == 0, 0,
                   2*(cm$byClass["Sensitivity"]*cm$byClass["Pos Pred Value"]) /
                     (cm$byClass["Sensitivity"]+cm$byClass["Pos Pred Value"]))
  
  pr_obj <- PRROC::pr.curve(scores.class0 = test_prob[test_true=="0"],
                            scores.class1 = test_prob[test_true=="1"], curve = FALSE)
  
  # Guardar métricas
  final_auc    <- c(final_auc, as.numeric(pROC::auc(pROC::roc(test_true, test_prob, quiet=TRUE))))
  final_aucpr  <- c(final_aucpr, pr_obj$auc.integral)
  final_f1     <- c(final_f1, f1_val)
  final_acc    <- c(final_acc, cm$overall["Accuracy"])
  final_balacc <- c(final_balacc, cm$byClass["Balanced Accuracy"])
  final_sens   <- c(final_sens, cm$byClass["Sensitivity"])
  final_spec   <- c(final_spec, cm$byClass["Specificity"])
  
  # Valores SHAP
  sv_train <- shapviz(xgb_model, X_pred = as.matrix(smote_data[, feature_cols]))
  sv_test  <- shapviz(xgb_model, X_pred = as.matrix(test_set[, feature_cols]))
  
  shap_train_list[[i]] <- sv_train
  shap_test_list[[i]]  <- sv_test
}

# ----------------------------
# Resumen métricas
# ----------------------------
metrics_test_summary <- data.frame(
  Metric = c("AUC","AUC_PR","F1","Accuracy","Balanced_Accuracy","Sensitivity","Specificity"),
  Mean   = c(mean(final_auc), mean(final_aucpr), mean(final_f1, na.rm=T),
             mean(final_acc), mean(final_balacc), mean(final_sens), mean(final_spec)),
  SD     = c(sd(final_auc), sd(final_aucpr), sd(final_f1, na.rm=T),
             sd(final_acc), sd(final_balacc), sd(final_sens), sd(final_spec)),
  Min    = c(min(final_auc), min(final_aucpr), min(final_f1, na.rm=T),
             min(final_acc), min(final_balacc), min(final_sens), min(final_spec)),
  Max    = c(max(final_auc), max(final_aucpr), max(final_f1, na.rm=T),
             max(final_acc), max(final_balacc), max(final_sens), max(final_spec))
)

print(metrics_test_summary)

# ----------------------------
# Importancia media de SHAP
# ----------------------------
feature_importance_train <- map_dfr(shap_train_list, function(sv){
  importance <- sv_importance(sv, kind = "bar")$data
  as.data.frame(importance)
}) %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(value)) %>%
  arrange(desc(mean_abs_shap))

print(feature_importance_train)
