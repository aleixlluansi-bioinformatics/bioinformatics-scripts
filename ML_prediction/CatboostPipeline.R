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


library(catboost)

# --- 1. Preparar datos ---

# Definir feature matrix y label
X <- smote_data[, core_common_all]
y <- smote_data$navm_confirm

# Crear Pool de CatBoost
train_pool <- catboost.load_pool(data = X, label = y)

# --- 2. Definir grid de hiperparámetros ---
grid <- expand.grid(
  depth = c(4, 5,6,7, 8,9,10),
  learning_rate = c(0.01, 0.03,0.05,0.07, 0.1,0.2),
  l2_leaf_reg = c(1,2 ,3, 5,7,9)
)

set.seed(123)

# --- 3. Hacer 5-fold CV con CatBoost ---
cv_results <- list()

for (i in 1:nrow(grid)) {
  params <- list(
    loss_function = 'Logloss',
    eval_metric = 'AUC',
    iterations = 200,
    depth = grid$depth[i],
    learning_rate = grid$learning_rate[i],
    l2_leaf_reg = grid$l2_leaf_reg[i],
    random_seed = 123,
    logging_level = 'Silent'
  )
  
  cv <- catboost.cv(
    pool = train_pool,
    params = params,
    fold_count = 5,
    type = 'Classical',
    partition_random_seed = 123,
    early_stopping_rounds = 20
  )
  
  best_iter <- which.max(cv$test.AUC.mean)
  best_auc <- max(cv$test.AUC.mean)
  
  cv_results[[i]] <- cbind(grid[i, ], best_iter, best_auc)
}

# --- 4. Combinar resultados ---
cv_results_df <- do.call(rbind, cv_results)
cv_results_df <- as.data.frame(cv_results_df)
cv_results_df <- cv_results_df[order(-cv_results_df$best_auc), ]

# --- 5. Mejor modelo ---
best_params <- cv_results_df[1, ]
print(best_params)






######## CV 5 x 10 ###########

set.seed(123)

repeat_cv <- createMultiFolds(smote_data$navm_confirm, k = 5, times = 10)

# --- 1. Preparar datos y variables ---
feature_cols <- core_common_all # todas las columnas excepto label
label_col    <- "navm_confirm"

# --- 2. Inicializar vectores de métricas ---
auc_train_results  <- c()
auc_val_results    <- c()
acc_val_results    <- c()
balacc_val_results <- c()
sens_val_results   <- c()
spec_val_results   <- c()

# --- 3. Lista para importancias ---
feature_importance_list <- list()

# --- 4. Loop de validación cruzada ---
for (i in seq_along(repeat_cv)) {
  
  # Dividir folds
  train_fold <- smote_data[repeat_cv[[i]], ]
  val_fold   <- smote_data[-repeat_cv[[i]], ]
  
  # --- Pools CatBoost ---
  train_pool <- catboost.load_pool(data = train_fold[, feature_cols], label = train_fold[[label_col]])
  val_pool   <- catboost.load_pool(data = val_fold[, feature_cols], label = val_fold[[label_col]])
  
  # --- Entrenar modelo ---
  model <- catboost.train(
    learn_pool = train_pool,
    test_pool  = val_pool,
    params = list(
      loss_function = 'Logloss',
      eval_metric   = 'AUC',
      iterations    = 39,
      depth         = 4,
      learning_rate = 0.2,
      l2_leaf_reg   = 5,
      random_seed   = 123,
      verbose       = 0
    )
  )
  
  # --- Predicciones ---
  train_pred_prob <- catboost.predict(model, train_pool, prediction_type = "Probability")
  val_pred_prob   <- catboost.predict(model, val_pool, prediction_type = "Probability")
  
  # --- AUC ---
  if(length(unique(train_fold[[label_col]])) == 2){
    roc_train <- pROC::roc(response = train_fold[[label_col]], predictor = train_pred_prob, quiet = TRUE)
    auc_train_results <- c(auc_train_results, as.numeric(pROC::auc(roc_train)))
  }
  
  if(length(unique(val_fold[[label_col]])) == 2){
    roc_val <- pROC::roc(response = val_fold[[label_col]], predictor = val_pred_prob, quiet = TRUE)
    auc_val_results <- c(auc_val_results, as.numeric(pROC::auc(roc_val)))
    
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
  
  # --- Importancias de variables ---
  fi <- catboost.get_feature_importance(model, pool = val_pool)
  feature_importance_list[[i]] <- fi
}

# --- 5. Resumen de métricas ---
metrics_summary <- data.frame(
  AUC_train   = mean(auc_train_results, na.rm = TRUE),
  AUC_val     = mean(auc_val_results, na.rm = TRUE),
  Accuracy    = mean(acc_val_results, na.rm = TRUE),
  Bal_Acc     = mean(balacc_val_results, na.rm = TRUE),
  Sensitivity = mean(sens_val_results, na.rm = TRUE),
  Specificity = mean(spec_val_results, na.rm = TRUE)
)

print(metrics_summary)

# --- 6. Calcular importancias promedio ---
fi_df <- do.call(cbind, feature_importance_list)
colnames(fi_df) <- paste0("Fold", 1:ncol(fi_df))
rownames(fi_df) <- feature_cols

mean_importance <- rowMeans(fi_df, na.rm = TRUE)

importance_df <- data.frame(
  Feature = feature_cols,
  Importance = mean_importance
)

importance_df <- importance_df[order(-importance_df$Importance), ]

# --- 7. Gráfico de importancias ---
ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Importancia promedio de predictores (CV 5x10 CatBoost)",
       x = "Predictor", y = "Importancia media")


#Este código realiza un CV 5×10 sobre un dataset balanceado con SMOTE, 
#entrenando un modelo CatBoost con hiperparámetros fijos. 
#Por cada fold se calculan predicciones probabilísticas, 
#se obtiene el threshold óptimo y se generan métricas de rendimiento como AUC, 
#Accuracy, Balanced Accuracy, Sensitivity y Specificity.
#Finalmente, se resumen las métricas promedio de todos los folds para evaluar
#el desempeño global del modelo.







###############################
# Modelo final 100 reps y test #
###############################

# --- 1. Preparar seeds y métricas ---
set.seed(123)
n_reps <- 100
seeds <- sample(1:10000, n_reps)

auc_test_results    <- c()
acc_test_results    <- c()
balacc_test_results <- c()
sens_test_results   <- c()
spec_test_results   <- c()

# --- 2. Loop de 100 modelos ---
for (s in seeds) {
  
  # --- Pools CatBoost ---
  train_pool <- catboost.load_pool(
    data  = smote_data[, core_common_all],
    label = as.numeric(as.character(smote_data$navm_confirm))
  )
  
  test_pool <- catboost.load_pool(
    data  = testing_set[, core_common_all],
    label = as.numeric(as.character(testing_set$navm_confirm))
  )
  
  # --- Entrenar modelo final ---
  final_model <- catboost.train(
    learn_pool = train_pool,
    test_pool  = test_pool,
    params = list(
      loss_function = 'Logloss',
      eval_metric   = 'Accuracy',
      iterations    = 39,
      depth         = 4,
      learning_rate = 0.2,
      l2_leaf_reg   = 5,
      random_seed   = s,
      verbose       = 0
    )
  )
  
  # --- Predicciones en test set ---
  test_pred_prob  <- catboost.predict(final_model, test_pool, prediction_type = "Probability")
  test_pred_class <- factor(ifelse(test_pred_prob > 0.48, 1, 0), levels = c(0,1))
  test_true       <- factor(testing_set$navm_confirm, levels = c(0,1))
  
  # --- Confusion matrix ---
  cm <- caret::confusionMatrix(test_pred_class, test_true, positive = "1")
  
  # --- Guardar métricas ---
  roc_test <- pROC::roc(test_true, test_pred_prob, quiet = TRUE)
  auc_test_results    <- c(auc_test_results, as.numeric(pROC::auc(roc_test)))
  acc_test_results    <- c(acc_test_results, cm$overall["Accuracy"])
  balacc_test_results <- c(balacc_test_results, cm$byClass["Balanced Accuracy"])
  sens_test_results   <- c(sens_test_results, cm$byClass["Sensitivity"])
  spec_test_results   <- c(spec_test_results, cm$byClass["Specificity"])
}

# --- 3. Resumen promedio ---
final_metrics <- data.frame(
  AUC               = mean(auc_test_results, na.rm = TRUE),
  Accuracy          = mean(acc_test_results, na.rm = TRUE),
  Balanced_Accuracy = mean(balacc_test_results, na.rm = TRUE),
  Sensitivity       = mean(sens_test_results, na.rm = TRUE),
  Specificity       = mean(spec_test_results, na.rm = TRUE)
)

print(final_metrics)








#########################################################################################################
# Entrenar modelo 100 reps y distintas semillas y testar sobre test_set (25%) con remplazo sistematico #
##########################################################################################################


set.seed(123)

# --- Core features originales ---
core_features <- core_common_all   # p.ej. c("GASV1","GASV2",...)

# --- Todas las features disponibles ---
all_features  <- setdiff(names(RF_input_generos_2), "navm_confirm")

# --- Features candidatas (las que no están en el core) ---
candidate_features <- setdiff(all_features, core_features)

training_set$navm_confirm <- as.numeric(as.factor(training_set$navm_confirm)) - 1

# BALANCEO DE TRAINING SET CON SMOTE

feature_cols <- core_common_all
training_set_reduced <- training_set[, c(feature_cols, "navm_confirm")]

smote_res <- SMOTE(
  X = training_set[, feature_cols],
  target = training_set$navm_confirm,
  K = 5
)
smote_data <- as.data.frame(smote_res$data)
colnames(smote_data)[ncol(smote_data)] <- "navm_confirm"
smote_data$navm_confirm <- factor(smote_data$navm_confirm, levels = c(0,1))



# ---- Evaluar reemplazo de una feature con N repeticiones ----
evaluate_replacement <- function(original_core, replaced_idx, new_feature, train_data, test_data, n_reps=100) {
  
  metrics_list <- list()
  
  for (s in 1:n_reps) {
    # Crear nueva combinación de predictores
    new_core <- original_core
    new_core[replaced_idx] <- new_feature
    
    train_pool <- catboost.load_pool(
      data  = train_data[, new_core],
      label = as.numeric(as.character(train_data$navm_confirm))
    )
    test_pool <- catboost.load_pool(
      data  = test_data[, new_core],
      label = as.numeric(as.character(test_data$navm_confirm))
    )
    
    # Entrenar modelo
    model <- catboost.train(
      learn_pool = train_pool,
      test_pool  = test_pool,
      params = list(
        loss_function = 'Logloss',
        eval_metric   = 'Accuracy',
        iterations    = 39,
        depth         = 4,
        learning_rate = 0.2,
        l2_leaf_reg   = 5,
        random_seed   = s,
        verbose       = 0
      )
    )
    
    # Predicciones
    pred_prob  <- catboost.predict(model, test_pool, prediction_type = "Probability")
    true_class <- as.numeric(as.character(test_data$navm_confirm))
    pred_class <- ifelse(pred_prob > 0.48, 1, 0)
    
    # AUC
    auc_val <- as.numeric(pROC::auc(pROC::roc(true_class, pred_prob)))
    
    # Accuracy
    acc_val <- mean(pred_class == true_class)
    
    # Sensitivity y Specificity manualmente
    tp <- sum(pred_class == 1 & true_class == 1)
    tn <- sum(pred_class == 0 & true_class == 0)
    fp <- sum(pred_class == 1 & true_class == 0)
    fn <- sum(pred_class == 0 & true_class == 1)
    
    sens_val <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
    spec_val <- ifelse((tn + fp) == 0, NA, tn / (tn + fp))
    balacc_val <- mean(c(sens_val, spec_val), na.rm = TRUE)
    
    metrics_list[[s]] <- c(Accuracy = acc_val,
                           Balanced_Accuracy = balacc_val,
                           Sensitivity = sens_val,
                           Specificity = spec_val,
                           AUC = auc_val)
  }
  
  metrics_df <- as.data.frame(do.call(rbind, metrics_list))
  metrics_mean <- colMeans(metrics_df, na.rm = TRUE)
  
  return(data.frame(
    replaced_idx = replaced_idx,
    replaced_feature = original_core[replaced_idx],
    new_feature = new_feature,
    Accuracy = metrics_mean["Accuracy"],
    Balanced_Accuracy = metrics_mean["Balanced_Accuracy"],
    Sensitivity = metrics_mean["Sensitivity"],
    Specificity = metrics_mean["Specificity"],
    AUC = metrics_mean["AUC"]
  ))
}



# --- 2. Evaluar todas las combinaciones de reemplazo ---
library(doParallel)

# Detectar cores y crear cluster
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# ---------------------------
# Loop paralelo: probar todos los reemplazos
# ---------------------------
results_replacement <- foreach(i = seq_along(core_features), .combine = rbind, 
                               .packages = c("catboost","pROC","caret")) %:%
  foreach(new_feat = candidate_features, .combine = rbind) %dopar% {
    evaluate_replacement(core_features, i, new_feat, training_set, testing_set)
  }

# Detener cluster
stopCluster(cl)


# --- 3. Revisar resultados ---
head(results_replacement)
summary(results_replacement)

stopImplicitCluster()
registerDoSEQ()
view(results_replacement)
