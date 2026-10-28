#################################################
# Metodo 1: combinacion de RFE con los 3 modelos #
##################################################
##################################################



# --- 1. Preparar datos ---
feature_cols <- setdiff(names(RF_input_generos_2), "navm_confirm")
label_col <- "navm_confirm"

X <- RF_input_generos_2[, feature_cols]
y <- RF_input_generos_2[[label_col]]

# --- 2. Parámetros de stability RFE ---
n_reps <- 20      # Número de repeticiones
top_features <- c() # Guardar features seleccionadas en cada repetición

# --- 3. Loop de RFE con distintos splits ---
for(i in 1:n_reps){
  
 
  
  # RFE control
  ctrl <- rfeControl(functions = rfFuncs,  # puedes cambiar rfFuncs por caretFuncs o treebagFuncs
                     method = "cv",
                     number = 5)
  
  # Ejecutar RFE
  results_rfe <- rfe(x = X, y = y,
                     sizes = seq(5, length(feature_cols), by = 1),
                     rfeControl = ctrl)
  
  # Guardar las features seleccionadas
  top_features <- c(top_features, predictors(results_rfe))
}

# --- 4. Calcular frecuencia de selección ---
feature_freq <- table(top_features) / n_reps
feature_freq_df <- data.frame(
  Feature = names(feature_freq),
  Frequency = as.numeric(feature_freq)
)

# --- 5. Filtrar features frecuentes (ej. ≥ 60% de las repeticiones) ---
core_features <- feature_freq_df$Feature[feature_freq_df$Frequency >= 0.95]
core_features
feature_freq_df[order(-feature_freq_df$Frequency), ]



library(catboost)
library(dplyr)

# --- 1. Preparar datos ---
feature_cols <- setdiff(names(RF_input_generos_2), "navm_confirm")
label_col <- "navm_confirm"
n_reps <- 20  # Número de repeticiones

# Matrices para guardar resultados
feature_selection_matrix <- matrix(0, nrow = n_reps, ncol = length(feature_cols))
feature_importance_matrix <- matrix(0, nrow = n_reps, ncol = length(feature_cols))

colnames(feature_selection_matrix) <- feature_cols
colnames(feature_importance_matrix) <- feature_cols

# --- 2. Repeticiones ---
set.seed(123)
for(i in 1:n_reps){
  
  # CatBoost Pool
  train_pool <- catboost.load_pool(
    data = RF_input_generos_2[, feature_cols],
    label = RF_input_generos_2[[label_col]]
  )
  
  # Entrenar CatBoost
  model <- catboost.train(
    learn_pool = train_pool,
    params = list(
      loss_function = 'Logloss',
      iterations = 2000,
      depth = 6,
      learning_rate = 0.1,
      l2_leaf_reg = 3,
      random_seed = i,
      verbose = 0
    )
  )
  
  # Importancias de features
  imp <- catboost.get_feature_importance(model, pool = train_pool, type = "FeatureImportance")
  
  # Guardar selección (binaria) e importancias reales
  feature_selection_matrix[i, ]  <- ifelse(imp > 0, 1, 0)
  feature_importance_matrix[i, ] <- imp
}

# --- 3. Calcular frecuencia de selección ---
feature_freq <- colMeans(feature_selection_matrix)

# --- 4. Calcular media de importancias ---
feature_mean_imp <- colMeans(feature_importance_matrix)

# --- 5. Crear dataframe combinado ---
feature_stats_df <- data.frame(
  Feature = feature_cols,
  Frequency = feature_freq,
  MeanImportance = feature_mean_imp
) %>% arrange(desc(MeanImportance))

print(feature_stats_df)


# --- 5. Filtrar features frecuentes (ej. ≥ 60% de las repeticiones) ---
core_features_cb <- feature_stats_df$Feature[feature_stats_df$MeanImportance >= 1]
core_features_cb

# --- Core set común entre RF y CatBoost ---
core_common <- intersect(core_features_cb, core_features)

# Mostrar resultados
print(core_common)




library(xgboost)
library(dplyr)

# --- Preparar datos ---
feature_cols <- setdiff(names(RF_input_generos_2), "navm_confirm")
label_col <- "navm_confirm"
X <- as.matrix(RF_input_generos_2[, feature_cols])
y <- as.numeric(RF_input_generos_2[[label_col]]) - 1  # convertir a 0/1

# --- Parámetros ---
n_reps <- 50   # número de repeticiones
top_threshold <- 0.039  # seleccionar features que aparezcan ≥60% de las repeticiones

# --- Matriz para guardar importancias ---
xgb_imp_matrix <- matrix(0, nrow = n_reps, ncol = length(feature_cols))
colnames(xgb_imp_matrix) <- feature_cols

# --- Loop de repeticiones ---
for(seed in 1:n_reps){
  
  set.seed(seed)
  
  # Entrenar modelo XGBoost
  xgb_model <- xgboost(
    data = X,
    label = y,
    nrounds = 52,
    max_depth = 6,
    eta = 0.1,
    objective = "binary:logistic",
    verbose = 0
  )
  
  # Obtener importancias por Gain
  imp <- xgb.importance(feature_names = feature_cols, model = xgb_model)
  
  # Mapear importancias a todas las features
  xgb_imp_matrix[seed, imp$Feature] <- imp$Gain
}

# --- Promediar importancias ---
xgb_mean_imp <- colMeans(xgb_imp_matrix)

# --- Normalizar para comparar ---
xgb_mean_imp <- xgb_mean_imp / sum(xgb_mean_imp)

# --- Ordenar y seleccionar features estables ---
xgb_importance_df <- data.frame(
  Feature = names(xgb_mean_imp),
  Importance = xgb_mean_imp
) %>% arrange(desc(Importance))

# Seleccionar core set según umbral
core_xgb_features <- xgb_importance_df$Feature[xgb_importance_df$Importance >= top_threshold]

# Mostrar resultados
print(core_xgb_features)




core_common_all <- Reduce(intersect, list(core_features, core_features_cb, core_xgb_features))
print(core_common_all)







#################################################
# Metodo 2: LASSO #
##################################################
##################################################

library(glmnet)

# Separar X y y
X <- as.matrix(RF_input_generos_2[, feature_cols])
y <- ifelse(RF_input_generos_2[[label_col]] == 2, 1, 0)

# LASSO con validación cruzada
set.seed(123)
cv_lasso <- cv.glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha = 1,          # LASSO
  nfolds = 10,
  type.measure = "auc"  # optimizar AUC
)

# Lambda óptimo
lambda_opt <- cv_lasso$lambda.min
cat("Lambda óptimo:", lambda_opt, "\n")

# Modelo final con lambda óptimo
lasso_model <- glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha = 1,
  lambda = lambda_opt
)

# Features seleccionadas
selected_features <- rownames(coef(lasso_model))[coef(lasso_model)[,1] != 0]
selected_features <- selected_features[selected_features != "(Intercept)"]
cat("Features seleccionadas:\n")
print(selected_features)

