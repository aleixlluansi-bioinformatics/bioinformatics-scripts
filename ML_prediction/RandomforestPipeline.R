X <- smote_data[, core_common_all]
y <- smote_data$navm_confirm

# --- 2. Definir grid de hiperparámetros ---
grid <- expand.grid(
  mtry = c(2, floor(sqrt(ncol(X))), ncol(X))
)

# --- 3. Configurar control de CV ---
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# --- 4. Entrenar con caret + RF ---
set.seed(123)
rf_cv <- train(
  x = X,
  y = y,
  method = "rf",
  metric = "ROC",       # optimizamos AUC
  trControl = control,
  tuneGrid = grid,
  ntree = 500
)

# --- 5. Resultados ---
print(rf_cv)
plot(rf_cv)

# --- 6. Mejor modelo ---
best_params <- rf_cv$bestTune
print(best_params)
