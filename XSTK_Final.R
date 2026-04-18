# ==============================================================================
# STEP 1: LOAD PACKAGES & SETUP ENVIRONMENT
# ==============================================================================
options(repos = c(CRAN = "https://cran.r-project.org"))

required_pkgs <- c("tidyverse", "ggcorrplot", "patchwork", "car", 
                   "lmtest", "broom", "sandwich", "caret", "glmnet", "Matrix")

missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(missing_pkgs)) install.packages(missing_pkgs, dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggcorrplot)
  library(patchwork)
  library(car)
  library(lmtest)
  library(broom)
  library(sandwich)
  library(caret)
  library(glmnet)
  library(Matrix)
})

# ==============================================================================
# STEP 2: IMPORT & PREPROCESS DATA
# ==============================================================================
raw_data <- read.csv("data.csv")
raw_data <- raw_data %>% drop_na()

if("tension_strenght" %in% names(raw_data)) {
  raw_data <- raw_data %>% rename(tension_strength = tension_strenght)
}

clean_data <- raw_data %>%
  mutate(
    infill_pattern_num = if_else(infill_pattern == "grid", 0, 1),
    material_num = if_else(material == "abs", 0, 1)
  ) %>% 
  select(-infill_pattern, -material) %>% 
  rename(
    infill_pattern = infill_pattern_num,
    material = material_num
  )

long_format_data <- clean_data %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  mutate(parameter = factor(parameter, levels = unique(parameter)))

data_model <- raw_data %>%
  mutate(
    infill_pattern = factor(infill_pattern) %>% relevel(ref = "grid"),
    material = factor(material) %>% relevel(ref = "abs")
  ) %>%
  select(-fan_speed) 

# ==============================================================================
# STEP 3: SUMMARY STATISTICS
# ==============================================================================
continous_function <- function(x) {
  c(
    n = length(x),
    Mean = mean(x, na.rm = TRUE),
    Sd = sd(x, na.rm = TRUE),
    `Q1.25%` = unname(quantile(x, 0.25, na.rm = TRUE)),
    `Q2.50%` = unname(quantile(x, 0.50, na.rm = TRUE)),
    Q2 = median(x, na.rm = TRUE),
    `Q3.75%` = unname(quantile(x, 0.75, na.rm = TRUE)),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

continous_data <- raw_data %>% select(-infill_pattern, -material)
continuous_stats <- apply(continous_data, 2, continous_function)

cat("\n--- CONTINUOUS VARIABLES STATISTICS ---\n")
print(continuous_stats)

cat("\n--- CATEGORICAL VARIABLES STATISTICS ---\n")
print(table(raw_data$infill_pattern))
print(table(raw_data$material))

# ==============================================================================
# STEP 4: DATA VISUALIZATION
# ==============================================================================
plot_hist <- ggplot(long_format_data, aes(x = value)) + 
  geom_histogram(fill = "green", color = "white", bins = 15, alpha = 0.8) +
  facet_wrap(~ parameter, scales = "free", nrow = 4) +
  theme_minimal()
print(plot_hist)

continuous_params_data <- long_format_data %>%
  filter(!parameter %in% c("infill_pattern", "material"))

plot_box <- ggplot(continuous_params_data, aes(y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(fill = "purple", outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(x = -0.75), width = 0.1, height = 0, color = "red", alpha = 0.4) +
  facet_wrap(~ parameter, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot_box)

indep_vars <- c("layer_height", "wall_thickness", "infill_density", 
                "nozzle_temperature", "bed_temperature", "print_speed", "fan_speed")

create_scatter_grid <- function(data, y_var) {
  plots <- lapply(indep_vars, function(x_var) {
    ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(shape = 1, color = "black", alpha = 0.7) + 
      labs(title = paste(x_var, "&", y_var), x = x_var, y = y_var) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  })
  wrap_plots(plots, ncol = 4)
}

print(create_scatter_grid(clean_data, "roughness"))
print(create_scatter_grid(clean_data, "tension_strength"))
print(create_scatter_grid(clean_data, "elongation"))

# ------------------------------------------------------------------------------
# 4.5. Correlation Matrix (Heatmap)
# ------------------------------------------------------------------------------
cat("\nGenerating Pearson Correlation Matrix...\n")
correlation_matrix <- cor(clean_data) %>% round(2)

plot_corr <- ggcorrplot(correlation_matrix,
                        method = "square", 
                        type = "lower", 
                        lab = TRUE,     
                        lab_size = 3.5,
                        colors = c("#F8766D", "white", "#00BFC4"), 
                        title = "Pearson Correlation Matrix",
                        legend.title = "Correlation\nCoefficient",
                        outline.color = "gray80") + 
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(b = 20)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "right" 
  )

print(plot_corr)

# ==============================================================================
# STEP 5: HELPER FUNCTIONS FOR MODELING
# ==============================================================================
check_model <- function(model, model_name){
  cat("\n====================================================\n")
  cat("MODEL:", model_name, "\n")
  cat("====================================================\n")
  print(summary(model))
  print(vif(model))
}

robust_check <- function(model, model_name = "MODEL"){
  cat("\n--- ROBUST HC3 STANDARD ERRORS FOR:", model_name, "---\n")
  print(coeftest(model, vcov = vcovHC(model, type = "HC3")))
}

nmae <- function(actual, pred){
  mean(abs(actual - pred)) / (max(actual) - min(actual))
}

safe_summary <- function(data, lev = NULL, model = NULL){
  rmse_val <- sqrt(mean((data$obs - data$pred)^2))
  mae_val  <- mean(abs(data$obs - data$pred))
  rsq_val <- if(sd(data$obs) == 0 || sd(data$pred) == 0) 0 else cor(data$obs, data$pred)^2
  if(is.na(rsq_val) || is.nan(rsq_val)) rsq_val <- 0
  c(RMSE = rmse_val, Rsquared = rsq_val, MAE = mae_val)
}

extract_en_coef <- function(en_model){
  coef_mat <- coef(en_model$finalModel, en_model$bestTune$lambda)
  data.frame(term = rownames(as.matrix(coef_mat)), estimate = as.numeric(coef_mat), row.names = NULL)
}

# ==============================================================================
# STEP 6: OLS MODELING
# ==============================================================================
rough_full <- lm(roughness ~ ., data = data_model %>% select(-tension_strength, -elongation))
rough_m2 <- lm(roughness ~ layer_height + nozzle_temperature + bed_temperature + print_speed + material, data = data_model)

tension_full <- lm(tension_strength ~ ., data = data_model %>% select(-roughness, -elongation))
tension_m2 <- lm(tension_strength ~ layer_height + wall_thickness + infill_density + nozzle_temperature + bed_temperature + material, data = data_model)

elong_full <- lm(elongation ~ ., data = data_model %>% select(-roughness, -tension_strength))
elong_m2 <- lm(elongation ~ layer_height + infill_density + nozzle_temperature + bed_temperature + material, data = data_model)
elong_m3 <- lm(elongation ~ layer_height + wall_thickness + infill_density + nozzle_temperature + bed_temperature + print_speed + material, data = data_model)

check_model(rough_m2, "ROUGHNESS - MODEL 2")
robust_check(rough_m2, "ROUGHNESS - MODEL 2")

check_model(tension_m2, "TENSION STRENGTH - MODEL 2")
check_model(elong_m3, "ELONGATION - MODEL 3")

par(mfrow = c(2,2))
plot(rough_m2)
plot(tension_m2)
plot(elong_m3)
par(mfrow = c(1,1))

# ==============================================================================
# STEP 7: ELASTIC NET MODELING
# ==============================================================================
set.seed(75)
train_rows <- sample(seq_len(nrow(data_model)), size = floor(0.8 * nrow(data_model)))
test_rows  <- setdiff(seq_len(nrow(data_model)), train_rows)

train_data <- data_model[train_rows, ]
test_data  <- data_model[test_rows, ]

control <- trainControl(method = "repeatedcv", number = 5, repeats = 5, summaryFunction = safe_summary)

run_en <- function(target_name){
  outputs <- c("roughness", "tension_strength", "elongation")
  outputs_to_drop <- setdiff(outputs, target_name)
  
  train_x <- train_data[, !colnames(train_data) %in% outputs_to_drop]
  test_x  <- test_data[, !colnames(test_data) %in% outputs_to_drop]
  
  set.seed(75)
  en_model <- train(
    as.formula(paste(target_name, "~ .")),
    data = train_x, method = "glmnet", trControl = control, metric = "RMSE",
    tuneGrid = expand.grid(alpha = 0.5, lambda = seq(0.001, 1, length = 100))
  )
  
  pred_train <- predict(en_model, newdata = train_x)
  pred_test  <- predict(en_model, newdata = test_x)
  
  actual_train <- train_x[[target_name]]
  actual_test  <- test_x[[target_name]]
  
  rsq_test <- 1 - sum((actual_test - pred_test)^2) / sum((actual_test - mean(actual_test))^2)
  
  return(list(
    model = en_model,
    nmae_train = nmae(actual_train, pred_train),
    nmae_test = nmae(actual_test, pred_test),
    rsq_test = rsq_test,
    coef_table = extract_en_coef(en_model)
  ))
}

en_rough   <- run_en("roughness")
en_tension <- run_en("tension_strength")
en_elong   <- run_en("elongation")

# ==============================================================================
# STEP 8: COMPARISON & EXPORT
# ==============================================================================
compare_table <- data.frame(
  Output = c("Roughness", "Tension Strength", "Elongation"),
  OLS_Selected_Model = c("rough_m2", "tension_m2", "elong_m3"),
  OLS_Adj_R2_InSample = c(summary(rough_m2)$adj.r.squared, summary(tension_m2)$adj.r.squared, summary(elong_m3)$adj.r.squared),
  EN_R2_Test_OutSample = c(en_rough$rsq_test, en_tension$rsq_test, en_elong$rsq_test),
  EN_NMAE_Train = c(en_rough$nmae_train, en_tension$nmae_train, en_elong$nmae_train),
  EN_NMAE_Test = c(en_rough$nmae_test, en_tension$nmae_test, en_elong$nmae_test)
)

cat("\n===== OLS VS ELASTIC NET =====\n")
print(compare_table)

write.csv(compare_table, "ols_vs_elastic_net.csv", row.names = FALSE)
write.csv(en_rough$coef_table, "en_coef_roughness.csv", row.names = FALSE)
write.csv(en_tension$coef_table, "en_coef_tension.csv", row.names = FALSE)
write.csv(en_elong$coef_table, "en_coef_elongation.csv", row.names = FALSE)

cat("\n===== DONE =====\n")