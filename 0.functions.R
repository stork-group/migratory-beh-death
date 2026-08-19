
# ------------------------------------------------------------------------------
# ---- List of all functions ---------------------------------------------------
# ------------------------------------------------------------------------------

# 1. Run the distance models ---------------------------------------------------
run_distance_models <- function(df, window_label_filter = "Window: 15 Days",
                                run_dharma = TRUE, save_csv = TRUE, 
                                output_dir = "model_outputs") {
  storks <- unique(df$ID)
  results <- list()
  models  <- list()
  dharma_res <- list()
  
  # Create output directory if it doesn't exist
  if (save_csv && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (stork in storks) {
    sub_df <- df |>
      dplyr::filter(ID == stork, window_label == window_label_filter)
    
    d_year <- death_years[stork]
    
    sub_df <- sub_df |>
      mutate(year_fact = factor(year),
             year_fact = relevel(year_fact, ref = as.character(d_year)))
    
    model_key <- paste(stork, window_label_filter)
    
    tryCatch({
      model <- glmmTMB(log(total_dist_km) ~ year_fact,
                       data = sub_df,
                       family = gaussian())
      
      cat("\n========================================================\n")
      cat("MODEL FOR:", stork, "| METRIC: total_dist_km | WINDOW:", window_label_filter, "\n")
      cat("Reference (Intercept): year of death (", d_year, ")\n")
      cat("========================================================\n")
      print(summary(model))
      
      # --- store the model object separately ---
      models[[model_key]] <- model
      
      # --- store only the tidy dataframe in results ---
      res <- tidy(model, effects = "fixed") |>
        mutate(ID = stork,
               metric = "total_dist_km",
               window = window_label_filter,
               death_year_ref = d_year) |>
        mutate(term = str_remove(term, "year_fact"))
      
      results[[model_key]] <- res
      
      # --- save clean, machine-readable CSV ---
      if (save_csv) {
        file_name <- file.path(output_dir, 
                               paste0(stork, "_total_dist_km_", 
                                      str_replace_all(window_label_filter, " |:", "_"), ".csv"))
        write_csv(res, file_name)
      }
      
      if (run_dharma) {
        sim <- simulateResiduals(model, n = 1000)
        dharma_res[[model_key]] <- sim
        plot(sim, title = paste0(stork, " — total_dist_km (", window_label_filter, ")"))
        
        disp_test <- testDispersion(sim, plot = FALSE)
      }
      
    }, error = function(e) {
      message(paste("Error in model:", stork, "for", window_label_filter))
    })
  }
  
  list(
    results = bind_rows(results),   # only ever dataframes go in here
    models  = models,
    dharma  = dharma_res
  )
}

# 2. Get the daily number of foraging, resting and flying locations ------------
prepare_daily_counts <- function(data, window_days) {
  data |>
    dplyr::filter(days_to_mig >= -window_days & days_to_mig <= -1) |>
    group_by(ID, year, date, days_to_mig) |> 
    summarise(
      total = n(),
      foraging_success = sum(predBeh == "Foraging", na.rm = TRUE),
      resting_success  = sum(predBeh == "Resting", na.rm = TRUE),
      flying_success   = sum(predBeh == "Flying", na.rm = TRUE),
      .groups = "drop"
    )
}

# 3. Run the behaviour models --------------------------------------------------
run_beh_models <- function(df, window_label, run_dharma = TRUE, 
                           save_csv = TRUE, output_dir = "model_outputs") {
  storks <- unique(df$ID)
  behaviors <- c("foraging", "resting", "flying")
  results <- list()
  models  <- list()
  dharma_res <- list()
  
  # Create output directory if it doesn't exist
  if (save_csv && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (stork in storks) {
    sub_df <- df |> dplyr::filter(ID == stork)
    d_year <- death_years[stork]
    
    sub_df <- sub_df |>
      mutate(year_fact = factor(year),
             year_fact = relevel(year_fact, ref = as.character(d_year)))
    
    for (beh in behaviors) {
      success_col <- paste0(beh, "_success")
      sub_df$success <- sub_df[[success_col]]
      sub_df$failure <- sub_df$total - sub_df$success
      
      model_key <- paste(stork, beh, window_label)
      
      tryCatch({
        model <- glmmTMB(cbind(success, failure) ~ year_fact,
                         data = sub_df,
                         family = betabinomial())
        
        cat("\n========================================================\n")
        cat("Name:", stork, "| Behaviour:", beh, "| Window:", window_label, "\n")
        cat("Intercept: Year of death (", d_year, ")\n")
        cat("========================================================\n")
        print(summary(model))
        
        # --- store the model object separately ---
        models[[model_key]] <- model
        
        # --- store only the tidy dataframe in results ---
        res <- tidy(model, effects = "fixed") |>
          mutate(ID = stork,
                 behavior = str_to_title(beh),
                 window = window_label,
                 death_year_ref = d_year) |>
          mutate(term = str_remove(term, "year_fact"))
        
        results[[model_key]] <- res
        
        # --- save individual CSV, with a title row identifying the model ---
        if (save_csv) {
          file_name <- file.path(output_dir, 
                                 paste0(stork, "_", beh, "_", 
                                        str_replace_all(window_label, " ", "_"), 
                                        ".csv"))
          write_csv(res, file_name)
        }
        
        if (run_dharma) {
          sim <- simulateResiduals(model, n = 1000)
          dharma_res[[model_key]] <- sim
          plot(sim, title = paste0(stork, " — ", str_to_title(beh), 
                                   " (", window_label, ")"))
        }
        
      }, error = function(e) {
        message(paste("Error in model:", stork, beh, "for the window", 
                      window_label))
      })
    }
  }
  
  list(
    results = bind_rows(results), 
    models  = models,
    dharma  = dharma_res
  )
}

# 4. Create a pseudo-year for plotting -----------------------------------------
make_plot_date <- function(date) {
  date <- as.Date(date)
  yr <- if_else(month(date) >= 5, 2000, 2001)
  make_date(yr, month(date), day(date))
}

# 5. Summarise the daily behaviour of the focal and non-focal inds -------------
summarise_daily_behaviour <- function(data) {
  data |>
    filter(!is.na(predBeh)) |>
    count(individual_local_identifier, focal, age_cat, date, predBeh, mcp_id,
          name = "n") |>
    pivot_wider(
      names_from = predBeh,
      values_from = n,
      values_fill = 0,
      names_prefix = "n_points_"
    ) |>
    rowwise() |>
    mutate(tot_n_points = sum(c_across(starts_with("n_points_")))) |>
    ungroup() |>
    mutate(
      n_points_flying = coalesce(n_points_Soaring, 0) + coalesce(n_points_Flapping, 0),
      pct_foraging = 100 * coalesce(n_points_Foraging, 0) / tot_n_points,
      pct_resting  = 100 * coalesce(n_points_Resting, 0)  / tot_n_points,
      pct_flying   = 100 * n_points_flying / tot_n_points
    )
}


# 6. Focal vs non-focal behaviour models ---------------------------------------
run_focal_behaviour_models <- function(data_list, min_points = 10, run_dharma = TRUE,
                                       skip_models = character(0)) {
  behaviour_cols <- c(
    foraging = "n_points_Foraging",
    resting  = "n_points_Resting",
    flying   = "n_points_flying"
  )
  
  results <- list()
  models  <- list()
  dharma_res <- list()
  
  for (id in names(data_list)) {
    
    sub_df_full <- data_list[[id]] |> filter(tot_n_points >= min_points)
    
    for (beh in names(behaviour_cols)) {
      model_key <- paste(id, beh, sep = " | ")
      
      # --- skip this specific ID/behaviour combination if requested ---
      if (model_key %in% skip_models) {
        message(paste("Skipping model:", model_key))
        next
      }
      
      beh_col <- behaviour_cols[[beh]]
      
      sub_df <- sub_df_full
      sub_df$success <- sub_df[[beh_col]]
      sub_df$failure <- sub_df$tot_n_points - sub_df$success
      
      tryCatch({
        model <- glmmTMB(cbind(success, failure) ~ focal + age_cat + 
                           (1 | individual_local_identifier),
                         data = sub_df,
                         family = betabinomial())
        
        cat("\n========================================================\n")
        cat("MODEL FOR:", id, "| BEHAVIOUR:", beh, "\n")
        cat("========================================================\n")
        print(summary(model))
        
        models[[model_key]] <- model
        
        res <- tidy(model, effects = "fixed", conf.int = TRUE) |>
          filter(term != "(Intercept)") |>
          mutate(ID = id,
                 behavior = str_to_title(beh))
        
        results[[model_key]] <- res
        
        if (run_dharma) {
          sim <- simulateResiduals(model, n = 1000)
          dharma_res[[model_key]] <- sim
          plot(sim, title = paste0(id, " — ", str_to_title(beh)))
          
          disp_test <- testDispersion(sim, plot = FALSE)
          cat("Dispersion test p-value:", round(disp_test$p.value, 3), "\n")
        }
        
      }, error = function(e) {
        message(paste("Error in model:", id, beh, "-", e$message))
      })
    }
  }
  
  list(
    results = bind_rows(results),
    models  = models,
    dharma  = dharma_res
  )
}


# 7. Extract predictions from each saved model ---------------------------------
extract_focal_predictions <- function(models_list) {
  map_dfr(names(models_list), function(key) {
    model <- models_list[[key]]
    
    # key format: "ID | behaviour"
    parts <- str_split(key, " \\| ", simplify = TRUE)
    id <- parts[1, 1]
    beh <- parts[1, 2]
    
    preds <- ggpredict(model, terms = "focal")
    
    as.data.frame(preds) |>
      transmute(
        ID = id,
        behavior = str_to_title(beh),
        focal = x,
        predicted_pct = predicted * 100,
        conf.low_pct  = conf.low  * 100,
        conf.high_pct = conf.high * 100
      )
  })
}

