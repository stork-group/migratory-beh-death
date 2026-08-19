
# ------------------------------------------------------------------------------
# -------- Changes behaviour during winter compared to others ------------------
# ------------------------------------------------------------------------------

# --- 1. Set up the R session --------------------------------------------------
source("Scripts_to_publish/0.setup.R")
source("Scripts_to_publish/0.functions.R")


# --- 2. Data upload -----------------------------------------------------------
all_data <- fread("Data/storks_behaviour.csv")

mig_dates <- read_xlsx("Data/migration_dates_v02.xlsx") |>
  mutate(across(matches("start|end"), as.Date))


# --- 3. Calculate the final year wintering MCPs -------------------------------
winter_periods <- mig_dates |>
  pivot_longer(
    cols = matches("^(start_spring|end_autumn)_\\d+$"),
    names_to = c(".value", "event_num"),
    names_pattern = "(start_spring|end_autumn)_(\\d+)"
  ) |>
  dplyr::filter(!is.na(start_spring), !is.na(end_autumn)) |>
  mutate(event_num = as.integer(event_num)) |>
  group_by(ID) |>
  slice_max(event_num, n = 1) |>
  ungroup() |>
  transmute(
    ID,
    period_start = as.Date(end_autumn) + 1,
    period_end   = as.Date(start_spring) - 1
  )

tracking_winter <- all_data |>
  inner_join(winter_periods, by = "ID") |>
  filter(timestamp > period_start & timestamp < period_end)

# 2. Convert to SpatialPointsDataFrame
sp_points <- tracking_winter |>
  dplyr::select(ID, long, lat) 

coordinates(sp_points) <- ~long + lat
proj4string(sp_points) <- CRS("+proj=longlat +datum=WGS84")

# 3. Compute 100% MCP per individual
mcp_100 <- mcp(sp_points[, "ID"], percent = 100)

# 4. Convert to sf object for easy plotting/export
mcp_100_sf <- st_as_sf(mcp_100)

# save(mcp_100_sf, file = "mcp_100.RData")


# --- 4. Select other stork's data inside the MCPS -----------------------------
# Here we have the tracking data of all other tracked storks inside Xarrama's and Malpartida's MCPs, at the same time that Xarrama and Malpartida were using it

other_data <- fread("Data/other_birds_mcp.csv")


# --- 5. Calculate the behaviour summaries -------------------------------------
focal <- c("1196", "1259")
ads <- c("O584_Barlavento", "1196", "1198_Ermida", "1259", "1317_Alvalade",
         "Pequenino_O103")

# 0. prepare dataset
other_data <- other_data |>
  mutate(
    date = as.Date(timestamp),
    focal = ifelse(individual_local_identifier %in% focal, "focal", "non-focal"),
    focal = as.factor(focal),
    age_cat = ifelse(individual_local_identifier %in% ads, "adult", "immature"),
    age_cat = as.factor(age_cat)
  )

# 1. Create daily datasets
all_daily_beh <- summarise_daily_behaviour(other_data)


# 2. Split into two separate dataframes, one per individual's MCP  
xarrama_daily_behaviour <- all_daily_beh |> 
  filter(mcp_id == "Xarrama")

malpartida_daily_behaviour <- all_daily_beh |> 
  filter(mcp_id == "Malpartida")



# --- 6. Model the behaviour of the focal vs non-focal storks ------------------
output_focal <- run_focal_behaviour_models(
  list(Malpartida = malpartida_daily_behaviour,
       Xarrama    = xarrama_daily_behaviour),
  skip_models = "Xarrama | flying"
)



# --- 7. Build figure 4 --------------------------------------------------------
# 0. Create significance labels
focal_comparison_results <- output_focal$results |>
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE ~ " "
  ))

# 1. Extract model predictions for plotting
focal_predictions <- extract_focal_predictions(output_focal$models)

# 2. Join significance codes
sig_labels <- focal_comparison_results |>
  filter(str_detect(term, "^focal")) |> 
  dplyr::select(ID, behavior, significance)

label_positions <- focal_predictions |>
  group_by(ID, behavior) |>
  summarise(y_pos = max(conf.high_pct, na.rm = TRUE) * 1.1, .groups = "drop") |>
  left_join(sig_labels, by = c("ID", "behavior"))

# 3. Figure 4
plot_focal_predictions <- ggplot(focal_predictions, 
                                 aes(x = behavior, y = predicted_pct, color = focal)) +
  geom_pointrange(aes(ymin = conf.low_pct, ymax = conf.high_pct),
                  position = position_dodge(width = 0.5), 
                  size = 1, linewidth = 1) +
  geom_text(data = label_positions, 
            aes(x = behavior, y = y_pos, label = significance),
            inherit.aes = FALSE,
            size = 6, fontface = "bold", colour = "black") +
  facet_wrap(~ID, scales = "free_x") +
  scale_color_manual(values = c("focal" = "#B88244FF", "non-focal" = "#B8B69EFF"),
                     name = "Individual",
                     labels = c("focal" = "Focal", "non-focal" = "Others")) +
  labs(x = "Behaviour", y = "Predicted % of time") +
  common_theme +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

print(plot_focal_predictions)


ggsave("Figures/figure_4.png", 
       plot = plot_focal_predictions, 
       width = 15, height = 8, 
       dpi = 300)










