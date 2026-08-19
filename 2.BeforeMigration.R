
# ------------------------------------------------------------------------------
# -------- Changes before migration --------------------------------------------
# ------------------------------------------------------------------------------

# --- 1. Set up the R session --------------------------------------------------
source("0.setup.R")
source("0.functions.R")

# --- 2. Data upload -----------------------------------------------------------
all_data <- fread("Data/storks_behaviour.csv")

mig_dates <- read_xlsx("Data/migration_dates_v02.xlsx") |>
  mutate(across(matches("start|end"), as.Date))


# 1. Create a window for 15-days before the start of autumn migration
migration_windows <- mig_dates |>
  pivot_longer(cols = contains("start_autumn"), 
               names_to = "mig_event", 
               values_to = "mig_start_date") |>
  dplyr::filter(!is.na(mig_start_date)) |>
  dplyr::select(ID, mig_start_date) |>
  crossing(days_before = c(15)) |> 
  mutate(
    window_start = as.POSIXct(mig_start_date - days(days_before)),
    window_end = as.POSIXct(mig_start_date) - 1,
    window_label = factor(paste0("Window: ", days_before, " Days"),
                          levels = c("Window: 15 Days"))
  )

# 2. Year of death for all individuals, to set as intercept of the models
death_years <- c("Bras" = 2022, "Malpartida" = 2024, "Xarrama" = 2024)



# --- 3. Changes in movement behaviour -----------------------------------------
# 1. Calculate the distance between consecutive positions
storks_distances <- all_data |>
  group_by(ID) |>
  arrange(ID, timestamp) |>
  mutate(
    dist_to_prev = distGeo(
      cbind(long, lat),
      cbind(lag(long), lag(lat))
    ),
    date = as.Date(timestamp)
  ) |>
  ungroup()


# 2. Calculate the total travelled distance per day and the net displacement
daily_metrics_full <- storks_distances |>
  left_join(migration_windows, by = "ID", relationship = "many-to-many") |>
  dplyr::filter(date >= window_start & date < window_end) |> 
  mutate(year = year(mig_start_date)) |>
  group_by(ID, year, window_label, date) |>
  summarise(
    total_dist_km = sum(dist_to_prev, na.rm = TRUE) / 1000,
    net_dist_km = distGeo(
      cbind(first(long), first(lat)),
      cbind(last(long), last(lat))
    ) / 1000,
    .groups = "drop"
  )

# 3. Run the models
out_dist_15d <- run_distance_models(daily_metrics_full, 
                                    window_label_filter = "Window: 15 Days")

distance_comparison_results <- out_dist_15d$results |>
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE ~ " "
  ))


# 4. Add age in the dataset
distance_comparison_results <- distance_comparison_results |>
  mutate(birth_year = 2018,
         age = as.numeric(term) - birth_year + 1)

# 5. Figure 2.A
plot_dist_comp <- ggplot(subset(distance_comparison_results, term != "(Intercept)"), 
                         aes(x = term, y = estimate, fill = factor(age))) +
  geom_hline(yintercept = 0, color = "black") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.2) +
  geom_text(aes(label = significance,
                vjust = ifelse(estimate >= 0, -0.6, 1.4)), 
            position = position_dodge(width = 0.9), 
            size = 3.2, fontface = "bold") +
  facet_wrap(~ID, scales = "free_x") + 
  scale_fill_manual(values = age_palette,
                    name = "Age\n(years)") + 
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    y = "Estimate (compared to year of death)", 
    x = "Year") +
  common_theme + 
  theme(legend.position = "right")
plot_dist_comp


# --- 4. Changes in behaviour --------------------------------------------------
# 1. Creating a "days-till-migration" dataset
beh_master_data <- all_data |>
  left_join(mig_dates |> 
              pivot_longer(cols = contains("start_autumn"), 
                           names_to = "event", 
                           values_to = "mig_date") |>
              dplyr::filter(!is.na(mig_date)), 
            by = "ID", relationship = "many-to-many") |>
  mutate(
    date = as.Date(timestamp),
    days_to_mig = as.numeric(date - as.Date(mig_date)),
    year = year(mig_date),
    predBeh = fct_collapse(as.factor(predBeh), 
                           Flying = c("Soaring", "Flapping"))
  )


# 3. Prepare the dataset
data_15d_model <- prepare_daily_counts(beh_master_data, 15)

# 4. Build the models
output_beh_15d <- run_beh_models(data_15d_model, "15 Days Before")

beh_comparison_results <- output_beh_15d$results |>
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.1   ~ ".",
    TRUE ~ " "
  ))

# 5. Figure 2.B
plot_beh_comp <- ggplot(subset(beh_comparison_results, term != "(Intercept)"), 
                          aes(x = term, y = estimate, fill = behavior)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = significance,
                vjust = ifelse(estimate >= 0, -0.6, 1.4)), 
            position = position_dodge(width = 0.9), 
            size = 3.2, fontface = "bold") +
  facet_wrap(~ID, scales = "free_x") + 
  scale_fill_manual(values = beh_palette,
                    name = "Behaviour") +
   scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    y = "Estimate (compared to year of death)", 
    x = "Year") +
  common_theme + 
  theme(legend.position = "right")

print(plot_beh_comp)


# --- 5. Build Figure 2 --------------------------------------------------------
plot_dist_comp_2 <- plot_dist_comp + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

plot_beh_comp_2 <- plot_beh_comp + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())


fig.2 <- (plot_dist_comp_2 / plot_beh_comp_2) + 
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.location = "panel",    
    plot.tag.position  = c(0.02, 0.93),  
    plot.tag = element_text(size = 20),
    legend.position = "right",
    legend.justification = "left",
    legend.box.just = "left") & 
  guides(fill = guide_legend(ncol = 1))

print(fig.2)

ggsave("Figures/figure_2.png", 
       plot = fig.2, 
       width = 15, height = 9, 
       dpi = 300)



