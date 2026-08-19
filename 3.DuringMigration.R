
# ------------------------------------------------------------------------------
# -------- Changes in migration phenology and performance ----------------------
# ------------------------------------------------------------------------------

# --- 1. Set up the R session --------------------------------------------------
source("Scripts_to_publish/0.setup.R")
source("Scripts_to_publish/0.functions.R")

# --- 2. Data upload -----------------------------------------------------------
all_data <- fread("Data/storks_behaviour.csv")

mig_dates <- read_xlsx("Data/migration_dates_v02.xlsx") |>
  mutate(across(matches("start|end"), as.Date))


# --- 3. Migration phenology ---------------------------------------------------
# 0. Add death date to the original wide df ----
mig_dates <- mig_dates |>
  mutate(death_date = as.Date(c("2025-01-25", "2024-11-25", "2022-10-11")))  
# order matches df$ID: Xarrama, Malpartida, Bras

# 1. Reshape to long format
mig_dates_long <- mig_dates |>
  pivot_longer(
    cols = -c(ID, year, age, death_date),
    names_to = c(".value", "tracked_year"),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  mutate(tracked_year = as.integer(tracked_year)) |>
  filter(!is.na(start_autumn))

# 2. Add the age of the individuals
mig_dates_long <- mig_dates_long |>
  group_by(ID) |>
  mutate(
    start_age = case_when(ID == "Bras" ~ 1, TRUE ~ 2),
    age_num   = start_age + (tracked_year - 1)
  ) |>
  ungroup()

# 3. Wrap dates onto pseudo-year (May -> April) 
mig_dates_long <- mig_dates_long |>
  mutate(
    start_autumn_plot = make_plot_date(start_autumn),
    end_spring_plot   = make_plot_date(end_spring),
    death_date_plot   = make_plot_date(death_date)
  )

# 4. Creating the death segments (end of spring migration -> death, final year only)
death_segments <- mig_dates_long |>
  group_by(ID) |>
  filter(tracked_year == max(tracked_year)) |>
  ungroup() |>
  filter(!is.na(death_date_plot))

# 5. Combine all point events (migration points + death) into one df
migration_points <- mig_dates_long |>
  dplyr::select(ID, age_num, start_autumn_plot, end_spring_plot) |>
  pivot_longer(
    cols = c(start_autumn_plot, end_spring_plot),
    names_to = "event", values_to = "date"
  ) |>
  mutate(event = recode(event,
                        start_autumn_plot = "Start autumn migration",
                        end_spring_plot   = "End spring migration"))

death_points <- death_segments |>
  transmute(ID, age_num, date = death_date_plot, event = "Death")

points_df <- bind_rows(migration_points, death_points) |>
  mutate(event = factor(event, levels = c("Start autumn migration",
                                          "End spring migration",
                                          "Death")))

# 6. Figure 3.A 
Sys.setlocale("LC_TIME", "English")

mig_dur_plot <- 
  ggplot() +
  geom_segment(
    data = mig_dates_long,
    aes(x = start_autumn_plot, xend = end_spring_plot, y = age_num, yend = age_num),
    color = "grey50", linewidth = 0.6
  ) +
  geom_segment(
    data = death_segments,
    aes(x = end_spring_plot, xend = death_date_plot, y = age_num, yend = age_num),
    color = "grey30", linewidth = 0.6, linetype = "dashed"
  ) +
  geom_point(
    data = points_df,
    aes(x = date, y = age_num, colour = event, shape = event),
    size = 4, stroke = 1.5
  ) +
  facet_wrap(~ID) +
  scale_x_date(
    date_breaks = "2 months", date_labels = "%b",
    limits = c(as.Date("2000-05-01"), as.Date("2001-04-30")),
    expand = expansion(mult = 0.02)
  ) +
  scale_colour_manual(values = c(
    "Start autumn migration" = "deepskyblue3",
    "End spring migration"   = "goldenrod1",
    "Death"                  = "firebrick"
  ), name = "Phenology") +
  scale_shape_manual(values = c(
    "Start autumn migration" = 16, 
    "End spring migration"   = 16, 
    "Death"                  = 4  
  ), name = "Phenology") +
  scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1)) +
  labs(x = NULL, y = "Age (years)", colour = NULL, shape = NULL) +
  common_theme +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right") 

print(mig_dur_plot)


# --- 4. Calculate delay in timing of migration -------------------------------- 
# 0. Reference point for day-counting: May 1st of the pseudo-year
ref_date <- as.Date("2000-05-01")

mig_dates_long <- mig_dates_long |>
  mutate(
    start_autumn_day = as.numeric(start_autumn_plot - ref_date),
    end_spring_day   = as.numeric(end_spring_plot - ref_date)
  )

# 1. Compute median "typical" timing per individual, excluding the final year
typical_timing <- mig_dates_long |>
  group_by(ID) |>
  dplyr::filter(tracked_year != max(tracked_year)) |> 
  summarise(
    median_start_autumn_day = median(start_autumn_day, na.rm = TRUE),
    median_end_spring_day   = median(end_spring_day, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Extract the last year's actual dates
last_year <- mig_dates_long |>
  group_by(ID) |>
  dplyr::filter(tracked_year == max(tracked_year)) |>
  ungroup() |>
  dplyr::select(ID, tracked_year, last_start_autumn_day = start_autumn_day,
         last_end_spring_day = end_spring_day)

# 3. Join and compute delay in days
delay_summary <- typical_timing |>
  left_join(last_year, by = "ID") |>
  mutate(
    delay_start_autumn = last_start_autumn_day - median_start_autumn_day,
    delay_end_spring    = last_end_spring_day - median_end_spring_day
  ) |>
  dplyr::select(ID, tracked_year, 
         median_start_autumn_day, last_start_autumn_day, delay_start_autumn,
         median_end_spring_day, last_end_spring_day, delay_end_spring)

# 4. Convert day-counts back to readable dates
delay_summary <- delay_summary |>
  mutate(
    median_start_autumn_date = format(ref_date + median_start_autumn_day, "%d %b"),
    last_start_autumn_date   = format(ref_date + last_start_autumn_day, "%d %b"),
    median_end_spring_date   = format(ref_date + median_end_spring_day, "%d %b"),
    last_end_spring_date     = format(ref_date + last_end_spring_day, "%d %b")
  )



# --- 4. Migration performance -------------------------------------------------
# 0. We first define the stop day and latitude of the terminal year
terminal_config <- tribble(
  ~ID,          ~term_year, ~manual_julian,
  "Bras",       2022,       252,
  "Malpartida", 2024,       247,
  "Xarrama",    2024,       239
)

# 1. And find when did the birds crossed the stop latitude in previous years
stop_latitudes <- all_data |>
  mutate(year = year(timestamp),
         julian_day = yday(timestamp)) |>
  inner_join(terminal_config, by = c("ID" = "ID", "year" = "term_year")) |>
  group_by(ID) |>
  filter(abs(julian_day - manual_julian) == min(abs(julian_day - manual_julian))) |>
  slice(1) |> 
  summarise(final_stop_lat = lat, .groups = "drop")

migration_performance <- mig_dates_long |>
  mutate(year = year(start_autumn)) |>
  inner_join(stop_latitudes, by = "ID") |>
  inner_join(all_data, by = "ID", relationship = "many-to-many") |>
  dplyr::select(-c(ID_raw, predBeh)) |>
  dplyr::filter(timestamp >= start_autumn, timestamp <= start_autumn + days(90)) |>
  arrange(ID, start_autumn, timestamp) |>
  group_by(ID, start_autumn) |>
  mutate(start_lat = first(lat),
         start_long = first(long)) |>
  dplyr::filter(
    ifelse(final_stop_lat < start_lat, lat <= final_stop_lat, lat >= final_stop_lat)
  ) |>
  slice(1) 

# 3. Correct Xarrama's data 
# In 2019, 2020 and 2021, Xarrama does a longitudinal migration first and only then migrates south, so the stop latitude did not work. So we manually corrected the timestamps and locations
migration_performance <- migration_performance |>
  mutate(
    timestamp = case_when(
      ID == "Xarrama" & year == 2019 ~ 
        as.POSIXct("2019-05-17 10:17:02", tz = "UTC"),
      ID == "Xarrama" & year == 2020 ~ 
        as.POSIXct("2020-06-12 14:34:46", tz = "UTC"),
      ID == "Xarrama" & year == 2021 ~ 
        as.POSIXct("2021-07-11 12:18:08", tz = "UTC"),
      TRUE ~ timestamp),
    lat = case_when(
      ID == "Xarrama" & year == 2019 ~ 37.1791183,
      ID == "Xarrama" & year == 2020 ~ 37.450965,
      ID == "Xarrama" & year == 2021 ~ 37.44885,
      TRUE ~ lat),
    long = case_when(
      ID == "Xarrama" & year == 2019 ~ -6.7206316,
      ID == "Xarrama" & year == 2020 ~ -6.638593,
      ID == "Xarrama" & year == 2019 ~ -6.642215,
      TRUE ~ long))

# 4. Calculate migration performance metrics
migration_performance <- migration_performance |>
  ungroup() |>
  mutate(
    days_to_reach = as.numeric(difftime(timestamp, start_autumn, units = "days"))) |>
  rename(arrival_date = timestamp,
         stop_lat = lat,
         stop_long = long
         ) |>
  filter(days_to_reach > 0 & days_to_reach < 70)

final_metrics <- migration_performance |>
  rowwise() |> 
  mutate(metrics = list({
    t_start <- start_autumn
    t_end   <- arrival_date
    curr_id <- ID
    
    segment <- all_data |>
      filter(ID == curr_id, timestamp >= t_start, timestamp <= t_end) |>
      arrange(timestamp)
    
    dist_total_km <- segment |>
      mutate(d = distGeo(cbind(long, lat), 
                         cbind(lag(long), lag(lat)))) |>
      summarise(total = sum(d, na.rm = TRUE) / 1000) |> pull(total)
    
    coords_start <- c(start_long, start_lat)
    coords_end   <- c(stop_long, stop_lat)
    dist_net_km  <- distGeo(coords_start, coords_end) / 1000
    
    data.frame(
      dist_net_km = dist_net_km,
      total_km = dist_total_km,
      tortuosity = dist_net_km / dist_total_km,
      mig_speed = dist_net_km / days_to_reach
    )
  })) |>
  unnest(metrics) |> 
  ungroup()

# 5. Build Figure 3.B
mig_speed_plot <- ggplot(final_metrics, 
                         aes(x = year, y = mig_speed, fill = factor(age_num))) +
  geom_path(aes(group = ID), linetype = "dashed", alpha = 0.3) + 
  geom_point(size = 4, shape = 21, color = "black") + 
  scale_fill_manual(values = age_palette, name = "Age (years)") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)),
                     breaks = seq(2018,2024,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Year",
    y = "Migration speed (Km / days)"
  ) +
  facet_wrap(~ID) +
  guides(fill = guide_legend(nrow = 1)) +
  common_theme + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_blank())


# ---- 5. Build Figure 3 ------------------------------------------------------- 
fig.3 <- mig_dur_plot / mig_speed_plot + 
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.location = "panel",
    plot.tag.position  = c(0.02, 0.93),
    plot.tag = element_text(size = 20),
    legend.position = "right",
    legend.justification = "left",
    legend.box.just = "left"
  ) & guides(fill = guide_legend(ncol = 1))

fig.3

ggsave("Figures/figure_4.png",
       plot = fig.3,
       width = 15, height = 9,
       dpi = 300)


# ---- 6. Build Supplementary Table S1 -----------------------------------------
# 1. Start building the table, with the migration dates
mig_dates_long <- mig_dates |>
  dplyr::select(ID, year, starts_with("start_autumn"), starts_with("end_autumn"),
         starts_with("start_spring"), starts_with("end_spring")) |>
  pivot_longer(
    cols = -c(ID, year),
    names_to = c(".value", "tracked_year"),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  dplyr::select(-year) |>
  mutate(tracked_year = as.integer(tracked_year)) |>
  rename(year = tracked_year) |>
  filter(!is.na(start_autumn)) |>
  mutate(year = year(start_autumn),
         age = year - 2018 + 1) |>
  arrange(ID) |>
  relocate(age, .after = year)

# 2. Select only the autumn migration period 
points_in_autumn <- all_data |>
  inner_join(mig_dates_long |> 
               dplyr::select(ID, year, start_autumn, end_autumn),
             by = "ID", relationship = "many-to-many") |>
  filter(timestamp >= start_autumn & timestamp <= end_autumn) |>
  arrange(ID, year, timestamp)

# 3. Calculate the autumn migration bee-line distance
bee_line_df <- points_in_autumn |>
  group_by(ID, year) |>
  arrange(timestamp) |>
  summarise(
    bee_line_km = round(distGeo(c(first(long), first(lat)),
                                c(last(long),  last(lat))) / 1000,1),
    .groups = "drop"
  )

# 4. Total travelled distance
total_dist_df <- points_in_autumn |>
  group_by(ID, year) |>
  mutate(
    step_dist_km = distGeo(cbind(long, lat), cbind(lag(long), lag(lat))) / 1000
  ) |>
  summarise(
    total_dist_km = round(sum(step_dist_km, na.rm = TRUE),1),
    .groups = "drop"
  )

# 5. Create Supplementary Table S1
migration_summary <- mig_dates_long |>
  left_join(bee_line_df, by = c("ID", "year")) |>
  left_join(total_dist_df, by = c("ID", "year"))

# 6. Add the final year's migration performance metrics
migration_summary <- migration_summary |>
  left_join(final_metrics[, c("ID", "year", "death_date", "days_to_reach",
                              "dist_net_km", "total_km", "mig_speed")],
            by = c("ID", "year")) |>
  rename(final_mig_dur = days_to_reach,
         final_net_dist_km = dist_net_km,
         final_tot_dist_km = total_km,
         final_mig_speed = mig_speed) |>
  relocate(death_date, .after = age) |>
  mutate(across(where(is.numeric), ~round(.x, 1)))


write.csv(migration_summary, "sup_table_s1.csv", row.names = FALSE)



