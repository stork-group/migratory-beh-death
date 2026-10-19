
# ------------------------------------------------------------------------
# - Sup. Material - Other storks' migration phenology and performance ----
# ------------------------------------------------------------------------

# --- 1. Set up the R session --------------------------------------------
source("Scripts_to_publish/0.setup.R")
source("Scripts_to_publish/0.functions.R")
Sys.setlocale("LC_TIME", "English")


# --- 2. Data upload ----------------------------------------------------
other_data <- fread("Data/other_storks.csv")

mig_dates <- read_xlsx("Data/migration_dates_v03.xlsx") |>
  mutate(across(matches("^(start|end)_"),
                ~ {
                  if (is.character(.x)) as.Date(as.numeric(.x), origin = "1899-12-30") else as.Date(.x)
                }))
mig_dates <- filter(mig_dates, focal == "no")


# --- 3. Map the migrations of other storks ------------------------------
# 1. Prepare data and calculate tracking year per individual
other_data_year <- other_data |>
  mutate(current_year = year(timestamp),
         current_month = month(timestamp)) |>
  group_by(ID) |>
  mutate(
    track_year = ifelse(ID == "Barlavento", 
                        (current_year - 2018) + 1, 
                        (current_year - 2019) + 1),
    track_year = ifelse(current_month %in% seq(1, 4, 1), track_year-1, track_year),
    track_year_fact = as.factor(track_year)
  ) |>
  ungroup()


# To keep the maps light, show only the months of main migratory activity 
map_data_year <- other_data_year |>
  filter(!(month(timestamp) %in% seq(2, 6, 1)))

# 2. Map Clipping (Standard Buffer)
buffer <- 1
world_clip <- map_data(map = "world") |>
  filter(long >= (min(map_data_year$long) - buffer),
         long <= (max(map_data_year$long) + buffer),
         lat >= (min(map_data_year$lat) - buffer),
         lat <= (max(map_data_year$lat) + buffer))

# 3. Create the maps 
plot_maps_clean <- ggplot() +
  geom_polygon(data = world_clip, aes(x = long, y = lat, group = group), 
               fill = "grey90", color = "grey80", linewidth = 0.1) +
  geom_path(data = map_data_year, 
            aes(x = long, y = lat, color = track_year_fact, 
                group = current_year),
            linewidth = 0.8, alpha = 0.8) +
  geom_point(data = map_data_year, 
             aes(x = long, y = lat, color = track_year_fact),
             size = 0.3, alpha = 0.4) +
  facet_wrap(~ID, ncol = 4) +
  scale_colour_manual(values = age_palette, name = "Age (Years)") +
  labs(y = "Latitude",
       x = "Longitude") +
  common_theme +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_blank()) + 
  guides(color = "none")

# 4. Latitude graphs
plot_lat_clean <- ggplot(other_data_year, aes(x = timestamp, y = lat, 
                                           color = track_year_fact)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~ID, scales = "free_x", ncol = 4) + 
  scale_colour_manual(values = age_palette, name = "Tracked year") +
  scale_x_datetime(breaks = date_breaks("12 months"), 
                   labels = date_format("%b %y")) +
  labs(y = "Latitude",
       x = "Date") +
  common_theme +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(t = -5)) +
  guides(colour = "none", fill = "none")


# ---- 4. Plot the migration phenology -----------------------------------
# 1. Reshape to long format
mig_dates_long <- mig_dates |>
  pivot_longer(
    cols = -c(ID, year),
    names_to = c(".value", "tracked_year"),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  mutate(tracked_year = as.integer(tracked_year)) |>
  filter(!is.na(start_autumn))

# 2. Wrap dates onto pseudo-year (May -> April) 
mig_dates_long <- mig_dates_long |>
  mutate(
    start_autumn_plot = make_plot_date(start_autumn),
    end_spring_plot   = make_plot_date(end_spring))


# 3. Combine all point events (migration points) into one df
migration_points <- mig_dates_long |>
  dplyr::select(ID, tracked_year, start_autumn_plot, end_spring_plot) |>
  pivot_longer(
    cols = c(start_autumn_plot, end_spring_plot),
    names_to = "event", values_to = "date"
  ) |>
  mutate(event = recode(event,
                        start_autumn_plot = "Start autumn migration",
                        end_spring_plot   = "End spring migration"))

# 4. Make the figure
mig_dur_plot <- 
  ggplot() +
  geom_segment(
    data = mig_dates_long,
    aes(x = start_autumn_plot, xend = end_spring_plot, y = tracked_year,
        yend = tracked_year),
    color = "grey50", linewidth = 0.6
  ) +
  geom_point(
    data = migration_points,
    aes(x = date, y = tracked_year, colour = event),
    size = 4, stroke = 1.5
  ) +
  facet_wrap(~ID, ncol = 4) +
  scale_x_date(
    date_breaks = "2 months", date_labels = "%b",
    limits = c(as.Date("2000-07-01"), as.Date("2001-02-15")),
    expand = expansion(mult = 0.02)
  ) +
  scale_colour_manual(values = c(
    "Start autumn migration" = "deepskyblue3",
    "End spring migration"   = "goldenrod1"
  ), name = "Phenology") +
  scale_shape_manual(values = c(
    "Start autumn migration" = 16, 
    "End spring migration"   = 16
  ), name = "Phenology") +
  scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]),
                                              by = 1)) +
  labs(x = "Date", y = "Tracked year", colour = NULL, shape = NULL) +
  common_theme +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom")

print(mig_dur_plot)


# --- 5. Migration performance -------------------------------------------
# 1. Calculate migration distance and duration
# Add a date column to the GPS data so it can be matched to the migration dates
other_data <- other_data %>%
  mutate(date = as.Date(timestamp, tz = "UTC"))

# Keep only the columns needed from the migration dates
autumn <- mig_dates_long %>%
  dplyr::select(ID, tracked_year, start_autumn, end_autumn)

# First fix on the start_autumn date
start_pts <- autumn %>%
  inner_join(other_data, by = c("ID", "start_autumn" = "date")) %>%
  group_by(ID, tracked_year) %>%
  slice_min(timestamp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(ID, tracked_year,
                start_time = timestamp, start_lat = lat, start_long = long)

# Last fix on the end_autumn date
end_pts <- autumn %>%
  inner_join(other_data, by = c("ID", "end_autumn" = "date")) %>%
  group_by(ID, tracked_year) %>%
  slice_max(timestamp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(ID, tracked_year,
                end_time = timestamp, end_lat = lat, end_long = long)

# Join start and end points and calculate the net distance in km
autumn_perf <- autumn %>%
  left_join(start_pts, by = c("ID", "tracked_year")) %>%
  left_join(end_pts,   by = c("ID", "tracked_year")) %>%
  mutate(dist_net_km = distGeo(cbind(start_long, start_lat),
                               cbind(end_long,   end_lat)) / 1000,
         duration = as.numeric(difftime(as.Date(end_autumn),
                                        as.Date(start_autumn), 
                                        units = "days")),
         mig_speed = dist_net_km / duration,
         year = year(start_autumn))

autumn_perf

# 2. Make the figure
mig_speed_plot <- ggplot(autumn_perf, 
                         aes(x = year, y = mig_speed, 
                             fill = factor(tracked_year))) +
  geom_path(aes(group = ID), linetype = "dashed", alpha = 0.3) + 
  geom_point(size = 4, shape = 21, color = "black") + 
  scale_fill_manual(values = age_palette, name = "Tracked year") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)),
                     breaks = seq(2018,2025,2)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Year",
    y = "Migration speed (Km / days)"
  ) +
  facet_wrap(~ID, ncol = 4) +
  guides(fill = guide_legend(nrow = 1)) +
  common_theme + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_blank())

print(mig_speed_plot)


# ---- 6. Figure Assembly ----------------------------------------------

sup_fig_s1 <- (plot_maps_clean / plot_lat_clean / mig_dur_plot / mig_speed_plot) +
  plot_layout(heights = c(2.6, 1.5, 1.5, 1.5), guides = "collect") &
  plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom",
        legend.box = "horizontal")

# Display and save
print(sup_fig_s1)

ggsave("Figures/sup_fig_s1.png", 
       sup_fig_s1, 
       width = 14, height = 18, 
       dpi = 300)


