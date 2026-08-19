
# ------------------------------------------------------------------------------
# -------- Describe the migrations ---------------------------------------------
# ------------------------------------------------------------------------------

# --- 1. Set up the R session --------------------------------------------------
source("0.setup.R")
source("0.functions.R")


# --- 2. Data upload -----------------------------------------------------------
all_data <- fread("Data/storks_behaviour.csv")

mig_dates <- read_xlsx("Data/migration_dates_v02.xlsx") |>
  mutate(across(matches("start|end"), as.Date))



# --- 3. Figure 1 --------------------------------------------------------------
Sys.setlocale("LC_TIME", "English")

# 1. Prepare data and calculate age per individual
# Note: age changes in May every year
all_data_age <- all_data |>
  mutate(current_year = year(timestamp),
         current_month = month(timestamp)) |>
  group_by(ID) |>
  mutate(
    age = ifelse(ID == "Bras", 
                 (current_year - 2018) + 1, 
                 (current_year - 2019) + 2),
    age = ifelse(current_month %in% seq(1, 4, 1), age-1, age),
    age_fact = as.factor(age)
  ) |>
  ungroup()

# To keep the maps light, show only the months of main migratory activity 
# (August-November)
map_data_age <- all_data_age |>
  filter(month(timestamp) >= 8 & month(timestamp) <= 11)

# 2. Map Clipping (Standard Buffer)
buffer <- 1
world_clip <- map_data(map = "world") |>
  filter(long >= (min(map_data_age$long) - buffer),
         long <= (max(map_data_age$long) + buffer),
         lat >= (min(map_data_age$lat) - buffer),
         lat <= (max(map_data_age$lat) + buffer))

# 3. Create the maps 
plot_maps_clean <- ggplot() +
  geom_polygon(data = world_clip, aes(x = long, y = lat, group = group), 
               fill = "grey90", color = "grey80", linewidth = 0.1) +
  geom_path(data = map_data_age, 
            aes(x = long, y = lat, color = age_fact, 
                group = current_year),
            linewidth = 0.8, alpha = 0.8) +
  geom_point(data = map_data_age, 
             aes(x = long, y = lat, color = age_fact),
             size = 0.3, alpha = 0.4) +
  facet_wrap(~ID, ncol = 3) +
  scale_colour_manual(values = age_palette, name = "Age (Years)") +
  labs(y = "Latitude",
       x = "Longitude") +
  common_theme +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_blank())

# 4. Latitude graphs
plot_lat_clean <- ggplot(all_data_age, aes(x = timestamp, y = lat, 
                                           color = age_fact)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~ID, scales = "free_x", ncol = 3) + 
  scale_colour_manual(values = age_palette, name = "Age (Years)") +
  scale_x_datetime(breaks = date_breaks("6 months"), 
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
  guides(color = guide_legend(nrow = 1, 
                              title.position = "left", 
                              title.vjust = 0.5,
                              label.position = "right", 
                              keywidth = unit(1, "lines"),
                              keyheight = unit(0.5, "lines"),
                              override.aes = list(linewidth = 5)))

# 5. Final figure assembly
final_figure <- (plot_maps_clean / plot_lat_clean) + 
  plot_layout(heights = c(2.6, 1)) & 
  theme(plot.margin = margin(5, 5, 5, 5))

# Display and save
print(final_figure)

ggsave("Figures/fig_1.png", 
       final_figure, 
       width = 14, height = 10.5, 
       dpi = 300)

