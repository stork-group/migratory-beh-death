
# ------------------------------------------------------------------------------
# ---- List of R setup options and libraries -----------------------------------
# ------------------------------------------------------------------------------

# ---- List of libraries -------------------------------------------------------
library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)
library(viridis)
library(patchwork)
library(scales)
library(geosphere)
library(glmmTMB)
library(DHARMa)
library(broom.mixed)
library(ggeffects)
library(sf)
library(adehabitatHR)


# ---- Plotting customization --------------------------------------------------
# Create a common theme for all ggplot graphs
common_theme <- theme_bw(base_size = 18) + 
  theme(strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey95"),
        panel.border = element_rect(color = "grey20", fill = NA, 
                                    linewidth = 0.8),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Create the palettes
beh_palette   <- c("Flying" = "#3F858CFF", "Foraging" = "#F2D43DFF", 
                   "Resting" = "#D9814EFF")

age_palette   <- c("1" = "#440154", "2" = "#443983", "3" = "#31688e",
                   "4" = "#21918c", "5" = "#35b779", "6" = "#90d743", 
                   "7" = "#fde725")


