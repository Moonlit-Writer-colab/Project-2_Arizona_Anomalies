# Helper functions for the Arizona Anomalies season analysis.
#
# The same functions serve both research questions:
#   - Question 1 (extreme heat): flag hot days, find the first and last hot day
#     in each calendar year, and plot them.
#   - Question 2 (freezes): flag freezing nights and run the same event-date
#     functions on the July-June cold-season year instead of the calendar year.
#
# Each function takes the daily data frame and a column name, so adding a new
# threshold or a new technique means adding a column, not rewriting the code.

library(tidyverse)

# ---- Settings shared by every chart ----------------------------------------

city_levels <- c("Tucson", "Phoenix", "Flagstaff")

city_lookup <- c(
  "Tucson International Airport" = "Tucson",
  "Phoenix Airport"              = "Phoenix",
  "Flagstaff Airport"            = "Flagstaff"
)

early_years  <- 1982:1991
recent_years <- 2015:2024
base_years   <- 1991:2020   # climate-normal period for city-relative thresholds

# Heat intensity: one sequential ramp, light yellow (warm) to deep red
# (extreme). The steps follow the thresholds in the proposal, so every tile
# tells you how hot the day was, not just whether it passed 100°F.
heat_levels <- c(
  "Below 80°F"    = "#ebeae6",
  "80–89°F"       = "#ffeda0",
  "90–99°F"       = "#fdbf5c",
  "100–104°F"     = "#f6803b",
  "105–109°F"     = "#d7301f",
  "110°F or more" = "#7f0000",
  "Missing reading" = "#ffffff"
)

heat_level <- function(tmax_f) {
  cut(tmax_f, breaks = c(-Inf, 80, 90, 100, 105, 110, Inf), right = FALSE,
      labels = names(heat_levels)[1:6]) |>
    fct_na_value_to_level("Missing reading") |>
    fct_expand(names(heat_levels))
}

# One fixed colour per city, used in every chart (checked for colour-blind
# separation). Blues are kept free for the freeze charts in question 2, and
# reds/oranges for heat intensity.
city_colors <- c(
  "Tucson"    = "#6250d6",  # violet
  "Phoenix"   = "#e87ba4",  # pink
  "Flagstaff" = "#008300"   # green
)

text_primary   <- "#0b0b0b"
text_secondary <- "#52514e"

theme_anomalies <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", color = text_primary),
      plot.subtitle      = element_text(color = text_secondary),
      plot.caption       = element_text(color = text_secondary, hjust = 0),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.text          = element_text(color = text_secondary),
      strip.text         = element_text(face = "bold", color = text_primary,
                                        hjust = 0, size = base_size),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(color = "#ecebe7", linewidth = 0.3),
      legend.position    = "top",
      legend.justification = "left",
      plot.background    = element_rect(fill = "white", color = NA)
    )
}

data_caption <- "Source: NOAA GHCN-Daily, airport stations, 1982–2024."

# Put any date on the leap-year 2000 calendar so dates align across years.
to_plot_date <- function(x) {
  as.Date(paste0("2000-", format(x, "%m-%d")), format = "%Y-%m-%d")
}

# ---- Data ------------------------------------------------------------------

# Read the cleaned daily file and add the columns every analysis needs.
# plot_date puts every date on a leap-year calendar (year 2000), so Feb 29
# keeps a fixed position and dates line up across years.
load_weather <- function(path = NULL) {
  # Works when run from the repo root or from scripts/ (where Quarto renders).
  path <- path %||% c("data/arizona_weather.csv", "../data/arizona_weather.csv") |>
    keep(file.exists) |> first()
  read_csv(path, show_col_types = FALSE) |>
    mutate(
      city        = factor(city_lookup[station_name], levels = city_levels),
      date        = as.Date(date),
      plot_date   = to_plot_date(date),
      cold_season = if_else(month >= 7, year, year - 1L)  # used in question 2
    )
}

# ---- Thresholds ------------------------------------------------------------

# City-relative hot-day threshold: the city's 90th-percentile June–August high
# over the base period (the TX90p idea, Zhang et al. 2011).
relative_heat_thresholds <- function(weather, prob = 0.90, years = base_years) {
  weather |>
    filter(month %in% 6:8, year %in% years) |>
    group_by(city) |>
    summarise(rel_hot_f = unname(quantile(tmax_f, prob, na.rm = TRUE)),
              .groups = "drop")
}

# Add one logical flag per hot-day definition. NA readings stay NA.
flag_hot_days <- function(weather, fixed_f = 100, thresholds = NULL) {
  thresholds <- thresholds %||% relative_heat_thresholds(weather)
  weather |>
    left_join(thresholds, by = "city") |>
    mutate(
      hot_fixed    = tmax_f >= fixed_f,
      hot_relative = tmax_f >= rel_hot_f
    )
}

# ---- Event dates -----------------------------------------------------------

# For each city and period, find the first and last day where `flag` is TRUE,
# the span between them (both ends counted), and how many days qualified.
# `period` is the grouping year: `year` for heat, `cold_season` for freezes.
# Periods with no qualifying day keep NA dates and n_days = 0.
event_dates <- function(weather, flag, period = year) {
  weather |>
    group_by(city, period = {{ period }}) |>
    summarise(
      n_days     = sum({{ flag }}, na.rm = TRUE),
      n_missing  = sum(is.na({{ flag }})),
      first_date = if (any({{ flag }}, na.rm = TRUE))
                     min(date[{{ flag }} %in% TRUE]) else as.Date(NA),
      last_date  = if (any({{ flag }}, na.rm = TRUE))
                     max(date[{{ flag }} %in% TRUE]) else as.Date(NA),
      .groups = "drop"
    ) |>
    mutate(
      span_days  = as.numeric(last_date - first_date) + 1,
      first_plot = to_plot_date(first_date),
      last_plot  = to_plot_date(last_date)
    )
}

# Label each year as early, recent, or neither.
add_period_group <- function(df, year_col = period) {
  df |>
    mutate(group = case_when(
      {{ year_col }} %in% early_years  ~ "1982–1991",
      {{ year_col }} %in% recent_years ~ "2015–2024",
      TRUE ~ NA_character_
    ))
}

# Linear trend per city, reported as change per decade.
trend_per_decade <- function(df, value, x = period) {
  df |>
    transmute(city, y = as.numeric({{ value }}), x = {{ x }}) |>
    filter(!is.na(y)) |>
    group_by(city) |>
    summarise(
      n_years = n(),
      slope_per_decade = if (n() > 2) 10 * unname(coef(lm(y ~ x))[2]) else NA_real_,
      .groups = "drop"
    )
}

# Centred running median that ignores missing years.
slider_median <- function(x, k = 9) {
  half <- k %/% 2
  n <- length(x)
  map_dbl(seq_len(n), function(i) {
    if (i - half < 1 || i + half > n) return(NA_real_)
    median(x[(i - half):(i + half)], na.rm = TRUE)
  })
}

# ---- Charts ----------------------------------------------------------------

# Season calendar: one row per year, one tile per date, shaded by how hot the
# day's high was. Thin lines trace the first and last hot day of each year.
plot_heat_calendar <- function(hot_weather, hot_events, cities = city_levels,
                               title = NULL, subtitle = NULL) {
  tiles <- hot_weather |>
    filter(city %in% cities) |>
    mutate(tile = heat_level(tmax_f))

  traces <- hot_events |>
    filter(city %in% cities) |>
    select(city, year = period, first_plot, last_plot, measure) |>
    pivot_longer(c(first_plot, last_plot), names_to = "event",
                 values_to = "plot_date") |>
    mutate(event = if_else(event == "first_plot", "First hot day", "Last hot day"))

  labels <- traces |>
    filter(!is.na(plot_date)) |>
    group_by(city, event) |>
    slice_max(year, n = 1) |>
    ungroup()

  ggplot(tiles, aes(x = plot_date, y = year)) +
    geom_tile(aes(fill = tile), height = 0.9, width = 1) +
    geom_path(data = traces, aes(group = event), color = text_primary,
              linewidth = 0.35, na.rm = TRUE) +
    geom_text(data = labels, aes(label = event),
              hjust = if_else(labels$event == "First hot day", 1.08, -0.08),
              vjust = 0.5, size = 3, color = text_primary) +
    facet_wrap(~ city, nrow = 1) +
    scale_fill_manual(values = heat_levels, drop = FALSE,
                      name = "Daily high") +
    scale_x_date(date_labels = "%b",
                 breaks = seq(as.Date("2000-01-01"), as.Date("2000-12-01"), by = "2 months"),
                 limits = as.Date(c("1999-12-31", "2001-01-01")),
                 expand = c(0, 0)) +
    scale_y_reverse(breaks = seq(1985, 2020, 5), expand = c(0, 0)) +
    guides(fill = guide_legend(nrow = 1, title.position = "left")) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle,
         caption = data_caption) +
    theme_anomalies() +
    theme(panel.grid = element_blank(),
          panel.spacing.x = unit(1.2, "lines"),
          legend.title = element_text(color = text_secondary, size = 10),
          legend.key = element_rect(color = "#d6d5d0", linewidth = 0.2))
}
