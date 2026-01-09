#!/usr/bin/env Rscript
#' BACI (Before-After-Control-Impact) Analysis for Upper Great Highway Closure
#'
#' This script implements a rigorous BACI experimental design to assess whether
#' the Upper Great Highway closure (March 14, 2025) affected traffic crash rates in
#' west side SF neighborhoods, controlling for citywide temporal variation.
#'
#' Design:
#' - Impact sites: Sunset area neighborhoods (Sunset/Parkside, Outer Sunset, Inner Sunset, Golden Gate Park)
#' - Control sites: Multiple east/south side neighborhoods for robustness
#' - Before period: Monthly data from April-November, 2019-2024
#' - After period: Monthly data from April-November, 2025
#' - Unit of analysis: Month for increased statistical power
#'
#' Statistical approach:
#' - GLM with negative binomial distribution for overdispersed count data
#' - Welch's t-test on difference-in-differences (non-parametric alternative)
#' - BACI effect = (Impact_after - Impact_before) - (Control_after - Control_before)
#' - Tests the site x period interaction term
#' - Multiple control sites for robustness checking

library(tidyverse)
library(MASS)
library(lubridate)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

DATA_PATH <- here::here("data", "raw", "Traffic_Crashes_Resulting_in_Injury_20260104.csv")
OUTPUT_DIR <- here::here("outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define neighborhood groups
# Impact: Sunset area neighborhoods most directly affected by closure
IMPACT_NEIGHBORHOODS <- c(
  "Sunset/Parkside",
  "Inner Sunset",
  "Outer Sunset",
  "Golden Gate Park"
)

# Multiple control neighborhoods for robustness
# Selected based on: east/south location, residential character, adequate crash volume
CONTROL_NEIGHBORHOODS <- list(
  "Bernal Heights" = "Bernal Heights",
  "Outer Mission" = "Outer Mission",
  "Excelsior" = "Excelsior",
  "Oceanview/Merced/Ingleside" = "Oceanview/Merced/Ingleside",
  "Portola" = "Portola",
  "West of Twin Peaks" = "West of Twin Peaks",
  "Bayview Hunters Point" = "Bayview Hunters Point",
  "Pacific Heights" = "Pacific Heights"
)

# Time periods - using full months for cleaner analysis
ANALYSIS_MONTHS <- 4:11  # April through November

# Years to include
BEFORE_YEARS <- c(2019, 2020, 2021, 2022, 2023, 2024)
AFTER_YEARS <- c(2025)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

load_and_prepare_data <- function() {
  df <- read_csv(DATA_PATH, show_col_types = FALSE)

  df <- df %>%
    mutate(
      collision_date = parse_date_time(collision_date, orders = "Y B d"),
      year = year(collision_date),
      month = month(collision_date),
      day = day(collision_date)
    )

  return(df)
}

get_days_in_month <- function(year, month) {
  days_in_month(ymd(paste(year, month, "01", sep = "-")))
}

create_monthly_baci_dataset <- function(df, control_name, control_neighborhoods) {
  all_years <- c(BEFORE_YEARS, AFTER_YEARS)

  records <- map_dfr(all_years, function(yr) {
    period <- ifelse(yr %in% AFTER_YEARS, "After", "Before")

    map_dfr(ANALYSIS_MONTHS, function(mo) {
      month_data <- df %>% filter(year == yr, month == mo)
      days <- get_days_in_month(yr, mo)

      # Impact counts
      impact_count <- month_data %>%
        filter(analysis_neighborhood %in% IMPACT_NEIGHBORHOODS) %>%
        nrow()

      # Control counts
      control_count <- month_data %>%
        filter(analysis_neighborhood %in% control_neighborhoods) %>%
        nrow()

      tibble(
        year = yr,
        month = mo,
        year_month = sprintf("%d-%02d", yr, mo),
        period = period,
        site = c("Impact", "Control"),
        crashes = c(impact_count, control_count),
        days = days
      )
    })
  })

  return(records)
}

create_pooled_control_dataset <- function(df) {
  all_control_neighborhoods <- unlist(CONTROL_NEIGHBORHOODS)
  all_years <- c(BEFORE_YEARS, AFTER_YEARS)

  records <- map_dfr(all_years, function(yr) {
    period <- ifelse(yr %in% AFTER_YEARS, "After", "Before")

    map_dfr(ANALYSIS_MONTHS, function(mo) {
      month_data <- df %>% filter(year == yr, month == mo)
      days <- get_days_in_month(yr, mo)

      # Impact counts
      impact_count <- month_data %>%
        filter(analysis_neighborhood %in% IMPACT_NEIGHBORHOODS) %>%
        nrow()

      # Pooled control counts
      control_count <- month_data %>%
        filter(analysis_neighborhood %in% all_control_neighborhoods) %>%
        nrow()

      tibble(
        year = yr,
        month = mo,
        year_month = sprintf("%d-%02d", yr, mo),
        period = period,
        site = c("Impact", "Control"),
        crashes = c(impact_count, control_count),
        days = days
      )
    })
  })

  return(records)
}

# =============================================================================
# BACI STATISTICAL ANALYSIS
# =============================================================================

calculate_simple_baci_effect <- function(baci_df) {
  impact_before <- baci_df %>%
    filter(site == "Impact", period == "Before") %>%
    pull(crashes) %>%
    mean()

  impact_after <- baci_df %>%
    filter(site == "Impact", period == "After") %>%
    pull(crashes) %>%
    mean()

  control_before <- baci_df %>%
    filter(site == "Control", period == "Before") %>%
    pull(crashes) %>%
    mean()

  control_after <- baci_df %>%
    filter(site == "Control", period == "After") %>%
    pull(crashes) %>%
    mean()

  impact_change <- impact_after - impact_before
  control_change <- control_after - control_before
  baci_effect <- impact_change - control_change

  impact_pct_change <- ifelse(impact_before > 0,
                               ((impact_after - impact_before) / impact_before) * 100,
                               0)
  control_pct_change <- ifelse(control_before > 0,
                                ((control_after - control_before) / control_before) * 100,
                                0)
  baci_effect_pct <- impact_pct_change - control_pct_change

  list(
    impact_before = impact_before,
    impact_after = impact_after,
    control_before = control_before,
    control_after = control_after,
    impact_change = impact_change,
    control_change = control_change,
    impact_pct_change = impact_pct_change,
    control_pct_change = control_pct_change,
    baci_effect = baci_effect,
    baci_effect_pct = baci_effect_pct
  )
}

fit_baci_glm <- function(baci_df) {
  df_model <- baci_df %>%
    mutate(
      is_impact = as.integer(site == "Impact"),
      is_after = as.integer(period == "After"),
      interaction = is_impact * is_after,
      log_days = log(days)
    )

  # Try negative binomial first, fall back to Poisson
  model <- tryCatch({
    glm.nb(crashes ~ is_impact + is_after + interaction + offset(log_days),
           data = df_model)
  }, error = function(e) {
    glm(crashes ~ is_impact + is_after + interaction + offset(log_days),
        data = df_model, family = poisson())
  })

  model_type <- if (inherits(model, "negbin")) "NegBin" else "Poisson"

  list(model = model, model_type = model_type)
}

welch_ttest_baci <- function(baci_df) {
  # Pivot to get Impact and Control side by side for each time point
  pivot <- baci_df %>%
    pivot_wider(
      id_cols = c(year, month, year_month, period),
      names_from = site,
      values_from = crashes
    ) %>%
    mutate(diff = Impact - Control)

  # Separate before and after periods
  diffs_before <- pivot %>% filter(period == "Before") %>% pull(diff)
  diffs_after <- pivot %>% filter(period == "After") %>% pull(diff)

  # Welch's t-test (unequal variances)
  ttest_result <- t.test(diffs_after, diffs_before, var.equal = FALSE)

  # Calculate effect size (Cohen's d)
  mean_before <- mean(diffs_before)
  mean_after <- mean(diffs_after)

  n_before <- length(diffs_before)
  n_after <- length(diffs_after)
  var_before <- var(diffs_before)
  var_after <- var(diffs_after)

  # Pooled standard deviation for Cohen's d
  pooled_sd <- sqrt(((n_before - 1) * var_before + (n_after - 1) * var_after) /
                      (n_before + n_after - 2))

  cohens_d <- if (pooled_sd > 0) (mean_after - mean_before) / pooled_sd else 0

  list(
    t_stat = ttest_result$statistic,
    p_value = ttest_result$p.value,
    mean_diff_before = mean_before,
    mean_diff_after = mean_after,
    baci_effect = mean_after - mean_before,
    cohens_d = cohens_d,
    ci_low = ttest_result$conf.int[1],
    ci_high = ttest_result$conf.int[2],
    n_before = n_before,
    n_after = n_after,
    df = ttest_result$parameter
  )
}

analyze_single_control <- function(df, control_name, control_neighborhoods) {
  baci_df <- create_monthly_baci_dataset(df, control_name, control_neighborhoods)
  effect <- calculate_simple_baci_effect(baci_df)
  glm_result <- fit_baci_glm(baci_df)
  ttest <- welch_ttest_baci(baci_df)

  model <- glm_result$model
  model_summary <- summary(model)
  coef_table <- coef(model_summary)

  irr <- exp(coef(model)["interaction"])
  p_value <- coef_table["interaction", "Pr(>|z|)"]
  ci <- confint(model, "interaction", level = 0.95)

  list(
    control_name = control_name,
    control_before = effect$control_before,
    control_after = effect$control_after,
    control_pct_change = effect$control_pct_change,
    baci_effect = effect$baci_effect,
    baci_effect_pct = effect$baci_effect_pct,
    # GLM results
    irr = irr,
    p_value_glm = p_value,
    ci_low = ci[1],
    ci_high = ci[2],
    model_type = glm_result$model_type,
    # Welch's t-test results
    t_stat = ttest$t_stat,
    p_value_ttest = ttest$p_value,
    cohens_d = ttest$cohens_d,
    ttest_ci_low = ttest$ci_low,
    ttest_ci_high = ttest$ci_high
  )
}

# =============================================================================
# VISUALIZATION
# =============================================================================

plot_all_controls_comparison <- function(results_df, impact_before, impact_after, output_path) {
  # Sort by BACI effect
  results_sorted <- results_df %>% arrange(baci_effect)

  # Plot 1: BACI effects
  p1 <- ggplot(results_sorted, aes(x = baci_effect,
                                    y = fct_reorder(control_name, baci_effect),
                                    fill = p_value_glm < 0.05)) +
    geom_col(alpha = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#1f77b4"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05")) +
    labs(x = "BACI Effect (crashes/month)",
         y = NULL,
         title = "BACI Effect by Control Neighborhood",
         fill = "Significance") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Plot 2: Before-After changes
  impact_change <- ((impact_after - impact_before) / impact_before) * 100

  changes_df <- bind_rows(
    tibble(name = "IMPACT (Sunset)", pct_change = impact_change, type = "Impact"),
    results_sorted %>%
      transmute(name = control_name, pct_change = control_pct_change, type = "Control")
  ) %>%
    mutate(name = fct_inorder(name))

  p2 <- ggplot(changes_df, aes(x = pct_change, y = name, fill = type)) +
    geom_col(alpha = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    scale_fill_manual(values = c("Impact" = "#d62728", "Control" = "#1f77b4")) +
    labs(x = "% Change (Before -> After)",
         y = NULL,
         title = "Percent Change by Site",
         fill = "Site Type") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Combine plots
  combined <- cowplot::plot_grid(p1, p2, ncol = 2)
  ggsave(output_path, combined, width = 14, height = 6, dpi = 150)

  output_path
}

plot_monthly_trends_multi <- function(df, output_path) {
  # Create dataset with ALL months (1-12), not just analysis months
  all_control_neighborhoods <- unlist(CONTROL_NEIGHBORHOODS)
  all_years <- c(BEFORE_YEARS, AFTER_YEARS)
  all_months <- 1:12

  baci_df <- map_dfr(all_years, function(yr) {
    period <- ifelse(yr %in% AFTER_YEARS, "After", "Before")

    map_dfr(all_months, function(mo) {
      month_data <- df %>% filter(year == yr, month == mo)
      days <- get_days_in_month(yr, mo)

      # Impact counts
      impact_count <- month_data %>%
        filter(analysis_neighborhood %in% IMPACT_NEIGHBORHOODS) %>%
        nrow()

      # Pooled control counts
      control_count <- month_data %>%
        filter(analysis_neighborhood %in% all_control_neighborhoods) %>%
        nrow()

      tibble(
        year = yr,
        month = mo,
        year_month = sprintf("%d-%02d", yr, mo),
        period = period,
        site = c("Impact", "Control"),
        crashes = c(impact_count, control_count),
        days = days
      )
    })
  })

  # Ensure proper ordering
  baci_df <- baci_df %>%
    filter(year_month < "2025-12") |> 
    arrange(year_month) %>%
    mutate(time_idx = as.integer(factor(year_month)))

  # Find intervention point (March 2025 = 2025-03)
  n_before <- baci_df %>%
    filter(site == "Impact", year_month < "2025-03") %>%
    nrow()

  p <- ggplot(baci_df, aes(x = time_idx, y = crashes, color = site)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    geom_vline(xintercept = n_before + 0.5, linetype = "dashed",
               color = "gray40", linewidth = 1) +
    scale_color_manual(values = c("Impact" = "#555599", "Control" = "#66BBBB"),
                       labels = c("Impact" = "Sunset Neighborhood",
                                  "Control" = "Pooled Controls (8 neighborhoods)")) +
    annotate("text", x = n_before + 0.5, y = max(baci_df$crashes),
             label = " UGH Closure  \n(March 2025)", hjust = -0.1, vjust = 1, color = "gray40") +
    labs(x = "", 
         y = "Crashes per Month",
         title = "Monthly Traffic Crashes for Select San Francisco Neighborhoods",
         subtitle = "Traffic Crashes Resulting in Injury Before and After the Closure of the Upper Great Highway",
         color = NULL) +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16), 
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          plot.title = element_text(face = "bold", size = 20),
          plot.subtitle = element_text(color = "gray40", size = 16),
          panel.grid.minor = element_blank()
        )

  # Custom x-axis labels (every 12th point = yearly)
  year_months <- baci_df %>%
    filter(site == "Impact") %>%
    arrange(time_idx) %>%
    pull(year_month)

  tick_positions <- seq(1, length(year_months), by = 12)
  tick_labels <- year_months[tick_positions]

  p <- p + scale_x_continuous(breaks = tick_positions, labels = tick_labels)

  ggsave(output_path, p, width = 14, height = 6, dpi = 300)
  output_path
}

plot_forest <- function(results_df, output_path) {
  results_sorted <- results_df %>%
    arrange(irr) %>%
    mutate(control_name = fct_inorder(control_name))

  p <- ggplot(results_sorted, aes(y = control_name)) +
    geom_errorbarh(aes(xmin = exp(ci_low), xmax = exp(ci_high),
                       color = p_value_glm < 0.05),
                   height = 0.3, linewidth = 1, alpha = 0.7) +
    geom_point(aes(x = irr, color = p_value_glm < 0.05), size = 4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#1f77b4"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05")) +
    scale_x_log10(limits = c(0.3, 3)) +
    labs(x = "Incidence Rate Ratio (IRR)",
         y = NULL,
         title = "Forest Plot: BACI Effect (IRR) by Control Neighborhood",
         subtitle = "IRR = 1 means no effect; green = p < 0.05",
         color = "Significance") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(output_path, p, width = 10, height = 8, dpi = 150)
  output_path
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

main <- function() {
  cat(strrep("=", 80), "\n")
  cat("BACI ANALYSIS: SUNSET DUNES CLOSURE IMPACT ON TRAFFIC CRASHES\n")
  cat("Multi-Control Robustness Analysis\n")
  cat(strrep("=", 80), "\n")

  # Load data
  cat("\nLoading crash data...\n")
  df <- load_and_prepare_data()
  cat(sprintf("Total records: %s\n", format(nrow(df), big.mark = ",")))
  cat(sprintf("Date range: %s to %s\n",
              min(df$collision_date, na.rm = TRUE),
              max(df$collision_date, na.rm = TRUE)))

  # Get impact area stats
  cat("\n", strrep("-", 80), "\n", sep = "")
  cat("IMPACT AREA: Sunset/Parkside, Inner Sunset, Golden Gate Park\n")
  cat(strrep("-", 80), "\n")

  impact_df <- create_monthly_baci_dataset(df, "dummy", character(0))
  impact_before <- impact_df %>%
    filter(site == "Impact", period == "Before") %>%
    pull(crashes) %>%
    mean()
  impact_after <- impact_df %>%
    filter(site == "Impact", period == "After") %>%
    pull(crashes) %>%
    mean()
  impact_pct <- ((impact_after - impact_before) / impact_before) * 100

  cat(sprintf("  Before (48 months): %.1f crashes/month\n", impact_before))
  cat(sprintf("  After (8 months):   %.1f crashes/month\n", impact_after))
  cat(sprintf("  Change:             %+.1f%%\n", impact_pct))

  # Analyze each control neighborhood
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("INDIVIDUAL CONTROL ANALYSIS\n")
  cat(strrep("=", 80), "\n")

  results <- map(names(CONTROL_NEIGHBORHOODS), function(control_name) {
    control_neighborhoods <- CONTROL_NEIGHBORHOODS[[control_name]]
    result <- analyze_single_control(df, control_name, control_neighborhoods)

    cat(sprintf("\n%s:\n", control_name))
    cat(sprintf("  Before: %.1f/mo, After: %.1f/mo, Change: %+.1f%%\n",
                result$control_before, result$control_after, result$control_pct_change))
    cat(sprintf("  BACI Effect: %+.2f/mo (%+.1f pp)\n",
                result$baci_effect, result$baci_effect_pct))
    cat(sprintf("  GLM:   IRR=%.3f, p=%.4f\n", result$irr, result$p_value_glm))
    cat(sprintf("  t-test: t=%.3f, p=%.4f, Cohen's d=%.3f\n",
                result$t_stat, result$p_value_ttest, result$cohens_d))

    result
  })

  results_df <- bind_rows(results)

  # Pooled control analysis
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("POOLED CONTROL ANALYSIS (All 8 neighborhoods combined)\n")
  cat(strrep("=", 80), "\n")

  pooled_df <- create_pooled_control_dataset(df)
  pooled_effect <- calculate_simple_baci_effect(pooled_df)
  pooled_glm <- fit_baci_glm(pooled_df)
  pooled_ttest <- welch_ttest_baci(pooled_df)

  pooled_model <- pooled_glm$model
  pooled_irr <- exp(coef(pooled_model)["interaction"])
  pooled_p_glm <- summary(pooled_model)$coefficients["interaction", "Pr(>|z|)"]
  pooled_p_ttest <- pooled_ttest$p_value

  cat(sprintf("\n  Pooled Control Before: %.1f/mo\n", pooled_effect$control_before))
  cat(sprintf("  Pooled Control After:  %.1f/mo\n", pooled_effect$control_after))
  cat(sprintf("  Pooled Control Change: %+.1f%%\n", pooled_effect$control_pct_change))
  cat(sprintf("\n  BACI Effect: %+.2f crashes/month\n", pooled_effect$baci_effect))
  cat(sprintf("  BACI Effect: %+.1f percentage points\n", pooled_effect$baci_effect_pct))
  cat(sprintf("\n  GLM:    IRR=%.3f, p=%.4f\n", pooled_irr, pooled_p_glm))
  cat(sprintf("  t-test: t=%.3f, p=%.4f, Cohen's d=%.3f\n",
              pooled_ttest$t_stat, pooled_p_ttest, pooled_ttest$cohens_d))
  cat(sprintf("          95%% CI for BACI effect: [%.2f, %.2f]\n",
              pooled_ttest$ci_low, pooled_ttest$ci_high))

  # Summary table
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("SUMMARY TABLE: ALL CONTROLS\n")
  cat(strrep("=", 80), "\n\n")

  cat(sprintf("%-28s %7s %7s %7s %6s %8s %8s %6s\n",
              "Control", "Before", "After", "BACI", "IRR", "p(GLM)", "p(t)", "d"))
  cat(strrep("-", 96), "\n")

  results_df %>%
    arrange(baci_effect) %>%
    pwalk(function(...) {
      row <- list(...)
      sig_glm <- if (row$p_value_glm < 0.05) "*" else ""
      sig_t <- if (row$p_value_ttest < 0.05) "*" else ""
      cat(sprintf("%-28s %7.1f %7.1f %+6.2f %6.2f %7.3f%s %7.3f%s %+5.2f\n",
                  row$control_name, row$control_before, row$control_after,
                  row$baci_effect, row$irr, row$p_value_glm, sig_glm,
                  row$p_value_ttest, sig_t, row$cohens_d))
    })

  cat(strrep("-", 96), "\n")
  cat(sprintf("%-28s %7.1f %7.1f %+6.2f %6.2f %7.3f  %7.3f  %+5.2f\n",
              "POOLED (all 8)", pooled_effect$control_before, pooled_effect$control_after,
              pooled_effect$baci_effect, pooled_irr, pooled_p_glm,
              pooled_p_ttest, pooled_ttest$cohens_d))
  cat(sprintf("%-28s %7.1f %7.1f %+5.1f%%\n",
              "IMPACT (Sunset area)", impact_before, impact_after, impact_pct))

  # Statistics across controls
  cat("\n", strrep("-", 80), "\n", sep = "")
  cat("ROBUSTNESS STATISTICS\n")
  cat(strrep("-", 80), "\n")

  mean_baci <- mean(results_df$baci_effect)
  std_baci <- sd(results_df$baci_effect)
  median_baci <- median(results_df$baci_effect)
  mean_irr <- mean(results_df$irr)
  median_irr <- median(results_df$irr)
  mean_cohens_d <- mean(results_df$cohens_d)
  n_sig_glm <- sum(results_df$p_value_glm < 0.05)
  n_sig_ttest <- sum(results_df$p_value_ttest < 0.05)
  n_positive <- sum(results_df$baci_effect > 0)
  n_negative <- sum(results_df$baci_effect < 0)

  cat(sprintf("\n  Number of control neighborhoods: %d\n", nrow(results_df)))
  cat(sprintf("  Mean BACI effect:   %+.2f crashes/month (SD: %.2f)\n", mean_baci, std_baci))
  cat(sprintf("  Median BACI effect: %+.2f crashes/month\n", median_baci))
  cat(sprintf("  Mean IRR:           %.3f\n", mean_irr))
  cat(sprintf("  Median IRR:         %.3f\n", median_irr))
  cat(sprintf("  Mean Cohen's d:     %+.3f\n", mean_cohens_d))
  cat(sprintf("\n  Controls showing increase (BACI > 0): %d/%d\n", n_positive, nrow(results_df)))
  cat(sprintf("  Controls showing decrease (BACI < 0): %d/%d\n", n_negative, nrow(results_df)))
  cat(sprintf("  Significant by GLM (p < 0.05):    %d/%d\n", n_sig_glm, nrow(results_df)))
  cat(sprintf("  Significant by t-test (p < 0.05): %d/%d\n", n_sig_ttest, nrow(results_df)))

  # Generate visualizations
  cat("\n", strrep("-", 80), "\n", sep = "")
  cat("GENERATING VISUALIZATIONS\n")
  cat(strrep("-", 80), "\n")

  comparison_path <- file.path(OUTPUT_DIR, "baci_multi_control_comparison_r.png")
  plot_all_controls_comparison(results_df, impact_before, impact_after, comparison_path)
  cat(sprintf("  Saved: %s\n", comparison_path))

  trends_path <- file.path(OUTPUT_DIR, "baci_monthly_trends_pooled_r.png")
  plot_monthly_trends_multi(df, trends_path)
  cat(sprintf("  Saved: %s\n", trends_path))

  forest_path <- file.path(OUTPUT_DIR, "baci_forest_plot_r.png")
  plot_forest(results_df, forest_path)
  cat(sprintf("  Saved: %s\n", forest_path))

  # Interpretation
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("INTERPRETATION\n")
  cat(strrep("=", 80), "\n")

  cat("
BACI Design Summary:
- Impact area: Sunset/Parkside, Inner Sunset, Golden Gate Park
- Control areas: 8 SF neighborhoods
- Before period: April-November 2019-2024 (48 months)
- After period: April-November 2025 (8 months)

Statistical Methods:
- Negative Binomial GLM with site x period interaction (tests IRR)
- Welch's t-test on difference-in-differences (tests mean change)

Key Findings:
")

  if (n_sig_glm == 0 && n_sig_ttest == 0) {
    cat("- NO SIGNIFICANT EFFECTS: Neither GLM nor t-test found significant\n")
    cat("  BACI effects in any of the 8 control comparisons (all p > 0.05).\n")
  } else if (n_sig_glm == nrow(results_df) || n_sig_ttest == nrow(results_df)) {
    cat(sprintf("- CONSISTENT SIGNIFICANT EFFECTS: Multiple control comparisons\n"))
    cat(sprintf("  show significant effects (GLM: %d/8, t-test: %d/8).\n", n_sig_glm, n_sig_ttest))
  } else {
    cat(sprintf("- GLM significant: %d/%d controls\n", n_sig_glm, nrow(results_df)))
    cat(sprintf("- t-test significant: %d/%d controls\n", n_sig_ttest, nrow(results_df)))
  }

  if (n_positive > n_negative) {
    cat(sprintf("\n- DIRECTION: %d/%d controls suggest the closure\n", n_positive, nrow(results_df)))
    cat("  was associated with INCREASED crashes in the impact area.\n")
  } else if (n_negative > n_positive) {
    cat(sprintf("\n- DIRECTION: %d/%d controls suggest the closure\n", n_negative, nrow(results_df)))
    cat("  was associated with DECREASED crashes in the impact area.\n")
  } else {
    cat(sprintf("\n- DIRECTION: Results are evenly split (%d increase vs %d decrease),\n",
                n_positive, n_negative))
    cat("  suggesting no consistent directional effect.\n")
  }

  cat(sprintf("
Effect Size:
- Mean BACI effect across controls: %+.2f crashes/month
- Mean IRR (GLM): %.3f (%+.1f%% relative change)
- Mean Cohen's d (t-test): %+.3f

Pooled Analysis (all 8 controls combined):
- GLM:    IRR = %.3f, p = %.4f
- t-test: d = %+.3f, p = %.4f
- 95%% CI for BACI effect: [%.2f, %.2f] crashes/month

Robustness:
- Results are %s across control neighborhoods (SD = %.2f)
- GLM and t-test %s on significance

Conclusion:
",
              mean_baci,
              mean_irr, (mean_irr - 1) * 100,
              mean_cohens_d,
              pooled_irr, pooled_p_glm,
              pooled_ttest$cohens_d, pooled_p_ttest,
              pooled_ttest$ci_low, pooled_ttest$ci_high,
              if (std_baci < 2) "consistent" else "variable", std_baci,
              if ((pooled_p_glm < 0.05) == (pooled_p_ttest < 0.05)) "agree" else "disagree"))

  both_nonsig <- pooled_p_glm >= 0.05 && pooled_p_ttest >= 0.05
  both_sig <- pooled_p_glm < 0.05 && pooled_p_ttest < 0.05

  if (both_nonsig && n_sig_glm == 0 && n_sig_ttest == 0) {
    cat("The Sunset Dunes closure does not appear to have caused a statistically\n")
    cat("detectable change in traffic crash rates in the Sunset area. Both GLM and\n")
    cat("Welch's t-test found no significant effects across 8 control neighborhoods.\n")
  } else if (both_sig) {
    direction <- if (pooled_irr > 1) "increase" else "decrease"
    cat(sprintf("Both statistical tests suggest the closure may be associated with a %s\n", direction))
    cat(sprintf("in crash rates (GLM: IRR=%.2f, p=%.4f;\n", pooled_irr, pooled_p_glm))
    cat(sprintf("t-test: p=%.4f).\n", pooled_p_ttest))
  } else {
    cat("Results are inconclusive. The two statistical methods show similar findings,\n")
    cat("but neither reaches conventional significance thresholds. The effect size\n")
    cat("is small and within the range of normal variation.\n")
  }

  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("ANALYSIS COMPLETE\n")
  cat(strrep("=", 80), "\n")

  invisible(list(results_df = results_df, pooled_effect = pooled_effect))
}

# Run main analysis
main()
