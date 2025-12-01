#### --- Steller Sea Lion Constraint Line Analysis --- ####
setwd("C:/Users/kpkolda/Desktop/columbus_analysis/Steller Sea Lions")
library(tidyverse)
library(readxl)

df_1 <- read_excel("tbl_All_Dives_with_Shape_with_Bout_with_Trip.xlsx")

start_age <- df_1 %>%
  group_by(ref) %>%
  summarize(start_age_months = min(AGE_MONTHS, na.rm = TRUE))

df_2 <- df_1 %>% 
  group_by(ref) %>% 
  mutate(start_age_months = min(AGE_MONTHS, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(cADL = case_when(
    start_age_months <= 10 ~ 108,
    start_age_months > 10 & start_age_months <= 17 ~ 138,
    start_age_months > 17 & start_age_months <= 22 ~ 162,
    start_age_months > 22 ~ 198,
    TRUE ~ NA_real_
  ))

####Moving Window####
per_sealion <- df_2 %>%
  group_by(ref, cADL) %>%
  reframe(
    max_dur_ind = max(DIVE_DUR, na.rm = TRUE),
    med_bADL_ind = if (all(DIVE_DUR <= cADL | is.na(DIVE_DUR))) NA_real_
    else median(DIVE_DUR[DIVE_DUR > cADL], na.rm = TRUE),
    diff_cADL_med_bADL = ifelse(is.na(med_bADL_ind), NA_real_,
                                med_bADL_ind - cADL),
    pct_beyond_cADL = sum(DIVE_DUR > cADL, na.rm = TRUE) / sum(!is.na(DIVE_DUR)) * 100,
    count_beyond_cADL = sum(DIVE_DUR > cADL, na.rm = TRUE)
  )

summary_groups <- per_sealion %>%
  group_by(cADL) %>%
  reframe(
    avg_max = mean(max_dur_ind, na.rm = TRUE),
    avg_median_bADL = mean(med_bADL_ind, na.rm = TRUE),
    moving_window = (avg_max - cADL) / (avg_median_bADL - cADL),
    n_animals = sum(!is.na(max_dur_ind)),
    # New stats
    avg_diff_cADL_med_bADL = mean(diff_cADL_med_bADL, na.rm = TRUE),
    sd_diff_cADL_med_bADL = sd(diff_cADL_med_bADL, na.rm = TRUE),
    avg_pct_beyond_cADL = mean(pct_beyond_cADL, na.rm = TRUE),
    sd_pct_beyond_cADL = sd(pct_beyond_cADL, na.rm = TRUE)
  )

summary_groups_unique <- summary_groups %>%
  group_by(cADL) %>%
  slice(1) %>%
  ungroup()

df_2 <- df_2 %>%
  left_join(
    summary_groups_unique %>% select(cADL, moving_window),
    by = "cADL"
  )


#round the moving window
df_2 <- df_2 %>% mutate(moving_window_rounded = round(moving_window, 0))

####Apply Moving Window####
#determine which dives are consecutive
df_2 <- df_2 %>% 
  mutate(DIVE_END = DIVE_START_DATE + DIVE_DUR)

df_2 <- df_2 %>%
  mutate(SURF_DUR = ifelse(is.na(SURF_DUR), 0, SURF_DUR)) %>%  #replace NAs with 0
  arrange(ref, DIVE_START_DATE) %>%
  group_by(ref) %>%
  mutate(
    PDI_END = DIVE_END + SURF_DUR,
    next_start = lead(DIVE_START_DATE),
    consecutive = !is.na(next_start) & !is.na(PDI_END) &
      (abs(as.numeric(difftime(next_start, PDI_END, units="secs"))) < 1),
    seq_group = cumsum(!consecutive | is.na(consecutive)) + 1
  ) %>%
  ungroup()


#moving window
library(RcppRoll)

df_2 <- df_2 %>%
  group_by(ref, seq_group) %>%
  mutate(
    #Only compute rolling sums if sequence is long enough
    moving_sum_dur = if(n() >= first(moving_window_rounded)) {
      roll_sum(
        x = DIVE_DUR,
        n = first(moving_window_rounded),
        align = "right",
        fill = NA
      )
    } else {
      rep(NA_real_, n())
    },
    moving_sum_surf = if(n() >= first(moving_window_rounded)) {
      roll_sum(
        x = SURF_DUR,
        n = first(moving_window_rounded),
        align = "right",
        fill = NA
      )
    } else {
      rep(NA_real_, n())
    }
  ) %>%
  ungroup()

summary_total <- df_2 %>%
  summarize(
    total_dives = n(),
    consecutive_dives = sum(consecutive, na.rm = TRUE),
    dives_in_moving_window = sum(!is.na(moving_sum_dur)),
  )

summary_total_animal <- df_2 %>%
  group_by(ref) %>%
  summarize(
    total_dives = n(),
    consecutive_dives = sum(consecutive, na.rm = TRUE),
    dives_in_moving_window = sum(!is.na(moving_sum_dur)),
    .groups = "drop"
  )

####Check Plots####
animals <- unique(df_2$ref)

for(a in animals) {
  
  #Filter data for this animal
  df_animal <- df_2[df_2$ref == a & !is.na(df_2$moving_sum_dur) & !is.na(df_2$moving_sum_surf), ]
  
  mw <- unique(na.omit(round(df_animal$moving_window_rounded)))
  
  #Skip if no data
  if(nrow(df_animal) == 0) next
  
  filename <- paste0(a, "_plot.png")
  
  #Open PNG device
  png(filename = filename, width = 800, height = 600)
  
  #Plot
  plot(df_animal$moving_sum_dur, df_animal$moving_sum_surf,
       main = paste(a),
       xlab = paste0(mw, " Successive Dive durations (s)"),
       ylab = paste0(mw, " Successive Post Dive Intervals (s)"),
       pch = 19, col = "black")
  
  #save PNG
  dev.off()
}

####Constraint Line Analysis####
library(quantreg)

animals_to_save <- c("st_Cliff_02", "st_Frazer_02", "st_Lindy_01")

#Save csv files from animals
for(a in animals_to_save) {
  df_animal <- df_2 %>% filter(ref == a)
  
  #Check if there is any data
  if(nrow(df_animal) == 0) next
  
  #Save as CSV
  write_csv(df_animal, paste0(a, "_mw.csv"))
}

prepare_data <- function(file) {
  df <- read_csv(file, show_col_types = FALSE) %>%
    rename(
      x = moving_sum_dur,
      y = moving_sum_surf
    ) %>%
    filter(!is.na(x) & !is.na(y))
  
  return(df)
}

compute_error_sum <- function(dataframe, tau) {
  if(nrow(dataframe) < 2) {
    return(list(tau = tau, k = 0, error_sum = NA_real_, slope = NA_real_, intercept = NA_real_))
  }
  
  mod1 <- tryCatch(
    rq(y ~ x, data = dataframe, tau = tau),
    error = function(e) return(NULL)
  )
  
  if(is.null(mod1)) return(list(tau = tau, k = 0, error_sum = NA_real_, slope = NA_real_, intercept = NA_real_))
  
  coefs <- coef(mod1)
  slope <- coefs[2]
  intercept <- coefs[1]
  
  dataframe$predicted_y <- intercept + slope * dataframe$x
  dataframe$residuals <- dataframe$y - dataframe$predicted_y
  below_line <- dataframe %>% filter(y < predicted_y)
  k <- nrow(below_line)
  
  if(k == 0) error_sum <- NA_real_ else error_sum <- sqrt(sum(below_line$residuals^2) / k)
  
  return(list(tau = tau, k = k, error_sum = error_sum, slope = slope, intercept = intercept))
}
#Find edge
find_edge <- function(df, min_k = 5, max_k = 100) {
  results <- data.frame(tau = numeric(),
                        k = numeric(),
                        error_sum = numeric(),
                        slope = numeric(),
                        intercept = numeric())
  
  tau <- 0.5
  while(tau >= 0.001) {
    res <- compute_error_sum(df, tau)
    if(!is.na(res$error_sum) && res$k >= min_k && res$k <= max_k) {
      results <- rbind(results, data.frame(
        tau = res$tau,
        k = res$k,
        error_sum = res$error_sum,
        slope = res$slope,
        intercept = res$intercept
      ))
    }
    tau <- tau - 0.001
  }
  
  #Best fit: either the one with lowest error_sum or a default if empty
  if(nrow(results) == 0) {
    return(data.frame(tau = NA_real_, error_sum = NA_real_, slope = NA_real_, intercept = NA_real_))
  } else {
    best <- results %>% arrange(error_sum) %>% slice(1)
    return(best)
  }
}

#Find Edges
find_edges <- function(file, split = 0.5, n_split = 100) {
  df <- prepare_data(file)
  if(nrow(df) == 0) return(data.frame())  # skip empty files
  
  increment <- (0.95 - split) / n_split
  all_results <- data.frame()
  
  split_val <- split
  while(split_val <= 0.95) {
    threshold <- tryCatch(
      quantile(df$x, split_val, na.rm = TRUE),
      error = function(e) NA_real_
    )
    if(is.na(threshold)) {
      split_val <- split_val + increment
      next
    }
    
    mat_aer <- df[df$x <= threshold, ]
    mat_anaer <- df[df$x > threshold, ]
    
    result_aer <- find_edge(mat_aer)
    result_anaer <- find_edge(mat_anaer)
    
    all_results <- rbind(all_results, data.frame(
      split = split_val,
      best_tau_aer = result_aer$tau,
      error_sum_aer = result_aer$error_sum,
      slope_aer = result_aer$slope,
      intercept_aer = result_aer$intercept,
      best_tau_anaer = result_anaer$tau,
      error_sum_anaer = result_anaer$error_sum,
      slope_anaer = result_anaer$slope,
      intercept_anaer = result_anaer$intercept
    ))
    
    split_val <- split_val + increment
  }
  
  return(all_results)
}

#Angle between slopes
angle_between_slopes <- function(m1, m2) {
  atan((m2 - m1) / (1 + m1 * m2)) * 180 / pi
}

#Calculate line intersection
x_intersection <- function(m1, b1, m2, b2){
  if(m1 == m2) {
    return(NA)
  }
  x <- (b2 - b1) / (m1 - m2)
  return(x)
}

y_intersection <- function(m1, b1, m2, b2, x){
  if(m1 == m2) {
    return(NA)
  }
  y <- m1 * x + b1
  return(y)
}

file_list <- list.files(pattern = "*_mw.csv")
combined_results <- data.frame()

for (file in file_list) {
  file_basename <- sub("_mw.csv$", "", basename(file))
  
  results <- find_edges(file, .5, 100)
  
  if(nrow(results) == 0) {
    best_fit <- data.frame(
      split = NA_real_,
      best_tau_aer = NA_real_,
      error_sum_aer = NA_real_,
      slope_aer = NA_real_,
      intercept_aer = NA_real_,
      best_tau_anaer = NA_real_,
      error_sum_anaer = NA_real_,
      slope_anaer = NA_real_,
      intercept_anaer = NA_real_,
      error_sum_cum = NA_real_
    )
  } else {
    results$error_sum_cum <- results$error_sum_aer + results$error_sum_anaer
    best_fit <- results %>% arrange(error_sum_cum) %>% slice(1)
  }
  
  best_fit$animal_id <- file_basename
  
  combined_results <- rbind(combined_results, best_fit)
}

#Add angle difference
combined_results$angle_diff <- mapply(
  angle_between_slopes,
  combined_results$slope_aer,
  combined_results$slope_anaer)

combined_results$x_intersect <- mapply(
  x_intersection,
  combined_results$slope_aer,
  combined_results$intercept_aer,
  combined_results$slope_anaer,
  combined_results$intercept_anaer)

combined_results$y_intersect <- mapply(
  y_intersection,
  combined_results$slope_aer,
  combined_results$intercept_aer,
  combined_results$slope_anaer,
  combined_results$intercept_anaer,
  combined_results$x_intersect)

combined_results$STI <- (combined_results$x_intersect / 5) / 60

####Continuity violations####

df_lindy <- df_2[df_2$ref == "st_Lindy_01", ]
df_cliff <- df_2[df_2$ref == "st_Cliff_02",]
df_frazer <- df_2[df_2$ref == "st_Frazer_02",]

analysis_results <- data.frame()


for (file in file_list){
  df <- prepare_data(file)
  file_basename <- sub("_mw.csv$", "", basename(file))
  
  if (file_basename %in% c("st_Cliff_02","st_Frazer_02")) {
    
    integration_factor <- 7
    cADL <- 108
    upper_limit <- 756
    
  } else {
    
    integration_factor <- 6
    cADL <- 138
    upper_limit <- 828
    
  }
  
  df_filter <- df[df$x < upper_limit, ]
  
  results <- find_edge(df_filter, 5, 100)
  results$animal_id <- file_basename
  
  png(filename = paste0(file_basename, "_final_plot.png"))
  plot(y ~ x, data = df,
       xlab = paste0(integration_factor, " Successive Dive Durations (s)"),
       ylab = paste0(integration_factor, " Successive Post Dive Intervals (s)"),
       main = file_basename)
  
  segments(x0 = min(df$x),
           y0 = results$intercept + results$slope * min(df$x),
           x1 = upper_limit,
           y1 = results$intercept + results$slope * upper_limit,
           col = "black", lwd = 2)
  
  if (!is.na(results$slope)) {
    segments(x0 = min(df$x, na.rm=TRUE),
             y0 = results$intercept + results$slope * min(df$x, na.rm=TRUE),
             x1 = upper_limit,
             y1 = results$intercept + results$slope * upper_limit,
             col = "black", lwd = 2)
  }
  
  dev.off()
}
