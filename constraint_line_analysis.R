### --- Habor Seal Analysis --- ###
setwd("C:/Users/kpkolda/Desktop/columbus_analysis/Harbor Seals")
library(tidyverse)
library(quantreg)
library(readxl)

####Filter and movign window summation####
df_capture <- read_csv("capture_data.csv")
df_capture <- df_capture[!is.na(df_capture$SerialNum_ID),]
colnames(df_capture) <- c("sn", "animalid", "mass", "sex", "age", "capture_location", "repro")

#filter dives > 4m and > 20
for (i in 1:nrow(df_capture)){
  sn <- df_capture$sn[i]
  file <- paste0("out-Archive_", sn, "_diveStats.csv")
  df <- read_csv(file)
  
  #filter for dives > 4m and > 20s
  df_filter <- df[df$descdist > 4 & df$divetim > 20,]
  
  #calculate dive times and post dive intervals
  df_calc <- df_filter %>% 
    mutate(
      end_time = begdesc + seconds(divetim),
      next_begdesc = lead(begdesc),
      pdi = as.numeric(difftime(next_begdesc, end_time, units = "secs"))
    )
  
  write.csv(df_calc, file = paste0(sn, "_filtered_divestats.csv"))
}


#Initialize columns with NA before the loop
df_capture$median_bADL <- NA_real_
df_capture$max_divetim <- NA_real_
df_capture$integration_factor <- NA_real_
df_capture$ndives_bADL <- NA_integer_
df_capture$percent_dives_beyond <- NA_real_

for (i in 1:nrow(df_capture)) {
  sn <- df_capture$sn[i]
  file <- paste0(sn, "_filtered_diveStats.csv")
  df <- read_csv(file)
  
  age <- df_capture$age[i]
  mass <- df_capture$mass[i]
  
  if (age == "yearling") {
    cADL <- 438
  } else {
    cADL <- 612
  }
  
  df_filtered <- df[df$divetim > cADL, ]
  
  median_bADL <- median(df_filtered$divetim, na.rm = TRUE)
  max_divetim <- max(df$divetim, na.rm = TRUE)
  ndives_bADL <- sum(df$divetim > cADL, na.rm = TRUE)
  total_dives <- nrow(df)
  
  #check if median_bADL is NA or if denominator is 0 before calculating integration_factor
  if (!is.na(median_bADL) && (median_bADL - cADL) != 0) {
    integration_factor <- (max_divetim - cADL) / (median_bADL - cADL)
  } else {
    integration_factor <- NA_real_
  }
  
  percent_dives_beyond <- ifelse(total_dives > 0, (ndives_bADL / total_dives) * 100, NA_real_)
  
  df_capture$median_bADL[i] <- median_bADL
  df_capture$max_divetim[i] <- max_divetim
  df_capture$integration_factor[i] <- integration_factor
  df_capture$ndives_bADL[i] <- ndives_bADL
  df_capture$percent_dives_beyond[i] <- percent_dives_beyond
}

if_adult <- median(df_capture$integration_factor[df_capture$age != "yearling"], na.rm = TRUE)
if_yearling <- median(df_capture$integration_factor[df_capture$age == "yearling"], na.rm = TRUE)

df_capture <- df_capture %>%
  mutate(if_group = case_when(
    age == "yearling" ~ if_yearling,
    TRUE              ~ if_adult
  ))


####Applying the moving sums####
library(tidyverse)
library(RcppRoll)

for (i in 1:nrow(df_capture)) {
  sn <- df_capture$sn[i]
  file <- paste0(sn, "_filtered_diveStats.csv")
  df <- read_csv(file)
  
  age <- df_capture$age[i]
  
  if (age == "yearling") {
    integration_factor <- if_yearling
  } else {
    integration_factor <- if_adult
  }
  
  df_mw <- df %>%  mutate(
    dd_mw = roll_sum(divetim, n = integration_factor, fill = NA, align = "left"),
    pdi_mw = roll_sum(postdive.dur, n = integration_factor, fill = NA, align = "left")
  )
  
  write_csv(df_mw, file = paste0(sn, "_mw.csv"))
}

#Plot
for (i in 1:nrow(df_capture)){
  sn <- df_capture$sn[i]
  file <- paste0(sn, "_mw.csv")
  df <- read_csv(file)
  
  if (age == "yearling") {
    integration_factor <- 13
  } else {
    integration_factor <- 11
  }
  
  plot(pdi_mw~dd_mw, data = df, ylim = c(0, 2000),
       xlab = paste0(integration_factor, " Successive Dive Durations (s)"),
       ylab = paste0(integration_factor, " Successive Post Dive Intervals (s)"),
       main = paste0(sn))
}

####Constraint Line analysis####
library(tidyverse)
library(quantreg)

### --- Functions --- ###

#Prepare data
prepare_data <- function(file) {
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  
  #Rename columns dd_mw -> x and pdi_mw -> y
  df <- df %>% rename(
    x = dd_mw,
    y = pdi_mw
  )
  
  return(df)
}

compute_error_sum <- function(dataframe, tau) {
  #Quantile regression
  mod1 <- rq(y ~ x, data = dataframe, tau = tau)
  coefficients <- coef(mod1)
  slope <- coefficients[2]
  intercept <- coefficients[1]
  
  #Find points below quantile regression line
  dataframe$predicted_y <- intercept + (slope * dataframe$x)
  dataframe$residuals <- dataframe$y - dataframe$predicted_y
  
  below_line <- dataframe %>% filter(y < predicted_y)
  k <- nrow(below_line)
  
  #Find error of points below regression line
  error_sum <- sqrt(sum(below_line$residuals^2) / k)
  
  return(list(tau = tau, 
              k = k, 
              error_sum = error_sum, 
              slope = slope, 
              intercept = intercept))
}

#Find edge
find_edge <- function(dataframe){
  
  #Prep results table to store quantile regression results
  results <- data.frame(tau = numeric(),
                        k = numeric(),
                        error_sum = numeric(),
                        slope = numeric(),
                        intercept = numeric())
  
  #Reduce tau from 0.5 by .001 increments
  tau <- 0.5
  while (tau >= 0.001) {
    result <- compute_error_sum(dataframe = dataframe, tau = tau)
    if (result$k >= 5 && result$k <= 40) {
      results <- rbind(results, data.frame(tau = result$tau,
                                           k = result$k,
                                           error_sum = result$error_sum,
                                           slope = result$slope,
                                           intercept = result$intercept))
    }
    tau <- tau - 0.001
  }
  
  best_tau <- results %>%
    filter(k >= 5 & k <= 40) %>%
    arrange(error_sum) %>%
    slice(1)
  
  if (nrow(best_tau) == 0) {
    best_tau <- data.frame(tau = NA, error_sum = NA)
  }
  
  return(data.frame(best_tau = best_tau$tau, 
                    error_sum = best_tau$error_sum,
                    slope = best_tau$slope,
                    intercept = best_tau$intercept))
}

#Find Edges
find_edges <- function(file, split, n_split){
  df <- prepare_data(file)
  increment <- (0.95 - split) / n_split
  
  #Create dataframe to store all results
  all_results <- data.frame(split = numeric(), 
                            best_tau_aer = numeric(), 
                            error_sum_aer = numeric(),
                            slope_aer = numeric(),
                            intercept_aer = numeric(),
                            best_tau_anaer = numeric(), 
                            error_sum_anaer = numeric(),
                            slope_anaer = numeric(),
                            intercept_anaer = numeric())
  
  #Split segments and perform quantreg on both splits
  split <- 0.5
  while (split <= 0.95) {
    threshold <- quantile(df$x, split)
    mat_aer <- df[df$x <= threshold, ]
    mat_anaer <- df[df$x > threshold, ]
    
    result_aer <- find_edge(mat_aer)
    result_anaer <- find_edge(mat_anaer)
    
    if (any(is.na(result_aer$best_tau), is.na(result_anaer$best_tau))) {
      split <- split + increment
      next
    }
    
    all_results <- rbind(all_results, data.frame(
      split = split,
      best_tau_aer = result_aer$best_tau,
      error_sum_aer = result_aer$error_sum,
      slope_aer = result_aer$slope,
      intercept_aer = result_aer$intercept,
      best_tau_anaer = result_anaer$best_tau,
      error_sum_anaer = result_anaer$error_sum,
      slope_anaer = result_anaer$slope,
      intercept_anaer = result_anaer$intercept))
    
    split <- split + increment
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

### --- Main Loop --- ###
file_list <- list.files(pattern = "*_mw.csv")
combined_results <- data.frame()

for (file in file_list) {
  file_basename <- sub("_mw.csv$", "", basename(file))
  results <- find_edges(file, .5, 100)
  
  results$error_sum_cum <- results$error_sum_aer + results$error_sum_anaer
  best_fit <- results %>% arrange(error_sum_cum) %>% slice(1)
  
  df <- prepare_data(file)
  
  threshold <- quantile(df$x, best_fit$split)
  mat_aer <- df[df$x <= threshold, ]
  mat_anaer <- df[df$x > threshold, ]
  
  
  aer_values <- list(
    slope = best_fit$slope_aer,
    intercept = best_fit$intercept_aer
  )
  
  anaer_values <- list(
    slope = best_fit$slope_anaer,
    intercept = best_fit$intercept_anaer
  )
  
  best_fit$animal_id <- file_basename
  
  integration_factor <- df_capture$integration_factor[df_capture$sn == file_basename]
  
  
  png(filename = paste0(file_basename, "_prelim_plot.png"))
  plot(y ~ x, data = df, ylim = c(0, 2000),
       xlab = paste0(integration_factor, " Successive Dive Durations (s)"),
       ylab = paste0(integration_factor, " Successive Post Dive Intervals (s)"),
       main = file_basename)
  
  segments(x0 = min(df$x),
           y0 = aer_values$intercept + aer_values$slope * min(df$x),
           x1 = threshold,
           y1 = aer_values$intercept + aer_values$slope * threshold,
           col = "black", lwd = 2)
  
  segments(x0 = threshold,
           y0 = anaer_values$intercept + anaer_values$slope * threshold,
           x1 = max(df$x),
           y1 = anaer_values$intercept + anaer_values$slope * max(df$x),
           col = "black", lwd = 2)
  
  abline(v = threshold, col = "black", lty = 4, lwd = 0.5)
  dev.off()
  
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

write_csv(combined_results, file = "constraint_line_analysis.csv")

####Slope and continuity violations####
library(tidyverse)

### --- Functions --- ###
prepare_data <- function(file) {
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  
  #Rename columns dd_mw -> x and pdi_mw -> y
  df <- df %>% rename(
    x = dd_mw,
    y = pdi_mw
  )
  
  return(df)
}

compute_error_sum <- function(dataframe, tau) {
  #Quantile regression
  mod1 <- rq(y ~ x, data = dataframe, tau = tau)
  coefficients <- coef(mod1)
  slope <- coefficients[2]
  intercept <- coefficients[1]
  
  #Find points below quantile regression line
  dataframe$predicted_y <- intercept + (slope * dataframe$x)
  dataframe$residuals <- dataframe$y - dataframe$predicted_y
  
  below_line <- dataframe %>% filter(y < predicted_y)
  k <- nrow(below_line)
  
  #Find error of points below regression line
  error_sum <- sqrt(sum(below_line$residuals^2) / k)
  
  return(list(tau = tau, 
              k = k, 
              error_sum = error_sum, 
              slope = slope, 
              intercept = intercept))
}

#Find edge
find_edge <- function(dataframe){
  
  #Prep results table to store quantile regression results
  results <- data.frame(tau = numeric(),
                        k = numeric(),
                        error_sum = numeric(),
                        slope = numeric(),
                        intercept = numeric())
  
  #Reduce tau from 0.5 by .001 increments
  tau <- 0.5
  while (tau >= 0.001) {
    result <- compute_error_sum(dataframe = dataframe, tau = tau)
    if (result$k >= 10 && result$k <= 40) {
      results <- rbind(results, data.frame(tau = result$tau,
                                           k = result$k,
                                           error_sum = result$error_sum,
                                           slope = result$slope,
                                           intercept = result$intercept))
    }
    tau <- tau - 0.001
  }
  
  best_tau <- results %>%
    filter(k >= 10 & k <= 40) %>%
    arrange(error_sum) %>%
    slice(1)
  
  if (nrow(best_tau) == 0) {
    best_tau <- data.frame(tau = NA, error_sum = NA)
  }
  
  return(data.frame(best_tau = best_tau$tau, 
                    error_sum = best_tau$error_sum,
                    slope = best_tau$slope,
                    intercept = best_tau$intercept))
}

### --- Main Loop --- ###

previous_results <- read_csv("constraint_line_analysis.csv")

analysis_results <- data.frame()

for (i in 1:nrow(previous_results)){
  id <- previous_results$animal_id[i]
  file <- paste0(id, "_mw.csv")
  
  
  if (df_capture$age[df_capture$sn == id] == "yearling"){integration_factor <- 12} else {
    integration_factor <- 11}
  
  if (df_capture$age[df_capture$sn == id] == "yearling"){upper_limit <- 5256} else {
    upper_limit <- 6732}
  
  #Pull row from constraint line results
  df_id <- previous_results[i, ]
  
  #Prepare data
  df_og <- prepare_data(file)
  
  if (df_id$angle_diff < 40){
    slope_violation <- "YES"
    continuity_violation <- "NA"
    
    #rerun QR on aerobic portion of data
    df_aer <- df_og %>% filter(x < upper_limit)
    best_fit <- find_edge(df_aer)
    
    #open png for plotting
    png(filename = paste0(id, "_plot.png"))
    
    #Plot new constraint line
    plot(y~x, data = df_og, 
         ylim = c(0,900),
         xlab = paste0(integration_factor, " Sucessive Dive Durations (s)"),
         ylab = paste0(integration_factor, " Successive Post Dive Intervals (s)"),
         main = paste0(id))
    
    segments(x0 = min(df_og$x), 
             y0 = best_fit$intercept + best_fit$slope * min(df_og$x),
             x1 = max(upper_limit), 
             y1 = best_fit$intercept + best_fit$slope * max(upper_limit),
             col = "grey0",
             lwd = 2)
    
    #close png
    dev.off()
    
  } else {
    slope_violation <- "NO"
    
    #Continuity violation check
    threshold <- quantile(df_og$x, df_id$split)
    x_intersection <- (df_id$intercept_anaer - df_id$intercept_aer) / (df_id$slope_aer - df_id$slope_anaer)
    
    #Check if x intersection falls within 10% of the split
    upper_limit <- threshold * .9
    lower_limit <- threshold * 1.1
    
    if (x_intersection > lower_limit & x_intersection < upper_limit){
      continuity_violation <- "NO"
      best_fit <- data.frame(slope = NA, intercept = NA)
    } else {
      continuity_violation <- "YES"
      
      best_fit <- data.frame(slope = NA, intercept = NA)
    }}
  
  #store results
  analysis_results <- rbind(analysis_results, data.frame(
    animal_id = id,
    s_violation = slope_violation,
    c_violation = continuity_violation,
    slope = best_fit$slope,
    intercept = best_fit$intercept))
}

write_csv(analysis_results, file = "violations.csv")
