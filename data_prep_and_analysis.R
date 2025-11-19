###Data Prep and Constraint Line analysis###
###Using quantreg package###
####Filter by 4m and 20s and add foraging trip ID ####
library(tidyverse)
setwd("C:/Users/kpkolda/Desktop/Thesis/NFS Data/Final_Run_1")

trips <- read_csv("trips.csv") #read in the foraging trips of all individuals
trips$Start_gmt <- as.POSIXct(trips$Start_gmt, format ="%m/%d/%Y %H:%M", tz = "GMT")  #turn into datetime
trips$End_gmt <- as.POSIXct(trips$End_gmt, format ="%m/%d/%Y %H:%M", tz = "GMT")  #turn into datetime


file_list <- list.files(pattern = '*_divestats.csv')
for (file in file_list){
  
  #load data
  df <- read_csv(file)
  base_name <- tools::file_path_sans_ext(basename(file))
  base_name <- sub("_CF_divestats$", "", base_name)
  
  #Pulls trip ID number
  ID_trips <- trips %>% 
    filter(Id == base_name)
  
  df_trip <- df %>% 
    rowwise() %>% 
    mutate(
      trip_no = ID_trips %>%
        filter(begdesc >= Start_gmt & begdesc <= End_gmt) %>%
        pull(Trip_no)
    ) %>%
    ungroup()
  
  #filter for dives > 4m and > 20s
  df_filter <- df_trip[df$maxdep > 4 & df$divetim > 20,]
  
  #calculate dive times and post dive intervals
  df_calc <- df_filter %>% 
    mutate(
      end_time = begdesc + seconds(divetim),
      next_begdesc = lead(begdesc),
      pdi = as.numeric(difftime(next_begdesc, end_time, units = "secs"))
    )
  
  write.csv(df_calc, file = paste0(base_name, "_filtered_divestats.csv"))
}

####Moving Sum Calculations####
library(tidyverse)
library(diveMove)
library(gdata)

#load in capture data and modify data frame
capture_data <- read.csv("Capture_data.csv")
capture_data$AnimalID <- sub("^FC", "", capture_data$AnimalID)

#Create empty results to store data from loop
results <- data.frame(AnimalID = character(),
                      mass = numeric(),
                      median_bADL = numeric(),
                      max = numeric(),
                      ndives_bADL = numeric(),
                      stringsAsFactors = FALSE)

file_list <- list.files(pattern = '*_filtered_divestats.csv')

for (file in file_list){
  df <- read_csv(file)
  file_basename <- basename(file)
  ID <- sub('*_filtered_divestats.csv', '', file_basename)
  ndive <- NROW(df)
  
  animal_mass <- capture_data[capture_data$AnimalID == ID, "Mass"]
  
  filtered_data <- df$divetim[df$divetim>156]
  median_bADL <- median(filtered_data, na.rm = TRUE)
  max <- max(df$divetim, na.rm = TRUE)
  ndives_bADL <- length(filtered_data)
  
  percent_dives_beyond <- (ndives_bADL/ndive) * 100
  
  temp_df <- data.frame(AnimalID = ID,
                        mass = animal_mass,
                        median_bADL = median_bADL,
                        max = max,
                        ndives_bADL = ndives_bADL,
                        ndive = ndive,
                        percent_dives_beyond = percent_dives_beyond,
                        stringsAsFactors = FALSE)
  results <- rbind(results, temp_df)
}
results_narm <- subset(results, AnimalID != "SP1406")

calculate_moving_sum <- function(df, ndives, cADL){
  
  df_temp <- subset(df, ndives_bADL > ndives)
  avg_mass <- mean(df_temp$mass)
  avg_max <- mean(df_temp$max)
  avg_median_bADL <- mean(df_temp$median_bADL, na.rm = TRUE)
  animal_ids <- unique(df_temp$AnimalID)
  N <- length(unique(df_temp$AnimalID))
  moving_window <- (avg_max - cADL)/(avg_median_bADL - cADL)
  p_anaerobic <- mean(df$percent_dives_beyond)
  sd <- sd(df$percent_dives_beyond)
  
  return(list(
    Integration_Factor = moving_window,
    Average_Mass = avg_mass,
    Average_Max_Dive_Duration = avg_max,
    Average_Median_bADL = avg_median_bADL,
    p_anaerobic = p_anaerobic,
    sd = sd,
    N = N))
}

calculate_moving_sum(results_narm, 100, 156)
calculate_moving_sum(results_narm, 50, 156)
calculate_moving_sum(results_narm, 30, 156)
calculate_moving_sum(results_narm, 10, 156)
calculate_moving_sum(results_narm, 0, 156)

####Applying the moving sums####
library(tidyverse)
library(RcppRoll)
library(zoo)


file_list <- list.files(pattern = '*filtered_divestats.csv')

for (file in file_list){
  df <- read_csv(file)
  df <- df %>% 
    rename(dive_no = ...2)  %>% 
    group_by(trip_no) %>%
    filter(row_number() < n()) %>% #removes final point of each foraging trip
    ungroup() %>%
    select(divetim, pdi, trip_no, dive_no, begdesc) %>%
    group_by(trip_no) %>%
    mutate(
      dd_mw = roll_sum(divetim, n = 5, fill = NA, align = "left"),
      pdi_mw = roll_sum(pdi, n = 5, fill = NA, align = "left"),
      first_dive_in_sum = ifelse(!is.na(dd_mw) & !is.na(pdi_mw), dive_no, NA),
      trip_no = ifelse(!is.na(dd_mw) & !is.na(pdi_mw), trip_no, NA),
      begdesc = ifelse(!is.na(dd_mw) & !is.na(pdi_mw), yday(begdesc), NA)
    ) %>% 
    ungroup %>% 
    select(first_dive_in_sum, dd_mw, pdi_mw, trip_no, begdesc) %>% 
    rename(jday = begdesc)
  df <- na.omit(df)
  
  base_name <- tools::file_path_sans_ext(basename(file))
  base_name <- sub("_filtered_divestats$", "", base_name)
  write_csv(df, file = paste0(base_name, "_mw.csv"))
}

#Plot
file_list <- list.files(pattern = '*mw.csv')
for (file in file_list){
  file_basename <- basename(file)
  file_basename <- sub('*_mw.csv', '', file_basename)
  df <- read_csv(file)
  names(df)[2] <- "x"
  names(df)[3] <- "y"
  df <- df %>% 
    na.omit(df) %>% 
    filter(y<900)
  
  plot(y~x, data = df,
       xlab = "5 Sucessive Dive Durations (s)",
       ylab = "5 Successive Post Dive Intervals (s)",
       main = paste0(file_basename))
}

####Constraint Line analysis####
library(tidyverse)
library(quantreg)

### --- Functions --- ###

#Prepare data
prepare_data <- function(file) {
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
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
  
  png(filename = paste0(file_basename, "_prelim_plot.png"))
  plot(y ~ x, data = df, ylim = c(0, 900),
       xlab = "5 Successive Dive Durations (s)",
       ylab = "5 Successive Post Dive Intervals (s)",
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
####Slope and continuity violations####
library(tidyverse)

### --- Functions --- ###

#Prepare data
prepare_data <- function(file) {
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
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
  
  #Pull row from constraint line results
  df_id <- previous_results[i, ]
  
  #Prepare data
  df_og <- prepare_data(file)
  
  if (df_id$angle_diff < 40){
    slope_violation <- "YES"
    continuity_violation <- "NA"
    
    #rerun QR on aerobic portion of data
    df_aer <- df_og %>% filter(x < 780)
    best_fit <- find_edge(df_aer)
    
    #open png for plotting
    png(filename = paste0(id, "_plot.png"))
    
    #Plot new constraint line
    plot(y~x, data = df_og, 
         ylim = c(0,900),
         xlab = "5 Sucessive Dive Durations (s)",
         ylab = "5 Successive Post Dive Intervals (s)",
         main = paste0(id))
    
    segments(x0 = min(df_og$x), 
             y0 = best_fit$intercept + best_fit$slope * min(df_og$x),
             x1 = max(780), 
             y1 = best_fit$intercept + best_fit$slope * max(780),
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

####Diving performance metrics####
library(tidyverse)

### --- Daily Values --- ###

combined_results <- read_csv("constraint_line_analysis.csv")
violation_results <- read_csv("violations.csv")

#Pull capture data
df_capture <- read_csv("Capture_data.csv")
df_capture$Date <- as.Date(df_capture$Date, format = "%m/%d/%Y")
df_capture$Re_date <- as.Date(df_capture$Re_date, format = "%m/%d/%Y")
df_capture <- df_capture[df_capture$AnimalID != "FCSP1406",]
df_capture$AnimalID <- sub("^FC", "", df_capture$AnimalID)
df_capture$mass_change <- (df_capture$Mass - df_capture$Re_mass)/df_capture$Mass * 100

residuals <- data.frame()

violation_results <- violation_results[violation_results$animal_id != "SP1406",]

for (i in 1:nrow(violation_results)){
  id <- violation_results$animal_id[i]
  file <- paste0(id, "_mw.csv")
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
  
  s_violation <- violation_results$s_violation[i]
  
  #pull capture data
  df_capture_ind <- df_capture[df_capture$AnimalID == id,]
  
  #load divestats
  file_2 <- paste0(id, "_filtered_divestats.csv")
  df_divestats <- read_csv(file_2)
  
  df_divestats <- df_divestats %>% 
    rename(dive_no = ...2)  %>% 
    group_by(trip_no) %>%
    filter(row_number() < n()) %>% # removes final point of each foraging trip
    ungroup() %>%
    select(dive_no, maxdep, begdesc) %>%
    mutate(begdesc = yday(begdesc)) %>% 
    rename(jday = begdesc) %>%
    na.omit()
  
  # Filter to only dives in df$first_dive_in_sum
  df_filtered <- df_divestats %>% 
    filter(dive_no %in% df$first_dive_in_sum)
  
  if (s_violation == "YES"){
    slope <- violation_results$slope[i]
    intercept <- violation_results$intercept[i]
    
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    
    # Aerobic boundary limit
    aerobic_limit <- 780
    
    # Filter df_2 for aerobic dives and residuals > 0
    df_2 <- df %>%
      filter(x < aerobic_limit, residuals > 0) %>%
      group_by(jday) %>%
      summarise(
        dres = median(residuals, na.rm = TRUE),
        prox = (median(x, na.rm = TRUE) / aerobic_limit) * 100
      ) %>%
      ungroup()
    
    # Find dives within aerobic boundary for median_max_depth calc
    aerobic_dive_nos <- df %>%
      filter(x < aerobic_limit) %>%
      pull(first_dive_in_sum)
    
    # Filter df_filtered to only those dives within aerobic limit
    df_filtered_aerobic <- df_filtered %>%
      filter(dive_no %in% aerobic_dive_nos)
    
    # Calculate median max depth for aerobic dives
    df_median_descdist <- df_filtered_aerobic %>%
      group_by(jday) %>%
      summarise(median_max_depth = median(maxdep, na.rm = TRUE)) %>%
      ungroup()
    
  } else {
    combined_results_id <- combined_results[combined_results$animal_id == id,]
    slope <- combined_results_id$slope_aer
    intercept <- combined_results_id$intercept_aer
    
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    
    x_max <- combined_results_id$x_intersect
    
    df_2 <- df %>%
      filter(x < x_max, residuals > 0) %>%
      group_by(jday) %>%
      summarise(
        dres = median(residuals, na.rm = TRUE),
        prox = (median(x, na.rm = TRUE) / x_max) * 100
      ) %>%
      ungroup()
    
    # Find dives within aerobic boundary for median_max_depth calc
    aerobic_limit <- x_max
    aerobic_dive_nos <- df %>%
      filter(x < aerobic_limit) %>%
      pull(first_dive_in_sum)
    
    df_filtered_aerobic <- df_filtered %>%
      filter(dive_no %in% aerobic_dive_nos)
    
    df_median_descdist <- df_filtered_aerobic %>%
      group_by(jday) %>%
      summarise(median_max_depth = median(maxdep, na.rm = TRUE)) %>%
      ungroup()
  }
  
  df_2$animalid <- id
  df_2$mass <- df_capture_ind$Mass
  df_2$year <- format(df_capture_ind$Date, "%Y")
  df_2$mass_change <- df_capture_ind$Mass_diff
  
  df_2 <- df_2 %>%
    left_join(df_median_descdist, by = "jday")
  
  residuals <- rbind(residuals, df_2)
}


write_csv(residuals, file = "daily_metrics.csv")







### --- Broad Dive Metrics --- ###

dive_metrics <- data.frame()

for (i in 1:nrow(violation_results)){
  id <- violation_results$animal_id[i]
  file <- paste0(id, "_mw.csv")
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
  
  s_violation <- violation_results$s_violation[i]
  
  if (s_violation == "YES"){
    slope <- violation_results$slope[i]
    intercept <- violation_results$intercept[i]
    
    #find aerobic median residual
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    median_residual <- median(df$residuals[df$x < 780 & df$residuals > 0])
    
    #Find % median dive duration of ADL
    valid_dives <- df$x[df$x < 780]
    median_duration <- median(valid_dives, na.rm = TRUE)
    prox <- (median_duration / 780) * 100
    
  } else {
    combined_results_id <- combined_results[combined_results$animal_id == id,]
    slope <- combined_results_id$slope_aer
    intercept <- combined_results_id$intercept_aer
    
    #find residuals and take daily averages
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    
    #find max aerobic x value
    x_max <- combined_results_id$x_intersect
    
    #find aerobic median residual
    median_residual <- median(df$residuals[df$x < x_max & df$residuals > 0])
    
    #Find % median dive duration of ADL
    valid_dives <- df$x[df$x < x_max]
    median_duration <- median(valid_dives, na.rm = TRUE)
    prox <- (median_duration / x_max) * 100
    
  }
  
  summary_row <- data.frame(
    animal_id = id,
    median_residual = median_residual,
    prox = prox)
  
  dive_metrics <- bind_rows(dive_metrics, summary_row)
}

#create dataframe with all metrics
names(df_capture)[1] <- "animal_id"
analysis_results <- left_join(df_capture, dive_metrics, by = "animal_id")
analysis_results <- left_join(analysis_results, combined_results, by = "animal_id")
analysis_results <- left_join(analysis_results, violation_results, by = "animal_id")

write_csv(analysis_results, file = "data.csv")







### --- Dive metrics during deployment overlap --- ###

df_1 <- read_csv("data.csv")
df_o_dive_metrics <- data.frame()

for (i in 1:nrow(df_1)){
  id <- df_1$animal_id[i]
  file <- paste0(id, "_mw.csv")
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
  
  s_violation <- df_1$s_violation[i]
  year <- df_1$Year[i]
  
  julian_day <- df[[5]]
  
  if (year %in% c(2010, 2014)) {
    df <- df[julian_day >= 225 & julian_day <= 272, ]
  } else if (year == 2016) {
    df <- df[julian_day >= 226 & julian_day <= 273, ]
  } else {
    next
  }
  
  if (nrow(df) == 0) next
  
  if (s_violation == "YES") {
    slope <- df_1$slope[i]
    intercept <- df_1$intercept[i]
    
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    median_residual <- median(df$residuals[df$x < 780 & df$residuals > 0], na.rm = TRUE)
    
    valid_dives <- df$x[df$x < 780]
    median_duration <- median(valid_dives, na.rm = TRUE)
    prox <- (median_duration / 780) * 100
    
  } else {
    slope <- df_1$slope_aer[i]
    intercept <- df_1$intercept_aer[i]
    x_max <- df_1$x_intersect[i]
    
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    median_residual <- median(df$residuals[df$x < x_max & df$residuals > 0], na.rm = TRUE)
    
    valid_dives <- df$x[df$x < x_max]
    median_duration <- median(valid_dives, na.rm = TRUE)
    prox <- (median_duration / x_max) * 100
  }
  
  summary_row <- data.frame(
    animal_id = id,
    overlap_median_residual = median_residual,
    overlap_prox = prox,
    year <- year
  )
  df_o_dive_metrics <- bind_rows(df_o_dive_metrics, summary_row)
}

df_o_dive_metrics$year....year <- as.factor(df_o_dive_metrics$year....year)
plot_overlap_mr <- boxplot(median~year, data = df_o_dive_metrics)

df_o_dive_metrics$year....year <- NULL

df_1 <- left_join(df_1, df_o_dive_metrics, by = "animal_id")

write_csv(df_1, file = "data.csv")






### --- Broad dive depths --- ###

summary_df <- data.frame()  # initialize empty summary data frame

for (i in 1:nrow(df_1)){
  id <- df_1$animal_id[i]
  file <- paste0(id, "_mw.csv")
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
  
  s_violation <- df_1$s_violation[i]
  
  file_2 <- paste0(id, "_filtered_divestats.csv")
  df_divestats <- read_csv(file_2)
  
  df_divestats <- df_divestats %>% 
    rename(dive_no = ...2)  %>% 
    group_by(trip_no) %>%
    filter(row_number() < n()) %>% 
    ungroup() %>%
    select(dive_no, maxdep, begdesc) %>%
    mutate(begdesc = yday(begdesc)) %>% 
    rename(jday = begdesc) %>%
    na.omit()
  
  df_filtered <- df_divestats %>% 
    filter(dive_no %in% df$first_dive_in_sum)
  
  if (s_violation == "YES"){
    slope <- df_1$slope[i]
    intercept <- df_1$intercept[i]
    aerobic_limit <- 780
    
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    
    aerobic_dive_nos <- df %>%
      filter(x < aerobic_limit) %>%
      pull(first_dive_in_sum)
    
    df_filtered_aerobic <- df_filtered %>%
      filter(dive_no %in% aerobic_dive_nos)
    
  } else {
    slope <- df_1$slope_aer[i]
    intercept <- df_1$intercept_aer[i]
    x_max <- df_1$x_intersect[i]
    
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    
    aerobic_limit <- x_max
    
    aerobic_dive_nos <- df %>%
      filter(x < aerobic_limit) %>%
      pull(first_dive_in_sum)
    
    df_filtered_aerobic <- df_filtered %>%
      filter(dive_no %in% aerobic_dive_nos)
  }
  
  year_num <- as.numeric(format(df_1$Date[i], "%Y"))
  
  #Median max depth across all aerobic dives
  median_max_depth_all <- median(df_filtered_aerobic$maxdep, na.rm = TRUE)
  
  #Filter by julian day depending on year
  if (year_num %in% c(2010, 2014)) {
    df_filtered_aerobic_filtered <- df_filtered_aerobic %>%
      filter(jday >= 225 & jday <= 272)
  } else if (year_num == 2016) {
    df_filtered_aerobic_filtered <- df_filtered_aerobic %>%
      filter(jday >= 226 & jday <= 273)
  } else {
    df_filtered_aerobic_filtered <- df_filtered_aerobic[0,]
  }
  
  #Median max depth filtered by julian day range
  median_max_depth_filtered <- median(df_filtered_aerobic_filtered$maxdep, na.rm = TRUE)
  
  #Create one summary row per individual
  summary_row <- data.frame(
    animal_id = id,
    median_max_depth_all = median_max_depth_all,
    median_max_depth_filtered = median_max_depth_filtered
  )
  
  summary_df <- bind_rows(summary_df, summary_row)
}


df_1 <- left_join(df_1, summary_df, by = "animal_id")

write_csv(df_1, file = "data.csv")






####Statistical Analyses####
library(tidyverse)
library(mgcv)
library(ggpubr)
library(FSA)

df_daily_metrics <- read_csv("daily_metrics.csv")
df_1 <- read_csv("data.csv")
df_1$Year <- as.factor(df_1$Year)

### --- Yearly comparisons --- ###
qqnorm(df_1$median_residual)
shapiro.test(df_1$median_residual)
ktest_mr <- kruskal.test(median_residual ~ Year, data = df_1)
dtest_mr <- dunnTest(median_residual ~ Year, data = df_1, method = "holm")
plot_mr <- boxplot(median_residual~Year, 
                   data = df_1, 
                   col = c("lightsteelblue1", "lightsteelblue3", "lightsteelblue4"),
                   xlab = "",
                   ylab = "",
                   main = "")

qqnorm(df_1$prox)
shapiro.test(df_1$prox)
aov_prox <- aov(prox ~Year, data = df_1)
summary(aov_prox)

qqnorm(df_1$median_max_depth_all)
shapiro.test(df_1$median_max_depth_all)
aov_depth <- aov(median_max_depth_all~Year, data = df_1)
summary(aov_depth)
plot_depth <- boxplot(df_1$median_max_depth_all~df_1$Year)

### --- Yearly Comparisons during Overlap --- ###
shapiro.test(df_1$overlap_median_residual)
ktest_mr_overlap <- kruskal.test(overlap_median_residual ~ Year, data = df_1)
dtest_mr_overlap <- dunnTest(overlap_median_residual ~ Year, data = df_1, method = "holm")
print(dtest_mr_overlap)

shapiro.test(df_1$overlap_prox)
aov_prox_overlap <- aov(overlap_prox ~Year, data = df_1)
summary(aov_prox_overlap)

qqnorm(df_1$median_max_depth_all)
shapiro.test(df_1$median_max_depth_all)
aov_depth_overlap <- aov(df_1$median_max_depth_filtered~df_1$Year)
summary(aov_depth_overlap)
boxplot(df_1$median_max_depth_filtered~df_1$Year)





### --- GAMM --- ###
df_daily_metrics <- read_csv("daily_metrics.csv")
df_daily_metrics$mc <- (df_daily_metrics$mass_change/df_daily_metrics$mass) * 100

plot <- plot(df_1$mass_change~df_1$Mass)
plot_2 <- plot(df_1$prox~df_1$Mass)
plot_3 <- plot(df_1$median_residual~df_1$Mass)
df_daily_metrics$animalid <- as.factor(df_daily_metrics$animalid)
df_daily_metrics$year <- as.factor(df_daily_metrics$year)

#check for normality: it is right skewed
plot_1 <- hist(df_daily_metrics$dres, breaks = 52)
qqnorm(df_daily_metrics$dres)
shapiro.test(df_daily_metrics$dres)


plot_2 <- hist(df_daily_metrics$median_max_depth, breaks = 100)
qqnorm(df_daily_metrics$median_max_depth)
shapiro.test(df_daily_metrics$median_max_depth)

#GAMM with log transform
df_daily_metrics$log_dres <- log(df_daily_metrics$dres + 1e-6)

mod1 <- gam(log_dres ~ year + 
              s(jday, by = year) + 
              s(mc, by = year, k = 40) +
              s(mass, by = year, k = 40) +
              s(median_max_depth, by = year, k = 30) +
              s(animalid, bs = "re"), 
            family = gaussian(), 
            data = df_daily_metrics,
            method = "REML")

summary(mod1)
gam.check(mod1)
plot(mod1)

mod2 <- gam(log_dres ~ year + 
              s(jday, by = year) + 
              s(mass, by = year, k = 40) +
              s(median_max_depth, by = year, k = 30) +
              s(animalid, bs = "re"), 
            family = gaussian(), 
            data = df_daily_metrics,
            method = "REML")

summary(mod2)
gam.check(mod2)
plot(mod2)

AIC(mod1, mod2)

mod3 <- gam(log_dres ~ year + 
              s(jday, by = year) +
              s(median_max_depth, by = year, k = 30) +
              s(animalid, bs = "re"), 
            family = gaussian(), 
            data = df_daily_metrics,
            method = "REML")

summary(mod3)
gam.check(mod3)
plot(mod3)

AIC(mod2, mod3)

mod_snull<- gam(log_dres ~ 1,
                family = gaussian(),
                data = df_daily_metrics,
                method = "REML")

mod_null <- gam(log_dres ~ s(animalid, bs = "re"),
                family = gaussian(),
                data = df_daily_metrics,
                method = "REML")

mod_no_interactions <- gam(log_dres ~ year + 
                             s(jday, k = 25) +
                             s(median_max_depth, k = 30) +
                             s(animalid, bs = "re"), 
                           family = gaussian(), 
                           data = df_daily_metrics,
                           method = "REML")

summary(mod_no_interactions)
gam.check(mod_no_interactions)
plot(mod_no_interactions)

AIC(mod1, mod2, mod3, mod_snull, mod_null, mod_no_interactions)

#Re-do all without interactions
mod1 <- gam(log_dres ~ year + 
              s(jday) + 
              s(mc, k = 40) +
              s(mass, k = 40) +
              s(median_max_depth, k = 30) +
              s(animalid, bs = "re"), 
            family = gaussian(), 
            data = df_daily_metrics,
            method = "REML")

summary(mod1)
gam.check(mod1)
plot(mod1)

mod2 <- gam(log_dres ~ year + 
              s(jday) +
              s(mass, k = 40) +
              s(median_max_depth, k = 30) +
              s(animalid, bs = "re"), 
            family = gaussian(), 
            data = df_daily_metrics,
            method = "REML")

summary(mod2)
gam.check(mod2)
plot(mod2)

mod3 <- gam(log_dres ~ year + 
              s(jday, k = 25) +
              s(median_max_depth, k = 30) +
              s(animalid, bs = "re"), 
            family = gaussian(), 
            data = df_daily_metrics,
            method = "REML")

summary(mod3)
gam.check(mod3)
plot(mod3)

AIC(mod1, mod2, mod3)

### --- Trips --- ###
df_trips <- read_csv("trips.csv")
plot_trip <- boxplot(df_trips$Trip_dur~df_trips$Year)

df_avg_trips <- df_trips %>% 
  group_by(Id) %>% 
  summarise(
    animalid = first(Id),
    avg_trip_dur = mean(Trip_dur),
    year = first(Year)
  )

plot_trip <- boxplot(df_avg_trips$avg_trip_dur~df_avg_trips$year)
hist(df_avg_trips$avg_trip_dur)
shapiro.test(df_avg_trips$avg_trip_dur)
ktest_trips <- kruskal.test(df_avg_trips$avg_trip_dur~df_avg_trips$year)

plot_weight <- boxplot(df_1$Mass~df_1$Year)
shapiro.test(df_1$Mass)
aov_mass <- aov(df_1$Mass~df_1$Year)
summary(aov_mass)

library(lme4)
library(MuMIn)
trips <- ggplot(df_trips, aes(x = start_julian, y = Trip_dur, color = Id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

df_trips$Year <- as.factor(df_trips$Year)
model_trips <- lmer(Trip_dur~start_julian * Year + (1 | Id), data = df_trips)
summary(model_trips)
r.squaredGLMM(model_trips)

ggplot(df_trips, aes(x = start_julian, y = Trip_dur, color = Year)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
       x = "Julian Day", y = "Trip Duration (days)") +
  theme_minimal()

#trips overlap
df_capture <- read_csv("Capture_data.csv")
df_capture$Date <- as.POSIXct(df_capture$Date, format = "%m/%d/%Y")
df_capture$Re_date <- as.POSIXct(df_capture$Re_date, format = "%m/%d/%Y")
df_trips$Start_gmt <- as.POSIXct(df_trips$Start_gmt, format = "%m/%d/%Y %H:%M")
df_trips$End_gmt <- as.POSIXct(df_trips$End_gmt, format = "%m/%d/%Y %H:%M")

df_trips$start_julian <- as.numeric(format(df_trips$Start_gmt, "%j"))
df_trips$end_julian   <- as.numeric(format(df_trips$End_gmt,   "%j"))

filter_trips <- function(df, year_num) {
  if (year_num %in% c(2010, 2014)) {
    df %>%
      filter(
        Year == year_num,
        start_julian >= 225 & start_julian <= 272,
        end_julian >= 225 & end_julian <= 272
      )
  } else if (year_num == 2016) {
    df %>%
      filter(
        Year == year_num,
        start_julian >= 226 & start_julian <= 273,
        end_julian >= 226 & end_julian <= 273
      )
  } else {
    df %>% filter(Year == year_num)
  }
}

years_to_check <- c(2010, 2014, 2016)


avg_trip_dur_by_id <- lapply(years_to_check, function(y) {
  filtered <- filter_trips(df_trips, y)
  filtered %>%
    group_by(Id) %>%
    summarise(
      Year = y,
      Avg_Trip_Duration = mean(Trip_dur, na.rm = TRUE),
      Trips_Count = n()
    ) %>%
    ungroup()
})

avg_trip_dur_df <- bind_rows(avg_trip_dur_by_id)

plot_trip <- boxplot(avg_trip_dur_df$Avg_Trip_Duration~avg_trip_dur_df$Year)
shapiro.test(avg_trip_dur_df$Avg_Trip_Duration)
aov_trips_overlap <- aov(avg_trip_dur_df$Avg_Trip_Duration~avg_trip_dur_df$Year)
summary(aov_trips_overlap)


#dataframe with overlapping trips
library(dplyr)

filter_trips <- function(df, year_num) {
  if (year_num %in% c(2010, 2014)) {
    df %>%
      filter(
        Year == year_num,
        start_julian >= 225 & start_julian <= 272,
        end_julian >= 225 & end_julian <= 272
      )
  } else if (year_num == 2016) {
    df %>%
      filter(
        Year == year_num,
        start_julian >= 226 & start_julian <= 273,
        end_julian >= 226 & end_julian <= 273
      )
  } else {
    df %>% filter(Year == year_num)
  }
}

years_to_include <- c(2010, 2014, 2016)

filtered_trips_all_years <- bind_rows(
  lapply(years_to_include, function(y) filter_trips(df_trips, y))
)


trips_overlap <- ggplot(filtered_trips_all_years, aes(x = start_julian, y = Trip_dur, color = Id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


# Combine all years into one dataframe
average_trip_durations_by_id_df <- do.call(rbind, average_trip_durations_by_id)

print(average_trip_durations_by_id_df)
#### Figures ####
library(tidyverse)
library(patchwork)
library(ggpubr)

df_1 <- read_csv("data.csv")
df_figure <- read_csv("NFSF2716_mw.csv")
names(df_figure)[2] <- "x"
names(df_figure)[3] <- "y"
data_figure <- df_1[df_1$animal_id == "NFSF2716",]

df_figure$predicted <- (df_figure$x * data_figure$slope) + data_figure$intercept
df_figure$residuals <- df_figure$y - df_figure$predicted
valid_dives <- df_figure$x[df_figure$x < 780]
median_duration <- median(valid_dives, na.rm = TRUE)
prox <- (median_duration / 780) * 100




### --- Figure 1 --- ###
#Define regression segment coordinates
x_start <- 115
x_end <- 780
y_start <- data_figure$intercept + data_figure$slope * 115
y_end   <- data_figure$intercept + data_figure$slope * x_end

#panel A
panel_A <- ggplot(df_figure, aes(x = x, y = y)) + 
  geom_point(shape = 1, color = "#222428") +
  
  
  annotate("segment",
           x = x_start, y = y_start,
           xend = x_end, yend = y_end,
           color = "#222428", linewidth = 1.25) +
  
  scale_x_continuous(limits = c(115, 790),
                     breaks = seq(100, 800, by = 300)) +
  scale_y_continuous(limits = c(0, 900),
                     breaks = seq(0, 900, by = 300)) +
  labs(x = NULL, y = NULL, title = NULL) +
  
  theme_classic() +
  theme(axis.text = element_text(size = 14))

#Print the plot
print(panel_A)

#Panel B
panel_b <- ggplot(df_figure, aes(x = x, y = y)) + 
  geom_point(shape = 1, color = "lightsteelblue3") + 
  
  annotate("segment",
           x = x_start, y = y_start,
           xend = x_end, yend = y_end,
           color = "lightsteelblue3", linewidth = 1.25) +
  
  annotate("point",
           x = 236.0448,
           y = 260.9280,
           size = 4,
           shape = 19,
           color = "#222428") + 
  
  annotate("segment",
           x = 236.0448, y = 260.9280,
           xend = 236.0448, yend = 79.00477,
           color = "#222428", linewidth = 1.5, 
           linetype = 3) +
  
  annotate("point",
           x = 354.1535,
           y = 540.8640,
           size = 4,
           shape = 19,
           color = "#222428") +
  
  annotate("segment",
           x = 354.1535, y = 540.8640,
           xend = 780, yend = 540.8640,
           color = "#191a19", linewidth = 1.5,
           linetype = 5) +
  
  scale_x_continuous(limits = c(115, 790),
                     breaks = seq(100, 800, by = 300)) +
  scale_y_continuous(limits = c(0, 900),
                     breaks = seq(0, 900, by = 300)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 14))

print(panel_b)

final_plot <- panel_A/panel_b

final_plot









### --- Figure 2 --- ###

ggplot(df_1, aes(x = as.factor(Year), y = median_residual)) +
  geom_boxplot(fill = c("lightsteelblue1", "lightsteelblue3", "lightsteelblue4")) +
  labs(x = "", y = "") +
  theme_classic(base_size = 14) +
  
  # Add significance comparisons
  stat_compare_means(comparisons = list(c("2010", "2014"), c("2010", "2016")),
                     method = "t.test", 
                     label = "p.signif")









### --- GAMMM --- ###
library(mgcv)
library(ggplot2)
library(dplyr)

library(mgcv)
library(ggplot2)
library(dplyr)

# Ensure 'year' is a factor with all three levels
df_daily_metrics$year <- factor(df_daily_metrics$year, levels = c("2010", "2014", "2016"))

# Reference animal ID and other constants
ref_animal <- levels(df_daily_metrics$animalid)[1]
median_depth <- median(df_daily_metrics$median_max_depth, na.rm = TRUE)
jday_vals <- sort(unique(df_daily_metrics$jday))
depth_vals <- sort(unique(df_daily_metrics$median_max_depth))
median_jday <- median(df_daily_metrics$jday, na.rm = TRUE)
unique_years <- levels(df_daily_metrics$year)
cap_val <- -10  # Cap for log underflow


newdata_jday <- expand.grid(
  jday = jday_vals,
  median_max_depth = median_depth,
  year = unique_years,
  animalid = ref_animal
)

newdata_jday$year <- factor(newdata_jday$year, levels = unique_years)
newdata_jday$animalid <- factor(newdata_jday$animalid, levels = levels(df_daily_metrics$animalid))

# Predict, excluding random effect of animal ID
pred_jday <- predict(mod3, newdata = newdata_jday, se.fit = TRUE, exclude = "s(animalid)")

newdata_jday <- newdata_jday %>%
  mutate(
    fit = exp(pmax(pred_jday$fit, cap_val)),
    upper = exp(pmax(pred_jday$fit + 1.96 * pred_jday$se.fit, cap_val)),
    lower = exp(pmax(pred_jday$fit - 1.96 * pred_jday$se.fit, cap_val))
  )

# Plot Julian Day effect
p_jday <- ggplot(newdata_jday, aes(x = jday, y = fit, color = year, fill = year)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(x = "Julian Day", y = "Predicted Median Residuals (s)") +
  scale_color_manual(values = c("2010" = "lightsteelblue2", "2014" = "lightsteelblue3", "2016" = "black")) +
  scale_fill_manual(values = c("2010" = "lightsteelblue2", "2014" = "lightsteelblue3", "2016" = "black")) +
  theme_classic(base_size = 14)

print(p_jday)


newdata_depth <- expand.grid(
  median_max_depth = depth_vals,
  jday = median_jday,
  year = unique_years,
  animalid = ref_animal
)

newdata_depth$year <- factor(newdata_depth$year, levels = unique_years)
newdata_depth$animalid <- factor(newdata_depth$animalid, levels = levels(df_daily_metrics$animalid))

# Predict, excluding random effect
pred_depth <- predict(mod3, newdata = newdata_depth, se.fit = TRUE, exclude = "s(animalid)")

newdata_depth <- newdata_depth %>%
  mutate(
    fit = exp(pmax(pred_depth$fit, cap_val)),
    upper = exp(pmax(pred_depth$fit + 1.96 * pred_depth$se.fit, cap_val)),
    lower = exp(pmax(pred_depth$fit - 1.96 * pred_depth$se.fit, cap_val))
  )

# Plot Median Depth effect
p_depth <- ggplot(newdata_depth, aes(x = median_max_depth, y = fit, color = year, fill = year)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(x = "Median Dive Depth (m)", y = "Predicted Median Residuals (s)") +
  scale_color_manual(values = c("2010" = "lightsteelblue2", "2014" = "lightsteelblue3", "2016" = "black")) +
  scale_fill_manual(values = c("2010" = "lightsteelblue2", "2014" = "lightsteelblue3", "2016" = "black")) +
  theme_classic(base_size = 14)

print(p_depth)

###---GAMM Fixed --- ###
# Create a grid of jday values across the observed range
jday_grid <- data.frame(
  jday = seq(min(df_daily_metrics$jday, na.rm = TRUE),
             max(df_daily_metrics$jday, na.rm = TRUE),
             length.out = 200),
  median_max_depth = median(df_daily_metrics$median_max_depth, na.rm = TRUE),  # hold constant
  year = as.integer(median(df_daily_metrics$year, na.rm = TRUE)),              # hold constant
  animalid = NA                                                                 # NA for random effect
)

####Percent beyond upper limit####
combined_results <- df_1
dive_metrics <- data.frame()

for (i in 1:nrow(violation_results)){
  id <- violation_results$animal_id[i]
  file <- paste0(id, "_mw.csv")
  df <- read_csv(file)
  df <- df %>% arrange(desc(dd_mw)) %>% na.omit()
  names(df)[2] <- "x"
  names(df)[3] <- "y"
  
  s_violation <- violation_results$s_violation[i]
  
  if (s_violation == "YES"){
    slope <- violation_results$slope[i]
    intercept <- violation_results$intercept[i]
    
    # Find aerobic median residual
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    median_residual <- median(df$residuals[df$x < 780 & df$residuals > 0])
    
    # Find % median dive duration of ADL
    valid_dives <- df$x[df$x < 780]
    median_duration <- median(valid_dives, na.rm = TRUE)
    prox <- (median_duration / 780) * 100
    
    # Calculate % dives beyond upper limit (780)
    percent_beyond_upper <- (sum(df$x > 780) / nrow(df)) * 100
    
  } else {
    combined_results_id <- combined_results[combined_results$animal_id == id,]
    slope <- combined_results_id$slope_aer
    intercept <- combined_results_id$intercept_aer
    
    # Find residuals and take daily averages
    df$predicted <- (df$x * slope) + intercept
    df$residuals <- df$y - df$predicted
    
    # Find max aerobic x value
    x_max <- combined_results_id$x_intersect
    
    # Find aerobic median residual
    median_residual <- median(df$residuals[df$x < x_max & df$residuals > 0])
    
    # Find % median dive duration of ADL
    valid_dives <- df$x[df$x < x_max]
    median_duration <- median(valid_dives, na.rm = TRUE)
    prox <- (median_duration / x_max) * 100
    
    # Calculate % dives beyond upper limit (x_max)
    percent_beyond_upper <- (sum(df$x > x_max) / nrow(df)) * 100
  }
  
  summary_row <- data.frame(
    animal_id = id,
    median_residual = median_residual,
    prox = prox,
    percent_beyond_upper = percent_beyond_upper
  )
  
  dive_metrics <- bind_rows(dive_metrics, summary_row)
}

