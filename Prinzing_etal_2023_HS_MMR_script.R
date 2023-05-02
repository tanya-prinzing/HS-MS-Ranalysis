
# Prinzing et al. 2023 Horn Shark metabolic rate and gill area scaling analysis
# Estimating MMR

# Load libraries
library("tidyverse")
library("rollRegres")
library("lubridate")
library("respR")


#################################### ---------------------- #################################

###### ---------------- ########## -------   FUNCTIONS   ------- ########### ------------------ ######

#################################### -------------------------- ######################################


# ---- respR ---Function to tidy Fibox 4 .cvs MMR files, remane columns, filter lag time, and convert %DO into mg/L oxygen --- #
F4_MMR_convert_DO <- function(df, lag, end, S, interval = 1) { # Enter the raw data frame, the amount of lag time in minutes, and the salinity in mbar
  
  keep_columns <- c("Date", "Time", "Value", "Temp", "patm") # Create object of columns to keep
  # Filter columns to keep, rename "Value" to "O2air", change Temp, patm and O2air into decimal values
  F4_trial_df <- dplyr::select(df, all_of(keep_columns)) %>% 
    mutate(Time = as.character(Time), patm = as.double(patm), 
           Value = as.double(Value), Temp = as.double(Temp)) %>%
    rename(O2air = Value) %>%
    slice(1:n()-1) %>%
    mutate(TimeMinutes = seq(0, by=interval/60, length.out=nrow(.))) %>% # Make minutes column for plotting
    filter(TimeMinutes >= lag,
           TimeMinutes <= end) %>% # Filter out lag time when no shark was in the chamber
    mutate(TimeMinutes = seq(0, by=interval/60, length.out=nrow(.)))
  # Convert %DO to mg/L oxygen with convert DO function
  converted_df <- F4_trial_df %>%
    rowwise() %>%
      mutate(O2mgL = convert_DO(x=O2air, from = "%Air", to = "mg/L", 
                S = S, t = Temp, P = (patm/1000)))
  
  return(converted_df)
  
  
}



####################################### ----------------------- ##############################################


###### ----- respR -- FUNCTION to tidy and convert FIBOX 3 MMR data ------ ##############
F3_MMR_convert_DO <- function(df, lag, end, S, patm) {
  # Tidy up the raw data set
  MO2_calculation_df <- df %>% 
    rename(Time = `time/hh:mm:ss`, LogTime = `logtime/min`,
           O2air = `oxygen/% airsatur.`,Temp = `temp/∞C`, Amp = amp, Phase = `phase/∞`)  %>% 
    dplyr::select(-ErrorMessage) %>%
    mutate(TimeMinutes = as.double(LogTime), O2air = as.double(O2air), Temp = as.double(Temp)) %>%
    filter(TimeMinutes >= lag, TimeMinutes <= end) %>%
    mutate(TimeMinutes = seq(0, by=1/60, length.out=nrow(.)))
  
  
  # Convert %DO to mg/L oxygen with convert DO function
  converted_df <- MO2_calculation_df %>%
    rowwise() %>%
      mutate(O2mgL = convert_DO(x=O2air, from = "%Air", to = "mg/L", 
                S = S, t = Temp, P = (patm/1000)))
  
  return(converted_df)
  
}

#################################### ------------------------------------ #########################################




### --- MMR rolling regression function incorporating the increasing window width functionality --- ###
# Metabolic rate as * mgO2/hr *
MMR_rolling_regression <- function(df, hs, vr, vf, sm, lg, by, rsq, maxMMR = 2000) { 
  
  widths_seq <- seq(from= (sm*60), to=(lg*60), by = by) # declare sequence of interval sizes for rolling regressions
  trial_df <- df
  MMR_MO2 <- c()
  IntervalLength <- c()
  R_squareds <- c()
  RollEnd <- c()

  
  for (i in widths_seq) {
    # Output of rolling regression
    #MMR_interval_df <- 
    
    MMR_coefs_df <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, trial_df,
                                              do_compute = c("sigmas", "r.squareds", "1_step_forecasts"), 
                                              width = i)) %>% 
      mutate(RollEnd = seq(0, by=1/60, length.out=nrow(.))) %>%
      arrange(coefs.TimeMinutes) %>%
      filter(r.squareds >= rsq,
             coefs.TimeMinutes < 0 ) %>%
      filter(!is.na(sigmas)) %>% 
      mutate(MMR_MO2 = ((((vr - vf) * coefs.TimeMinutes)/vf)* -60),
             Interval = (i / 60),
             MO2_sd = sigmas) %>%
      filter(MMR_MO2 <= maxMMR) %>%
      slice(which.max(MMR_MO2))
    
    
    MMR_MO2 <- c(MMR_MO2, MMR_coefs_df$MMR_MO2)
    IntervalLength <- c(IntervalLength, MMR_coefs_df$Interval[1])
    R_squareds <- c(R_squareds, MMR_coefs_df$r.squareds)
    RollEnd <- c(RollEnd, MMR_coefs_df$RollEnd)

  }
  
  HSID <- rep(hs, length(IntervalLength)) # The specimen ID will be coded as a column on the output data frame here
  # Combine all estimates into a data frame
  MMR_combined_df <- data.frame(HSID, MMR_MO2, IntervalLength, R_squareds, RollEnd)
  
  MMR_RR_out_df <- MMR_combined_df %>%
    mutate(MassKG = vf,
           ChamVol = vr,
           Window_mid = (RollEnd - (IntervalLength/2))) # Calculate mid-point of the MMR window
  
  return(MMR_RR_out_df)
}


####################################### ----------------------- ##############################################



### --- MMR rolling regression function incorporating the increasing window width functionality --- ###
# Function for oxygen consumption mistakenly measured using a 5 second interval instead of 1 second interval
# Metabolic rate as * mgO2/hr *
MMR_5sec_rolling_regression <- function(df, hs, vr, vf, sm, lg, by, rsq) { 
  
  widths_seq <- seq(from= (sm*12), to=(lg*12), by = by) # declare sequence of interval sizes for rolling regressions
  trial_df <- df
  MMR_MO2 <- c()
  IntervalLength <- c()
  R_squareds <- c()
  RollEnd <- c()
  
  
  for (i in widths_seq) {
    # Output of rolling regression
    #MMR_interval_df <- 
    
    MMR_coefs_df <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, trial_df,
                                              do_compute = c("sigmas", "r.squareds", "1_step_forecasts"), 
                                              width = i)) %>% 
      mutate(RollEnd = seq(0, by=5/60, length.out=nrow(.))) %>%
      arrange(coefs.TimeMinutes) %>%
      filter(r.squareds >= rsq,
             coefs.TimeMinutes < 0 ) %>%
      filter(!is.na(sigmas)) %>% 
      mutate(MMR_MO2 = ((((vr - vf) * coefs.TimeMinutes)/vf)* -60),
             Interval = (i / 12),
             MO2_sd = sigmas)  %>%
      slice(which.max(MMR_MO2))
    
    
    MMR_MO2 <- c(MMR_MO2, MMR_coefs_df$MMR_MO2)
    IntervalLength <- c(IntervalLength, MMR_coefs_df$Interval[1])
    R_squareds <- c(R_squareds, MMR_coefs_df$r.squareds)
    RollEnd <- c(RollEnd, MMR_coefs_df$RollEnd)
    
  }
  
  HSID <- rep(hs, length(IntervalLength)) # The specimen ID will be coded as a column on the output data frame here
  # Combine all estimates into a data frame
  MMR_combined_df <- data.frame(HSID, MMR_MO2, IntervalLength, R_squareds, RollEnd)
  
  MMR_RR_out_df <- MMR_combined_df %>%
    mutate(MassKG = vf,
           ChamVol = vr,
           Window_mid = (RollEnd - (IntervalLength/2))) # Calculate mid-point of the MMR window
  
  return(MMR_RR_out_df)
}


####################################### ----------------------- ##############################################


#---- Function to plot basic MAXIMUM MR O2mgL/TimeMinutes for individual sharks
basic_MMR_plot <- function(df, hs, scaleyby, scalexby) {
  ggplot(data = df, aes(x = TimeMinutes, y = O2mgL)) +
    geom_point() + 
    ggtitle(hs)+
    scale_y_continuous(breaks = seq(0, 12, scaleyby)) + 
    scale_x_continuous(breaks = seq(0, 300, scalexby)) + theme_bw()
}


####################################### ----------------------- ##############################################


####################################### ----------------------- ##############################################

# Define salinity
HS_S <- 33.5
# Set amount of time to increase RR window width with each roll
RRstep <- 1
# Smallest interval width to test for HORN SHARKS Rolling Regression
smHS <- 1
lgHS <- 6


####################################### ----------------------- ##############################################



####################################### ----------------------- ##############################################

### ------------------------ Horn shark 36 ---------------------- #### 
vfHS36chase <- 0.0432 
vrHS36 <- 2.605
raw_HS36chase_df <- read_csv2("RawHornSharkCSV/BABYHORN3_MMR_CHASE_ONLY.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS36chase_trial_df <- F4_MMR_convert_DO(raw_HS36chase_df, lag = 8, end = 40, S = HS_S, interval = 5)

# Plot to check things look right
#basic_MMR_plot(HS36chase_trial_df, "HS36", 0.2, 2)

# Calculate MMR 
HS36_MMR_rollout_chase_df <- MMR_5sec_rolling_regression(HS36chase_trial_df, hs = "HS36", 
                                                         vr = vrHS36, vf = vfHS36chase, 
                                                         sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------ Horn shark 35 ---------------------- #### 
vfHS35chase <- 0.058
vrHS35 <- 2.585
raw_HS35chase_df <- read_csv2("RawHornSharkCSV/HARRIETA_2_MMR_CHASE_ONLY.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
                           
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS35chase_trial_df <- F4_MMR_convert_DO(raw_HS35chase_df, lag = 0, end = 40, S = HS_S)

# Plot to check things look right
 #basic_MMR_plot(HS35chase_trial_df, "HS35", 0.2, 5)

# Calculate MMR
HS35_MMR_rollout_chase_df <- MMR_rolling_regression(HS35chase_trial_df, hs = "HS35", 
                                                         vr = vrHS35, vf = vfHS35chase, 
                                                         sm = smHS, lg = lgHS, by = RRstep, rsq = 0.8) 



### ------------------------ Horn shark 34 ---------------------- #### 
vfHS34 <- 0.0387
vrHS34 <- 2.573
raw_HS34chase_df <- read_csv2("RawHornSharkCSV/HENRY_1.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Remove the RMR data from the top of the raw df
HS34chase_MMR_df <- tail(raw_HS34chase_df, -15280)

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS34chase_trial_df <- F4_MMR_convert_DO(HS34chase_MMR_df, lag = 0, end = 44, S = HS_S, interval = 5)

# Plot to check things look right
#basic_MMR_plot(HS34chase_trial_df, "HS34", 0.2, 1)

# Calculate MMR
HS34_MMR_rollout_chase_df <- MMR_5sec_rolling_regression(HS34chase_trial_df, hs = "HS34", 
                                                      vr = vrHS34, vf = vfHS34, 
                                                      sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------ Horn shark 33 ---------------------- #### 
vHS33 <- 3.38
vrHS33 <- 52.5
raw_HS33chase_df <- read_csv2("RawHornSharkCSV/MMR_HS33_1211_1CHASE.csv", skip = 1,
                         col_types = cols(Value = col_character(), patm = col_character(),
                                          Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS33chase_trial_df <- F4_MMR_convert_DO(raw_HS33chase_df, lag = 3.25, end = 17, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS33chase_trial_df, "HS33", 0.2, 1)

# Run rolling regression function to calculate MMR
HS33_MMR_rollout_chase_df <- MMR_rolling_regression(HS33chase_trial_df, hs = "HS33", 
                                              vr = vrHS33, vf = vHS33, 
                                              sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### --------------------- Horn shark 32  -----------------#####
vHS32 <- 1.58
vrHS32 <- 30.2
raw_HS32chase_df <- read_csv2("RawHornSharkCSV/MMR_HS32_1212_2CHASE.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS32chase_trial_df <- F4_MMR_convert_DO(raw_HS32chase_df, lag = 3, end = 16, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS32chase_trial_df, "HS32", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS32_MMR_rollout_chase_df <- MMR_rolling_regression(HS32chase_trial_df, hs = "HS32", 
                                              vr = vrHS32, vf = vHS32, 
                                              sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)



### ------------------------- Horn shark 31  -------------------------#####
vHS31 <- 3.82
vrHS31 <- 52.5
raw_HS31chase_df <- read_csv2("RawHornSharkCSV/MMR_HS31_1213_2CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS31chase_trial_df <- F4_MMR_convert_DO(raw_HS31chase_df, lag = 3.25, end = 17, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS31chase_trial_df, "HS31", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS31_MMR_rollout_chase_df <- MMR_rolling_regression(HS31chase_trial_df, hs = "HS31", 
                                                    vr = vrHS31, vf = vHS31, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)


### ------------------------- Horn shark 30  -------------------------#####
vHS30 <- 1.7
vrHS30 <- 30.2
raw_HS30chase_df <- read_csv2("RawHornSharkCSV/MMR_HS30_1203_1CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS30chase_trial_df <- F4_MMR_convert_DO(raw_HS30chase_df, lag = 4, end = 17, S = HS_S)

# Plot to check things look right
# basic_MMR_plot(HS30chase_trial_df, "HS30", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS30_MMR_rollout_chase_df <- MMR_rolling_regression(HS30chase_trial_df, hs = "HS30", 
                                                    vr = vrHS30, vf = vHS30, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




### ------------------------- Horn shark 27  -------------------------#####
vHS27 <- 2.3
vrHS27 <- 40
raw_HS27chase_df <- read_csv2("RawHornSharkCSV/MMR_HS27_627_2CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS27chase_trial_df <- F4_MMR_convert_DO(raw_HS27chase_df, lag = 7, end = 20, S = HS_S)

# Plot to check things look right
# basic_MMR_plot(HS27chase_trial_df, "HS27", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS27_MMR_rollout_chase_df <- MMR_rolling_regression(HS27chase_trial_df, hs = "HS27", 
                                                    vr = vrHS27, vf = vHS27, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 26  -------------------------#####
vHS26 <- 4.44
vrHS26 <- 52.5
raw_HS26chase_df <- read_csv2("RawHornSharkCSV/MMR_HS26_702_1CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS26chase_trial_df <- F4_MMR_convert_DO(raw_HS26chase_df, lag = 4.5, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS26chase_trial_df, "HS26", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS26_MMR_rollout_chase_df <- MMR_rolling_regression(HS26chase_trial_df, hs = "HS26", 
                                                    vr = vrHS26, vf = vHS26, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 22  -------------------------#####
vHS22 <- 0.697
vrHS22 <- 8.3
raw_HS22chase_df <- read_csv2("RawHornSharkCSV/MMR_190630_HS22_1chase.csv", skip = 57,
                              col_types = cols(`oxygen/% airsatur.` = col_character(),`temp/∞C` = col_character()))

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS22chase_trial_df <- F3_MMR_convert_DO(raw_HS22chase_df, lag = 4.5, end = 18, S = HS_S, patm = 1015)

# Plot to check things look right
#basic_MMR_plot(HS22chase_trial_df, "HS22", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS22_MMR_rollout_chase_df <- MMR_rolling_regression(HS22chase_trial_df, hs = "HS22", 
                                                    vr = vrHS22, vf = vHS22, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 21  -------------------------#####
vHS21 <- 0.439
vrHS21 <- 8.3
raw_HS21chase_df <- read_csv2("RawHornSharkCSV/MMR_HS21_628_2CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 

# 2. Filter out column we want to use, rename "Value" to "O2chase", change Temp and O2chase into decimal values
HS21chase_trial_df <- F4_MMR_convert_DO(raw_HS21chase_df, lag = 4.5, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS21chase_trial_df, "HS21", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS21_MMR_rollout_chase_df <- MMR_rolling_regression(HS21chase_trial_df, hs = "HS21", 
                                                    vr = vrHS21, vf = vHS21, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 20  -------------------------#####
vHS20 <- 0.715
vrHS20 <- 8.3
raw_HS20chase_df <- read_csv2("RawHornSharkCSV/MMR_190624_HS20_1chase.csv", skip = 57,
                              col_types = cols(`oxygen/% airsatur.` = col_character(),`temp/∞C` = col_character()))

# 2. Filter out column we want to use, rename "Value" to "O2chase", change Temp and O2chase into decimal values
HS20chase_trial_df <- F3_MMR_convert_DO(raw_HS20chase_df, lag = 5, end = 18, S = HS_S, patm = 1013)

# Plot to check things look right
#basic_MMR_plot(HS20chase_trial_df, "HS20", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS20_MMR_rollout_chase_df <- MMR_rolling_regression(HS20chase_trial_df, hs = "HS20", 
                                                    vr = vrHS20, vf = vHS20, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 


### ------------------------- Horn shark 19  -------------------------#####
vHS19 <- 0.738
vrHS19 <- 8.3
raw_HS19chase_df <- read_csv2("RawHornSharkCSV/MMR_HS19_627_2CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 

# 2. Filter out column we want to use, rename "Value" to "O2chase", change Temp and O2chase into decimal values
HS19chase_trial_df <- F4_MMR_convert_DO(raw_HS19chase_df, lag = 4.5, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS19chase_trial_df, "HS19", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS19_MMR_rollout_chase_df <- MMR_rolling_regression(HS19chase_trial_df, hs = "HS19", 
                                                    vr = vrHS19, vf = vHS19, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 


### ------------------------- Horn shark 18  -------------------------#####
vHS18 <- 3.02
vrHS18 <- 40
raw_HS18chase_df <- read_csv2("RawHornSharkCSV/MMR_190627_HS18_1chase.csv", skip = 57,
                              col_types = cols(`oxygen/% airsatur.` = col_character(),`temp/∞C` = col_character()))

# 2. Filter out column we want to use, rename "Value" to "O2chase", change Temp and O2chase into decimal values
HS18chase_trial_df <- F3_MMR_convert_DO(raw_HS18chase_df, lag = 4.5, end = 18, S = HS_S, patm = 1014)

# Plot to check things look right
#basic_MMR_plot(HS18chase_trial_df, "HS18", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS18_MMR_rollout_chase_df <- MMR_rolling_regression(HS18chase_trial_df, hs = "HS18", 
                                              vr = vrHS18, vf = vHS18, 
                                              sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)


### --------------------- Horn shark 17 ----------------------#####
vHS17 <- 2.46
vrHS17 <- 52.5
raw_HS17chase_df <- read_csv2("RawHornSharkCSV/MMR_HS17_704_3CHASE.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS17chase_trial_df <- F4_MMR_convert_DO(raw_HS17chase_df, lag = 4.5, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS17chase_trial_df, "HS17", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS17_MMR_rollout_chase_df <- MMR_rolling_regression(HS17chase_trial_df, hs = "HS17", 
                                                    vr = vrHS17, vf = vHS17, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)



### ------------------------- Horn shark 16  -------------------------#####
vHS16 <- 2.46
vrHS16 <- 30.2
raw_HS16chase_df <- read_csv2("RawHornSharkCSV/MMR_HS16_630_2CHASE.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS16chase_trial_df <- F4_MMR_convert_DO(raw_HS16chase_df, lag = 4.4, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS16chase_trial_df, "HS16", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS16_MMR_rollout_chase_df <- MMR_rolling_regression(HS16chase_trial_df, hs = "HS16", 
                                              vr = vrHS16, vf = vHS16, 
                                              sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 




### ------------------------- Horn shark 15  -------------------------#####
vHS15 <- 3.55
vrHS15 <- 40
raw_HS15chase_df <- read_csv2("RawHornSharkCSV/MMR_HS15_624_1.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS15chase_trial_df <- F4_MMR_convert_DO(raw_HS15chase_df, lag = 5, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS15chase_trial_df, "HS15", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS15_MMR_rollout_chase_df <- MMR_rolling_regression(HS15chase_trial_df, hs = "HS15", 
                                              vr = vrHS15, vf = vHS15, 
                                              sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




### ------------------------- Horn shark 14  -------------------------#####
vHS14 <- 0.893
vrHS14 <- 8.3
raw_HS14chase_df <- read_csv2("RawHornSharkCSV/MMR_190624_HS14_2chase.csv", skip = 57,
                            col_types = cols(`oxygen/% airsatur.` = col_character(),`temp/∞C` = col_character()))

# 2. Filter out column we want to use, rename "Value" to "O2air", change Temp and O2air into decimal values
HS14chase_trial_df <- F3_MMR_convert_DO(raw_HS14chase_df, lag = 5, end = 18, S = HS_S, patm = 1013)

# Plot to check things look right
#basic_MMR_plot(HS14chase_trial_df, "HS14", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS14_MMR_rollout_chase_df <- MMR_rolling_regression(HS14chase_trial_df, hs = "HS14", 
                                                  vr = vrHS14, vf = vHS14, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)



### ------------------------- Horn shark 13  -------------------------#####
vHS13 <- 0.469
vrHS13<- 8.3
raw_HS13chase_df <- read_csv2("RawHornSharkCSV/MMR_HS13_617_1.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS13chase_trial_df <- F4_MMR_convert_DO(raw_HS13chase_df, lag = 6, end = 19, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS13chase_trial_df, "HS13", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS13_MMR_rollout_chase_df <- MMR_rolling_regression(HS13chase_trial_df, hs = "HS13", 
                                                    vr = vrHS13, vf = vHS13, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 12  -------------------------#####
vHS12 <- 0.203
vrHS12 <- 5.825
raw_HS12chase_df <- read_csv2("RawHornSharkCSV/MMR_HS12_617_2.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS12chase_trial_df <- F4_MMR_convert_DO(raw_HS12chase_df, lag = 7, end = 20, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS12chase_trial_df, "HS12", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS12_MMR_rollout_chase_df <- MMR_rolling_regression(HS12chase_trial_df, hs = "HS12", 
                                                    vr = vrHS12, vf = vHS12, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




#####################################################################################
#####################################################################################

# -------------------------- CHASE + AIR MMR DATA SETS ---------------------------- #

#####################################################################################
#####################################################################################

### ------------------------ Horn shark 36 ---------------------- #### 
vfHS36air <- 0.0432 
vrHS36air <- 2.59
raw_HS36air_df <- read_csv2("RawHornSharkCSV/BABYHORN3_MMR_AIR.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS36air_trial_df <- F4_MMR_convert_DO(raw_HS36air_df, lag = 1, end = 40, S = HS_S, interval = 1)

# Plot to check things look right
#basic_MMR_plot(HS36air_trial_df, "HS36", 0.2, 2)

# Calculate MMR
HS36_MMR_rollout_air_df <- MMR_rolling_regression(HS36air_trial_df, hs = "HS36", 
                                                         vr = vrHS36air, vf = vfHS36air, 
                                                         sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 




### ------------------------ Horn shark 35 ---------------------- #### 
vfHS35 <- 0.058
vrHS35 <- 2.585
raw_HS35air_df <- read_csv2("RawHornSharkCSV/HARRIETA_MMR.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS35air_trial_df <- F4_MMR_convert_DO(raw_HS35air_df, lag = 0, end = 100, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS35air_trial_df, "HS35", 0.2, 5)

# Calculate MMR 
HS35_MMR_rollout_air_df <- MMR_rolling_regression(HS35air_trial_df, hs = "HS35", 
                                                    vr = vrHS35, vf = vfHS35, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.8) 



### ------------------------ Horn shark 34 ---------------------- #### 
vfHS34air <- 0.0387
vrHS34air <- 2.637
raw_HS34air_df <- read_csv2("RawHornSharkCSV/HENRY_3_MMR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS34air_trial_df <- F4_MMR_convert_DO(raw_HS34air_df, lag = 1, end = 40, S = HS_S)

# Plot to check things look right
# basic_MMR_plot(HS34air_trial_df, "HS34", 0.2, 2)

# Calculate MMR 
HS34_MMR_rollout_air_df <- MMR_rolling_regression(HS34air_trial_df, hs = "HS34", 
                                                  vr = vrHS34air, vf = vfHS34air, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------ Horn shark 33 ---------------------- #### 
vHS33 <- 3.38
vrHS33 <- 52.5
raw_HS33air_df <- read_csv2("RawHornSharkCSV/MMR_HS33_1213_2AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS33air_trial_df <- F4_MMR_convert_DO(raw_HS33air_df, lag = 1, end = 20, S = 3.5)

# Plot to check things look right
#basic_MMR_plot(HS33air_trial_df, "HS33", 0.2, 1)

# Run rolling regression function to calculate MMR
HS33_MMR_rollout_air_df <- MMR_rolling_regression(HS33air_trial_df, hs = "HS33", 
                                                    vr = vrHS33, vf = vHS33, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### --------------------- Horn shark 32  -----------------#####
vHS32 <- 1.58
vrHS32 <- 30.2
raw_HS32air_df <- read_csv2("RawHornSharkCSV/MMR_HS32_1210_1AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS32air_trial_df <- F4_MMR_convert_DO(raw_HS32air_df, lag = 1, end = 16, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS32air_trial_df, "HS32", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS32_MMR_rollout_air_df <- MMR_rolling_regression(HS32air_trial_df, hs = "HS32", 
                                                    vr = vrHS32, vf = vHS32, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




### ------------------------- Horn shark 31  -------------------------#####
vHS31 <- 3.82
vrHS31 <- 52.5
raw_HS31air_df <- read_csv2("RawHornSharkCSV/MMR_HS31_1216_3AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS31air_trial_df <- F4_MMR_convert_DO(raw_HS31air_df, lag = 1, end = 17, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS31air_trial_df, "HS31", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS31_MMR_rollout_air_df <- MMR_rolling_regression(HS31air_trial_df, hs = "HS31", 
                                                  vr = vrHS31, vf = vHS31, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




### ------------------------- Horn shark 30  -------------------------#####
vHS30 <- 1.7
vrHS30 <- 30.2
raw_HS30air_df <- read_csv2("RawHornSharkCSV/MMR_HS30_1209_3AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS30air_trial_df <- F4_MMR_convert_DO(raw_HS30air_df, lag = 1, end = 17, S = HS_S)

# Plot to check things look right
# basic_MMR_plot(HS30air_trial_df, "HS30", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS30_MMR_rollout_air_df <- MMR_rolling_regression(HS30air_trial_df, hs = "HS30", 
                                                    vr = vrHS30, vf = vHS30, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)


### ------------------------- Horn shark 27  -------------------------#####
vHS27 <- 2.3
vrHS27 <- 40
raw_HS27air_df <- read_csv2("RawHornSharkCSV/MMR_HS27_625_1AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS27air_trial_df <- F4_MMR_convert_DO(raw_HS27air_df, lag = 1, end = 20, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS27air_trial_df, "HS17", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS27_MMR_rollout_air_df <- MMR_rolling_regression(HS27air_trial_df, hs = "HS27", 
                                                  vr = vrHS27, vf = vHS27, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 




### ------------------------- Horn shark 26  -------------------------#####
vHS26 <- 4.44
vrHS26 <- 52.5
raw_HS26air_df <- read_csv2("RawHornSharkCSV/MMR_HS26_704_2AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS26air_trial_df <- F4_MMR_convert_DO(raw_HS26air_df, lag = 1, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS26air_trial_df, "HS26", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS26_MMR_rollout_air_df <- MMR_rolling_regression(HS26air_trial_df, hs = "HS26", 
                                                    vr = vrHS26, vf = vHS26, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 


### ------------------------- Horn shark 22  -------------------------#####
vHS22 <- 0.697
vrHS22 <- 8.3
raw_HS22air_df <- read_csv2("RawHornSharkCSV/MMR_HS22_702_2AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS22air_trial_df <- F4_MMR_convert_DO(raw_HS22air_df, lag = 1, end = 20, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS22air_trial_df, "HS22", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS22_MMR_rollout_air_df <- MMR_rolling_regression(HS22air_trial_df, hs = "HS22", 
                                                  vr = vrHS22, vf = vHS22, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 21  -------------------------#####
vHS21 <- 0.439
vrHS21 <- 8.3
raw_HS21air_df <- read_csv2("RawHornSharkCSV/MMR_190626_HS21_1air.csv", skip = 57,
                            col_types = cols(`oxygen/% airsatur.` = col_character(),
                                             `temp/∞C` = col_character()))

# 2. Filter out column we want to use, rename "Value" to "O2chase", change Temp and O2chase into decimal values
HS21air_trial_df <- F3_MMR_convert_DO(raw_HS21air_df, lag = 0, end = 40, S = HS_S, patm = 1014)

# Plot to check things look right
#basic_MMR_plot(HS21air_trial_df, "HS21", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS21_MMR_rollout_air_df <- MMR_rolling_regression(HS21air_trial_df, hs = "HS21", 
                                                  vr = vrHS21, vf = vHS21, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 




### ------------------------- Horn shark 20  -------------------------#####
vHS20 <- 0.715
vrHS20 <- 8.3
raw_HS20air_df <- read_csv2("RawHornSharkCSV/MMR_190626_HS20_2air.csv", skip = 57,
                            col_types = cols(`oxygen/% airsatur.` = col_character(),
                                             `temp/∞C` = col_character()))

# 2. Filter out column we want to use, rename "Value" to "O2chase", change Temp and O2chase into decimal values
HS20air_trial_df <- F3_MMR_convert_DO(raw_HS20air_df, lag = 1, end = 18, S = HS_S, patm = 1013)

# Plot to check things look right
#basic_MMR_plot(HS20air_trial_df, "HS20", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS20_MMR_rollout_air_df <- MMR_rolling_regression(HS20air_trial_df, hs = "HS20", 
                                                  vr = vrHS20, vf = vHS20, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 



### ------------------------- Horn shark 19  -------------------------#####
vHS19 <- 0.738
vrHS19 <- 8.3
raw_HS19air_df <- read_csv2("RawHornSharkCSV/MMR_190625_HS19_1air.csv", skip = 57,
                            col_types = cols(`oxygen/% airsatur.` = col_character(),
                                             `temp/∞C` = col_character()))
# 2. Filter out column we want to use, rename "Value" to "O2air", change Temp and O2air into decimal values
HS19air_trial_df <- F3_MMR_convert_DO(raw_HS19air_df, lag = 1, end = 18, S = HS_S, patm = 1013)

# Plot to check things look right
#basic_MMR_plot(HS19air_trial_df, "HS19", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS19_MMR_rollout_air_df <- MMR_rolling_regression(HS19air_trial_df, hs = "HS19", 
                                                  vr = vrHS19, vf = vHS19, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 




### ------------------------- Horn shark 18  -------------------------#####
vHS18 <- 3.02
vrHS18 <- 40
raw_HS18air_df <- read_csv2("RawHornSharkCSV/MMR_HS18_629_2AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 

# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS18air_trial_df <- F4_MMR_convert_DO(raw_HS18air_df, lag = 1, end = 20, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS18air_trial_df, "HS18", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS18_MMR_rollout_air_df <- MMR_rolling_regression(HS18air_trial_df, hs = "HS18", 
                                                  vr = vrHS18, vf = vHS18, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)



### --------------------- Horn shark 17 ----------------------#####
vHS17 <- 2.46
vrHS17air <- 40
raw_HS17air_df <- read_csv2("RawHornSharkCSV/MMR_HS17_629_1AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS17air_trial_df <- F4_MMR_convert_DO(raw_HS17air_df, lag = 1, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS17air_trial_df, "HS17", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS17_MMR_rollout_air_df <- MMR_rolling_regression(HS17air_trial_df, hs = "HS17", 
                                                  vr = vrHS17air, vf = vHS17, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)



### ------------------------- Horn shark 16  -------------------------#####
vHS16 <- 2.46
vrHS16 <- 30.2
raw_HS16air_df <- read_csv2("RawHornSharkCSV/MMR_HS16_628_1AIR.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS16air_trial_df <- F4_MMR_convert_DO(raw_HS16air_df, lag = 1, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS16air_trial_df, "HS16", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS16_MMR_rollout_air_df <- MMR_rolling_regression(HS16air_trial_df, hs = "HS16", 
                                                    vr = vrHS16, vf = vHS16, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 




### ------------------------- Horn shark 15  -------------------------#####
vHS15 <- 3.55
vrHS15air <- 52.5
raw_HS15air_df <- read_csv2("RawHornSharkCSV/MMR_HS15_626_2.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS15air_trial_df <- F4_MMR_convert_DO(raw_HS15air_df, lag = 1, end = 18, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS15air_trial_df, "HS15", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS15_MMR_rollout_air_df <- MMR_rolling_regression(HS15air_trial_df, hs = "HS15", 
                                                    vr = vrHS15air, vf = vHS15, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




### ------------------------- Horn shark 14  -------------------------#####
vHS14 <- 0.893
vrHS14air <- 16.4
raw_HS14air_df <- read_csv2("RawHornSharkCSV/MMR_190622_HS14_1.csv", skip = 57,
                            col_types = cols(`oxygen/% airsatur.` = col_character(),
                                             `temp/∞C` = col_character()))

# 2. Filter out column we want to use, rename "Value" to "O2air", change Temp and O2air into decimal values
HS14air_trial_df <- F3_MMR_convert_DO(raw_HS14air_df, lag = 1, end = 18, S = HS_S, patm = 1013)

# Plot to check things look right
#basic_MMR_plot(HS14air_trial_df, "HS14", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS14_MMR_rollout_air_df <- MMR_rolling_regression(HS14air_trial_df, hs = "HS14", 
                                                    vr = vrHS14air, vf = vHS14, 
                                                    sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)






### ------------------------- Horn shark 13  -------------------------#####
vHS13 <- 0.469
vrHS13<- 8.3
raw_HS13air_df <- read_csv2("RawHornSharkCSV/MMR_HS13_619_2.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS13air_trial_df <- F4_MMR_convert_DO(raw_HS13air_df, lag = 1, end = 19, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS13air_trial_df, "HS13", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS13_MMR_rollout_air_df <- MMR_rolling_regression(HS13air_trial_df, hs = "HS13", 
                                                  vr = vrHS13, vf = vHS13, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7) 


### ------------------------- Horn shark 12  -------------------------#####
vHS12 <- 0.203
vrHS12 <- 5.825
raw_HS12air_df <- read_csv2("RawHornSharkCSV/MMRHS12_613_1.csv", skip = 1,
                            col_types = cols(Value = col_character(), patm = col_character(),
                                             Temp = col_character())) 
# Clean the raw data, rename columns and convert %DO to mg/L oxygen
HS12air_trial_df <- F4_MMR_convert_DO(raw_HS12air_df, lag = 1, end = 20, S = HS_S)

# Plot to check things look right
#basic_MMR_plot(HS12air_trial_df, "HS12", 0.2, 1)

# Run rolling regression function to calculate MMR at 0.5-3 minutes, increasing by 10 seconds each time
HS12_MMR_rollout_air_df <- MMR_rolling_regression(HS12air_trial_df, hs = "HS12", 
                                                  vr = vrHS12, vf = vHS12, 
                                                  sm = smHS, lg = lgHS, by = RRstep, rsq = 0.7)




#######################################################################################################################################
######################  ----------------------------------------------------------------------------- ################################
#######################################################################################################################################

##### --------- Bind individual ** CHASE ** outputs into one df  ---------- ###
HS_CHASE_joined_df <- rbind(HS12_MMR_rollout_chase_df, HS13_MMR_rollout_chase_df, HS22_MMR_rollout_chase_df,
                            HS32_MMR_rollout_chase_df, HS30_MMR_rollout_chase_df, HS27_MMR_rollout_chase_df,
                            HS16_MMR_rollout_chase_df, HS17_MMR_rollout_chase_df, HS18_MMR_rollout_chase_df, 
                            HS33_MMR_rollout_chase_df, HS31_MMR_rollout_chase_df, HS26_MMR_rollout_chase_df, 
                            HS15_MMR_rollout_chase_df, HS19_MMR_rollout_chase_df, HS14_MMR_rollout_chase_df,
                            HS21_MMR_rollout_chase_df, HS20_MMR_rollout_chase_df, HS34_MMR_rollout_chase_df,
                            HS35_MMR_rollout_chase_df, HS36_MMR_rollout_chase_df) %>%
          mutate(Model = "MMR_chase",
                 MassG_log = (log10(MassKG * 1000)),
                 MassG = MassKG*1000,
                 MassKG_log = log10(MassKG),
                 WholeMMR_KG = (MMR_MO2 * MassKG),
                 WholeMMR_log_KG = log10(MMR_MO2 * MassKG),
                 Cham_fish_ratio = ChamVol/MassKG) %>%
         filter(IntervalLength == 2)



#######################################################################################################################################
#######################################################################################################################################

#######################################################################################################################################

##### --------- Bind individual ** AIR ** outputs into one df ---------- ###
HS_AIR_joined_df <- rbind(HS12_MMR_rollout_air_df, HS13_MMR_rollout_air_df, HS22_MMR_rollout_air_df,
                          HS32_MMR_rollout_air_df, HS30_MMR_rollout_air_df, HS27_MMR_rollout_air_df,
                          HS16_MMR_rollout_air_df, HS17_MMR_rollout_air_df, HS18_MMR_rollout_air_df, 
                          HS33_MMR_rollout_air_df, HS31_MMR_rollout_air_df, HS26_MMR_rollout_air_df, 
                          HS15_MMR_rollout_air_df, HS19_MMR_rollout_air_df, HS14_MMR_rollout_air_df,
                          HS21_MMR_rollout_air_df, HS20_MMR_rollout_air_df, HS34_MMR_rollout_air_df, 
                          HS35_MMR_rollout_air_df, HS36_MMR_rollout_air_df) %>%
  mutate(Model = "MMR_air",
         MassG_log = (log10(MassKG * 1000)),
         MassG = MassKG*1000,
         MassKG_log = log10(MassKG),
         WholeMMR_KG = (MMR_MO2 * MassKG),
         WholeMMR_log_KG = log10(MMR_MO2 * MassKG),
         Cham_fish_ratio = ChamVol/MassKG)%>%
  filter(IntervalLength == 2)









