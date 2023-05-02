

# Prinzing et al. 2023 Horn Shark metabolic rate and gill area scaling analysis
# Estimating RMR

# Load libraries
library("tidyverse")
library("rollRegres")
library("lubridate")
library("respR")


###########################################################################################################################
####################################### ----- FUNCTIONS FOR CALCULATING RMR  ------ ########################################
##############################################################################################################################



# respR -- FUNCTION to tidy raw FIBOX 4 df and convert %O2 to O2mgL, rename columns, 
#Filter columns to keep, add column of time in a more managable format, make cols into "double" numbers where needed
F4_RMR_trial_cleanup <- function(df, lag, end, S) { # Enter the raw data frame, the amount of lag time in minutes, and the salinity in mbar
  
  # Columns to keep
  keep_columns <- c("Date", "Time", "Value", "Temp", "patm") # Create object of columns to keep
  
  # Filter columns to keep, rename "Value" to "O2air", change Temp, patm and O2air into decimal values
  F4_trial_df <- dplyr::select(df, all_of(keep_columns)) %>% 
    mutate(Time = as.character(Time), patm = as.double(patm), 
           Value = as.double(Value), Temp = as.double(Temp)) %>%
    rename(O2air = Value) %>%
    slice(1:n()-1) %>%
    mutate(TimeMinutes = seq(0, by=5/60, length.out=nrow(.))) %>% # Make minutes column for plotting
    filter(TimeMinutes >= lag,
           TimeMinutes <= end) %>% # Filter out lag time when no shark was in the chamber
    mutate(TimeMinutes = seq(0, by=5/60, length.out=nrow(.)))
  
  # Convert %DO to mg/L oxygen with convert DO function
    converted_df <- F4_trial_df %>%
      rowwise() %>%
      mutate(O2mgL = convert_DO(x=O2air, from = "%Air", to = "mg/L", 
                S = S, t = Temp, P = (patm/1000)))
  
  # converted_df <- F4_trial_df %>%
  #   rowwise() %>%
  #   dplyr::mutate(O2mgL = ((respR::convert_DO(x=100, from = "%Air", to = "mg/L", S = S, 
  #                                      t = Temp, P = (patm/1000))[[2]]) * (O2air/100)))
  
  return(converted_df)
  
  
}



###### ----- respR -- FUNCTION to tidy and convert FIBOX 3 MMR data ------ ##############
F3_RMR_trial_cleanup <- function(df, lag, end, S, patm) {
  # Tidy up the raw data set
  MO2_calculation_df <- df %>% 
    rename(Time = `time/hh:mm:ss`, LogTime = `logtime/min`,
           O2air = `oxygen/% airsatur.`,Temp = `temp/∞C`, Amp = amp, Phase = `phase/∞`)  %>% 
    dplyr::select(-ErrorMessage) %>%
    mutate(TimeMinutes = as.double(LogTime), O2air = as.double(O2air), Temp = as.double(Temp)) %>%
    filter(TimeMinutes >= lag, TimeMinutes <= end) %>%
    mutate(TimeMinutes = seq(0, by=5/60, length.out=nrow(.)))
  
  
  # Convert %DO to mg/L oxygen with convert DO function
  converted_df <- MO2_calculation_df %>%
    rowwise() %>%
      mutate(O2mgL = convert_DO(x=O2air, from = "%Air", to = "mg/L", 
                S = S, t = Temp, P = (patm/1000)))
  
  return(converted_df)
  
}


#---- Function to plot basic RESTING MR O2mgL/TimeMinutes for individual sharks
basic_RMR_plot <- function(df, hs, scaleyby, scalexby) {
  
  plot <- ggplot(data = df, aes(x = TimeMinutes, y = O2mgL)) +
    geom_point() + 
    ggtitle(hs)+
    scale_y_continuous(breaks = seq(0, 12, scaleyby)) + 
    scale_x_continuous(breaks = seq(0, 2000, scalexby)) + theme_bw()
  
  return(plot)
}



# --- FUNCTION TO CALCULATE RMR ---- #
# Rolling regression with the width of the measurement window as the width of the regressions
HS_RMR_mgL_hour <- function(df, flush, close, wait, measure, vf, vr, hs, rsqrd, 
                            low_n = 2, high_n = 10, acclimation = 0) {
  
  # Total time of one measurement cycle
  ClosePeriod <- (flush * 12L) + (close *12L) # should be called MeasePeriod or something but too late
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts <- (close*12L) - (wait * 12L)
  # Measurement period
  windowwidth <- (measure * 12L)
  mean_loop_vec <- seq(low_n, high_n, 1)
  MeanOfN <- c()
  MeanSlope <- c()
  
  # --- Rolling regression with the width of the measurement window as the width of the regressions
  roll_df <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df, width = windowwidth, 
                                       do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df <- roll_df %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df <- RMR_logtime_df %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod)) %>%
    slice_tail(n = closedEsts) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)%>%
    filter(LogTime >= (acclimation*3600)) %>%
    arrange(desc(coefs.TimeMinutes))  # Arrange slopes so lowest RMR slope is first
  
  # Take mean of lowest n slope values, loop to get a sample of what that mean is using different numbers of slopes
  for(i in mean_loop_vec) {
    
    rmr_loop_df <- rmr_rr_close_df %>%
      head(i) # Take lowest n estimates
    
    MeanSlopeVal <- mean(rmr_loop_df$coefs.TimeMinutes)
    
    # Mean low slope of the slopes kept from above
    MeanSlope <- c(MeanSlope, MeanSlopeVal)
    # Number of slopes used to calulate the mean
    MeanOfN <- c(MeanOfN, i)
    
  }
  
  HSID <- rep(hs, length(MeanSlope)) # The specimen ID will be coded as a column on the output data frame here
  # Calculate MO2 per kg per hour
  RMR_kghr_df <- data.frame(HSID) %>%
    mutate(MassKG = vf,
           RMR = (((vr - vf) * MeanSlope) / vf) * -60,# Negative b/c want a positive RMR estimate
           MeanOfN = MeanOfN, 
           MeanSlope = MeanSlope,
           ChamVol = vr)
  
  # Return the slope estimate for that shark
  return(RMR_kghr_df)
  
}



# --- FUNCTION TO CALCULATE RMR ** TWO DATA FRAMES ---- #
# Rolling regression with the width of the measurement window as the width of the regressions
HS_RMR_2df_mgL_hour <- function(df1, df2, df3, flush1, close1, flush2, close2,
                                wait, measure, vf, vr, hs, rsqrd, low_n = 2, high_n = 10) {
  
  # Total time of one measurement cycle
  ClosePeriod1 <- (flush1 * 12L) + (close1 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts1 <- (close1*12L) - (wait * 12L)
  # Total time of one measurement cycle
  ClosePeriod2 <- (flush2 * 12L) + (close2 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts2 <- (close2*12L) - (wait * 12L)
  # Measurement period
  windowwidth <- (measure * 12L)
  mean_loop_vec <- seq(low_n, high_n, 1)
  MeanOfN <- c()
  MeanSlope <- c()
  
  # - ONE -- Rolling regression with the width of the measurement window as the width of the regressions
  roll_df1 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df1, width = windowwidth, 
                                        do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df1 <- roll_df1 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df1) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df1 <- RMR_logtime_df1 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod1)) %>%
    slice_tail(n = closedEsts1) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  # -- TWO - Rolling regression with the width of the measurement window as the width of the regressions
  roll_df2 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df2, width = windowwidth, 
                                        do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df2 <- roll_df2 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df2) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df2 <- RMR_logtime_df2 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod2)) %>%
    slice_tail(n = closedEsts2) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  
  # Join all the slope estimates from all the rolling regression outputs
  RR_out_slopes_all_df <- rbind(rmr_rr_close_df1, rmr_rr_close_df2) 
  RR_out_slopes_arranged_df <- RR_out_slopes_all_df %>%
    arrange(desc(coefs.TimeMinutes))  # Arrange slopes so lowest RMR slope is first
  
  ################### --------------------- ####################
  
  # Take mean of lowest n slope values, loop to get a sample of what that mean is using different numbers of slopes
  for(i in mean_loop_vec) {
    
    rmr_loop_df <- RR_out_slopes_arranged_df %>%
      head(i) # Take lowest n estimates
    # Mean of lowest n estimates
    MeanSlopeVal <- mean(rmr_loop_df$coefs.TimeMinutes)
    # Mean low slope of the slopes kept from above, added to growing vector
    MeanSlope <- c(MeanSlope, MeanSlopeVal)
    # Number of slopes used to calulate the mean
    MeanOfN <- c(MeanOfN, i)
    
  }
  
  HSID <- rep(hs, length(MeanSlope)) # The specimen ID will be coded as a column on the output data frame here
  # Calculate MO2 
  RMR_kghr_df <- data.frame(HSID) %>%
    mutate(MassKG = vf,
           RMR = (((vr - vf) * MeanSlope) / vf) * -60,# Negative b/c want a positive RMR estimate
           MeanOfN = MeanOfN, 
           MeanSlope = MeanSlope,
           ChamVol = vr)
  
  # Return the slope estimate for that shark
  return(RMR_kghr_df)
  
}



# --- FUNCTION TO CALCULATE RMR ** THREE DATA FRAMES ---- #
# Rolling regression with the width of the measurement window as the width of the regressions
HS_RMR_3df_mgL_hour <- function(df1, df2, df3, flush1, close1, flush2, close2, flush3, close3 ,
                                wait, measure, vf, vr, hs, rsqrd, 
                                low_n = 2, high_n = 10) {
  
  # Total time of one measurement cycle
  ClosePeriod1 <- (flush1 * 12L) + (close1 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts1 <- (close1*12L) - (wait * 12L)
  # Total time of one measurement cycle
  ClosePeriod2 <- (flush2 * 12L) + (close2 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts2 <- (close2*12L) - (wait * 12L)
  # Total time of one measurement cycle
  ClosePeriod3 <- (flush3 * 12L) + (close3 *12L)
  # Time of the closed cycle minus the "wait" or lag period
  closedEsts3 <- (close3*12L) - (wait * 12L)
  # Measurement period
  windowwidth <- (measure * 12L)
  mean_loop_vec <- seq(low_n, high_n, 1)
  MeanOfN <- c()
  MeanSlope <- c()
  
  # - ONE -- Rolling regression with the width of the measurement window as the width of the regressions
  roll_df1 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df1, width = windowwidth, 
                                       do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df1 <- roll_df1 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df1) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df1 <- RMR_logtime_df1 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod1)) %>%
    slice_tail(n = closedEsts1) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  # -- TWO - Rolling regression with the width of the measurement window as the width of the regressions
  roll_df2 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df2, width = windowwidth, 
                                       do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df2 <- roll_df2 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df2) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df2 <- RMR_logtime_df2 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod2)) %>%
    slice_tail(n = closedEsts2) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  # -- THREE - Rolling regression with the width of the measurement window as the width of the regressions
  roll_df3 <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, df3, width = windowwidth, 
                                        do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))
  
  # Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
  # --> LogTime indicates the end time of the regression window
  RMR_logtime_df3 <- roll_df3 %>%
    mutate(RowNum = row_number()) %>%
    mutate(LogTime = seq(0, ((nrow(roll_df3) * 5) - 1), 5))
  
  # Filter the slopes df to keep the estimates at the desired time point within each closed period
  rmr_rr_close_df3 <- RMR_logtime_df3 %>%
    group_by(grp = rep(row_number(), length.out = n(), each = ClosePeriod3)) %>%
    slice_tail(n = closedEsts3) %>% # Keep the closed period estimates minus the lag period
    slice_head(n = windowwidth) %>% # Keep the first n estimates from what was just selected
    slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at fluse+wait+measure
    filter(coefs.TimeMinutes < 0,
           r.squareds >= rsqrd)
  
  
  # Join all the slope estimates from all the rolling regression outputs
  RR_out_slopes_all_df <- rbind(rmr_rr_close_df1, rmr_rr_close_df2, rmr_rr_close_df3) 
  RR_out_slopes_arranged_df <- RR_out_slopes_all_df %>%
                        arrange(desc(coefs.TimeMinutes))  # Arrange slopes so lowest RMR slope is first

  ################### --------------------- ####################
  
  # Take mean of lowest n slope values, loop to get a sample of what that mean is using different numbers of slopes
  for(i in mean_loop_vec) {
    
    rmr_loop_df <- RR_out_slopes_arranged_df %>%
      head(i) # Take lowest n estimates
    # Mean of lowest n slope values
    MeanSlopeVal <- mean(rmr_loop_df$coefs.TimeMinutes)
    # Mean low slope of the slopes kept from above, added to growing vector
    MeanSlope <- c(MeanSlope, MeanSlopeVal)
    # Number of slopes used to calulate the mean
    MeanOfN <- c(MeanOfN, i)
    
  }
  
  HSID <- rep(hs, length(MeanSlope)) # The specimen ID will be coded as a column on the output data frame here
  # Calculate MO2 
  RMR_kghr_df <- data.frame(HSID) %>%
    mutate(MassKG = vf,
           RMR = (((vr - vf) * MeanSlope) / vf) * -60,# Negative b/c want a positive RMR estimate
           MeanOfN = MeanOfN, 
           MeanSlope = MeanSlope,
           ChamVol = vr)
  
  # Return the slope estimate for that shark
  return(RMR_kghr_df)
  
}




# --- FUNCTION to tidy Fibox 4 .cvs MMR files, remane columns, filter lag time, and convert %DO into mg/L oxygen --- #
F4_br_convert_DO <- function(df, lag, end, S, interval = 5) { # Enter the raw data frame, the amount of lag time in minutes, and the salinity in mbar
  
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


##############################################################################################################################
##############################################################################################################################




#### --- Pre-sets for all sharks -- #####
# Define salinity
HS_S <- 33.5
# Wait period (lag period) in minutes
LagPeriod <- 3
# R-Squared for filtering slopes
HS_RSqrd <- 0.85
# Measurement window length
HSmeasure <- 7
# Lowest and highest number of measurement periods to take mean RMR from
HSlow_n = 2
HShigh_n = 10


##############################################################################################################################



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 12 ----- ##


# Chamber volume
vrHS12 <- 5.825
# Fish mass (mass = volume)
vfHS12 <- 0.203
# FLush cycle length 
FlushHS12_1 <- 10
# Close cycle length
CloseHS12_1 <- 11
# FLush cycle length 
FlushHS12_2 <- 10
# Close cycle length
CloseHS12_2 <- 20
# FLush cycle length 
FlushHS12_3 <- 8
# Close cycle length
CloseHS12_3 <- 30

# Read in MMR df 1
raw_HS12_1_RMR_df <- read_csv2("RawHornSharkCSV/RMRHS12_613_INTER1.csv", skip = 1,
                       col_types = cols(Value = col_character(), patm = col_character(),
                                        Temp = col_character())) 
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS12_1_converted_df <- F4_RMR_trial_cleanup(df = raw_HS12_1_RMR_df, lag = 0, end = 47, S = HS_S) 

# Read in MMR df 2
raw_HS12_2_RMR_df<- read_csv2("RawHornSharkCSV/RMRHS12_613_INTER2.csv", skip = 1,
             col_types = cols(Value = col_character(), patm = col_character(),
                              Temp = col_character())) 
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS12_2_converted_df <- F4_RMR_trial_cleanup(df = raw_HS12_2_RMR_df, lag = 0, end = 180, S = HS_S)

# Read in MMR df 3
raw_HS12_3_RMR_df <- read_csv2("RawHornSharkCSV/RMRHS12_613_INTER3.csv", skip = 1,
                        col_types = cols(Value = col_character(), patm = col_character(),
                                         Temp = col_character())) 
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS12_3_converted_df <- F4_RMR_trial_cleanup(df = raw_HS12_3_RMR_df, lag = 0, end = 76, S = HS_S)

# Plot the converted df
# basic_RMR_plot(df = HS12_2_converted_df, hs = "HS12", scaleyby = 0.1, scalexby = (10))

# Calculate RMR - mean of n lowest slope estimates
HS12_RMR_est_df <- HS_RMR_3df_mgL_hour(df1 = HS12_1_converted_df, flush1 = FlushHS12_1, close1 = CloseHS12_1, 
                                       df2 = HS12_2_converted_df, flush2 = FlushHS12_2, close2 = CloseHS12_2, 
                                       df3 = HS12_3_converted_df, flush3 = FlushHS12_3, close3 = CloseHS12_3, 
                                       wait = LagPeriod, measure = HSmeasure,  
                                       vf = vfHS12, vr = vrHS12, hs = "HS12", rsqrd = HS_RSqrd, 
                                       low_n = HSlow_n, high_n = HShigh_n) # Already a 12-hour acclimation in there




########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 13 ----- ##

# Chamber volume
vrHS13 <- 8.3
# Fish mass (mass = volume)
vfHS13 <- 0.469
# FLush cycle length 
FlushHS13_1 <- 20
# Close cycle length
CloseHS13_1 <- 15
# FLush cycle length 
FlushHS13_2 <- 20
# Close cycle length
CloseHS13_2 <- 20

# Read in MMR df 1
raw_HS13_1_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS13_616_INTER1.csv", skip = 1,
                               col_types = cols(Value = col_character(), patm = col_character(),
                                                Temp = col_character())) 
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS13_1_converted_df <- F4_RMR_trial_cleanup(df = raw_HS13_1_RMR_df, lag = 665, end = 1225, S = HS_S) 

# Read in MMR df 2
raw_HS13_2_RMR_df<- read_csv2("RawHornSharkCSV/RMR_HS13_617_INTER2.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS13_2_converted_df <- F4_RMR_trial_cleanup(df = raw_HS13_2_RMR_df, lag = 0, end = 120, S = HS_S)

# Plot the converted df
basic_RMR_plot(df = HS13_2_converted_df, hs = "HS13", scaleyby = 0.1, scalexby = (40))

# Calculate RMR - mean of n lowest slope estimates
HS13_RMR_est_df <- HS_RMR_2df_mgL_hour(df1 = HS13_1_converted_df, flush1 = FlushHS13_1, close1 = CloseHS13_1, 
                                   df2 = HS13_2_converted_df, flush2 = FlushHS13_2, close2 = CloseHS13_2,
                                   wait = LagPeriod, measure = HSmeasure,
                                   vf = vfHS13, vr = vrHS13, hs = "HS13", rsqrd = HS_RSqrd, 
                                   low_n = HSlow_n, high_n = HShigh_n)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 14 ----- ##

# Air Pressure - patm
HS14pressure <- 1013
# Chamber volume
vrHS14 <- 16.4
# Fish mass (mass = volume)
vfHS14 <- 0.893
# FLush cycle length 
FlushHS14 <- 15
# Close cycle length
CloseHS14 <- 19

# Read in intermittent flow resprometry data data file 
raw_HS14_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190621_HS14_inter1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS14_converted_df <- F3_RMR_trial_cleanup(df = raw_HS14_RMR_df, lag = 0, end = 1020, patm =  HS14pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS14_converted_df, hs = "HS14", scaleyby = 0.1, scalexby = (FlushHS14 +CloseHS14))

# Calculate RMR - mean of n lowest slope estimates
HS14_RMR_est_df <- HS_RMR_mgL_hour(df = HS14_converted_df, flush = FlushHS14, close = CloseHS14, wait = LagPeriod, 
                                   measure = HSmeasure, vf = vfHS14, vr = vrHS14, 
                                   hs = "HS14", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 7)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 15 ----- ##


# Chamber volume
vrHS15 <- 40
# Fish mass (mass = volume)
vfHS15 <- 3.55
# FLush cycle length 
FlushHS15 <- 25
# Close cycle length 
CloseHS15 <- 12

# Read in intermittent flow resprometry data data file 
raw_HS15_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS15_623_INTER1.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS15_converted_df <- F4_RMR_trial_cleanup(df = raw_HS15_RMR_df, lag = 0, end = 999, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS15_converted_df, hs = "HS15", scaleyby = 0.1, scalexby = (FlushHS15 +CloseHS15))

# Calculate RMR - mean of n lowest slope estimates
HS15_RMR_est_df <- HS_RMR_mgL_hour(df = HS15_converted_df, flush = FlushHS15, close = CloseHS15, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS15, vr = vrHS15, hs = "HS15", 
                                   rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 8.5)


########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 16 ----- ##

# Air Pressure - patm
HS16pressure <- 1015
# Chamber volume
vrHS16 <- 30.2
# Fish mass (mass = volume)
vfHS16 <- 2.46
# FLush cycle length 
FlushHS16 <- 15
# Close cycle length
CloseHS16 <- 15

# Read in intermittent flow resprometry data data file 
raw_HS16_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190627_HS16_inter2.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS16_converted_df <- F3_RMR_trial_cleanup(df = raw_HS16_RMR_df, lag = 0, end = 1080, patm =  HS16pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS16_converted_df, hs = "HS16", scaleyby = 0.1, scalexby = (FlushHS16 + CloseHS16))

# Calculate RMR - mean of n lowest slope estimates
HS16_RMR_est_df <- HS_RMR_mgL_hour(df = HS16_converted_df, flush = FlushHS16, close = CloseHS16, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS16, vr = vrHS16, 
                                   hs = "HS16", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 10)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 17 ----- ##


# Chamber volume
vrHS17 <- 40
# Fish mass (mass = volume)
vfHS17 <- 2.46
# FLush cycle length 
FlushHS17 <- 15
# Close cycle length 
CloseHS17 <- 18

# Read in intermittent flow resprometry data data file 
raw_HS17_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS17_628_INTER1.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS17_converted_df <- F4_RMR_trial_cleanup(df = raw_HS17_RMR_df, lag = 0, end = 1023, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS17_converted_df, hs = "HS17", scaleyby = 0.1, scalexby = (FlushHS17 +CloseHS17))

# Calculate RMR - mean of n lowest slope estimates
HS17_RMR_est_df <- HS_RMR_mgL_hour(df = HS17_converted_df, flush = FlushHS17, close = CloseHS17, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS17, vr = vrHS17, 
                                   hs = "HS17", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 11)






########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 18 ----- ##

# Air Pressure - patm
HS18pressure <- 1013
# Chamber volume
vrHS18 <- 40
# Fish mass (mass = volume)
vfHS18 <- 3.02
# FLush cycle length 
FlushHS18 <- 15
# Close cycle length
CloseHS18 <- 25

# Read in intermittent flow resprometry data data file 
raw_HS18_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190626_HS18_inter1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS18_converted_df <- F3_RMR_trial_cleanup(df = raw_HS18_RMR_df, lag = 0, end = 1040, patm =  HS18pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS18_converted_df, hs = "HS18", scaleyby = 0.1, scalexby = (FlushHS18 +CloseHS18))

# Calculate RMR - mean of n lowest slope estimates
HS18_RMR_est_df <- HS_RMR_mgL_hour(df = HS18_converted_df, flush = FlushHS18, close = CloseHS18, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS18, vr = vrHS18, hs = "HS18", 
                                   rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 10)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 19 ----- ##

# Air Pressure - patm
HS19pressure <- 1013
# Chamber volume
vrHS19 <- 8.3
# Fish mass (mass = volume)
vfHS19 <- 0.738
# FLush cycle length 
FlushHS19 <- 15
# Close cycle length
CloseHS19 <- 20

# Read in intermittent flow resprometry data data file 
raw_HS19_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190624_HS19_inter1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS19_converted_df <- F3_RMR_trial_cleanup(df = raw_HS19_RMR_df, lag = 0, end = 980, patm =  HS19pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS19_converted_df, hs = "HS19", scaleyby = 0.1, scalexby = (FlushHS19 +CloseHS19))

# Calculate RMR - mean of n lowest slope estimates
HS19_RMR_est_df <- HS_RMR_mgL_hour(df = HS19_converted_df, flush = FlushHS19, close = CloseHS19, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS19, vr = vrHS19, 
                                   hs = "HS19", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 5)


########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 20 ----- ##

# Air Pressure - patm
HS20pressure <- 1015
# Chamber volume
vrHS20 <- 8.3
# Fish mass (mass = volume)
vfHS20 <- 0.715
# FLush cycle length 
FlushHS20 <- 10
# Close cycle length
CloseHS20 <- 10

# Read in intermittent flow resprometry data data file 
raw_HS20_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190623_HS20_inter1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS20_converted_df <- F3_RMR_trial_cleanup(df = raw_HS20_RMR_df, lag = 0, end = 1000, patm =  HS20pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS20_converted_df, hs = "HS20", scaleyby = 0.1, scalexby = (FlushHS20 +CloseHS20))

# Calculate RMR - mean of n lowest slope estimates
HS20_RMR_est_df <- HS_RMR_mgL_hour(df = HS20_converted_df, flush = FlushHS20, close = CloseHS20, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS20, vr = vrHS20, hs = "HS20", 
                                   rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 12)





########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 21 ----- ##

# Air Pressure - patm
HS21pressure <- 1012
# Chamber volume
vrHS21 <- 8.3
# Fish mass (mass = volume)
vfHS21 <- 0.439
# FLush cycle length 
FlushHS21 <- 10
# Close cycle length
CloseHS21 <- 20

# Read in intermittent flow resprometry data data file 
raw_HS21_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190625_HS21_intermit1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS21_converted_df <- F3_RMR_trial_cleanup(df = raw_HS21_RMR_df, lag = 28, end = 1018, patm =  HS21pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS21_converted_df, hs = "HS21", scaleyby = 0.1, scalexby = (FlushHS21 +CloseHS21))

# Calculate RMR - mean of n lowest slope estimates
HS21_RMR_est_df <- HS_RMR_mgL_hour(df = HS21_converted_df, flush = FlushHS21, close = CloseHS21, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS21, vr = vrHS21, hs = "HS21", 
                                   rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 10)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 22 ----- ##

# Air Pressure - patm
HS22pressure <- 1015
# Chamber volume
vrHS22 <- 8.3
# Fish mass (mass = volume)
vfHS22 <- 0.697
# FLush cycle length 
FlushHS22 <- 15
# Close cycle length
CloseHS22 <- 15

# Read in intermittent flow resprometry data data file 
raw_HS22_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190629_HS22_inter1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS22_converted_df <- F3_RMR_trial_cleanup(df = raw_HS22_RMR_df, lag = 0, end = 1170, patm =  HS22pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS22_converted_df, hs = "HS22", scaleyby = 0.1, scalexby = (FlushHS22 +CloseHS22))

# Calculate RMR - mean of n lowest slope estimates
HS22_RMR_est_df <- HS_RMR_mgL_hour(df = HS22_converted_df, flush = FlushHS22, close = CloseHS22, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS22, vr = vrHS22, hs = "HS22", 
                                   rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 12)




########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 26 ----- ##

# Air Pressure - patm
HS26pressure <- 1013
# Chamber volume
vrHS26 <- 52.5
# Fish mass (mass = volume)
vfHS26 <- 4.44
# FLush cycle length 
FlushHS26_1 <- 15
# Close cycle length
CloseHS26_1 <- 20
# FLush cycle length 
FlushHS26_2 <- 20
# Close cycle length
CloseHS26_2 <- 15

# Read in MMR df 1
raw_HS26_1_RMR_df <- read_csv2("RawHornSharkCSV/RMR_190701_HS26_intermit1.csv", skip = 57,
                               col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                                `temp/∞C` = col_character()))
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS26_1_converted_df <- F3_RMR_trial_cleanup(df = raw_HS26_1_RMR_df, lag = 595, end = 735, patm =  HS26pressure, S = HS_S) 

# Read in MMR df 2
raw_HS26_2_RMR_df<- read_csv2("RawHornSharkCSV/RMR_190702_HS26_intermit2.csv",skip = 57,
                              col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                               `temp/∞C` = col_character()))
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS26_2_converted_df <- F3_RMR_trial_cleanup(df = raw_HS26_2_RMR_df, lag = 0, end = 245, patm =  HS26pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS26_2_converted_df, hs = "HS26", scaleyby = 0.1, scalexby = (35))

# Calculate RMR - mean of n lowest slope estimates
HS26_RMR_est_df <- HS_RMR_2df_mgL_hour(df1 = HS26_1_converted_df, flush1 = FlushHS26_1, close1 = CloseHS26_1, 
                                       df2 = HS26_2_converted_df, flush2 = FlushHS26_2, close2 = CloseHS26_2,
                                       wait = LagPeriod, measure = HSmeasure, vf = vfHS26, 
                                       vr = vrHS26, hs = "HS26", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n)


########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 27 ----- ##

# Chamber volume
vrHS27 <- 40
# Fish mass (mass = volume)
vfHS27 <- 2.3
# FLush cycle length 
FlushHS27_1 <- 15
# Close cycle length
CloseHS27_1 <- 20
# FLush cycle length 
FlushHS27_2 <- 10
# Close cycle length
CloseHS27_2 <- 25

# "test" is saved as an updated .csv format to match the function 
# The original format CSV refused to run for some reason
raw_HS27_1_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS27_624_INTER1_test.csv", skip = 1,
                               col_types = cols(Value = col_character(), patm = col_character(),
                                                Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS27_1_converted_df <- F4_RMR_trial_cleanup(df = raw_HS27_1_RMR_df, lag = 525, end = 910, S = HS_S) 

# Read in MMR df 2
raw_HS27_2_RMR_df<- read_csv2("RawHornSharkCSV/RMR_HS27_625_INTER2.csv", skip = 1,
                              col_types = cols(Value = col_character(), patm = col_character(),
                                               Temp = col_character())) 
# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS27_2_converted_df <- F4_RMR_trial_cleanup(df = raw_HS27_2_RMR_df, lag = 0, end = 315, S = HS_S)

# Plot the converted df
#basic_RMR_plot(df = HS27_2_converted_df, hs = "HS27", scaleyby = 0.1, scalexby = (35))

# Calculate RMR - mean of n lowest slope estimates
HS27_RMR_est_df <- HS_RMR_2df_mgL_hour(df1 = HS27_1_converted_df, flush1 = FlushHS27_1, close1 = CloseHS27_1, 
                                       df2 = HS27_2_converted_df, flush2 = FlushHS27_2, close2 = CloseHS27_2,
                                       wait = LagPeriod, measure = HSmeasure,  
                                       vf = vfHS27, vr = vrHS27, hs = "HS27", rsqrd = HS_RSqrd, 
                                       low_n = HSlow_n, high_n = HShigh_n)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 30 ----- ##

# Air Pressure - patm
HS30pressure <- 1018
# Chamber volume
vrHS30 <- 30.2
# Fish mass (mass = volume)
vfHS30 <- 1.7
# FLush cycle length 
FlushHS30 <- 20
# Close cycle length
CloseHS30 <- 10

# Read in intermittent flow resprometry data data file 
raw_HS30_RMR_df <- read_csv2("RawHornSharkCSV/RMR_191202_HS30_intermit1.csv", skip = 57,
                             col_types = cols(`oxygen/% airsatur.` = col_character(), 
                                              `temp/∞C` = col_character()))

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS30_converted_df <- F3_RMR_trial_cleanup(df = raw_HS30_RMR_df, lag = 690, end = 990, patm =  HS30pressure, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS30_converted_df, hs = "HS30", scaleyby = 0.1, scalexby = 10)

# Calculate RMR - mean of n lowest slope estimates
HS30_RMR_est_df <- HS_RMR_mgL_hour(df = HS30_converted_df, flush = FlushHS30, close = CloseHS30, wait = LagPeriod, 
                                   measure = HSmeasure, vf = vfHS30,  vr = vrHS30, hs = "HS30", 
                                   rsqrd = 0.1, low_n = HSlow_n, high_n = HShigh_n)



########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 31 ----- ##

# Chamber volume
vrHS31 <- 52.5
# Fish mass (mass = volume)
vfHS31 <- 3.82
# FLush cycle length 
FlushHS31 <- 25
# Close cycle length 
CloseHS31 <- 10

# Read in intermittent flow resprometry data data file 
raw_HS31_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS31_1212_2.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS31_converted_df <- F4_RMR_trial_cleanup(df = raw_HS31_RMR_df, lag = 0, end = 980, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS31_converted_df, hs = "HS31", scaleyby = 0.1, scalexby = (FlushHS31 +CloseHS31))

# Calculate RMR - mean of n lowest slope estimates
HS31_RMR_est_df <- HS_RMR_mgL_hour(df = HS31_converted_df, flush = FlushHS31, close = CloseHS31, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS31, vr = vrHS31, 
                                   hs = "HS31", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 11)






########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 32 ----- ##

# Chamber volume
vrHS32 <- 30.2
# Fish mass (mass = volume)
vfHS32 <- 1.58
# FLush cycle length 
FlushHS32 <- 20
# Close cycle length 
CloseHS32 <- 15

# Read in intermittent flow resprometry data data file 
raw_HS32_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS32_1209_1.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS32_converted_df <- F4_RMR_trial_cleanup(df = raw_HS32_RMR_df, lag = 0, end = 1155, S = HS_S) 

# Plot the converted df
 #basic_RMR_plot(df = HS32_converted_df, hs = "HS32", scaleyby = 0.1, scalexby = (FlushHS32 +CloseHS32))

# Calculate RMR - mean of n lowest slope estimates
HS32_RMR_est_df <- HS_RMR_mgL_hour(df = HS32_converted_df, flush = FlushHS32, close = CloseHS32, wait = LagPeriod, 
                                   measure = HSmeasure,  
                                   vf = vfHS32, vr = vrHS32, hs = "HS32", rsqrd = HS_RSqrd, 
                                   low_n = HSlow_n, high_n = HShigh_n, acclimation = 11)




########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 33 ----- ##

# Chamber volume
vrHS33 <- 52.5
# Fish mass (mass = volume)
vfHS33 <- 3.38
# FLush cycle length 
FlushHS33 <- 20
# Close cycle length 
CloseHS33 <- 10

# Read in intermittent flow resprometry data data file 
raw_HS33_RMR_df <- read_csv2("RawHornSharkCSV/RMR_HS33_1210_2.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS33_converted_df <- F4_RMR_trial_cleanup(df = raw_HS33_RMR_df, lag = 0, end = 960, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS33_converted_df, hs = "HS33", scaleyby = 0.1, scalexby = (FlushHS33 +CloseHS33))

# Calculate RMR - mean of n lowest slope estimates
HS33_RMR_est_df <- HS_RMR_mgL_hour(df = HS33_converted_df, flush = FlushHS33, close = CloseHS33, wait = LagPeriod, 
                                   measure = HSmeasure,  
                                   vf = vfHS33, vr = vrHS33, hs = "HS33", rsqrd = HS_RSqrd, 
                                   low_n = HSlow_n, high_n = HShigh_n, acclimation = 10)





########## ----------------------------------------------------------- ###########
vrHS34bg <- 2.573
## ------- HORN SHARK 34 ----- ##

# Chamber volume
vrHS34 <- 2.573
# Fish mass (mass = volume)
vfHS34 <- 0.0387
# FLush cycle length 
FlushHS34 <- 5 
# Close cycle length 
CloseHS34 <- 20 

# Read in intermittent flow resprometry data data file 
raw_HS34_RMR_df <- read_csv2("RawHornSharkCSV/HENRY_1.csv", skip = 1,
                       col_types = cols(Value = col_character(), patm = col_character(),
                                        Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS34_converted_df <- F4_RMR_trial_cleanup(df = raw_HS34_RMR_df, lag = 20, end = 1145, S = HS_S) 

# Plot the converted df
# basic_RMR_plot(df = HS34_converted_df, hs = "HS34", scaleyby = 0.1, scalexby = 25)

### HS34 background respiration was higher than 3% of the shark's total RMR so here we calculate
# what the BR would be at the time each RMR slope was estimated, assuming a linear increase in BR rate
# with time, then correct the RMR estimate by removing background respiration  ------ ###

# Take O2 measurements from after shark was removed from the chamber
HS34_BG_trial_df <- F4_br_convert_DO(raw_HS34_RMR_df, lag = 1330, end = 2000, S = HS_S)
# Calculate background MR
HS34_BG_MR_lm <- lm(O2mgL ~ TimeMinutes, data = HS34_BG_trial_df)
# Slope value
HS34_BG_slope <- ((summary(HS34_BG_MR_lm))$coefficients)[2]
# BG MR Estimate - mg O2 consumed within the respirometer per minute
HS34_BG_MR_est <- vrHS34 * HS34_BG_slope * -60

# Slope of change in background respiration rate (mg O2 /Hr in the chamber) while HS34 was in the chamber
HS34_BG_slope_perHr <- HS34_BG_MR_est/24
# --- Rolling regression with the width of the measurement window as the width of the regressions
HS34_RMR_roll_df <- as.data.frame(roll_regres(O2mgL ~ TimeMinutes, HS34_converted_df, width = HSmeasure*12L, 
                                     do_compute = c("sigmas", "r.squareds", "1_step_forecasts")))

# Add a column LogTime so we can check that the slope values are being pulled from the right timepoint
# --> LogTime indicates the end time of the regression window
HS34_RMR_logtime_df <- HS34_RMR_roll_df %>%
  mutate(RowNum = row_number()) %>%
  mutate(LogTime = seq(0, ((nrow(HS34_RMR_roll_df) * 5) - 1), 5))

# Filter the slopes df to keep the estimates at the desired time point within each closed period
HS34_rmr_rr_close_df <- HS34_RMR_logtime_df %>%
  group_by(grp = rep(row_number(), length.out = n(), each = (5 * 12L) + (20 *12L))) %>%
  slice_tail(n = (20*12L) - (3 * 12L)) %>% # Keep the closed period estimates minus the lag period
  slice_head(n = 7*12L) %>% # Keep the first n estimates from what was just selected
  slice_tail(n = 1) %>% # Keep the last estimate from above - corresponds to the regression ending at flush+wait+measure
  filter(coefs.TimeMinutes < 0,
         r.squareds >= 0.7)%>%
  filter(LogTime >= (12*3600)) %>%
  arrange(desc(coefs.TimeMinutes)) %>%  # Arrange slopes so lowest RMR slope is first
  head(3) %>% # mean of lowest 3 estimates
  mutate(TimeHrs = LogTime/60/60) %>%
  mutate(RMR = (((vrHS34 - vfHS34) * coefs.TimeMinutes)) * -60) %>% # **Calculate whole organism RMR/Hr ***
  mutate(WholeOrg_RMR_corrected = RMR - (HS34_BG_slope_perHr*TimeHrs)) # Calculate what time the 3 lowest RMR slopes were taken from

# -- Re-calculated RMR for HS34 -- ##
HS34_RMR_corrected_WO_est <- mean(HS34_rmr_rr_close_df$WholeOrg_RMR_corrected)
HS34_RMR_est_df <- data.frame(HSID = "HS34", MassKG = vfHS34, 
                                        RMR = (HS34_RMR_corrected_WO_est /vfHS34),
                                        MeanOfN = 3, 
                                        MeanSlope = (mean(HS34_rmr_rr_close_df$coefs.TimeMinutes)),
                                        ChamVol = vrHS34)
                                        
# UNCORRECTED calculation of RMR for HS 34
# HS34_RMR_est_df_UNCORRECTED <- HS_RMR_mgL_hour(df = HS34_converted_df, flush = FlushHS34, close = CloseHS34, wait = LagPeriod, 
#                                    measure = HSmeasure,  vf = vfHS34, vr = vrHS34, hs = "HS34", 
#                                    rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 12)


########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 35 ----- ##

# Chamber volume
vrHS35 <- 2.585
# Fish mass (mass = volume)
vfHS35 <- 0.058
# FLush cycle length 
FlushHS35 <- 5 
# Close cycle length 
CloseHS35 <- 15

# "test" is saved as an updated .csv format to match the function 
# Read in intermittent flow resprometry data data file 
raw_HS35_RMR_df <- read_csv2("RawHornSharkCSV/HARRIETA_RMR_test.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS35_converted_df <- F4_RMR_trial_cleanup(df = raw_HS35_RMR_df, lag = 18, end = 1278, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS35_converted_df, hs = "HS35", scaleyby = 0.1, scalexby = 20)

# Calculate RMR - mean of n lowest slope estimates
HS35_RMR_est_df <- HS_RMR_mgL_hour(df = HS35_converted_df, flush = FlushHS35, close = CloseHS35, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS35,  vr = vrHS35, 
                                   hs = "HS35", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 12)






########## ----------------------------------------------------------- ###########

## ------- HORN SHARK 36 ----- ##

# Chamber volume
vrHS36 <- 2.605
# Fish mass (mass = volume)
vfHS36 <- 0.0432 # Matches final weight as that's what was done for all other sharks
# FLush cycle length 
FlushHS36 <- 5 
# Close cycle length 
CloseHS36 <- 15

# Read in intermittent flow resprometry data data file 
raw_HS36_RMR_df <- read_csv2("RawHornSharkCSV/BABYHORN3_RMR.csv", skip = 1,
                             col_types = cols(Value = col_character(), patm = col_character(),
                                              Temp = col_character())) 

# Convert O2% into O2mgL, rename and tidy columns, add column of time in 5/60
HS36_converted_df <- F4_RMR_trial_cleanup(df = raw_HS36_RMR_df, lag = 18, end = 1278, S = HS_S) 

# Plot the converted df
#basic_RMR_plot(df = HS36_converted_df, hs = "HS36", scaleyby = 0.1, scalexby = 60)

# Calculate RMR - mean of n lowest slope estimates
HS36_RMR_est_df <- HS_RMR_mgL_hour(df = HS36_converted_df, flush = FlushHS36, close = CloseHS36, wait = LagPeriod, 
                                   measure = HSmeasure,  vf = vfHS36,  vr = vrHS36, 
                                   hs = "HS36", rsqrd = HS_RSqrd, low_n = HSlow_n, high_n = HShigh_n, acclimation = 12)








############################################################################################################
###################### ----------------------------------------------------------- #########################
############################################################################################################
############################################################################################################

# --------- Join all RMR data frames together for further analysis ------- #
# RMR is oxygen consumption rate per kg body mass
HS_RMR_all_ests_df <- rbind(HS12_RMR_est_df, HS13_RMR_est_df, HS14_RMR_est_df,HS15_RMR_est_df, 
                            HS16_RMR_est_df, HS17_RMR_est_df, HS18_RMR_est_df, HS19_RMR_est_df,
                            HS20_RMR_est_df, HS21_RMR_est_df, HS22_RMR_est_df, HS26_RMR_est_df,
                            HS27_RMR_est_df,HS30_RMR_est_df, HS31_RMR_est_df, 
                            HS32_RMR_est_df, HS33_RMR_est_df, HS34_RMR_est_df, HS35_RMR_est_df,
                            HS36_RMR_est_df) %>%
  mutate(Model = "RMR",
         MassG_log = (log10(MassKG * 1000)),
         MassG = MassKG*1000,
         MassKG_log = log10(MassKG),
         WholeMMR_KG = (RMR * MassKG), # Whole animal oxygen consumption per hour
         WholeMMR_log_KG = log10(RMR * MassKG),
         Cham_fish_ratio = ChamVol/MassKG) 






