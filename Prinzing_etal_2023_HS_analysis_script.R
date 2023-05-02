
# Prinzing et al. 2023 Horn Shark metabolic rate and gill area scaling analysis

library("readxl")
library("tidyverse")
library("segmented")
library("emmeans")
library("lme4") # lmer()
library("lattice") #qqmath()
library("multcomp") # Comparing MMR methods, cld

#library("lmtest") # Comparing MMR methods?



## ---------------- LOAD IN TRAIT ESTIMATE  DATA SETS -------------------- ##

# MMR CHASE df
MMR_CHASE_df <- read_csv("HS_CHASE_MMR_df.csv") %>%
  arrange(FishID) %>%
  rename(MeasureType = Model) %>%
  mutate(Estimate = MMR_MO2, 
         WholeEstimate_log_KG = WholeMMR_log_KG) %>%
  mutate(WholeMMR_KG = Estimate*MassKG) %>%
  dplyr::select(FishID, MassKG, MeasureType, Estimate, 
                WholeEstimate_log_KG, MassKG_log, Cham_fish_ratio, WholeMMR_KG)

###################------------------------ ########################

# MMR AIR df
MMR_AIR_df <- read_csv("HS_AIR_MMR_df.csv") %>%
  arrange(FishID) %>%
  rename(MeasureType = Model) %>%
  mutate(Estimate = MMR_MO2, 
         WholeEstimate_log_KG = WholeMMR_log_KG) %>%
  mutate(WholeMMR_KG = Estimate*MassKG) %>%
  dplyr::select(FishID, MassKG, MeasureType, Estimate,  
                WholeEstimate_log_KG, MassKG_log, Cham_fish_ratio, WholeMMR_KG)

###################------------------------ ########################

#RMR df 
RMR_df <- read_csv("HS_RMR_df.csv") %>%
  filter(MeanOfN == 3) %>%
  arrange(MassKG) %>%
  rename(MeasureType = Model) %>%
  mutate(Estimate = RMR, 
         WholeEstimate_log_KG = WholeMMR_log_KG,
         FishID = seq(1, 20, 1)) %>%
  mutate(WholeMMR_KG = Estimate*MassKG) %>%
  dplyr::select(FishID, MassKG, MeasureType, Estimate, 
                WholeEstimate_log_KG, MassKG_log, Cham_fish_ratio, WholeMMR_KG)


###################------------------------ ########################

# GSA df
GSA_raw_df <- read_excel("HS_GSA_df.xlsx") %>%
  mutate(Lam_SA_log = log10(BilateralLamellarSA),
         Lam_Freq_log = log10(LamellarFrequency),
         Fil_Length_log = log10(TotalFilamentLength2Sides),
         MassKG_log = log10(MassKG)) %>%
  mutate(GSAmm2 = GSAcm2*100)
# Tidy
GSA_df <- GSA_raw_df %>%
  mutate(Estimate = as.numeric(GSAcm2), 
         WholeEstimate_log_KG = log10(GSAcm2),
         MassKG_log = log10(MassKG), 
         MeasureType = "GSA", 
         Cham_fish_ratio = NA)%>%
  mutate(WholeMMR_KG = "NA") %>%
  arrange(MassKG) %>%
  dplyr::select(FishID, MassKG, MeasureType, Estimate, 
                WholeEstimate_log_KG, MassKG_log, Cham_fish_ratio, WholeMMR_KG)


########################################################################
########################################################################

# -- Bind RMR, MMR and GSA dfs together, stacked -- #
HornShark_all_df <- rbind(RMR_df, MMR_CHASE_df, GSA_df) %>%
  mutate(ScaledMassKG_log = scale(MassKG_log, center = TRUE, scale = FALSE))
# # Just MR data
HornShark_justMR_df <- HornShark_all_df %>%
  filter(MeasureType != "GSA",
         MeasureType != "MMR_air")


########################################################################
###################### ------------------------ ########################
########################################################################

# 1. Linear regression for each of RMR, MMR, and GSA across body mass after log10
# transformation of the raw estimates. MR estimates are mgO2/hr

# 2. Broken regression for GSA ~ body mass

# 3. Difference between slopes of RMR and MMR 

# 4. GSA and MR ratio analysis - * >0.203kg only

# 5. Gill components analysis - fillament length, lamellar frequency, lamellar surface area


# --- Function to extract coefs from linear models of each trait --- #
shark_coef_extract_func <- function(trait_lm, MeasureType, Analysis) {
  # coefs
  lm_coefs <- summary(trait_lm)$coefficients
  trait_slope <- lm_coefs[2]
  trait_intercept <- lm_coefs[1]
  slope_SE <- lm_coefs[2,2]
  slope_p <- lm_coefs[2,4]
  slope_confint <- trait_slope - ((confint(trait_lm))[2])
  intercept_confint <- trait_intercept - ((confint(trait_lm))[1])
  intercept_SE <- lm_coefs[1,2]
  intercept_p <- lm_coefs[1,4]
  lm_adj_rsqrd <- summary(trait_lm)$adj.r.squared
  coefs_df <- data.frame(Intercept = trait_intercept, Intercept_CI = intercept_confint,
                         Slope = trait_slope, Slope_CI = slope_confint, Slope_pval = slope_p,
                         Intercept_pval = intercept_p, Adj_rsqrd = lm_adj_rsqrd, 
                         Slope_SE = slope_SE, Intercept_SE = intercept_SE) %>%
    mutate(MeasureType = MeasureType,
           Analysis = Analysis) %>%
    dplyr::select(MeasureType, Intercept, Intercept_CI, Slope, Slope_CI, 
                  Slope_pval, Intercept_pval, Adj_rsqrd, Slope_SE, 
                  Intercept_SE, Analysis) %>%
    mutate_at(vars(Slope, Slope_CI, Intercept, Intercept_CI, Adj_rsqrd,
                   Slope_SE, Intercept_SE), funs(round(., 3)))
  
  return(coefs_df) 
}


# Trait over body mass basic plot function
trait_mass_plot <- function(df, trait) {
  
  ggplot(data = df, aes(x = MassKG_log, y =  WholeEstimate_log_KG)) +
    geom_point() +
    scale_y_continuous(name = trait)+
    theme_bw()
  
}


################# ------------------------------------------- ####################
################# ------------------------------------------- ####################
################# ------------------------------------------- ####################


### ---- >>>> LINEAR MODELS FOR EACH TRAIT - <<<----- ###

# ** RMR lm ** #
HS_RMR_lm <- lm(WholeEstimate_log_KG ~ MassKG_log, data = RMR_df)
#plot(HS_RMR_lm)
# Get coef values
HS_RMR_SCL <- shark_coef_extract_func(trait_lm = HS_RMR_lm, MeasureType = "RMR",
                                      Analysis = "FullN")

# # --- Test for normal distribution -- # 
#shapiro.test(RMR_df$WholeEstimate_log_KG) # Test if values are normally distributed

################# --------------------- ####################

# ** MMR CHASE lm ** #
HS_CHASE_lm <- lm(WholeEstimate_log_KG ~ MassKG_log, data = MMR_CHASE_df)
#plot(HS_CHASE_lm)
# Get coef values
HS_CHASE_SCL <- shark_coef_extract_func(trait_lm = HS_CHASE_lm, MeasureType = "MMR_chase",
                                        Analysis = "FullN")

# --- Test for normal distribution -- #
#shapiro.test(MMR_CHASE_df$WholeEstimate_log_KG) # Test if values are normally distributed

################# --------------------- ####################

# ** MMR AIR lm ** #
HS_AIR_lm <- lm(WholeEstimate_log_KG ~ MassKG_log, data = MMR_AIR_df)
#summary(HS_AIR_lm)
# Get coef values
HS_AIR_SCL <- shark_coef_extract_func(trait_lm = HS_AIR_lm, MeasureType = "MMR_air",
                                      Analysis = "FullN")

#  --- Test for normal distribution -- # 
#shapiro.test(MMR_AIR_df$WholeEstimate_log_KG) 
#plot(HS_AIR_lm) 

################# --------------------- ####################

#* GSA lm ** #
HS_GSA_lm <- lm(WholeEstimate_log_KG ~ MassKG_log, data = GSA_df)
#plot(HS_GSA_lm)
#summary(HS_GSA_lm)
# Get coef values
HS_GSA_SCL <- shark_coef_extract_func(trait_lm = HS_GSA_lm, MeasureType = "GSA",
                                      Analysis = "FullN")

# # --- Test for normal distribution -- #
#shapiro.test(GSA_df$WholeEstimate_log_KG) # Test if values are normally distributed



################### -------------------------------------- ########################
# ------------------------- GSA BROKEN REGRESSION ------------------------------  # 
################### -------------------------------------- ########################

# Segmented regression function that estimates whether there should be a break point!!
HS_GSA_segmod <- selgmented(HS_GSA_lm, seg.Z = ~MassKG_log, type = "aic",
                                 return.fit = TRUE)
GSA_seg_summary <- summary.segmented(HS_GSA_segmod)
GSA_seg_breakpoint <- GSA_seg_summary$psi # Where on the x-axis the break falls
GSA_brkpnt_confint <- 
  confint.segmented(HS_GSA_segmod)[1] - confint.segmented(HS_GSA_segmod)[2]
# Plot
# plot.segmented(HS_GSA_segmod, conf.level = .95, shade = TRUE, res = TRUE)
# Predict confidence intervals at each point in the model
HS_GSA_CI_fitted_df <- data.frame(predict.segmented(HS_GSA_segmod, interval = "confidence")) %>%
  mutate(FishID = seq(1, 20, 1))
# Pull out slope values
HS_GSA_slopes <- unname(slope(HS_GSA_segmod)$MassKG_log[ ,1])
HS_GSA_seg_confints<- slope(HS_GSA_segmod)

# Gather coefs into df - PLOTTING
HS_GSA_BR_SCL <- 
  data.frame(MeasureType = c("GSA <0.203kg", "GSA >0.203kg"), 
             Intercept = c((intercept(HS_GSA_segmod)$MassKG_log)[1], 
                           (intercept(HS_GSA_segmod)$MassKG_log)[2]),
             Intercept_CI = 0, Slope_pval = 0, Intercept_pval = 0, Intercept_SE = 0,
             Slope = HS_GSA_slopes,
             Slope_CI = c((HS_GSA_slopes[1] - slope(HS_GSA_segmod)$MassKG_log[1,4]),
                          (HS_GSA_slopes[2] - slope(HS_GSA_segmod)$MassKG_log[2,4])),
             Slope_SE = c(slope(HS_GSA_segmod)$MassKG_log[1,2], 
                          slope(HS_GSA_segmod)$MassKG_log[2,2]),
             Adj_rsqrd = GSA_seg_summary$adj.r.squared,
             Analysis = "FullN") %>%
  arrange(MeasureType)


# ---- BROKEN REGRESSION POINTS EXTRACTION PLOTTING FUNCTION ------ #
BR_plot_points_func <- function(br_segmod){
# Slope values
  slopeJuv <- unname(slope(br_segmod)$MassKG_log[1,1])
  slopeAdult <- unname(slope(br_segmod)$MassKG_log[2,1])
# Intercepts
  intcpJuv <- intercept(br_segmod)$MassKG_log[1]
  intcpAdult <- intercept(br_segmod)$MassKG_log[2]
# Calculate the start and end points of each regression line using body mass and 
# the slope and intercept points
  plot_df <- data.frame(
  Y_axis_spots = c(((slopeJuv * log10(0.037)) + intcpJuv),
                   ((slopeJuv * (confint.segmented(br_segmod)[1])) + intcpJuv),
                   ((slopeAdult * log10(4.4)) + intcpAdult)),
  X_axis_spots = c(log10(0.037), (confint.segmented(br_segmod)[1]), log10(4.4)))

return(plot_df)
  
}

# Pull out broken regression start, middle and end points for plotting
HS_GSA_BR_lines_plot_df <- BR_plot_points_func(HS_GSA_segmod)


# Gather coefs into df TABLE FOR MS - made this extra one which doesn't match with the 
# ones for the other regressions so can't be layered nicely, needed to add breakpoints etc. 
HS_GSA_BR_coefs_table <- 
  data.frame(MeasureType = c("GSA Youth", "GSA Adult"), Slope = HS_GSA_slopes, 
             Slope_SE = c(slope(HS_GSA_segmod)$MassKG_log[1,2], 
                          slope(HS_GSA_segmod)$MassKG_log[2,2]),
             Slope_CI = c((HS_GSA_slopes[1] - slope(HS_GSA_segmod)$MassKG_log[1,4]),
                         (HS_GSA_slopes[2] - slope(HS_GSA_segmod)$MassKG_log[2,4])),
             Intercept = c((intercept(HS_GSA_segmod)$MassKG_log)[1], 
                           (intercept(HS_GSA_segmod)$MassKG_log)[2]), 
             Adj_rsqrd = GSA_seg_summary$adj.r.squared, 
             BrkPoint = GSA_seg_breakpoint[2], BrkPoint_SE = GSA_seg_breakpoint[3],
             BrkPoint_CI = GSA_brkpnt_confint, Intercept_CI = 0) %>%
  arrange(MeasureType)


########################## --------------------------------------- ########################


########################## --------------------------------------- ########################
########################## --------------------------------------- ########################
########################## --------------------------------------- ########################


########################## --------------------------------------- ########################


##### -- TEST "DO SLOPES DIFFER ACROSS ESTIMATE TYPE?" (MMR and RMR) -- #####

#  lmer for MR measurements 
HS_MR_slopes_lme <- lmer(WholeEstimate_log_KG ~ MassKG_log*MeasureType + (1|FishID),
                         data = HornShark_justMR_df)
summary(HS_MR_slopes_lme)
# Compare slopes across all MeasureTypes 
HS_MR_slopes_lme_trend <- emtrends(HS_MR_slopes_lme, pairwise ~ MeasureType, var = "MassKG_log")
test(HS_MR_slopes_lme_trend)

# Levene's test
HornShark_tests_df <- HornShark_justMR_df %>%
  mutate(Residuals = (abs(residuals(HS_MR_slopes_lme))^2))
HS_Levene_test <- lm(Residuals ~ MeasureType, data=HornShark_tests_df)
HS_Leven_ANOVA <-  anova(HS_Levene_test) 
#plot(HS_MR_slopes_lme)
#qqmath(HS_MR_slopes_lme, id = 0.05) 



########################## --------------------------------------- ########################
########################## --------------------------------------- ########################
########################## --------------------------------------- ########################



## ---------------------------- RATIO ANALYSIS ------------------------------ ####

################### -------------------------------------- ########################

HS_calculate_ratios_df <- data.frame(FishID = GSA_df$FishID, 
                                     MassKG_log = GSA_df$MassKG_log, 
                                     MassKG = GSA_df$MassKG, GSA_wholeLog = GSA_df$WholeEstimate_log_KG,
                           MMR_wholeLog = MMR_CHASE_df$WholeEstimate_log_KG, 
                           RMR_wholeLog = RMR_df$WholeEstimate_log_KG, 
                           GSA_whole = GSA_df$Estimate,
                           MMR_whole = MMR_CHASE_df$WholeMMR_KG, 
                           RMR_whole = RMR_df$WholeMMR_KG) %>%
  mutate(GSA_RMR_ratio = (GSA_whole/RMR_whole), 
         GSA_MMR_ratio = (GSA_whole/MMR_whole)) %>% 
  mutate(GSA_RMR_ratio_log10 = log10(GSA_RMR_ratio), 
         GSA_MMR_ratio_log10 = log10(GSA_MMR_ratio)) 

# Ratios DF of just > 0.1 kg
HS_calculate_ratios_adult_df <- HS_calculate_ratios_df %>% 
  filter(MassKG > 0.1)


######  GSA:RMR ratio ~ Mass OLS ----------------------- #######
GSA_RMR_ratio_adult_lm <- lm(GSA_RMR_ratio_log10 ~ MassKG_log, 
                             data = HS_calculate_ratios_adult_df)
#summary(GSA_RMR_ratio_adult_lm)
# Extract coef values
GSA_RMR_ratio_SCL_ads <- shark_coef_extract_func(trait_lm =GSA_RMR_ratio_adult_lm,
                                               MeasureType = "GSA-RMR", Analysis = ">0.203kg")

# ------ GSA:MMR ratio ~ Mass OLS ---------- #
GSA_MMR_ratio_adult_lm <- lm(GSA_MMR_ratio_log10 ~ MassKG_log, 
                             data = HS_calculate_ratios_adult_df)
#summary(GSA_MMR_ratio_adult_lm)
# Extract coef values
GSA_MMR_ratio_SCL_ads <- shark_coef_extract_func(trait_lm =GSA_MMR_ratio_adult_lm,
                                                 MeasureType = "GSA-MMR", Analysis = ">0.203kg")




########################## --------------------------------------- ########################
########################## --------------------------------------- ########################

# ------------------------------- Gill COMPONENT ANALYSIS ------------------------------- #

########################## --------------------------------------- ########################
########################## --------------------------------------- ########################

# ---- Fillament length --- #
FilLength_lm <- lm(log10(TotalFilamentLength2Sides) ~ MassKG_log, data = GSA_raw_df)
#summary(FilLength_lm)
# Extract coef values
FilLength_SCL <- shark_coef_extract_func(trait_lm =FilLength_lm,  
                                             MeasureType = "Fil Length", Analysis = "FullN")

######## --------------------------------------- #########
# --- Lamellar frequency --- #
Lam_FREQ_lm <- lm(log10(LamellarFrequency) ~ MassKG_log, data = GSA_raw_df)
#summary(Lam_FREQ_lm)
# Extract coef values
LamFreq_SCL <- shark_coef_extract_func(trait_lm =Lam_FREQ_lm,  
                                         MeasureType = "Lam Freq", Analysis = "FullN")


######## --------------------------------------- #########
#  --- Lamellar SURFACE AREA --- ##
Lam_SA_lm <- lm(log10(BilateralLamellarSA) ~ MassKG_log, data = GSA_raw_df)
#summary(Lam_SA_lm)
# Extract coef values
LamSA_SCL <- shark_coef_extract_func(trait_lm =Lam_SA_lm,  
                                       MeasureType = "Lam SA", Analysis = "FullN")

# Lam SA broken regression 
HS_LamSA_segmod <- selgmented(Lam_SA_lm, seg.Z = ~MassKG_log, Kmax = 1,
                                      type = "aic", return.fit = TRUE)
#plot.segmented(HS_LamSA_segmod, conf.level = .95, shade = TRUE, res = TRUE)

LamSA_seg_summary <- summary.segmented(HS_LamSA_segmod)
LamSA_seg_breakpoint <- LamSA_seg_summary$psi # Where on the x-axis the break falls
LamSA_brkpnt_confint <- 
  confint.segmented(HS_LamSA_segmod)[1] - confint.segmented(HS_LamSA_segmod)[2]
LamSA_brkpnt_SE <- LamSA_seg_summary$psi[3]

# Pull out slope values
HS_LamArea_slopes <- unname(slope(HS_LamSA_segmod)$MassKG_log[ ,1])
HS_LamArea_seg_confints <- slope(HS_LamSA_segmod)
# Pull out slope values
HS_LamArea_BR_SCL <- 
  data.frame(MeasureType = c("LamArea <0.203kg", "LamArea >0.203kg"), 
             Intercept = c((intercept(HS_LamSA_segmod)$MassKG_log)[1], 
                           (intercept(HS_LamSA_segmod)$MassKG_log)[2]),
             Intercept_CI = 0, Slope_pval = 0, Intercept_pval = 0, Intercept_SE = 0,
             Slope = HS_LamArea_slopes, 
             Slope_CI = c((HS_LamArea_slopes[1] - slope(HS_LamSA_segmod)$MassKG_log[1,4]),
                         (HS_LamArea_slopes[2] - slope(HS_LamSA_segmod)$MassKG_log[2,4])),
             Slope_SE = c(slope(HS_LamSA_segmod)$MassKG_log[1,2], 
                          slope(HS_LamSA_segmod)$MassKG_log[2,2]),
             Adj_rsqrd = LamSA_seg_summary$adj.r.squared,
             Analysis = "FullN") %>%
  arrange(MeasureType)

# Pull out broken regression start, middle and end points for plotting
HS_LamSA_BR_lines_plot_df <- BR_plot_points_func(HS_LamSA_segmod)



###########################################################################################
########################## --------------------------------------- ########################
###########################################################################################

####----------------------- MODEL OUTPUTS FOR TABLES AND PLOTTING -------------- #####

###########################################################################################

# DF of slope and intercept estimates
HornShark_coeficient_df_ms <- rbind(HS_RMR_SCL, HS_CHASE_SCL, HS_GSA_SCL, HS_GSA_BR_SCL,
                                    GSA_RMR_ratio_SCL_ads, GSA_MMR_ratio_SCL_ads,
                                 FilLength_SCL, LamFreq_SCL, HS_LamArea_BR_SCL, LamSA_SCL) %>%
  mutate(MeasureType = factor(MeasureType, levels = MeasureType)) %>%
  mutate(Conf_up = Slope + Slope_CI, 
         Conf_low = Slope - Slope_CI,
         SE_up = Slope + Slope_SE, 
         SE_low = Slope - Slope_SE,
         ColorCol = case_when(MeasureType == "GSA >0.203kg" ~ "b",
                              MeasureType == "GSA <0.203kg" ~ "b",
                              MeasureType != "GSA >0.203kg" ~ "a"))%>%
  mutate_at(vars(Slope, Slope_CI, Intercept, Intercept_CI,
                 Slope_SE, Intercept_SE), funs(round(., 3))) %>%
  dplyr::select(MeasureType,  Intercept, Intercept_SE, 
                Slope, Slope_CI, Slope_SE, Conf_up, Conf_low, Adj_rsqrd,
                SE_up, SE_low, ColorCol, Intercept_CI, Analysis)



########################## --------------------------------------- ########################
########################## --------------------------------------- ########################
########################## --------------------------------------- ########################

#### ----- SUPPLEMENTARY MATERIAL - DOES MMR METHOD EFFECT MMR ESTIMATE?? ---- #####

# Mean mass, slopes and intercepts 
HS_mean_massKG <- mean(MMR_CHASE_df$MassKG)
HS_CHASE_intercept <- HS_CHASE_SCL$Intercept
HS_CHASE_slope <- HS_CHASE_SCL$Slope

HS_AIR_intercept <- HS_AIR_SCL$Intercept
HS_AIR_slope <- HS_AIR_SCL$Slope
########################## --------------------------------------- ########################

# MMR CHASE DF with standardized MMR and MASS estimates for MMR METHOD COMPARISON
HS_CHASE_comp_df <- MMR_CHASE_df %>%
  mutate(Pred_Whole_MMR_KG = (10^(HS_CHASE_intercept) * (MassKG^HS_CHASE_slope))) %>%
  mutate(Pred_MMR_diff =  (WholeMMR_KG - Pred_Whole_MMR_KG),
         HS_MMR_meanKG = (10^(HS_CHASE_intercept) * (HS_mean_massKG^HS_CHASE_slope))) %>% 
  mutate(MMR_scaled_xKG = HS_MMR_meanKG + Pred_MMR_diff)%>%
  dplyr::select(FishID, MeasureType, Estimate, MassKG, MassKG_log, MMR_scaled_xKG, 
                Pred_Whole_MMR_KG, Pred_MMR_diff, HS_MMR_meanKG)

# Smaller df for figuring out how estimates vary between individual sharks
HS_sup_chase_comp_df <- HS_CHASE_comp_df %>%
  mutate(EstimateChase = MMR_scaled_xKG) %>%
  dplyr::select(FishID, EstimateChase)

# MMR AIR DF with standardized MMR and MASS estimates for MMR METHOD COMPARISON
HS_AIR_comp_df <- MMR_AIR_df %>%
  mutate(Pred_Whole_MMR_KG = (10^(HS_AIR_intercept) * (MassKG^HS_AIR_slope))) %>%
  mutate(Pred_MMR_diff =  (WholeMMR_KG - Pred_Whole_MMR_KG),
         HS_MMR_meanKG = (10^(HS_AIR_intercept) * (HS_mean_massKG^HS_AIR_slope))) %>% 
  mutate(MMR_scaled_xKG = HS_MMR_meanKG + Pred_MMR_diff)%>%
  dplyr::select(FishID, MeasureType, Estimate, MassKG, MassKG_log, MMR_scaled_xKG, 
                Pred_Whole_MMR_KG, Pred_MMR_diff, HS_MMR_meanKG)

# Smaller df for figuring out how estimates vary between individual sharks
HS_sup_air_comp_df <- HS_AIR_comp_df %>%
  mutate(EstimateAir = MMR_scaled_xKG) %>%
  dplyr::select(FishID, MassKG, EstimateAir)

# Bind smaller dfs by FishID 
HS_sup_MMR_bound_df <- left_join(HS_sup_chase_comp_df, HS_sup_air_comp_df, 
                                 by = "FishID") %>%
  mutate(MMR_calc = EstimateChase - EstimateAir)

# Mean MMR for each method
HS_mean_scaled_MMR_chase <- mean(HS_CHASE_comp_df$MMR_scaled_xKG)
HS_mean_scaled_MMR_air <- mean(HS_AIR_comp_df$MMR_scaled_xKG)
# Percent difference between means for each method
HS_mean_scaled_MMR_chase/HS_mean_scaled_MMR_air

################################# ------------------------ ##################################

# Joind the dfs of MMR CHASE and AIR for comparison
HS_MMR_chaseVSair_df <- rbind(HS_CHASE_comp_df, HS_AIR_comp_df) %>%
  mutate(emmean = MMR_scaled_xKG,
         WholeEstimate_log_KG = log10(MassKG * Estimate)) 

## --- Compare between models with fish as a random effect ---- ###
HS_MMR_chaseVSair_lmer <- lmer(MMR_scaled_xKG ~ MeasureType + (1|FishID), data = HS_MMR_chaseVSair_df,
                               REML=TRUE)
# TEST - "DOES METHOD EFFECT MMR ESTIMATE?"
HS_MMR_chaseVSair_glht <- glht(HS_MMR_chaseVSair_lmer, linfct = mcp(MeasureType = "Tukey"))
summary(HS_MMR_chaseVSair_glht) # ** ALMOST 
HS_MMR_cVa_glht_summ <- (summary(HS_MMR_chaseVSair_glht))$test
# GLHT coefficients
HS_MMR_cVa_glht_coefs_df <- data.frame(HS_MMR_cVa_glht_summ$coefficients)
HS_MMR_cVa_glht_coefs <- HS_MMR_cVa_glht_coefs_df$HS_MMR_cVa_glht_summ.coefficients

# ANOVA to produce F-Value
HS_MMR_cVa_anova <- anova(HS_MMR_chaseVSair_lmer)

# F-test to compare two variances
#var.test(HS_CHASE_comp_df$MMR_scaled_xKG, HS_AIR_comp_df$MMR_scaled_xKG)
#plot(HS_MMR_chaseVSair_lmer)
#qqmath(HS_MMR_chaseVSair_lmer, id = 0.05)





################################# ------------------------ ##################################




