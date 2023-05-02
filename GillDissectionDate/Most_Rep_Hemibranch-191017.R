
# Tanya Prinzing Oct 15 2019
# Estimate surface area of each hemibranch to determine most representative for Horn Sharks
# Dissected and photographed lamellae from all hemibrnahcs on HS 13

library(dplyr)
library(knitr)
library(readxl)

# Load HS 13 data frame
GSALams_raw_df <-read_excel("HetFran-TPSD-190604-m13_GSAData.xlsx", 
                                    sheet="Filaments", col_names = TRUE)

# Calculate total lamellar surface area for each bin on each hemibranch (BOTH sides of the filament)
GSALams_df <- GSALams_raw_df %>%
  rowwise() %>%
  mutate(Arch_Hemibranch = paste(ArchNum, ArchCode),
         LamFreqBase = (LamCountBase/LFBaseDistanceMM), #Lam frequency at fillament base, at each spot
         LamFreqMid = (LamCountMid/LFMidDistanceMM), # Lam frequency at fillament mid, at each spot
         LamFreqTip = (LamCountTip/LFTipDistanceMM)) %>% # Lam frequency at fillament tip, at each spot
  mutate(AvgLamFreqMM = mean(c(LamFreqBase, LamFreqMid, LamFreqTip), #Total average Lam FREQUENCY for the shark
                             na.rm = TRUE), 
         AvLamAreaMM2 = mean(c(LamAreaBaseMM2, LamAreaMidMM2, # Toal average lam AREA for the shark
                               LamAreaTipMM2), na.rm = TRUE)) %>%
  mutate(LamFreqBin = (AvgLamFreqMM * MedFilLengthMM * NumFilaments)
         ) %>% # Total lamellar count for the shark using total average lamellar frequency 
  mutate(LamAreaBin = (AvLamAreaMM2 * LamFreqBin * 2) # Total SA of lamellae for the shark using total count
                                                      # of lamellae and average Lamellar area for the shark
         ) # Doubled for both sides of the filament


##### ----- Calculate which hemibranch has average SA closest to mean SA for all hemibranchs
GSAAvLam_df <- GSALams_df %>%
  group_by(Arch_Hemibranch) %>%
  summarise(AvLamHemi = mean(c(LamAreaBaseMM2, # Average lamellar SA on each hemibranch
                               LamAreaMidMM2, LamAreaTipMM2), na.rm = TRUE)) 

# Mean SA for all hemibranchs
MeanLamSA <- mean(GSAAvLam_df$AvLamHemi)
# Calculate difference between mean SA and SA of each hemibranch to find which hemibranch is 
# closest to the mean
SmallestLamDif_df <- GSAAvLam_df %>%
  mutate(Difference = (AvLamHemi - MeanLamSA))%>% 
  mutate(RelativeDif = (((abs(AvLamHemi - MeanLamSA))/((MeanLamSA + AvLamHemi)/2)) * 100))
# Gives Posterior 2 as the most representative





