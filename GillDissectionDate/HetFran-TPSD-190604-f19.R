
# HetFran-TPSD-190604-f19 GSA Caclculation Script 
# Estimate total length of all filaments on all hemibranchs on both sides of head. Median filament is representative of all filaments in that bin
# Estimate mean lamellar frequency (lamellar frequency per unit length on one side of the filament, multiplied by 2 to account for both sides of the filament)
# Estimate mean bilateral surface area of a lamellae
# Nov 13, 2019, Tanya Prinzing

# Load packages
library(dplyr)
library(knitr)
library(readxl)

# Read in filament length, lamellar frequency, and lamellar surface area data created using ImageJ
GSA_Raw_df <- read_excel("HetFran-TPSD-190604-f19-GSAData.xlsx", sheet="Filaments", col_names = TRUE)

# Calculate necessary parameters for final GSA calculation
GSA_df <- GSA_Raw_df %>%
  mutate(Arch_Hemibranch = paste0(ArchNum, ArchCode), # Create a column for specific arch position
         BinFilLenMM = (NumFilaments * MedFilLengthMM)) # Calculate total length of filaments in each bin in mm
  
# Calculate lamellar frequency and surface area from arch Posterior 2
GSALamellae_df <- GSA_df %>%
  rowwise() %>%
  filter(Arch_Hemibranch == "2Posterior") %>%
  mutate(LamFreqBase = (LamCountBase / LFBaseDistanceMM), # Lamellar frequency at base, middle and tip per mm
         LamFreqMid = (LamCountMid / LFMidDistanceMM),
         LamFreqTip = (LamCountTip / LFTipDistanceMM)) %>%
  mutate(AvgLamFreqBinMM = mean(c(LamFreqBase, LamFreqMid, LamFreqTip), na.rm = TRUE), # Average lamellar frequency on the median filament
         AvLamAreaBinMM2 = 2*(mean(c(LamAreaBaseMM2, LamAreaMidMM2, LamAreaTipMM2), na.rm = TRUE))) %>% # Average lamellar area on the median fillament, multiplied by 2 for both sides of the lamellae
  mutate(LamFreqTotalBin = AvgLamFreqBinMM * (MedFilLengthMM * NumFilaments)) %>% # Number of lamellae in that bin, based on how long the filaments are
  mutate(LamAreaTotalBin = LamFreqTotalBin * AvLamAreaBinMM2) # Total area of lamellae in that bin, bases on how many lamellae there are (from step above)

# Calculate the hemibranch filament length (length of all filaments on each hemibranch). 
SumFilLenHemi_df <- GSA_df %>% group_by(Arch_Hemibranch) %>% 
  summarise(sumTotalFilLengMM = sum(BinFilLenMM))

# Calculate average lamellae area 
AvLamAreaMM2 <- ((sum(GSALamellae_df$LamAreaTotalBin))/(sum(GSALamellae_df$LamFreqTotalBin)))
# Calculate average lamellar frequency /mm on both sides of the filament
AvLamFreq <- (sum(GSALamellae_df$LamFreqTotalBin))/(sum(GSALamellae_df$BinFilLenMM)) 
# Calculate total filament length on one side of the head
FilLeng2Side <- (sum(SumFilLenHemi_df$sumTotalFilLengMM)) * 2

# Calculate Gill Surface Area in cm2
GSA <- (FilLeng2Side * (2*AvLamFreq) * AvLamAreaMM2) * 0.01

















