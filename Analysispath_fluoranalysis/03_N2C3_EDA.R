#!bin/bash
library(lubridate)
library(tidyverse)

############# Exploratory Data Analysis #############
dir.create("03_N2C3_EDA")
# Notice that I name it 02. The numbers in front of the code/directories reminds me which order I need to run the code in
# This is "03" because it depends on output from "02" (compile_data). 
# If I have outputs from the EDA that I want to use downstream, then I will name them "03_"

dat <- read.delim("02_compile_data/dat_n2c3_vs_n2c3.txt")
dat_noplant <- read.delim("02_compile_data/dat_n2c3_vs_n2c3_noplant.txt")

dat
gg_n2c3_vs_cellcounts <- dat %>%
  filter(experiment == "2023-08-08_fluor") %>%
  mutate(StrainMix2 = gsub("N2C3","N2C3", gsub("WCS365","LacZ",StrainMix))) %>%
  mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2C3","LacZ-MOCK","LacZ-N2C3"))) %>%
  mutate(ratio = factor(ratio, levels=c("0-0","1-1","1-0","100-1","10000-1","0-1","1-100","1-10000"))) %>%
  ggplot(aes(x=NN_FC_CFU, y=NN_FC_fluor_offset)) +
  geom_point(aes(col=StrainMix2, pch=ratio))+
  ylab("Crim:Neon Fluorescence-estimated\nLog2Fold Change in cell count\n(Crim>0, Neon<0)")+
  xlab("Crim:Neon LacZ CFU counts\nLog2Fold Change\n(Crim>0, Neon<0)") +
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  facet_grid(.~fluor_protect, labeller = labeller(fluor_protect=c(crim="LacZ is crim", neon="LacZ is neon")))+
  geom_hline(aes(yintercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  labs(shape="LacZ:N2C3 ratios", col="Strain Mix") +
  scale_colour_manual(values=c("MOCK-MOCK"="black","MOCK-N2C3"="lightblue","LacZ-MOCK"="blue","LacZ-N2C3"="darkred")) 
gg_n2c3_vs_cellcounts
ggsave("03_N2C3_EDA/gg_n2c3_vs_cellcounts.png", gg_n2c3_vs_cellcounts, width=8, height=4)

gg_n2c3_vs_cellcounts_nottransformed <- dat %>%
  filter(experiment == "2023-08-08_fluor") %>%
  mutate(StrainMix2 = gsub("N2C3","N2C3", gsub("WCS365","LacZ",StrainMix))) %>%
  mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2C3","LacZ-MOCK","LacZ-N2C3"))) %>%
  mutate(ratio = factor(ratio, levels=c("0-0","1-1","1-0","100-1","10000-1","0-1","1-100","1-10000"))) %>%
  ggplot(aes(x=NN_FC_CFU, y=NN_FC_fluor, )) +
  geom_point(aes(col=StrainMix2, pch=ratio))+
  ylab("Crim:Neon Fluorescence-estimated\nLog2Fold Change in cell count\n(Crim>0, Neon<0)")+
  xlab("Crim:Neon LacZ CFU counts\nLog2Fold Change\n(Crim>0, Neon<0)") +
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  facet_grid(.~fluor_protect, labeller = labeller(fluor_protect=c(crim="LacZ is crim", neon="LacZ is neon")))+
  geom_hline(aes(yintercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  labs(shape="LacZ:N2C3 ratios", col="Strain Mix") +
  scale_colour_manual(values=c("MOCK-MOCK"="black","MOCK-N2C3"="lightblue","LacZ-MOCK"="blue","LacZ-N2C3"="darkred")) 
gg_n2c3_vs_cellcounts_nottransformed
ggsave("03_N2C3_EDA/gg_n2c3_vs_cellcounts_nottransformed.png", gg_n2c3_vs_cellcounts_nottransformed, width=8, height=4)


gg_n2c3_ratios_fluor <- dat %>%
  filter(is.finite(NN_FC_fluor_offset)) %>%
  mutate(StrainMix2 = gsub("N2C3","N2", gsub("WCS365","N1",StrainMix))) %>%
  mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2","N1-MOCK","N1-N2"))) %>%
  # filter(protect =="WCS365") %>%
  ggplot(aes(x=ratio, y=NN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=fluor_protect), width=0.1, height=0, alpha=0.5) +
  facet_grid(.~StrainMix2, drop=TRUE, scales="free_x", space = "free")+
  xlab("Inoculation Ratio\nN2C3 (Strain 1):(Strain2)") +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values=c("orange","darkgreen")) +
  labs(col="Fluor of N2C3\nStrain1") +
  ylab("Log2Fold change (N1 : N2)")
gg_n2c3_ratios_fluor
ggsave("03_N2C3_EDA/gg_n2c3_ratios_fluor.png", gg_n2c3_ratios_fluor, width=8, height=4)



