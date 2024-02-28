#!bin/bash
library(lubridate)
library(tidyverse)
library(brms)
library(cowplot)

############# Load #############
dir.create("04_WCSvsN2C3_EDA")

dat <- read.delim("02_compile_data/dat_wcs_vs_n2c3.txt") %>%
  mutate(Alive = (Healthiness_hsv>400)) %>%
  mutate(WN_FC_inoc = log10(protect_cells/path_cells)) 
N2C3_labeller <- c(`0`="No N2C3", `2`="N2C3 ~ 10^2 cells", `3`="N2C3 ~ 10^3 cells", `4`="N2C3 ~ 10^4 cells", `5`="N2C3 ~ 10^5 cells",`6`="N2C3 ~ 10^6 cells", `7`="N2C3 ~ 10^7 cells", `8`="N2C3 ~ 10^8 cells")

##### Plate maps #########

gg_platemap_fluor <- dat %>%
  mutate(row = factor(row, levels=rev(LETTERS[1:8]))
         , col=factor(col, levels=seq(1,12))) %>%
  mutate(CFU_ratios_withNAs = ifelse((!is.na(neon_CFU)|!is.na(crim_CFU)), "CFU COUNTED", NA)) %>%
  mutate(protect_cells_log = ifelse(protect=="MOCK", max(protect_cells_log), protect_cells_log)) %>%
  mutate(path_cells_log = ifelse(path=="MOCK", max(path_cells_log), path_cells_log)) %>%
  ggplot()+
  geom_tile(aes(x=col, y=row, fill=protect, alpha=protect_cells_log)) +
  geom_point(aes(x=col, y=row, cex=CFU_ratios_withNAs)) +
  geom_point(aes(x=col, y=row, col=path, alpha=path_cells_log)) +
  scale_fill_manual(values=c(MOCK="grey", WCS365 = "darkgreen"))  +
  scale_color_manual(values=c(MOCK="white", N2C3 = "darkorange"))  +
  facet_wrap(experiment~plate) +
  labs(col="Pathogen", fill="Commensal", alpha = "Cells (log10)", size="CFU counted") +
  ylab("Row") + xlab("Column")
gg_platemap_fluor
ggsave(filename = "04_WCSvsN2C3_EDA/gg_platemap_fluor.png", gg_platemap_fluor, height=8, width=10)
# Count the counted wells
dat %>%
  mutate(row = factor(row, levels=rev(LETTERS[1:8]))
         , col=factor(col, levels=seq(1,12))) %>%
  mutate(CFU_ratios_withNAs = ifelse((!is.na(neon_CFU)|!is.na(crim_CFU)), "CFU COUNTED", NA)) %>%
  select(CFU_ratios_withNAs) %>%table()

# Checking fluorescence
gg_fluor_sanitycheck <- dat %>%
  mutate(row = factor(row, levels=rev(LETTERS[1:8]))
         , col=factor(col, levels=seq(1,12))) %>%
  mutate(protect_cells_log = ifelse(protect=="MOCK", max(protect_cells_log), protect_cells_log)) %>%
  mutate(path_cells_log = ifelse(path=="MOCK", max(path_cells_log), path_cells_log)) %>%
  ggplot()+
  geom_tile(aes(x=col, y=row, fill=crim_raw)) +
  geom_point(aes(x=col, y=row, col=neon_raw)) +
  scale_fill_gradient(low="grey", high="darkred")  +
  scale_color_gradient(low="white", high="green")  +
  facet_wrap(experiment~plate) +
  # labs(col="Pathogen", fill="Commensal", alpha = "Cells (log10)") +
  ylab("Row") + xlab("Column")
gg_fluor_sanitycheck
ggsave(filename="04_WCSvsN2C3_EDA/gg_fluor_sanitycheck.png", gg_fluor_sanitycheck, height=8, width=10)



# Checking plant health

gg_plant_sanitycheck <- dat %>%
  mutate(row = factor(row, levels=rev(LETTERS[1:8]))
         , col=factor(col, levels=seq(1,12))) %>%
  mutate(fluor_ratios_withNAs = ifelse(!is.na(WN_FC_fluor_offset),WN_FC_fluor_offset,
                                     ifelse((!is.na(neon_CFU)|!is.na(crim_CFU)), 0, NA))) %>%
  ggplot()+
  geom_tile(aes(x=col, y=row, fill=Healthiness_hsv)) +
  geom_point(aes(x=col, y=row, col=fluor_ratios_withNAs))+
  scale_fill_gradient(low="grey", high="darkgreen")  +
  scale_color_gradient2(low="darkorange", mid="white", high="green", na.value = "black")  +
  facet_wrap(experiment~plate) +
  labs(col="Log2 FC\nfluor", fill="Health Score") +
  ylab("Row") + xlab("Column")
gg_plant_sanitycheck
ggsave(filename="04_WCSvsN2C3_EDA/gg_plant_sanitycheck.png", gg_plant_sanitycheck, height=8, width=10)


############# Fluor vs cell counts #############

gg_wcs_vs_cellcounts <- dat %>%
  # mutate(StrainMix2 = gsub("N2C3","N2C3", gsub("WCS365","LacZ",StrainMix))) %>%
  # mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2C3","LacZ-MOCK","LacZ-N2C3"))) %>%
  mutate(ratio = factor(ratio, levels=c("0-0","1-1","1-0","100-1","10000-1","0-1","1-100","1-10000"))) %>%
  ggplot(aes(x=WN_FC_CFU, y=WN_FC_fluor_offset )) +
  geom_point(aes(col=StrainMix, pch=ratio))+
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  ylab("WCS365:N2C3 Fluorescence-estimated\nLog2Fold Change\n(WCS>0, N2C3<0)")+
  xlab("WCS365:N2C3 CFU counts\nLog2Fold Change\n(WCS>0, N2C3<0)") +
  geom_hline(aes(yintercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  geom_abline(aes(slope=1, intercept=0), lty=2, col="grey")+
  labs(shape="WCS365:N2C3 ratios", col="Strain Mix") +
  # ylim(-1.5, 1.5) +
  scale_colour_manual(values=c("MOCK-MOCK"="black","MOCK-N2C3"="darkred","WCS365-MOCK"="blue","WCS365-N2C3"="purple"))
gg_wcs_vs_cellcounts
ggsave("04_WCSvsN2C3_EDA/gg_wcs_vs_cellcounts.png", gg_wcs_vs_cellcounts, width=6, height=4)

gg_wcs_vs_cellcounts_nottransformed <- dat %>%
  # mutate(StrainMix2 = gsub("N2C3","N2C3", gsub("WCS365","LacZ",StrainMix))) %>%
  # mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2C3","LacZ-MOCK","LacZ-N2C3"))) %>%
  mutate(ratio = factor(ratio, levels=c("0-0","1-1","1-0","100-1","10000-1","0-1","1-100","1-10000"))) %>%
  ggplot(aes(x=WN_FC_CFU, y=WN_FC_fluor )) +
  geom_point(aes(col=StrainMix, pch=ratio))+
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  ylab("WCS365:N2C3 Fluorescence-estimated\nLog2Fold Change\n(WCS>0, N2C3<0)")+
  xlab("WCS365:N2C3 CFU counts\nLog2Fold Change\n(WCS>0, N2C3<0)") +
  geom_hline(aes(yintercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  geom_abline(aes(slope=1, intercept=0), lty=2, col="grey")+
  labs(shape="WCS365:N2C3 ratios", col="Strain Mix") +
  # ylim(-1.5, 1.5) +
  scale_colour_manual(values=c("MOCK-MOCK"="black","MOCK-N2C3"="darkred","WCS365-MOCK"="blue","WCS365-N2C3"="purple"))
gg_wcs_vs_cellcounts
ggsave("04_WCSvsN2C3_EDA/gg_wcs_vs_cellcounts_nottransformed.png", gg_wcs_vs_cellcounts_nottransformed, width=6, height=4)


#### panelled by exp #######

gg_wcs_vs_cellcounts_byexp <- dat %>%
  filter(experiment!="2023-11-28_fluor") %>%
  # mutate(StrainMix2 = gsub("N2C3","N2C3", gsub("WCS365","LacZ",StrainMix))) %>%
  # mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2C3","LacZ-MOCK","LacZ-N2C3"))) %>%
  mutate(ratio = factor(ratio, levels=c("0-0","1-1","1-0","100-1","10000-1","0-1","1-100","1-10000"))) %>%
  ggplot(aes(x=WN_FC_CFU, y=WN_FC_fluor_offset )) +
  geom_point(aes(col=StrainMix, pch=ratio))+
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  ylab("WCS365:N2C3 Fluorescence-estimated\nLog2Fold Change\n(WCS>0, N2C3<0)")+
  xlab("WCS365:N2C3 CFU counts\nLog2Fold Change\n(WCS>0, N2C3<0)") +
  geom_hline(aes(yintercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  geom_abline(aes(slope=1, intercept=0), lty=2, col="grey")+
  labs(shape="WCS365:N2C3 ratios", col="Strain Mix") +
  facet_grid(.~experiment) +
  # ylim(-1.5, 1.5) +
  scale_colour_manual(values=c("MOCK-MOCK"="black","MOCK-N2C3"="darkred","WCS365-MOCK"="blue","WCS365-N2C3"="purple"))
gg_wcs_vs_cellcounts_byexp
ggsave("04_WCSvsN2C3_EDA/gg_wcs_vs_cellcounts_byexp.png", gg_wcs_vs_cellcounts_byexp, width=8, height=4)

gg_wcs_vs_cellcounts_byexp_nottransformed <- dat %>%
  filter(experiment!="2023-11-28_fluor") %>%
  # mutate(StrainMix2 = gsub("N2C3","N2C3", gsub("WCS365","LacZ",StrainMix))) %>%
  # mutate(StrainMix2 = factor(StrainMix2, levels=c("MOCK-MOCK","MOCK-N2C3","LacZ-MOCK","LacZ-N2C3"))) %>%
  mutate(ratio = factor(ratio, levels=c("0-0","1-1","1-0","100-1","10000-1","0-1","1-100","1-10000"))) %>%
  ggplot(aes(x=WN_FC_CFU, y=WN_FC_fluor )) +
  geom_point(aes(col=StrainMix, pch=ratio))+
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  ylab("WCS365:N2C3 Fluorescence-estimated\nLog2Fold Change\n(WCS>0, N2C3<0)")+
  xlab("WCS365:N2C3 CFU counts\nLog2Fold Change\n(WCS>0, N2C3<0)") +
  geom_hline(aes(yintercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  geom_abline(aes(slope=1, intercept=0), lty=2, col="grey")+
  labs(shape="WCS365:N2C3 ratios", col="Strain Mix") +
  facet_grid(.~experiment) +
  scale_colour_manual(values=c("MOCK-MOCK"="black","MOCK-N2C3"="darkred","WCS365-MOCK"="blue","WCS365-N2C3"="purple"))
gg_wcs_vs_cellcounts_byexp_nottransformed
ggsave("04_WCSvsN2C3_EDA/gg_wcs_vs_cellcounts_byexp_nottransformed.png", gg_wcs_vs_cellcounts_byexp_nottransformed, width=8, height=4)


############# Health vs log2foldratios #############

gg_ratios_fluor  <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  # filter(protect =="WCS365") %>%
  ggplot(aes(x=ratio, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=fluor_protect), width=0.1, height=0, alpha=0.5) +
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x", space = "free") +
  xlab("Inoculation Ratio\n WCS365:N2C3") +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values=c("orange","darkgreen")) +
  labs(col="Fluor of\nWCS365") +
  ylab("Log2Fold change (WCS365 : N2C3)")
gg_ratios_fluor
ggsave("04_WCSvsN2C3_EDA/gg_ratios_fluor.png", gg_ratios_fluor, width=10, height=4)


gg_ratios_fluor_wplanthealth  <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  # filter(protect =="WCS365") %>%
  ggplot(aes(x=ratio, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv), width=0.1, height=0, alpha=0.5) +
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x", space = "free") +
  xlab("Inoculation Ratio\n WCS365:N2C3") +
  geom_hline(aes(yintercept=0)) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  labs(col="Health score") +
  ylab("Log2Fold change (WCS365 : N2C3)")
gg_ratios_fluor_wplanthealth
ggsave("04_WCSvsN2C3_EDA/gg_ratios_fluor_wplanthealth.png", gg_ratios_fluor_wplanthealth, width=10, height=4)


gg_ratios_fluor_wplantbin  <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  # mutate(Alive = Healthiness_hsv>400) %>%
  # filter(protect =="WCS365") %>%
  ggplot(aes(x=ratio, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive), width=0.1, height=0, alpha=0.5) +
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x", space = "free") +
  xlab("Inoculation Ratio\n WCS365:N2C3") +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values=c("orange","darkgreen")) +
  labs(col="Healthy plant\n(>400)") +
  ylab("Log2Fold change (WCS365 : N2C3)")
gg_ratios_fluor_wplantbin
ggsave("04_WCSvsN2C3_EDA/gg_ratios_fluor_wplantbin.png", gg_ratios_fluor_wplantbin, width=10, height=4)

### Plots for pub ####
gg_ratios_fluor_mixtreats  <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  filter(StrainMix == "WCS365-N2C3") %>%
  # mutate(Alive = Healthiness_hsv>400) %>%
  # filter(protect =="WCS365") %>%
  ggplot(aes(x=ratio, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive), width=0.1, height=0, alpha=0.5) +
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x", space = "free", 
             labeller=labeller(StrainMix=c("WCS365-N2C3"="Mixed inoculation treatments"))) +
  xlab("Inoculation Ratio\n WCS365:N2C3") +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values=c("orange","darkgreen")) +
  labs(col="Healthy plant\n(>400)") +
  ylab("Log2Fold change (WCS365 : N2C3)") +
  ylim(-15, 15) # make even on top and bottom
gg_ratios_fluor_mixtreats

gg_ratios_fluor_monotreats  <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  filter(StrainMix != "WCS365-N2C3") %>%
  mutate(StrainMix2 = ifelse(StrainMix=="MOCK-MOCK", "MOCK",
                             ifelse(StrainMix == "MOCK-N2C3", "N2C3","WCS365"))) %>%
  mutate(NewCat = "Mock/Monoculture\ntreatments") %>%
  mutate(absFluor = (ifelse(WN_FC_fluor_offset>0,"More WCS365", "More N2C3"))) %>%
  ggplot(aes(x=StrainMix2, y=absFluor, col=Alive)) +
  geom_jitter(height=0.1, width=0.3, alpha=0.5, show.legend = FALSE) +
  facet_grid(.~NewCat, drop=TRUE, scales="free_x", space = "free") +
  # xlab("Inoculation Ratio\n WCS365:N2C3") +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values=c("orange","darkgreen")) +
  labs(col="Healthy plant\n(>400)") +
  ylab("Dominant strain\naccording to fluorescence ratio") + xlab("Treatment")
gg_ratios_fluor_monotreats

gg_ratios_flour_forpub <- plot_grid(gg_ratios_fluor_monotreats, gg_ratios_fluor_mixtreats, align = "h"
                                    , rel_widths = c(3,5))
gg_ratios_flour_forpub
ggsave("04_WCSvsN2C3_EDA/gg_ratios_flour_forpub.png", gg_ratios_flour_forpub, width=8, height=4)


### Other plots ####
gg_allFluor_simple <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  mutate(Controls=ifelse(protect=="MOCK" | path=="MOCK", "No strains / single strains", "Both strains")) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Healthiness_hsv, group=ratio, col=ratio)) +
  geom_point(aes(pch=ratio)) +
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  scale_color_manual(values=c("0-0"="black", "0-1"="red", "1-0"="blue","1-1"="darkgrey","1-100"="orange","1-10000"="brown","100-1"="lightblue","10000-1"="darkblue"))+
  # geom_smooth(method="lm", se=FALSE)+
  geom_hline(aes(yintercept=400), lty=2)+geom_vline(aes(xintercept=0), lty=2) +
  facet_grid(.~Controls)+
  ylab("Health Score") + xlab("Log2 Fold-change (WCS365:N2C3)")+
  labs(col="Ratio (WCS365:N2C3)",pch="Ratio (WCS365:N2C3)")
gg_allFluor_simple
ggsave("04_WCSvsN2C3_EDA/gg_allFluor_simple.png", gg_allFluor_simple, width=8, height=4)

gg_allFluor_simple_bin <- dat %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  mutate(Aliven=as.numeric(Alive)) %>%
  mutate(Controls=ifelse(protect=="MOCK" | path=="MOCK", "No strains / single strains", "Both strains")) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Aliven)) +
  geom_point(aes(pch=ratio, group=ratio, col=ratio)) +
  geom_smooth(method="glm", method.args=c(family="binomial"), col="black") +
  scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  scale_color_manual(values=c("0-0"="black", "0-1"="red", "1-0"="blue","1-1"="darkgrey","1-100"="orange","1-10000"="brown","100-1"="lightblue","10000-1"="darkblue"))+
  # geom_smooth(method="lm", se=FALSE)+
  # geom_hline(aes(yintercept=400), lty=2)+
  geom_vline(aes(xintercept=0), lty=2) +
  facet_grid(.~Controls)+
  ylab("Probability of healthy plant") + xlab("Log2 Fold-change (WCS365:N2C3)")+
  labs(col="Ratio (WCS365:N2C3)",pch="Ratio (WCS365:N2C3)")
gg_allFluor_simple_bin
ggsave("04_WCSvsN2C3_EDA/gg_allFluor_simple_bin.png", gg_allFluor_simple_bin, width=8, height=4)


gg_allFluor_data_cells <- dat %>%
  # filter(StrainMix !="MOCK-MOCK") %>%
  mutate(inoc_ratios = log10(protect_od/path_od)) %>%
  # select(inoc_ratios) %>% unique()
  mutate(inoc_ratios = ifelse(is.na(inoc_ratios), "No bacteria"
                              , ifelse(inoc_ratios == -Inf, "No WCS365"
                                       , ifelse(inoc_ratios==Inf, "No N2C3", inoc_ratios)))) %>%
  mutate(inoc_ratios = factor(inoc_ratios, levels=c("No bacteria","No WCS365", "-4","-2","0","2","4","No N2C3"))) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Healthiness_hsv)) +
  geom_point(aes(fill=factor(inoc_ratios)),pch=21, col="black") +
  geom_smooth(aes(group=factor(inoc_ratios), col=factor(inoc_ratios)), method="lm", show.legend = FALSE, se=FALSE)+
  geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=400))  +
  facet_grid(.~round(path_cells_log),  labeller = labeller(`round(path_cells_log)` = N2C3_labeller)) +
  ylab("Health Score") +
  xlab("WCS365:N2C3 Log2 Fold Change\n(positive = WCS365 more abundant)\nat END of experiment") +
  labs(fill="Log10 Fold difference\n of Inoculation ratio\nWCS365:N2C3\nat START of experiment") +
  scale_fill_manual(values=c("darkgrey","brown","darkorange","yellow","black","lightblue","blue","darkblue"))+
  scale_color_manual(values=c("darkgrey","brown","darkorange","yellow","black","lightblue","blue","darkblue"))
gg_allFluor_data_cells
ggsave("04_WCSvsN2C3_EDA/gg_allFluor_data_cells.png", gg_allFluor_data_cells, width=15, height=4)

gg_allFluor_data_cells_bin <- dat %>%
  mutate(Aliven=as.numeric(Alive)) %>%
  # filter(StrainMix !="MOCK-MOCK") %>%
  mutate(inoc_ratios = log10(protect_od/path_od)) %>%
  # select(inoc_ratios) %>% unique()
  mutate(inoc_ratios = ifelse(is.na(inoc_ratios), "No bacteria"
                              , ifelse(inoc_ratios == -Inf, "No WCS365"
                                       , ifelse(inoc_ratios==Inf, "No N2C3", inoc_ratios)))) %>%
  mutate(inoc_ratios = factor(inoc_ratios, levels=c("No bacteria","No WCS365", "-4","-2","0","2","4","No N2C3"))) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Aliven)) +
  geom_jitter(aes(fill=factor(inoc_ratios)),pch=21, col="black", height=0.1, width=0) +
  geom_smooth( method="glm", show.legend = FALSE, method.args=c(family="binomial"), col="black")+
  geom_vline(aes(xintercept=0))+
  # geom_hline(aes(yintercept=400))  +
  facet_grid(.~round(path_cells_log),  labeller = labeller(`round(path_cells_log)` = N2C3_labeller)) +
  ylab("Probability of healthy plant") +
  xlab("WCS365:N2C3 Log2 Fold Change\n(positive = WCS365 more abundant)\nat END of experiment") +
  labs(fill="Log10 Fold difference\n of Inoculation ratio\nWCS365:N2C3\nat START of experiment") +
  scale_fill_manual(values=c("darkgrey","brown","darkorange","yellow","black","lightblue","blue","darkblue"))+
  scale_color_manual(values=c("darkgrey","brown","darkorange","yellow","black","lightblue","blue","darkblue"))
gg_allFluor_data_cells_bin
ggsave("04_WCSvsN2C3_EDA/gg_allFluor_data_cells_bin.png", gg_allFluor_data_cells_bin, width=15, height=4)


gg_11ratio_only <- dat %>%
  filter(StrainMix == "WCS365-N2C3") %>%
  filter(ratio=="1-1") %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Healthiness_hsv)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(.~path_od, labeller = labeller(path_od=c(`1e-06`="Each strain inoculation\nAbsorbance=1e-06", `1e-04`="Each strain inoculation\nAbsorbance=1e-04", `0.01`="Each strain inoculation\nAbsorbance=1e-02"))) +
  xlab("Log2 fold-change (WCS365:N2C3)") + ylab("Health Score") +
  geom_hline(aes(yintercept=400), lty=2) + geom_vline(aes(xintercept=0), lty=2)+
  labs(title="Total inoculation absorbance (WCS365+N2C3)")
gg_11ratio_only
ggsave("04_WCSvsN2C3_EDA/gg_11ratio_only.png", gg_11ratio_only, height=4, width=6)

gg_11ratio_only_bin <- dat %>%
  mutate(Aliven=as.numeric(Alive)) %>%
  filter(StrainMix == "WCS365-N2C3") %>%
  filter(ratio=="1-1") %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Aliven)) +
  geom_point() +
  geom_smooth(method="glm", method.args=c(family="binomial")) +
  facet_grid(.~path_od, labeller = labeller(path_od=c(`1e-06`="Each strain inoculation\nAbsorbance=1e-06", `1e-04`="Each strain inoculation\nAbsorbance=1e-04", `0.01`="Each strain inoculation\nAbsorbance=1e-02"))) +
  xlab("Log2 fold-change (WCS365:N2C3)") + ylab("Probability of healthy plant") +
  # geom_hline(aes(yintercept=400), lty=2) +
  geom_vline(aes(xintercept=0), lty=2)+
  labs(title="Total inoculation absorbance (WCS365+N2C3)")
gg_11ratio_only_bin
ggsave("04_WCSvsN2C3_EDA/gg_11ratio_only_bin.png", gg_11ratio_only_bin, height=4, width=6)


gg_11ratio_only_bin_collapsed <- dat %>%
  mutate(Aliven=as.numeric(Alive)) %>%
  filter(StrainMix == "WCS365-N2C3") %>%
  filter(ratio=="1-1") %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Aliven)) +
  geom_jitter(aes(group=factor(round(path_cells_log)), col=factor(round(path_cells_log))), height=0.1, width=0) +
  geom_smooth(method="glm", method.args=c(family="binomial"), col="black") +
  geom_smooth(aes(group=factor(round(path_cells_log)), col=factor(round(path_cells_log))), method="glm", method.args=c(family="binomial"), se=FALSE) +
  # facet_grid(.~path_od, labeller = labeller(path_od=c(`1e-06`="Each strain inoculation\nAbsorbance=1e-06", `1e-04`="Each strain inoculation\nAbsorbance=1e-04", `0.01`="Each strain inoculation\nAbsorbance=1e-02"))) +
  xlab("Log2 fold-change (WCS365:N2C3)") + ylab("Probability of healthy plant") +
  # geom_hline(aes(yintercept=400), lty=2) +
  geom_vline(aes(xintercept=0), lty=2)+
  ylim(0,1) + labs(col="Number of cells (each)\n(log10)")+
  scale_color_manual(values=c("lightgreen","forestgreen","darkgreen"))
gg_11ratio_only_bin_collapsed
ggsave("04_WCSvsN2C3_EDA/gg_11ratio_only_bin_collapsed.png", gg_11ratio_only_bin_collapsed, height=4, width=6)

############# BRMS #############

# Does inoculation ratio or dominance better explain N2C3..?
dat_noNAs <- dat %>% select(Alive, WN_FC_inoc, WN_FC_fluor_offset, total_cells_log, plate, experiment, ratio, StrainMix) %>%
  mutate(WN_FC_inoc = ifelse(is.finite(WN_FC_inoc),WN_FC_inoc,NA)
         , WN_FC_fluor_offset = ifelse(is.finite(WN_FC_fluor_offset),WN_FC_fluor_offset,NA)) %>%
  drop_na() 
dat_noNAs %>% select(ratio, StrainMix) %>% table() # Confirmed that there are mixed treatments only

brm_fluorinoc_bin <- brm(Alive ~ WN_FC_fluor_offset + WN_FC_inoc+total_cells_log, data=dat_noNAs
                     , family="bernoulli"
                     , file="04_WCSvsN2C3_EDA/brm_fluorinoc_bin"
                     , seed=54398)

brm_fluor_bin <- brm(Alive ~ WN_FC_fluor_offset + total_cells_log, data=dat_noNAs
                     , family="bernoulli"
                     , file="04_WCSvsN2C3_EDA/brm_fluor_bin"
                     , seed=54398)
fixef(brm_fluor_bin)

brm_inoc_bin <- brm(Alive ~ WN_FC_inoc + total_cells_log, data=dat_noNAs
                    , family="bernoulli"
                    , file="04_WCSvsN2C3_EDA/brm_inoc_bin"
                    , seed=54398)
fixef(brm_inoc_bin)

# According to the BRMS documentation:
#The WAIC is an approximation of LOO that is faster and easier to compute. However, according to
# Vehtari et al. (2015), LOO may be the preferred method to perform model comparisons.
loo_compare <- LOO(brm_fluorinoc_bin, brm_fluor_bin, brm_inoc_bin)
loo_compare
# fluorescence model is best

############# add predictions to plots #############
# Add predictions to plots
dat_wfluorbrm <- posterior_linpred(brm_fluor_bin, transform=TRUE) %>%
  apply(MARGIN=2, FUN=function(x) c(mean = mean(x), sd = sd(x), quantile(x,prob=0.5),quantile(x, prob=0.05), quantile(x, prob=0.95)
                                    , quantile(x, prob=0.025), quantile(x, prob=0.975))) %>%
  t() %>% as.data.frame() %>%
  rename(pred_mean = mean
         , pred_sd = sd
         , Q50 = `50%`
         , Q5 = `5%`
         , Q95 = `95%`
         , Q2.5 = `2.5%`
         , Q97.5 = `97.5%`) %>%
  bind_cols(brm_fluor_bin$data) %>%
  left_join(dat)

dat_wfluorbrm %>%
  select(protect_od, path_od) %>% table()
dat_wfluorbrm %>%
  select(WN_FC_inoc,total_cells_log) %>% table()


dat_wfluorbrm %>%
  mutate(Aliven = as.numeric(Alive)) %>%
  mutate(total_cells_log_round = factor(round(total_cells_log,1))) %>%
  unite(total_cells_log_round, ratio, col="ribbongroups", remove=FALSE) %>%
  select(ribbongroups)%>% unique()
gg_brmsfit_fluor_health <- dat_wfluorbrm %>%
  mutate(Aliven = as.numeric(Alive)) %>%
  mutate(total_cells_log_round = factor(round(total_cells_log,1))) %>%
  unite(total_cells_log_round, ratio, col="ribbongroups", remove=FALSE) %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  mutate(ratio = factor(ratio, levels=c("10000-1","100-1","1-1","1-100","1-10000"))) %>%
  # mutate(Controls=ifelse(protect=="MOCK" | path=="MOCK", "No strains / single strains", "Both strains")) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Aliven)) +
  geom_jitter(aes(group=ratio, col=ratio), width=0, height=0.1) +
  geom_line(aes(x=WN_FC_fluor_offset, y=pred_mean, col = ratio, lty=total_cells_log_round)) +
  geom_vline(aes(xintercept=0), lty=2, col="black") +
  geom_ribbon(aes(x=WN_FC_fluor_offset, ymin=Q2.5, ymax=Q97.5, fill=ribbongroups, group=ribbongroups), alpha=0.25, show.legend = FALSE) +
  # scale_shape_manual(values=c("0-0"=3, "0-1"=16, "1-0"=21,"1-1"=8,"1-100"=17,"1-10000"=15,"100-1"=24,"10000-1"=22))+
  # scale_shape_manual(values=c(`3.5`=10, `5.2`=19, `5.5`=21, `7.2`=4, `7.5`=8)) +
  scale_color_manual(values=c("0-0"="black", "1-1"="darkgrey","1-100"="orange","1-10000"="darkorange4","100-1"="green","10000-1"="darkgreen"))+
  scale_fill_manual(values=c("5.5_1-1"="darkgrey","3.5_1-1"="darkgrey","7.5_1-1"="darkgrey"
                             ,"7.2_1-100"="orange","5.2_1-100"="orange"
                             ,"7.2_1-10000"="darkorange4","5.2_1-10000"="darkorange4"
                             ,"7.2_100-1"="green","5.2_100-1"="green"
                             ,"7.2_10000-1"="darkgreen","5.2_10000-1"="darkgreen"))+
  scale_linetype_manual(values=c(`3.5`=3, `5.2`=4, `5.5`=6, `7.2`=5, `7.5`=1)) +
    # facet_grid(WN_FC_inoc~.)+
  ylab("Probability of healthy plant") + xlab("Log2 Fold-change fluorescence (WCS365:N2C3)")+
  labs(col="Ratio (WCS365:N2C3)", lty = "Total cells inoculated\n(log10)") +
  ylim(0,1) +xlim(-15,15)
gg_brmsfit_fluor_health
ggsave(filename="04_WCSvsN2C3_EDA/gg_brmsfit_fluor_health.png", gg_brmsfit_fluor_health, height=4, width=7)


##### Numbers for paper ######
sink("04_WCSvsN2C3_EDA/numbers_for_paper.txt")
print("Inoculation ratio (dom) and rate of plant healthiness")
dat %>%
  filter(StrainMix!="MOCK-MOCK") %>%
  mutate(dom = ifelse(path_od>protect_od, "PATH", ifelse(protect_od>path_od, "PROTECT", "EQUAL"))) %>%
  select(dom, WN_FC_fluor_offset, Alive) %>%
  group_by(dom) %>%
  summarise(mean_WN_FC_fluor_offset = mean(WN_FC_fluor_offset, na.rm = TRUE)
            , min_diff = min(WN_FC_fluor_offset, na.rm = TRUE)
            , max_diff = max(WN_FC_fluor_offset, na.rm = TRUE)
            , healthy = sum(Alive)
            , rate = sum(Alive)/n()
            , raten = 1-rate)
print("ONLY MONOCULTURE TREATMENTS") 
dat %>%
  # filter(path=="MOCK" | protect == "MOCK") %>%
  # mutate(dom = ifelse(path_od>protect_od, "PATH", ifelse(protect_od>path_od, "PROTECT", "EQUAL"))) %>%
  select(StrainMix, WN_FC_fluor_offset, Alive) %>%
  group_by(StrainMix) %>%
  summarise(mean_WN_FC_fluor_offset = mean(WN_FC_fluor_offset, na.rm = TRUE)
            , min_diff = min(WN_FC_fluor_offset, na.rm = TRUE)
            , max_diff = max(WN_FC_fluor_offset, na.rm = TRUE)
            , healthy = sum(Alive)
            , rate = sum(Alive)/n()
            , raten = 1-rate)

dat %>%
  # select(ratio,WN_FC_CFU) %>%
  filter(is.finite(WN_FC_CFU)) %>%
  mutate(dom = ifelse(path_od>protect_od, "PATH", ifelse(protect_od>path_od, "PROTECT", "EQUAL"))) %>%
  select(dom, WN_FC_CFU) %>%
  drop_na() %>%
  group_by(dom) %>%
  summarise(mean_WN_FC_CFU = mean(WN_FC_CFU, na.rm = TRUE))
sink()
