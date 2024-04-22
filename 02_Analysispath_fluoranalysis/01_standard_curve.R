library(gridExtra)
library(tidyverse)
library(ggpubr)
dir.create("01_standard_curve")

dat_all <- read.csv("00_raw_data/dat/standard_curve_data.csv") %>%
  filter(!(is.na(read_fluorcell) | is.na(read_otherfluor_fluorcell)))
dat_mix <- read.csv("00_raw_data/dat/standard_curve_mixeddata_wide.csv") 

# 
blankonly <- dat_all %>% filter(OD==0) %>%
  group_by(volume, culture_stage, fluor_cell) %>%
  summarise(Abs600_blanks = mean(Abs600), read_fluorcell_blanks = mean(read_fluorcell)) %>% ungroup()
blankonly_both <- blankonly %>% mutate(fluor_cell=ifelse(fluor_cell=="crim","neon","crim")) %>%
  rename(read_otherfluor_fluorcell_blanks = read_fluorcell_blanks) %>% select(-Abs600_blanks) %>%
  right_join(blankonly)
dat_adj <- left_join(dat_all, blankonly_both) %>%
  mutate(Abs600_blanked = Abs600 - Abs600_blanks, read_fluorcell_blanked = read_fluorcell-read_fluorcell_blanks
         , read_otherfluor_fluorcell_blanked = read_otherfluor_fluorcell-read_otherfluor_fluorcell_blanks
  ) %>%
  unite(culture_stage, volume, col=experiment, remove=FALSE, sep="_") %>%
  mutate(CFU_fromOD = OD*5E8*volume)


#### Absorbance 600 to OD conversion ###########
gg_ODvsAbs <- dat_adj %>%
  unite(culture_stage, volume, col="experiment", remove=FALSE) %>%
  ggplot() +
  geom_point(aes(x=OD, y=Abs600_blanked, col=experiment))+
  # geom_vline(aes(xintercept=0.2)) +
  # geom_hline(aes(yintercept=0.05)) +
  scale_x_log10() + scale_y_log10() +
  labs(col="CultureStage_volume") +
  xlab("OD (blanked, 1cm cuvette)") + ylab("Abs600 (blanked, variable volumes)")
gg_ODvsAbs
ggsave("01_standard_curve/gg_ODvsAbs.png", gg_ODvsAbs, width=6, height=4)

# The exponential culture with 275uL doesn't have absorbance measurements, so let's estimate them using this
# Use the other 275 one
dat_on275_forcurve <- dat_adj %>%
  filter(volume==275, culture_stage == "overnight", OD>0.01)
lm_on275_forcurve <- lm(Abs600_blanked ~ OD, data=dat_on275_forcurve)
dat_adj <- lm_on275_forcurve$coefficients %>% t() %>%as.data.frame() %>%
  rename(abs_conv_intercept = `(Intercept)`, abs_conv_slope = OD) %>%
  cbind(dat_adj) %>%
  mutate(Abs600_blanked = ifelse(is.na(Abs600_blanked), OD*abs_conv_slope+abs_conv_intercept, Abs600_blanked)) %>%
  select(-abs_conv_intercept, -abs_conv_slope)

#### Converting between fluor and OD/Abs ####

# OD vs raw fluorescence
gg_summary_stcurve_ONandExp <- dat_adj %>%
  select(OD, read_fluorcell, read_otherfluor_fluorcell, fluor_cell, experiment) %>%
  pivot_longer(-c(OD, fluor_cell,experiment), names_to="ReadType", values_to="read") %>%
  mutate(
    read_range = ifelse(ReadType=="read_fluorcell", ifelse(fluor_cell=="crim","crim","neon"),
                        ifelse(ReadType=="read_otherfluor_fluorcell", ifelse(fluor_cell=="crim","neon","crim"), NA))
  ) %>%
  mutate(read_range_fluor = ifelse(read_range=="crim","Crimson spectra", "Neon spectra")) %>%
  mutate(cell_type_fluor = ifelse(fluor_cell=="crim", "Crimson Cells", "Neon Cells")) %>%
  unite(read_range_fluor, cell_type_fluor, remove=FALSE, col="Category", sep=" reading\n ") %>%
  # select(Category) %>% unique()
  mutate(Category = factor(Category, levels=c("Crimson spectra reading\n Crimson Cells","Neon spectra reading\n Neon Cells"
                                              ,"Crimson spectra reading\n Neon Cells","Neon spectra reading\n Crimson Cells"))) %>%
  ggplot(aes(x=OD, y=read)) +
  geom_point(aes(fill=fluor_cell), pch=21, cex=0.5,col="NA") +
  geom_point(aes(fill=fluor_cell, pch=experiment), col="NA", show.legend = FALSE) +
  geom_smooth(aes(col=read_range, lty=experiment), method="lm") +
  facet_wrap(.~Category, scales="free", strip.position = "right")+
  ylab("Fluorescence (raw)") + xlab("OD") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  scale_fill_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(fill="Cell strain", col="Reading em/exc\nrange", lty="CultureStage_volume") +
  scale_shape_manual(values=c(21, 22, 23, 24))+
  scale_linetype_manual(values=c(1,2,3,6))
gg_summary_stcurve_ONandExp
ggsave(filename = "01_standard_curve/gg_summary_stcurve_ONandExp.png", gg_summary_stcurve_ONandExp, height=6, width=10)

# Absorbance vs raw fluorescence
gg_summary_stcurve_ONandExp_abs <- dat_adj %>%
  select(OD, Abs600_blanked, read_fluorcell, read_otherfluor_fluorcell, fluor_cell, experiment) %>%
  pivot_longer(-c(OD, Abs600_blanked, fluor_cell,experiment), names_to="ReadType", values_to="read") %>%
  mutate(
    read_range = ifelse(ReadType=="read_fluorcell", ifelse(fluor_cell=="crim","crim","neon"),
                        ifelse(ReadType=="read_otherfluor_fluorcell", ifelse(fluor_cell=="crim","neon","crim"), NA))
  ) %>%
  mutate(read_range_fluor = ifelse(read_range=="crim","Crimson spectra", "Neon spectra")) %>%
  mutate(cell_type_fluor = ifelse(fluor_cell=="crim", "Crimson Cells", "Neon Cells")) %>%
  unite(read_range_fluor, cell_type_fluor, remove=FALSE, col="Category", sep=" reading\n ") %>%
  mutate(Category = factor(Category, levels=c("Crimson spectra reading\n Crimson Cells","Neon spectra reading\n Neon Cells"
                                              ,"Crimson spectra reading\n Neon Cells","Neon spectra reading\n Crimson Cells"))) %>%
  ggplot(aes(x=Abs600_blanked, y=read)) +
  geom_point(aes(fill=fluor_cell), pch=21, cex=0.5,col="NA") +
  geom_point(aes(fill=fluor_cell, pch=experiment), col="NA", show.legend = FALSE) +
  geom_smooth(aes(col=read_range, lty=experiment), method="lm") +
  facet_wrap(.~Category, scales="free", strip.position = "right")+
  ylab("Fluorescence (raw)") + xlab("Absorbance 600nm") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  scale_fill_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(fill="Cell strain", col="Reading em/exc\nrange", lty="CultureStage_volume") +
  scale_shape_manual(values=c(21, 22, 23, 24))+
  scale_linetype_manual(values=c(1,2,3,6))
gg_summary_stcurve_ONandExp_abs

ggsave(filename = "01_standard_curve/gg_summary_stcurve_ONandExp_abs.png", gg_summary_stcurve_ONandExp_abs, height=6, width=10)

# Raw fluorescence vs raw CFUs
gg_summary_stcurve_ONandExp_bycells <- dat_adj %>%
  select(OD, CFU_fromOD, read_fluorcell, read_otherfluor_fluorcell, fluor_cell, experiment) %>%
  pivot_longer(-c(OD, CFU_fromOD, fluor_cell,experiment), names_to="ReadType", values_to="read") %>%
  mutate(
    read_range = ifelse(ReadType=="read_fluorcell", ifelse(fluor_cell=="crim","crim","neon"),
                        ifelse(ReadType=="read_otherfluor_fluorcell", ifelse(fluor_cell=="crim","neon","crim"), NA))
  ) %>%
  mutate(read_range_fluor = ifelse(read_range=="crim","Crimson spectra", "Neon spectra")) %>%
  mutate(cell_type_fluor = ifelse(fluor_cell=="crim", "Crimson Cells", "Neon Cells")) %>%
  unite(read_range_fluor, cell_type_fluor, remove=FALSE, col="Category", sep=" reading\n ") %>%
  mutate(Category = factor(Category, levels=c("Crimson spectra reading\n Crimson Cells","Neon spectra reading\n Neon Cells"
                                              ,"Crimson spectra reading\n Neon Cells","Neon spectra reading\n Crimson Cells"))) %>%
  ggplot(aes(x=CFU_fromOD, y=read)) +
  geom_point(aes(fill=fluor_cell), pch=21, cex=0.5,col="NA") +
  geom_point(aes(fill=fluor_cell, pch=experiment), col="NA", show.legend = FALSE) +
  geom_smooth(aes(col=read_range, lty=experiment), method="lm") +
  facet_wrap(.~Category, scales="free", strip.position = "right")+
  ylab("Fluorescence (raw)") + xlab("Approx. CFU") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  scale_fill_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(fill="Cell strain", col="Reading em/exc\nrange", lty="CultureStage_volume") +
  scale_shape_manual(values=c(21, 22, 23, 24))+
  scale_linetype_manual(values=c(1,2,3,6))
gg_summary_stcurve_ONandExp_bycells
ggsave(filename = "01_standard_curve/gg_summary_stcurve_ONandExp_bycells.png", gg_summary_stcurve_ONandExp_bycells, height=6, width=10)


# Get boundaries for standard curves; manually selected
stcurve_bounds <- data.frame(odabs_type = c("OD","Abs600_blanked","CFU_fromOD")
                             , lwrBound = c(0.025, 0.005, 1E9)
                             , uprBound = c(0.75, 0.3, 1E11))

# Log the values
# pivot longer makes it so you get the same fluor values, repeated across each 'x' value (OD, Abs600, cell counts)
gg_summary_stcurve_loglog_all <- dat_adj %>%
  select(OD, Abs600_blanked, CFU_fromOD, read_fluorcell_blanked, fluor_cell, experiment) %>%
  pivot_longer(-c(read_fluorcell_blanked, fluor_cell, experiment), names_to="odabs_type", values_to = "odabs_val") %>%
  full_join(stcurve_bounds) %>%
  mutate(Reading = ifelse(fluor_cell=="crim","Crimson cells\nCrimson ex/em", "Neon cells\nNeon ex/em")) %>%
  mutate(odabs_type_adj = ifelse(odabs_type=="OD", "OD600 blanked (cuvette)"
                                 , ifelse(odabs_type=="Abs600_blanked", "Abs600 blanked (96 well plate)\nvolume variable", 
                                          ifelse(odabs_type=="CFU_fromOD", "Estimated CFU (from OD)\nOD*5E8*volume", NA)) )) %>%
  ggplot(aes(x=odabs_val, y=(read_fluorcell_blanked))) +
  geom_point(aes(col=fluor_cell, pch=experiment)) +
  geom_vline(aes(xintercept=lwrBound), lty=2) +
  geom_vline(aes(xintercept=uprBound), lty=2) +
  facet_grid(Reading~odabs_type_adj, scales="free", switch = "x")+
  ylab("Blanked Fluorescence (log scale)") + xlab("Known cell density (log scale)") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(col="Cell strain", pch="CultureStage_volume") +
  scale_shape_manual(values=c(15,22,19,21))+
  scale_x_log10() +
  scale_y_log10() 
gg_summary_stcurve_loglog_all
ggsave("01_standard_curve/gg_summary_stcurve_loglog_all.png", gg_summary_stcurve_loglog_all, height=5, width=10)


gg_summary_stcurve_loglog_simple <- dat_adj %>%
  select(OD, Abs600_blanked, CFU_fromOD, read_fluorcell_blanked, fluor_cell, experiment) %>%
  pivot_longer(-c(read_fluorcell_blanked, fluor_cell, experiment), names_to="odabs_type", values_to = "odabs_val") %>%
  full_join(stcurve_bounds) %>%
  mutate(Reading = ifelse(fluor_cell=="crim","Crimson cells\nCrimson ex/em", "Neon cells\nNeon ex/em")) %>%
  mutate(odabs_type_adj = ifelse(odabs_type=="OD", "OD600 blanked (cuvette)"
                                 , ifelse(odabs_type=="Abs600_blanked", "Abs600 blanked (96 well plate)\nvolume variable", 
                                          ifelse(odabs_type=="CFU_fromOD", "Estimated CFU (from OD)\nOD*5E8*volume", NA)) )) %>%
  filter(odabs_type == "Abs600_blanked") %>% #### Use only ABs600 here because it is technically the "most accurate"
  ggplot(aes(x=odabs_val, y=(read_fluorcell_blanked))) +
  geom_point(aes(col=fluor_cell, pch=experiment)) +
  geom_vline(aes(xintercept=lwrBound), lty=2) +
  geom_vline(aes(xintercept=uprBound), lty=2) +
  facet_grid(Reading~., scales="free", switch = "x")+
  ylab("Blanked Fluorescence (log scale)") + xlab("Known cell density\n(Absorbance 600, log scale)") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(col="Cell strain", pch="Replicate attempt") +
  scale_shape_manual(values=c(15,22,19,21))+
  scale_x_log10() +
  scale_y_log10() 
gg_summary_stcurve_loglog_simple
ggsave("01_standard_curve/gg_summary_stcurve_loglog_simple.png", gg_summary_stcurve_loglog_simple, height=5, width=6)


gg_summary_stcurve_loglog_simple2 <- dat_adj %>%
  select(OD, Abs600_blanked, CFU_fromOD, read_fluorcell_blanked, fluor_cell, experiment) %>%
  pivot_longer(-c(read_fluorcell_blanked, fluor_cell, experiment), names_to="odabs_type", values_to = "odabs_val") %>%
  full_join(stcurve_bounds) %>%
  mutate(Reading = ifelse(fluor_cell=="crim","Crimson cells\nCrimson ex/em", "Neon cells\nNeon ex/em")) %>%
  mutate(odabs_type_adj = ifelse(odabs_type=="OD", "OD600 blanked (cuvette)"
                                 , ifelse(odabs_type=="Abs600_blanked", "Abs600 blanked (96 well plate)\nvolume variable", 
                                          ifelse(odabs_type=="CFU_fromOD", "Estimated CFU (from OD)\nOD*5E8*volume", NA)) )) %>%
  filter(odabs_type == "Abs600_blanked") %>%
  ggplot(aes(x=odabs_val, y=(read_fluorcell_blanked))) +
  geom_point(aes(col=fluor_cell, pch=experiment)) +
  geom_smooth(aes(col=fluor_cell), se=FALSE, method="lm") +
  # facet_grid(Reading~., scales="free", switch = "x")+
  ylab("Blanked Fluorescence (log scale)") + xlab("Known cell density\n(Absorbance 600, log scale)") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(col="Cell strain", pch="Replicate attempt") +
  scale_shape_manual(values=c(15,22,19,21))+
  scale_x_log10() +
  scale_y_log10() 
gg_summary_stcurve_loglog_simple2
ggsave("01_standard_curve/gg_summary_stcurve_loglog_simple2.png", gg_summary_stcurve_loglog_simple2, height=3, width=5)


gg_summary_stcurve_loglog_all_BG <- dat_adj %>%
  select(OD, Abs600_blanked, CFU_fromOD, read_otherfluor_fluorcell_blanked, fluor_cell, experiment) %>%
  pivot_longer(-c(read_otherfluor_fluorcell_blanked, fluor_cell, experiment), names_to="odabs_type", values_to = "odabs_val") %>%
  mutate(Reading = ifelse(fluor_cell=="neon","Neon cells\nCrimson ex/em", "Crimson cells\nNeon ex/em")) %>%
  mutate(odabs_type_adj = ifelse(odabs_type=="OD", "OD600 blanked (cuvette)"
                                 , ifelse(odabs_type=="Abs600_blanked", "Abs600 blanked (96 well plate)\nvolume variable", 
                                          ifelse(odabs_type=="CFU_fromOD", "Estimated CFU (from OD)\nOD*5E8*volume", NA)) )) %>%
  ggplot(aes(x=odabs_val, y=(read_otherfluor_fluorcell_blanked))) +
  geom_point(aes(col=fluor_cell, pch=experiment)) +
  facet_grid(Reading~odabs_type_adj, scales="free", switch = "x")+
  ylab("Blanked Fluorescence (log scale)") + xlab("Known cell density (log scale)") +
  scale_color_manual(values=c(crim="darkorange", neon="limegreen"))+
  labs(col="Cell strain", pch="CultureStage_volume") +
  scale_shape_manual(values=c(15,22,19,21))+
  scale_x_log10() +
  scale_y_log10() 
gg_summary_stcurve_loglog_all_BG
ggsave("01_standard_curve/gg_summary_stcurve_loglog_all_BG.png", gg_summary_stcurve_loglog_all_BG, height=5, width=10)



dat_bgonly <- dat_adj %>%
  select(OD, experiment, rep, read_otherfluor_fluorcell_blanked, fluor_cell) %>%
  mutate(fluor_cell = ifelse(fluor_cell=="crim","neon","crim")) %>%
  rename(read_fluor_otherfluorcell_blanked = read_otherfluor_fluorcell_blanked) 
dat_adj <- dat_adj %>%
  full_join(dat_bgonly) %>%
  filter(!(Abs600_blanked>0.01 & (read_fluorcell_blanked)<=1.2E5 & fluor_cell=="neon")) ### REMOVING THAT RANDOM POINT- if you plot there's a weird outlier
  

gg_summary_stcurve_loglog_BGfluor <- dat_adj %>%
  # mutate(fluor_read = ifelse(fluor_cell=="crim","Neon ex/em spectra reading\nof crimson cells", "Crimson ex/em spectra reading\nof neon cells")) %>%
  mutate(fluor_reading = ifelse(fluor_cell=="crim","Crimson ex/em spectra reading", "Neon ex/em spectra reading")) %>%
  ggplot(aes(x=read_fluorcell_blanked, y=read_fluor_otherfluorcell_blanked)) +
  geom_point(aes(pch=experiment)) +
  # geom_hline(aes(yintercept=lwrBoundy), lty=2) +
  # geom_vline(aes(xintercept=lwrBound), lty=2) +
  facet_grid(.~fluor_reading, scales="free")+
  ylab("Fluorescence reading of non-matching cell (log scale)") + xlab("Fluoresecence reading of matching cell (log scale)") +
  labs(col="Cell strain", pch="CultureStage_volume") +
  geom_hline(aes(yintercept=1.5E4)) +
  scale_shape_manual(values=c(15,22,19,21))+
  scale_x_log10() +
  scale_y_log10() 
gg_summary_stcurve_loglog_BGfluor
ggsave("01_standard_curve/gg_summary_stcurve_loglog_BGfluor.png", gg_summary_stcurve_loglog_BGfluor, height=5, width=10)

######## Get linear fit ##########
dat_neon <- dat_adj %>%
  # filter((Abs600_blanked>0.01 & (read_fluorcell_blanked)<=1.2E5 & fluor_cell=="neon")) %>% # This is for checking that random point is gone
  select(OD, Abs600_blanked, CFU_fromOD, read_fluorcell_blanked, fluor_cell, experiment) %>%
  pivot_longer(-c(read_fluorcell_blanked, fluor_cell, experiment), names_to="odabs_type", values_to = "odabs_val") %>%
  full_join(stcurve_bounds) %>%
  filter(odabs_val> lwrBound,odabs_val<uprBound ) %>%
  filter(fluor_cell=="neon") %>%
  filter(odabs_type=="Abs600_blanked") %>%
  mutate(read_fluorcell_blanked_log10 = log10(read_fluorcell_blanked), odabs_val_log10 = log10(odabs_val))
dat_crim <- dat_adj %>%
  select(OD, Abs600_blanked, CFU_fromOD, read_fluorcell_blanked, fluor_cell, experiment) %>%
  pivot_longer(-c(read_fluorcell_blanked, fluor_cell, experiment), names_to="odabs_type", values_to = "odabs_val") %>%
  full_join(stcurve_bounds) %>%
  filter(odabs_val> lwrBound,odabs_val<uprBound ) %>%
  filter(fluor_cell=="crim")%>%
  filter(odabs_type=="Abs600_blanked")%>%
  mutate(read_fluorcell_blanked_log10 = log10(read_fluorcell_blanked), odabs_val_log10 = log10(odabs_val))

## Force through zero if blanked; do a separate curve for every experiment
dat_neon1 <- dat_neon %>% filter(experiment=="overnight_275") 
dat_neon2 <- dat_neon %>% filter(experiment=="exponential_275")
dat_neon3 <- dat_neon %>% filter(experiment=="overnight_200")
dat_neon4 <- dat_neon %>% filter(experiment=="exponential_200")
neon_fit1 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_neon1)
neon_fit2 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_neon2)
neon_fit3 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_neon3)
neon_fit4 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_neon4)
# neon_fitlog10 <- lm(read_fluorcell_blanked_log10 ~ (-1+odabs_val_log10), data=dat_neon)
neon_fit1_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_neon1)
neon_fit2_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_neon2)
neon_fit3_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_neon3)
neon_fit4_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_neon4)
neon_fitlog10_winter <- lm(read_fluorcell_blanked_log10 ~ odabs_val_log10, data=dat_neon)
# Look at fit


dat_crim1 <- dat_crim %>% filter(experiment=="overnight_275")
dat_crim2 <- dat_crim %>% filter(experiment=="exponential_275")
dat_crim3 <- dat_crim %>% filter(experiment=="overnight_200")
dat_crim4 <- dat_crim %>% filter(experiment=="exponential_200")
crim_fit1 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_crim1)
crim_fit2 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_crim2)
crim_fit3 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_crim3)
crim_fit4 <- lm(read_fluorcell_blanked ~ (-1+odabs_val), data=dat_crim4)
# crim_fitlog10 <- lm(read_fluorcell_blanked_log10 ~ (-1+odabs_val_log10), data=dat_crim)
crim_fit1_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_crim1)
crim_fit2_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_crim2)
crim_fit3_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_crim3)
crim_fit4_winter <- lm(read_fluorcell_blanked ~ odabs_val, data=dat_crim4)
crim_fitlog10_winter <- lm(read_fluorcell_blanked_log10 ~ odabs_val_log10, data=dat_crim)

# make into df
neon_df <- t(data.frame(neon_fit1$coefficients,
                        neon_fit2$coefficients,
                        neon_fit3$coefficients,
                        neon_fit4$coefficients)) %>%
  as.data.frame() %>%
  mutate(experiment = c("overnight_275","exponential_275", "overnight_200","exponential_200")) %>%
  mutate(fluor="neon") %>%
  # rename(od_conv_slope=read_fluorcell_blanked, od_conv_intercept = `(Intercept)`) 
rename(od_conv_slope=odabs_val) 
crim_df <- t(data.frame(crim_fit1$coefficients,
                        crim_fit2$coefficients,
                        crim_fit3$coefficients,
                        crim_fit4$coefficients)) %>%
  as.data.frame() %>%
  mutate(experiment = c("overnight_275","exponential_275", "overnight_200","exponential_200")) %>%
  mutate(fluor="crim") %>%
  # rename(od_conv_slope=read_fluorcell_blanked, od_conv_intercept = `(Intercept)`) 
rename(od_conv_slope=odabs_val) 


# make into df
neon_df_winter <- t(data.frame(neon_fit1_winter$coefficients,
                               neon_fit2_winter$coefficients,
                               neon_fit3_winter$coefficients,
                               neon_fit4_winter$coefficients,
                               neon_fitlog10_winter$coefficients)) %>%
  as.data.frame() %>%
  mutate(experiment = c("overnight_275_winter","exponential_275_winter", "overnight_200_winter","exponential_200_winter", "allLog10_winter")) %>%
  mutate(fluor="neon") %>%
  rename(od_conv_slope=odabs_val, od_conv_intercept = `(Intercept)`)
  # rename(od_conv_slope=OD) 
crim_df_winter <- t(data.frame(crim_fit1_winter$coefficients,
                               crim_fit2_winter$coefficients,
                               crim_fit3_winter$coefficients,
                               crim_fit4_winter$coefficients,
                               crim_fitlog10_winter$coefficients)) %>%
  as.data.frame() %>%
  mutate(experiment = c("overnight_275_winter","exponential_275_winter", "overnight_200_winter","exponential_200_winter","allLog10_winter")) %>%
  mutate(fluor="crim") %>%
  rename(od_conv_slope=odabs_val, od_conv_intercept = `(Intercept)`)
  # rename(od_conv_slope=OD) 

allConv <-full_join(neon_df, crim_df) %>% full_join(neon_df_winter) %>% full_join(crim_df_winter) %>%
  mutate(od_conv_intercept = ifelse(is.na(od_conv_intercept), 0, od_conv_intercept))
write.table(allConv, file="01_standard_curve/conversion_models.txt", quote=FALSE, row.names = FALSE, sep="\t")
save(allConv, file = "01_standard_curve/allConv.RData")

####### Get for autofluorescence ##########

stcurve_fluorbounds <- data.frame(fluor_cell=c("crim","neon")
                                  , lwrBoundy = c(1.5E4,1.5E4))
gg_background_fluor <- dat_adj %>%
  filter(OD!=0) %>%
  left_join(stcurve_fluorbounds) %>%
  mutate(fluor_cell = ifelse(fluor_cell=="crim", "Crimson cells", "Neon cells")) %>%
  ggplot() +
  geom_point(aes(x=(read_fluorcell_blanked), y=(read_otherfluor_fluorcell_blanked))) +
  facet_grid(.~fluor_cell, scales="free_x") +
  geom_hline(aes(yintercept=lwrBoundy)) +
  # xlab("Fluorescence of NEON\n from neon-cells (blanked)") +ylab("Fluorescence of CRIMSON\n from neon-cells (blanked)")+
  xlab("Blanked fluorescence reading\nthat matches fluorescent cell type\n(log scale)") +
  ylab("Fluoresence reading of\nopposite fluorophore\n(log scale)")+
  scale_x_log10() + scale_y_log10()
gg_background_fluor
ggsave(filename = "01_standard_curve/gg_background_fluor.png", gg_background_fluor, height=4, width=6)


dat_crimcells_neonbackground <-  dat_adj %>%
  left_join(stcurve_fluorbounds) %>%
  filter(read_otherfluor_fluorcell_blanked>lwrBoundy) %>%
  filter(fluor_cell=="crim") %>%
  mutate(read_otherfluor_fluorcell_blanked_log10 = log10(read_otherfluor_fluorcell_blanked)
         , read_fluorcell_blanked_log10 = log10(read_fluorcell_blanked))
dat_neoncells_crimbackground <-  dat_adj %>%
  left_join(stcurve_fluorbounds) %>%
  filter(read_otherfluor_fluorcell_blanked>lwrBoundy) %>%
  filter(fluor_cell=="neon")%>%
  mutate(read_otherfluor_fluorcell_blanked_log10 = log10(read_otherfluor_fluorcell_blanked)
         , read_fluorcell_blanked_log10 = log10(read_fluorcell_blanked))

## Okay, so we can subtract baseline excitation from crimson measurements using a simple blank.
# But for neon measurements, we need to subtract a non-flat amount based on how many crimson cells there are
# Calculate crimson first. Then calculate neon.

# See if there is background crimson when there are a lot of neon cells
# bgcrimfluor_fromneoncells<- lm(mean_bg_fluor_raw_blanked ~ mean_raw_blanked, data=dat_neoncells_crimbackground_nozero)
# bgcrimfluor_fromneoncells<- lm(read_otherfluor_fluorcell_blanked ~ 1, data=dat_neoncells_crimbackground)
# ggplot(dat_neoncells_crimbackground) +
  # geom_point(aes(x=mean_read_fluorcell_blanked, y=mean_read_bgfluorcell_blanked)) + 
  # geom_smooth(aes(x=mean_read_fluorcell_blanked, y=mean_read_bgfluorcell_blanked), method="lm")
# NO EFFECT BECAUSE NO POINTS ABOVE DETECTION THRESHOLD

# So, that means we can calculate crimson cells directly from crimson fluor

# Then, calculate hypothetical background fluorescence when there are two types of cells
bgneonfluor_fromcrimcells<- lm(read_otherfluor_fluorcell_blanked ~ -1 + read_fluorcell_blanked, data=dat_crimcells_neonbackground)
bgneonfluor_fromcrimcells_winter<- lm(read_otherfluor_fluorcell_blanked ~ 1 + read_fluorcell_blanked, data=dat_crimcells_neonbackground)
bgneonfluor_fromcrimcells
ggplot(dat_crimcells_neonbackground) +
  geom_point(aes(x=read_fluorcell_blanked, y=read_otherfluor_fluorcell_blanked)) +
  geom_abline(aes(slope=0.07782, intercept=814.98225), col="red") +
  geom_abline(aes(slope=0.07841, intercept=0), col="blue") 

# Finally, calcuate the actual neon ~ od rule for neon cells
# NOTE: here, there are no crimson cells which means we can just use direct measurements
#THIS IS BELOW

### Make bg coefficients into tables
# 
# bgcrimfluor_fromneoncells_fit_df <- t(bgcrimfluor_fromneoncells$coefficients) %>%
#   as.data.frame() %>%
#   mutate(fluor="crim") %>% ###### THIS IS THE BACKGROUND FOR OPPOSITE
#   rename(bg_conv_intercept = `(Intercept)`) %>%
#   mutate(bg_conv_slope=0)
# rename(bg_conv_intercept = `(Intercept)`, bg_conv_slope=mean_read_fluorcell_blanked)

# Make mock for crim
bgcrimfluor_fromneoncells_fit_df <- data.frame(bg_conv_intercept=0, bg_conv_slope=1, version=c("winter", "intercept0"),fluor="crim")

bgneonfluor_fromcrimcells_fit_df <- t(data.frame(c(bgneonfluor_fromcrimcells_winter$coefficients, version="winter"),
                                                 c(0,bgneonfluor_fromcrimcells$coefficients, version="intercept0"))) %>%
  as.data.frame() %>%
  mutate(fluor="neon") %>%
  as_tibble() %>%
  rename(bg_conv_intercept = `(Intercept)`, bg_conv_slope=read_fluorcell_blanked) %>%
  mutate(bg_conv_intercept = as.numeric(bg_conv_intercept), bg_conv_slope = as.numeric(bg_conv_slope))


allConv_bg <-full_join(bgcrimfluor_fromneoncells_fit_df, bgneonfluor_fromcrimcells_fit_df)
# allConv_bg <-full_join(bgneonfluor_fromcrimcells_fit_df,)
write.table(allConv_bg, file="01_standard_curve/conversion_models_bg.txt", quote=FALSE, row.names = FALSE, sep="\t")
save(allConv_bg, file = "01_standard_curve/allConv_bg.RData")

#### Check standard curve and bg fluor fits ####

allConv_long <- allConv %>%
  separate(experiment, into=c("experiment","inter"), remove=FALSE, sep="_win") %>%
  mutate(inter = ifelse(is.na(inter),"No intercept","Variable intercept")) %>%
  # filter(experiment !="allLog10") %>%
  rename(fluor_cell=fluor)
# Make a make set so I can draw on a log10 scale below
mock_allLog10_inter <- data.frame(fluor_cell=c("neon","neon","crim","crim")
                                  , Abs600_blanked = c(0.01, 0.1, 0.01, 0.1)) %>%
  full_join(allConv_long, relationship="many-to-many") %>%
  filter(experiment =="allLog10") %>%
  mutate(estFluorlog = log10(Abs600_blanked)*od_conv_slope + od_conv_intercept
         , estFluor = 10^estFluorlog) %>%
  select(fluor_cell, Abs600_blanked, estFluor)
  
mock_dat <- data.frame(fluor_cell=c("neon","neon","crim","crim")
           , Abs600_blanked = c(0.01, 0.1, 0.01, 0.1)) %>%
  full_join(allConv_long, relationship="many-to-many") %>%
  filter(experiment !="allLog10") %>%
  mutate(estFluor = Abs600_blanked*od_conv_slope + od_conv_intercept) %>%
  select(fluor_cell, Abs600_blanked, estFluor, experiment, inter)

gg_stcurve_fits <- dat_adj %>%
  # full_join(allConv_long, relationship = "many-to-many") %>%
  ggplot(aes(x=Abs600_blanked, y=read_fluorcell_blanked)) +
  geom_point()+
  geom_smooth(data=mock_dat, aes(x=Abs600_blanked, y=estFluor, lty=inter, col=inter), alpha=0.5, method="lm", fullrange=TRUE)+
  geom_smooth(data=mock_allLog10_inter, aes(x=Abs600_blanked, y=estFluor, col="Pooled all data", lty="Pooled all data"), alpha=0.5, method="lm", fullrange=TRUE)+
  facet_grid(fluor_cell~experiment, scales="free") +
  scale_x_log10() +
  scale_y_log10()+
  ylab("Fluorescence (log scale)") + xlab("Abs600 blanked (log scale)")+
  labs(col="Model parameters",lty="Model parameters")
gg_stcurve_fits
ggsave("01_standard_curve/gg_stcurve_fits.png",gg_stcurve_fits,width=10, height=6)




allConv_slopeconstants <- allConv %>%
  pivot_wider(id_cols=c(experiment), names_from=fluor, values_from=od_conv_slope) %>%
  mutate(slope_ratio = (crim/neon)
        , FC_slope_ratio = log2(crim/neon))

#### Mixed populations ####
# Get blanks for mixed pop standards
blanks_datmix <- dat_mix %>%
  filter(prop_crim==0 & prop_neon==0) %>%
  group_by(experiment) %>%
  summarise(read_crim_blanks = mean(read_crim), read_neon_blanks = mean(read_neon))
allConv_wide_2 <- allConv %>%
  pivot_wider(id_cols="experiment", names_from="fluor", values_from=c("od_conv_intercept", "od_conv_slope")) %>%
  rename(stcurve=experiment) %>%
  mutate(joinyby="yes") %>%
  filter(stcurve !="allLog10_winter")
allConv_bg_wide_2 <- allConv_bg%>%
  full_join(data.frame(bg_conv_intercept=c(0,0)
                       , bg_conv_slope = c(1,1)
                       , version = c("Nobgadj","Nobgadj")
                       , fluor=c("crim","neon"))) %>%
  filter(version!="intercept0") %>%
  pivot_wider(id_cols="version", names_from="fluor", values_from=c("bg_conv_intercept","bg_conv_slope")) %>%
  mutate(joinyby="yes") 

dat_mix_adj <- dat_mix %>%
  mutate(ratio_real = (prop_crim*100+1)/(prop_neon*100+1)
         , ratio_od = OD_crim/OD_neon
         , FC_ratio_real = log2(ratio_real)
         , FC_ratio_od = log2(ratio_od)) %>%
  full_join(blanks_datmix) %>%
  mutate(read_crim_blanked = read_crim-read_crim_blanks, read_neon_blanked = read_neon-read_neon_blanks
) %>%
  mutate(joinyby="yes") %>%
  full_join(allConv_bg_wide_2, relationship="many-to-many") %>%
  # rowwise() %>%
  mutate(read_neon_blanked_bgadj = read_neon_blanked- (read_crim_blanked*bg_conv_slope_neon + bg_conv_intercept_neon)) %>%
  # ungroup() %>%
  mutate(read_crim_blanked_zeroed = ifelse(read_crim_blanked<0, 0, read_crim_blanked)
         , read_neon_blanked_bgadj_zeroed = ifelse(read_neon_blanked_bgadj<0, 0, read_neon_blanked_bgadj)
         , ratio_cn = (read_crim_blanked_zeroed/read_neon_blanked_bgadj_zeroed)
         , FC_cn = log2(ratio_cn)
         ) %>%
  # arrange(experiment, rep, OD_crim, OD_neon) %>% View()
  # select(FC_ratio_real, FC_cn, experiment, version) %>% distinct() %>% View()
  full_join(allConv_wide_2, relationship="many-to-many") %>%
  mutate(crim_estAbs = (read_crim_blanked_zeroed-od_conv_intercept_crim)/od_conv_slope_crim,
         neon_estAbs = (read_neon_blanked_bgadj_zeroed-od_conv_intercept_neon)/od_conv_slope_neon  ) %>%
  mutate(
    FC_estAbs = log2(crim_estAbs/neon_estAbs) ) %>%
  # filter(!(prop_crim==0 & prop_neon==0)) %>%
  mutate(OD_both = OD_crim + OD_neon) %>%
  mutate(OD_both = ifelse(is.na(OD_both), 0.1, OD_both))  

# dat_mix_adj %>%
#   filter(experiment=="test3", prop_crim==0.5, prop_neon==0.5, stcurve=="overnight_275", rep==1) %>% 
#   View()
### Check fit of standard curves ###
gg_mixed_cultures_w_stcurves_mixtest2 <- dat_mix_adj %>% 
  filter(experiment =="test2", !(prop_neon==0 & prop_crim==0)) %>%
  ggplot() +
  geom_point(aes(x=prop_crim, y=crim_estAbs, col="crim")) +
  geom_point(aes(x=1-prop_neon, y=neon_estAbs, col="neon")) +
  facet_wrap(version~stcurve, scales="free", nrow=4
             , labeller = labeller(stcurve = c(exponential_200="St. curve made with\n Exponential, 200uL, intercept=0", exponential_275="St. curve made with\n Exponential, 275uL, intercept=0"
                                                 ,overnight_200="St. curve made with\n Overnight, 200uL, intercept=0", overnight_275="St. curve made with\n Overnight, 275uL, intercept=0"
                                                 ,exponential_200_winter="St. curve made with\n Exponential, 200uL, variable intercept", exponential_275_winter="St. curve made with\n Exponential, 275uL, variable intercept"
                                                 ,overnight_200_winter="St. curve made with\n Overnight, 200uL, variable intercept", overnight_275_winter="St. curve made with\n Overnight, 275uL, variable intercept")
                                   , version = c(Nobgadj = "Did not subtract\nbackground fluorescence", winter = "Subtracted background fluorescence"))) +
  xlab("Proportion Crimson cells") +ylab("Estimated Absorbance 600")+
  scale_color_manual(values=c(crim="darkorange", neon="green"))+
  labs(col="Strain", title="Mixed cultures test 2\n200uL, total OD = 0.3 or Abs600 ~ 0.05")
gg_mixed_cultures_w_stcurves_mixtest2
ggsave("01_standard_curve/gg_mixed_cultures_w_stcurves_mixtest2.png", gg_mixed_cultures_w_stcurves_mixtest2, height=10, width=14)

gg_mixed_cultures_w_stcurves_mixtest3 <- dat_mix_adj %>% 
  filter(experiment =="test3", !(prop_neon==0 & prop_crim==0)) %>%
  ggplot() +
  geom_point(aes(x=prop_crim, y=crim_estAbs, col="crim")) +
  geom_point(aes(x=1-prop_neon, y=neon_estAbs, col="neon")) +
  facet_wrap(version~stcurve, scales="free", nrow=4
             , labeller = labeller(stcurve = c(exponential_200="St. curve made with\n Exponential, 200uL, intercept=0", exponential_275="St. curve made with\n Exponential, 275uL, intercept=0"
                                               ,overnight_200="St. curve made with\n Overnight, 200uL, intercept=0", overnight_275="St. curve made with\n Overnight, 275uL, intercept=0"
                                               ,exponential_200_winter="St. curve made with\n Exponential, 200uL, variable intercept", exponential_275_winter="St. curve made with\n Exponential, 275uL, variable intercept"
                                               ,overnight_200_winter="St. curve made with\n Overnight, 200uL, variable intercept", overnight_275_winter="St. curve made with\n Overnight, 275uL, variable intercept")
                                   , version = c(Nobgadj = "Did not subtract\nbackground fluorescence", winter = "Subtracted background fluorescence"))) +
  xlab("Proportion Crimson cells") +ylab("Estimated Absorbance 600")+
  scale_color_manual(values=c(crim="darkorange", neon="green"))+
  labs(col="Strain", title="Mixed cultures test 3\n275uL, total Abs600 ~ 0.01")
gg_mixed_cultures_w_stcurves_mixtest3
ggsave("01_standard_curve/gg_mixed_cultures_w_stcurves_mixtest3.png", gg_mixed_cultures_w_stcurves_mixtest3, height=10, width=14)




gg_mixed_cultures_simple <- dat_mix_adj %>% 
  filter(experiment =="test3", !(prop_neon==0 & prop_crim==0)) %>%
  filter(stcurve=="exponential_275", version=="winter") %>%
  ggplot() +
  geom_point(aes(x=prop_crim, y=crim_estAbs, col="crim")) +
  geom_point(aes(x=1-prop_neon, y=neon_estAbs, col="neon")) +
  xlab("Proportion Crimson cells") +ylab("Estimated Absorbance 600")+
  scale_color_manual(values=c(crim="darkorange", neon="green"))+
  labs(col="Strain")
gg_mixed_cultures_simple
ggsave("01_standard_curve/gg_mixed_cultures_simple.png", gg_mixed_cultures_simple, height=3, width=5)


dat_mix_adj %>%
  select(read_crim_blanked, read_neon_blanked, read_neon_blanked_bgadj, read_neon_blanked_bgadj_zeroed, version) %>%
  ggplot() + 
  geom_point(aes(x=read_crim_blanked, y=read_neon_blanked_bgadj, col=version)) 
  
gg_ratio_standard <- dat_mix_adj %>%
  filter(!(prop_crim==0 & prop_neon==0)) %>%
  filter(OD_both>0.001) %>% #filter(experiment=="test2") %>%View()
  ggplot(aes(x=FC_ratio_real, y=FC_cn, col=experiment, pch=version)) +
  geom_point() +
  geom_abline(aes(slope=1, intercept=0)) +
  xlab("Fold2-change of OD of\n cell proportion ratios (Crim:Neon)") +
  ylab("Fold2-change of fluorescence (Crim:Neon)")
gg_ratio_standard
ggsave("01_standard_curve/gg_ratio_standard.png", gg_ratio_standard, height=4, width=6)

allConv_slopeconstants 
allConv_slopeconstants_adj <- allConv_slopeconstants %>%
  select(experiment, FC_slope_ratio) %>%
  filter(experiment!="allLog10_winter") %>%
  rename(stcurve = experiment)
gg_ratio_standard_adjsted <- dat_mix_adj %>%
  filter(!(prop_crim==0 & prop_neon==0)) %>%
  filter(OD_both>0.001) %>%
  full_join(allConv_slopeconstants_adj) %>%
  ggplot(aes(x=FC_ratio_real, y=FC_cn - FC_slope_ratio,  col=experiment)) +
  geom_point() +
  geom_abline(aes(slope=1, intercept=0))+
  xlab("Fold2-change of \ncell proportion ratios (Crim:Neon)") +
  ylab("Fold2-change of Abs600 (Crim:Neon)\n plus log2 slope ratio constant calculated by standard curve") +
  labs(title = "Removed total OD<0.001; and OD==0 for either strain") +
  facet_grid(version~stcurve
             , labeller = labeller(version=c(Nobgadj="No background\nneon adjustment", winter="Neon autofluorescence\nadjustment")
                                   ,stcurve=c(exponential_200="Standard curve made with\nExponential culture, 200uL\nIntercept forced through 0"
                                             ,exponential_200_winter="Standard curve made with\nExponential culture, 200uL\nFree intercept"
                                             ,exponential_275="Standard curve made with\nExponential culture, 275uL\nIntercept forced through 0"
                                             ,exponential_275_winter="Standard curve made with\nExponential culture, 275uL\nFree intercept"
                                             ,overnight_200="Standard curve made with\nOvernight culture, 200uL\nIntercept forced through 0"
                                             ,overnight_200_winter="Standard curve made with\nOvernight culture, 200uL\nFree intercept"
                                             ,overnight_275="Standard curve made with\nOvernight culture, 275uL\nIntercept forced through 0"
                                             ,overnight_275_winter="Standard curve made with\nOvernight culture, 275uL\nFree intercept"
                                   )))
gg_ratio_standard_adjsted
ggsave("01_standard_curve/gg_ratio_standard_adjsted.png", gg_ratio_standard_adjsted, height=6, width=16)

# Using a constant is the SAME as converting to "pure cells" and then calculating ratios.
gg_ratio_using_standard_curve <- dat_mix_adj %>%
  filter(!(prop_crim==0 & prop_neon==0)) %>%
  filter(OD_both>0.001) %>%
  ggplot(aes(x=FC_ratio_real, y=FC_estAbs, col=experiment)) +
  geom_point() +
  # geom_smooth(method="lm", se = FALSE)+
  geom_abline(aes(slope=1, intercept=0))+
  xlab("Fold2-change of \ncell proportion ratios (Crim:Neon)") +
  ylab("Fold2-change of estimated Abs600 (Crim:Neon)\n using different standard curves")+
  labs(title = "Removed total OD<0.001; and OD==0 for either strain") +
  facet_grid(version~stcurve
             , labeller = labeller(version=c(Nobgadj="No background\nneon adjustment", winter="Neon autofluorescence\nadjustment")
                                   ,stcurve=c(exponential_200="Standard curve made with\nExponential culture, 200uL\nIntercept forced through 0"
                                              ,exponential_200_winter="Standard curve made with\nExponential culture, 200uL\nFree intercept"
                                              ,exponential_275="Standard curve made with\nExponential culture, 275uL\nIntercept forced through 0"
                                              ,exponential_275_winter="Standard curve made with\nExponential culture, 275uL\nFree intercept"
                                              ,overnight_200="Standard curve made with\nOvernight culture, 200uL\nIntercept forced through 0"
                                              ,overnight_200_winter="Standard curve made with\nOvernight culture, 200uL\nFree intercept"
                                              ,overnight_275="Standard curve made with\nOvernight culture, 275uL\nIntercept forced through 0"
                                              ,overnight_275_winter="Standard curve made with\nOvernight culture, 275uL\nFree intercept"
                                   )))

gg_ratio_using_standard_curve
ggsave("01_standard_curve/gg_ratio_using_standard_curve.png", gg_ratio_using_standard_curve, height=6, width=16)
##
# In order to calculate "absolute" values first, you need 2 coefficients: slope of neon and slope of crim
# To get slope, you need 2 (probably ideally 3) counts for each; thus you need to count 4-6 wells worth of cells
# Alternatively: you can just calculate the offset factor which is the conversion of 5:5 to 5:5
# You theoretically only need ONE value for this
dat_mix_adj
gg_ratio_stcurve_vs_constant_compare <- dat_mix_adj %>%
  rowwise() %>%
  mutate(minOD = min(c(OD_crim,OD_neon))) %>%
  ungroup() %>%
  full_join(allConv_slopeconstants_adj) %>%
  mutate(FC_cn_constantadj = FC_cn - FC_slope_ratio) %>%
  filter(!(prop_crim==0 & prop_neon==0)) %>%
  filter(OD_both>0.001) %>%
  select(FC_ratio_real, FC_estAbs, FC_cn_constantadj, experiment, version, stcurve,minOD) %>%
  filter(version=="winter") %>% select(-version) %>%
  pivot_longer(c(FC_estAbs, FC_cn_constantadj), names_to="type", values_to="FC_cn") %>%
  ggplot(aes(x=FC_ratio_real, y=FC_cn, pch=experiment, col=(minOD))) +
  geom_point(position = position_dodge()) +
  # geom_smooth(method="lm", se = FALSE)+
  geom_abline(aes(slope=1, intercept=0))+
  xlab("Fold2-change of \ncell proportion ratios (Crim:Neon)") +
  ylab("Fold2-change of estimated Abs600 (Crim:Neon)")+
  labs(title = "Removed total OD<0.001; and OD==0 for either strain", col="Minimum OD\nin mix", shape="Experiment") +
  scale_shape_manual(values=c(15,17,19)) +
  facet_grid(type~stcurve
             ,labeller = labeller(type=c(FC_cn_constantadj="Log2-ratios of fluorescence plus\nlog2-slope-ratio-constant",
                                         FC_estAbs ="Log2-ratios of estimated absorbance\n(using st. curves)")
                                  , stcurve=c(exponential_200="Standard curve made with\nExponential culture, 200uL\nIntercept forced through 0"
                                              ,exponential_200_winter="Standard curve made with\nExponential culture, 200uL\nFree intercept"
                                              ,exponential_275="Standard curve made with\nExponential culture, 275uL\nIntercept forced through 0"
                                              ,exponential_275_winter="Standard curve made with\nExponential culture, 275uL\nFree intercept"
                                              ,overnight_200="Standard curve made with\nOvernight culture, 200uL\nIntercept forced through 0"
                                              ,overnight_200_winter="Standard curve made with\nOvernight culture, 200uL\nFree intercept"
                                              ,overnight_275="Standard curve made with\nOvernight culture, 275uL\nIntercept forced through 0"
                                              ,overnight_275_winter="Standard curve made with\nOvernight culture, 275uL\nFree intercept"
                                  )
                                  ))
  

gg_ratio_stcurve_vs_constant_compare
ggsave("01_standard_curve/gg_ratio_stcurve_vs_constant_compare.png", gg_ratio_stcurve_vs_constant_compare, height=7, width=16)


#### For simplified visual explanation ######
dat_mix_adj$stcurve
gg_simple_curve_fit <- dat_mix_adj %>%
  full_join(allConv_slopeconstants_adj) %>%
  mutate(FC_cn_constantadj = FC_cn - FC_slope_ratio) %>%
  filter(!(prop_crim==0 & prop_neon==0), experiment=="test3") %>%
  filter(OD_both>0.001) %>%
  filter(version=="winter", stcurve == "exponential_275") %>% select(-version, -stcurve) %>%
  select(FC_ratio_real, FC_cn_constantadj, FC_cn, experiment) %>%
  pivot_longer(c(FC_cn_constantadj, FC_cn), names_to="type", values_to="FC_fluor") %>%
  ggplot(aes(x=FC_ratio_real, y=FC_fluor)) +
  geom_point() +
  # geom_smooth(method="lm", se = FALSE)+
  geom_abline(aes(slope=1, intercept=0))+
  xlab("Fold2-change of \nculture proportion ratios\n(log2(CrimOD:NeonOD))") +
  ylab("Fold2-change of fluorescence ratios\n(log2(CrimFluor/NeonFluor))")+
  labs(title = "Removed total OD<0.001; when OD==0 for either strain") +
  facet_grid(.~type, labeller=labeller(type=c(FC_cn="Original plot\n(blanked fluoresence values)", FC_cn_constantadj="Y-transformed using\nlog2(m1/m2) constant")))+
  ylim(-5,5) +xlim(-5,5)
gg_simple_curve_fit
ggsave("01_standard_curve/gg_simple_curve_fit.png", gg_simple_curve_fit, height=5, width=8)

gg_simple_log2fc_example <- dat_mix_adj %>%
  full_join(allConv_slopeconstants_adj) %>%
  mutate(FC_cn_constantadj = FC_cn - FC_slope_ratio) %>%
  filter(!(prop_crim==0 & prop_neon==0), experiment=="test2") %>%
  # filter(OD_both>0.001) %>%
  filter(version=="winter", stcurve == "exponential_275") %>% select(-version, -stcurve) %>%
  select(FC_ratio_real, FC_cn, experiment) %>%
  # pivot_longer(c(FC_cn_constantadj, FC_cn), names_to="type", values_to="FC_fluor") %>%
  ggplot(aes(y=FC_ratio_real, x=FC_cn)) +
  geom_point() +
  ylab("Fold-change of \nculture proportion ratios\nlog2(CrimOD:NeonOD)") +
  xlab("Fold-change of fluorescence ratios\nlog2(CrimFluor/NeonFluor)")+
  # labs(title = "Removed total OD<0.001; when OD==0 for either strain") +
  # facet_grid(.~type, labeller=labeller(type=c(FC_cn="Original plot\n(blanked fluoresence values)", FC_cn_constantadj="Y-transformed using\nlog2(m1/m2) constant")))+
  ylim(-4,4) +xlim(-4,4)
gg_simple_log2fc_example
ggsave("01_standard_curve/gg_simple_log2fc_example.png", gg_simple_log2fc_example, height=3, width=4)

