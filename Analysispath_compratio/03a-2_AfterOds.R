library(tidyverse)
library(brms)
load("03_EDA_and_adj/allDat_final.RData")
load("03_EDA_and_adj/plant_colours.RData")
dir.create("03a-2_AfterOds")

od <- read.csv("00_raw_data/dat/abs_data/ODs_beforeafter_longformat.csv")
abs <- read.csv("00_raw_data/dat/abs_data/abs_to_od_data.csv")

allDat <- od %>% 
  mutate(col=as.factor(col)) %>%
  left_join(allDat_final) %>%
  filter(od_after<0.4) %>%
  mutate(log10_totalOD = log10(totalOD)) %>%
  mutate(log10_totalOD = ifelse(!is.finite(log10_totalOD),0,log10_totalOD)) %>%
  filter(!is.na(StrainMix))

allDat %>%
  filter(plate=="N2C3_WCS365_5day001_b") %>%
  select(path_cells_log, protect_cells_log) %>% table()

gg_beforeafterod <- allDat %>% 
  filter(plate=="N2C3_WCS365_5day001_b") %>%
  filter(col!=6) %>% # this col is contaminated I think
  ggplot(aes(x=total_cells_log, y=od_after)) +
  # geom_boxplot() +
  # geom_jitter(height=0) +
  geom_point() +
  # geom_smooth(method="lm")+
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x") +
  scale_x_continuous(breaks=c(0,3,4,5,6,7,8)) +
  xlab("Initial number of cells at inoculation") +
  ylab("Absorbance (600nm) of well\nafter 1 week growth with plant")
gg_beforeafterod
ggsave(filename = "03a-2_AfterOds/gg_beforeafterod.png", gg_beforeafterod, height=4, width=6)


dat_for_brms_cellload <-  allDat %>% 
  filter(plate=="N2C3_WCS365_5day001_b") %>%
  filter(col!=6) %>%
  filter(StrainMix %in% c("MOCK-N2C3","WCS365-MOCK"))
brm_cellload <- brm(od_after ~ total_cells_log*StrainMix, data=dat_for_brms_cellload
                    , seed=793
                    , iter=4000
                    , file="03a-2_AfterOds/brm_celload")
brm_cellload
fixef(brm_cellload)

dat_for_brms_cellload %>%
  select(experiment) %>% table()
# Make a table for this
gg_platemap_concentration <- dat_for_brms_cellload %>%
  mutate(row = factor(row, levels=rev(LETTERS[1:8]))) %>%
  mutate(Strain = gsub("-MOCK","",gsub("MOCK-","",StrainMix))) %>%
  ggplot()+
  geom_tile(aes(x=col, y=row, fill=total_cells_log, col=Strain), lwd=2) +
  scale_color_manual(values=c("orange","darkgreen")) +
  scale_fill_gradient(low="white", high="black")+
  labs(fill="Cell abundance\nat inoculation") +
  ylab("Row") + xlab("Column")
gg_platemap_concentration
ggsave("03a-2_AfterOds/gg_platemap_concentration.png", gg_platemap_concentration
       , height=5, width=10)

dat_for_brms_cellload %>%
  ggplot() +
  geom_point(aes(x=StrainMix, y=Healthiness_hsv, col=total_cells_log))
  

cells_per_Abs1 = (3.166*1-0.268)*5E8
lm_abs_to_od <- lm(OD600~ABS600, data=abs)
  ## Abs to OD

gg_abs_vs_od <- abs %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(x=ABS600, y=OD600, col=Isolate)) +
  geom_abline(aes(slope=lm_abs_to_od$coefficients[2], intercept = lm_abs_to_od$coefficients[1])) +
  geom_text(aes(x=0.2, y=0.75, label=paste0("y = ",signif(lm_abs_to_od$coefficients[2],3),"*x + ",signif(lm_abs_to_od$coefficients[1],3))))
gg_abs_vs_od

gg_abs_vs_od
ggsave("03a-2_AfterOds/gg_abs_vs_od.png", gg_abs_vs_od, height=4, width=5)

abs_filt <- abs %>% filter(ABS600>0.125)
lm_abs_to_od_filt <- lm(OD600~ABS600, data=abs_filt)
gg_abs_vs_od_filt <- abs %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(x=ABS600, y=OD600, col=Isolate)) +
  geom_abline(aes(slope=lm_abs_to_od_filt$coefficients[2], intercept = lm_abs_to_od_filt$coefficients[1])) +
  geom_text(aes(x=0.2, y=0.75, label=paste0("y = ",signif(lm_abs_to_od_filt$coefficients[2],3),"*x + ",signif(lm_abs_to_od_filt$coefficients[1],3)))) +
  labs(title="ABS600 > 0.125")
gg_abs_vs_od_filt
ggsave("03a-2_AfterOds/gg_abs_vs_od_filt.png", gg_abs_vs_od_filt, height=4, width=5)



