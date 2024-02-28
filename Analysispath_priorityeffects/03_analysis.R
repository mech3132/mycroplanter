#!bin/bash
library(brms)
library(cowplot)
library(tidyverse)
dir.create("03_analysis")

##### Load data #######
load("01_compile_data/dat_plant.RData")

#### Make platemaps ####

### Platecheck 
plateconv <-  dat_plant %>% select(plate, experiment) %>% distinct() %>%
  group_by(experiment) %>% mutate(platen=paste0("Plate",rank(plate, ties.method = "min")))  %>%
  ungroup() 
gg_platemap_priority <- dat_plant %>%
  left_join(plateconv) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2, pch=factor(transfer_h))) +
  facet_grid(experiment~platen)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", N2C3="orange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkgreen", N2C3="orange"), na.value = "black") +
  scale_shape_manual(values=c(`0`=21, `0.01`=19, `3`=17, `6`=15, `24`=4, `48`=8))+
  labs(shape="Time delay (h)", fill="Strain1", col="Strain2")+
  ylab("row") + xlab("column")
gg_platemap_priority
ggsave("03_analysis/gg_platemap_priority.png", gg_platemap_priority, width=16, height=6)

#### Make figure ####

### Look at raw neon and crimson for N2C3 alone only to see if fluorescence is noticably dimming over time
gg_effectdelay_onfluorloss <- dat_plant %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  select(experiment, StrainMix, Alive, strain1, strain2, fluor_strain1, fluor_strain2, crim_raw_blanked, neon_raw_blanked, transfer_h) %>%
  filter(!StrainMix %in% c("N2C3-MOCK","WCS365-MOCK", "MOCK-MOCK")) %>%
  filter(!is.na(fluor_strain2)) %>%
  filter(transfer_h>0) %>%
  mutate(fluor_val = ifelse(fluor_strain2=="crim", crim_raw_blanked, 
                            ifelse(fluor_strain2=="neon", neon_raw_blanked, NA))) %>%
  # unite(fluor_strain1, fluor_strain2, col="FluorMix", sep="-", remove=FALSE) %>%
  ggplot() +
  geom_jitter(aes(x=factor(transfer_h), y=log10(fluor_val+1), col=strain1), width=0.1, height=0) +
  facet_grid(fluor_strain2 ~ strain2,
             labeller=labeller(strain2=c(WCS365="WCS365 as second strain", N2C3="N2C3 as second strain")
                               , fluor_strain2 = c(neon="Second strain\nwas mNeon", crim="Second strain\nwas mCrimson"))) +
  ylab("Raw fluorescence (blanked)\n of SECOND strain at end of experiment") +
  xlab("Delay in inoculation of second strain (h)")+
  labs(col="First strain") +
  scale_color_manual(values=c("MOCK"="black", WCS365="darkgreen", N2C3="darkorange"))
gg_effectdelay_onfluorloss
ggsave("03_analysis/gg_effectdelay_onfluorloss.png", gg_effectdelay_onfluorloss, height=4, width=8)
  

## only controls (mono)
gg_priority_controls <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400, domStrain = ifelse(WN_FC_fluor_offset>0, "WCS365 more\nabundant", "N2C3 more\nabundant")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  ggplot(aes(x=StrainMix, y=domStrain)) +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0.2) +
  facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dipped in first", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "darkorange")) +
  geom_hline(aes(yintercept=0)) +
  ylab("Which strain dominates\naccording to fluorescence")+xlab("First Strain - Second Strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(col="Plant healthy\n(>400)", shape="Experiment") 
gg_priority_controls
ggsave("03_analysis/gg_priority_controls.png", gg_priority_controls, height=4, width=8)

gg_priority_controls_v2 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400, domStrain = ifelse(WN_FC_fluor_offset>0, "WCS365 more\nabundant", "N2C3 more\nabundant")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first\nonly", paste0(transfer_h," hours"))) %>%
  mutate(Delay = factor(Delay, levels=c("Dipped in first\nonly", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=domStrain)) +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0.2) +
  facet_grid(.~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-MOCK"="WCS365, then mock", "MOCK-WCS365"="Mock, then WCS365", "N2C3-MOCK"="N2C3, then mock", "MOCK-N2C3"="Mock, then N2C3"))) +
  scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "darkorange")) +
  geom_hline(aes(yintercept=0)) +
  ylab("Which strain dominates\naccording to fluorescence")+xlab("Time between first and second strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(col="Plant healthy\n(>400)", shape="Experiment") 
gg_priority_controls_v2
ggsave("03_analysis/gg_priority_controls_v2.png", gg_priority_controls_v2, height=4, width=8)



gg_priority_controls_v3 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(fluor_strain1=ifelse(strain1=="MOCK", NA, fluor_strain1), fluor_strain2=ifelse(strain2=="MOCK", NA, fluor_strain2)) %>%
  mutate(Alive = Healthiness_hsv>400, domStrain = ifelse(WN_FC_fluor_offset>0, "WCS365 more\nabundant", "N2C3 more\nabundant")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first\nonly", paste0(transfer_h," hours"))) %>%
  mutate(Delay = factor(Delay, levels=c("Dipped in first\nonly", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  mutate(fluor_of_strain = ifelse(!is.na(fluor_strain1), fluor_strain1, 
                                  ifelse(!is.na(fluor_strain2), fluor_strain2, 
                                         NA))) %>%
  mutate(fluor_of_oppositestrain = ifelse(fluor_of_strain=="crim", "neon", 
                                  ifelse(fluor_of_strain=="neon", "crim", 
                                         NA))) %>%
  mutate(fluor_val_of_strain = ifelse(fluor_of_strain=="crim", crim_raw_blanked, 
                                      ifelse(fluor_of_strain=="neon", neon_raw_blanked,
                                             NA))) %>%
  mutate(fluor_val_of_oppositestrain = ifelse(fluor_of_strain=="crim", neon_raw_blanked, 
                                      ifelse(fluor_of_strain=="neon", crim_raw_blanked,
                                             NA))) %>%
  ggplot(aes(x=Delay)) +
  geom_col(aes(y=(fluor_val_of_strain), fill=fluor_of_strain), position=position_dodge2()) +
  geom_col(aes(y=-(fluor_val_of_oppositestrain), fill=fluor_of_oppositestrain), position=position_dodge2()) +
  facet_grid(experiment~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-MOCK"="WCS365, then mock", "MOCK-WCS365"="Mock, then WCS365", "N2C3-MOCK"="N2C3, then mock", "MOCK-N2C3"="Mock, then N2C3"))) +
  scale_fill_manual(values=c(neon="green", crim = "darkred")) +
  geom_hline(aes(yintercept=0)) +
  ylab("Positive bars = Raw fluorescence (blanked) of strain\nNegative bars = Raw fluorescence (blanked) of other ex/em spectra")+xlab("Time between first and second strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(fill="Fluorescence\nof strain") 
gg_priority_controls_v3
ggsave("03_analysis/gg_priority_controls_v3.png", gg_priority_controls_v3, height=6, width=10)

gg_priority_mixedonly <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dipped in first", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "darkorange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  # scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First Strain - Second Strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(col="Plant healthy\n(>400)", shape="Experiment") +
  ylim(-10,10)
gg_priority_mixedonly
ggsave("03_analysis/gg_priority_mixedonly.png", gg_priority_mixedonly, height=6, width=10)

gg_priority_mixedonlyV2 <- dat_plant %>%
  mutate(experiment = ifelse(experiment == "2023-12_priorityeffects", "2023-11-08", 
                             ifelse(experiment == "2023-12_priorityeffects2", "2023-12-06", NA))) %>%
  filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first\nonly", paste0(transfer_h," hours"))) %>%
  mutate(Delay = factor(Delay, levels=c("Dipped in first\nonly", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  facet_grid(.~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-N2C3"="WCS365, then N2C3", "N2C3-WCS365"="N2C3, then WCS365"))) +
  scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "darkorange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  # scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Delay between first and second strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(col="Plant healthy\n(>400)", shape="Experiment") +
  ylim(-10,10)
gg_priority_mixedonlyV2
ggsave("03_analysis/gg_priority_mixedonlyV2.png", gg_priority_mixedonlyV2, height=4, width=6)

gg_priority_mixedonlyV3 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first\nonly", paste0(transfer_h," hours"))) %>%
  mutate(Delay = factor(Delay, levels=c("Dipped in first\nonly", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  facet_grid(experiment~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-N2C3"="WCS365, then N2C3", "N2C3-WCS365"="N2C3, then WCS365"))) +
  scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "darkorange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  # scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Delay between first and second strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(col="Plant healthy\n(>400)", shape="Experiment") +
  ylim(-10,10)
gg_priority_mixedonlyV3
ggsave("03_analysis/gg_priority_mixedonlyV3.png", gg_priority_mixedonlyV3, height=6, width=8)


dat_plant %>%
  filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
  geom_jitter(aes(col=factor(transfer_h)), width=0, height=0.1) +
  geom_smooth(aes(col=factor(transfer_h)),method="glm", method.args=c(family="binomial"), se=FALSE) +
  facet_grid(strain1~.) +
  scale_color_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))


##### BRMS ####
dat_plant %>%
  filter(StrainMix == "MOCK-MOCK")
dat_plant_forbrms <- dat_plant %>%
  filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("WCS365-N2C3", "N2C3-WCS365"))) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  # mutate(FirstStrain = ifelse(transfer_h==0, "Both", as.character(strain1))) %>%
  select(Alive, StrainMix, WN_FC_fluor_offset, transfer_h, experiment) %>%
  filter(transfer_h!=0)  %>%
  mutate(WN_FC_fluor_offset = ifelse(StrainMix=="MOCK-MOCK",0,WN_FC_fluor_offset)) %>%
  drop_na()
dat_plant_forbrms %>% select(transfer_h, StrainMix) %>% table()
# dat_plant_forbrms_exp1 <-dat_plant_forbrms %>%
#   filter(experiment=="2023-12_priorityeffects")
# dat_plant_forbrms_exp2 <-dat_plant_forbrms %>%
#   filter(experiment=="2023-12_priorityeffects2")

brms_mixedonly <- brm(Alive ~ WN_FC_fluor_offset + StrainMix + StrainMix:transfer_h, data=dat_plant_forbrms
    , seed=52093
    , iter=4000
    , family="bernoulli"
    , file="03_analysis/brms_mixedonly"
    )
brms_mixedonly

brms_mixedonly_nostrainmix <- brm(Alive ~ WN_FC_fluor_offset, data=dat_plant_forbrms
                           , seed=52093
                           , family="bernoulli"
                           , file="03_analysis/brms_mixedonly_nostrainmix"
)
brms_mixedonly_nostrainmix

brms_mixedonly_nofluor <- brm(Alive ~ StrainMix + StrainMix:transfer_h, data=dat_plant_forbrms
                           , seed=52093
                           , family="bernoulli"
                           , file="03_analysis/brms_mixedonly_nofluor"
)
brms_mixedonly_nofluor
# 
# brms_mixedonly_exp2 <- brm(Alive ~ WN_FC_fluor_offset + StrainMix + transfer_h, data=dat_plant_forbrms_exp2
#                            , seed=52093
#                            , family="bernoulli"
#                            , file="03_analysis/brms_mixedonly_exp2"
# )
# brms_mixedonly_exp2
# brms_mixedonly_exp2_nostrainmix <- brm(Alive ~ WN_FC_fluor_offset + transfer_h, data=dat_plant_forbrms_exp2
#                            , seed=52093
#                            , file="03_analysis/brms_mixedonly_exp2_nostrainmix"
# )
# brms_mixedonly_exp2_nostrainmix
# brms_mixedonly_exp2_nofluor <- brm(Alive ~ StrainMix + transfer_h, data=dat_plant_forbrms_exp2
#                                        , seed=52093
#                                    , family="bernoulli"
#                                    , file="03_analysis/brms_mixedonly_exp2_nofluor"
# )
# brms_mixedonly_exp2_nofluor


LOO_compare <- LOO(brms_mixedonly, brms_mixedonly_nostrainmix,brms_mixedonly_nofluor)
LOO_compare$diffs

LOO_compare2 <- LOO(brms_mixedonly_nostrainmix,brms_mixedonly_nofluor)
LOO_compare2$diffs

# LOO_compare_exp2 <- LOO(brms_mixedonly_exp2, brms_mixedonly_exp2_nostrainmix,brms_mixedonly_exp2_nofluor)
# LOO_compare_exp2$diffs

sink("03_analysis/model_loo_compare.txt")
print("both experimetnts 1")
LOO_compare$diffs

print("only single models")
LOO_compare2$diffs

print("both model") 
brms_mixedonly
print("fluor only")
brms_mixedonly_nostrainmix
print("Strainmix only")
brms_mixedonly_nofluor

sink()

#### Plot with brms ####

dat_mixedonly_exp2_withpredictions <- posterior_linpred(brms_mixedonly, transform=TRUE) %>%
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
  bind_cols(brms_mixedonly$data) %>%
  left_join(dat_plant)

gg_brmfit_exp2_all <- dat_mixedonly_exp2_withpredictions %>%
  filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
  geom_jitter(aes(col=factor(transfer_h)), width=0, height=0.1) +
  geom_ribbon(aes(fill=factor(transfer_h), ymin=Q2.5, ymax=Q97.5), alpha=0.1, show.legend = FALSE)+
  geom_line(aes(col=factor(transfer_h), y=Q50))+
  # geom_smooth(aes(col=factor(transfer_h)),method="glm", method.args=c(family="binomial"), se=FALSE) +
  facet_grid(strain1~.,
             labeller=labeller(strain1=c(WCS365="WCS365 as first strain", N2C3="N2C3 as first strain"))
             ) +
  scale_color_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))+
  scale_fill_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))+
  ylab("Probability of healthy plant") +
  xlab("WCS365 : N2C3 (log2-fold difference)")+ 
  labs(col="Delay in second strain (h)", fill="Delay in\nsecond strain (h)")
gg_brmfit_exp2_all
ggsave("03_analysis/gg_brmfit_exp2_all.png", gg_brmfit_exp2_all, height=4, width=6)


### second experiment


dat_plant %>%
  ggplot(aes(x=strain1, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv, pch=factor(transfer_h)), width=0.1, height=0) +
  facet_grid(.~strain2) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") 

dat_plant %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  filter(experiment=="2023-12_priorityeffects2") %>%
  ggplot(aes(x=strain1, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv), width=0.1, height=0) +
  facet_grid(transfer_h~strain2) +
  # scale_shape_manual(values=c(`TRUE` = 19, `FALSE` =21))+
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") 


dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = Healthiness_hsv>400) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  filter(experiment=="2023-12_priorityeffects2") %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv), width=0.1, height=0) +
  facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x") +
  # scale_shape_manual(values=c(`TRUE` = 19, `FALSE` =21))+
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))



gg_ratio_planthealth_exp2 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = Healthiness_hsv>400) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  # filter(experiment=="2023-12_priorityeffects2") %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(15,19))+
  scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  # scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_ratio_planthealth_exp2
ggsave("01_compile_data/gg_ratio_planthealth_exp2.png", gg_ratio_planthealth_exp2, height=5, width=8)


gg_ratio_planthealth_exp3 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = Healthiness_hsv>400) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  # filter(experiment=="2023-12_priorityeffects2") %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv, pch=experiment), width=0.1, height=0) +
  facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(15,19))+
  # scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_ratio_planthealth_exp3
ggsave("01_compile_data/gg_ratio_planthealth_exp3.png", gg_ratio_planthealth_exp3, height=5, width=8)

gg_ratio_planthealth_exp4 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = Healthiness_hsv>400) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  # filter(experiment=="2023-12_priorityeffects2") %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv, pch=experiment), width=0.1, height=0) +
  facet_grid(experiment~transfer_h, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(15,19))+
  # scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_ratio_planthealth_exp4
ggsave("01_compile_data/gg_ratio_planthealth_exp4.png", gg_ratio_planthealth_exp4, height=6, width=10)


gg_ratio_planthealth_exp5 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = Healthiness_hsv>400) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  # filter(experiment=="2023-12_priorityeffects2") %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
  ggplot(aes(x=factor(transfer_h), y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv, pch=experiment), width=0.1, height=0) +
  facet_grid(experiment~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(15,19))+
  # scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
  # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_ratio_planthealth_exp5
ggsave("01_compile_data/gg_ratio_planthealth_exp5.png", gg_ratio_planthealth_exp4, height=6, width=10)



dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  filter(experiment=="2023-12_priorityeffects2") %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365"), transfer_h>0) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
  geom_jitter(aes(col=StrainMix), width=0, height=0.2)+
  geom_smooth(method="glm",
              method.args=list(family="binomial"), col="black")+
  # facet_grid(transfer_h~StrainMix)
  facet_grid(.~transfer_h)

dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  filter(experiment=="2023-12_priorityeffects2") %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Alive, col=StrainMix)) +
  geom_point()+
  geom_smooth(method="glm",
              method.args=list(family="binomial"))+
  facet_grid(.~transfer_h)


## Look at ones were strains are the same
dat_plant %>%
  mutate(strain1 = as.character(strain1), strain2 = as.character(strain2)) %>%
  filter(strain1==strain2) %>% 
  select(strain1, strain2, crim_top, neon_top) %>%
  pivot_longer(-c(strain1, strain2), names_to="fluortype", values_to="read") %>%
  ggplot() +
  geom_point(aes(x=fluortype, y=log10(read))) +
  facet_grid(strain1 ~ strain2)

dat_plant %>%
  filter(strain2 !="No second strain") %>%
  ggplot() +
  geom_point(aes(x=WN_FC_fluor, y=Healthiness_hsv, col=strain1)) +
  facet_grid(strain2~transfer_h)


