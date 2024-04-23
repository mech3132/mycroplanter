#!bin/bash
library(brms)
library(cowplot)
library(tidyverse)
dir.create("03_analysis")

##### Load data #######
load("02_compile_data/dat_plant.RData")
load("02_compile_data/dat_plant_dd.RData")
load("02_compile_data/allCFU.RData")

# Fix experiment
dat_plant_fullset <- dat_plant
dat_plant <- dat_plant %>%
  filter(experiment !="2024-02-14", experiment !="2024-03-20") 
  
# dat_plant %>%
  # filter(Type=="MIX") %>% View()
#### Filter priority into different platemaps ####

dat_plant_peffect <- dat_plant 
# %>%
  # filter(strain2_od%in%c(0.001,0))
# dat_plant_peffect %>% select(strain2_od, strain1_od) %>% table()

## Figure out concentration stuff later...
# dat_plant_conc <- dat_plant %>%
#   filter(strain2_od !=0.001, strain2_od>0)
# dat_plant_conc %>%
#   select(strain2_od, transfer_h) %>% table()

#### Make platemaps ####

### Platecheck 
plateconv <-  full_join(dat_plant_fullset, dat_plant_dd) %>% select(plate, experiment) %>% distinct() %>%
  group_by(experiment) %>% mutate(platen=paste0("Plate",rank(plate, ties.method = "min")))  %>%
  ungroup() 
gg_platemap_priority <- full_join(dat_plant_fullset, dat_plant_dd) %>%
  left_join(plateconv) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2, pch=factor(transfer_h))) +
  facet_grid(experiment~platen)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkseagreen4", N2C3="darkorange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkseagreen4", N2C3="darkorange"), na.value = "black") +
  scale_shape_manual(values=c(`0`=21, `0.01`=19, `3`=17, `6`=15, `24`=4, `48`=8))+
  labs(shape="Time delay (h)", fill="Strain1", col="Strain2")+
  # theme_bw() +
  ylab("row") + xlab("column")
gg_platemap_priority
ggsave("03_analysis/gg_platemap_priority.png", gg_platemap_priority, width=16, height=6)


plateconv <-  full_join(dat_plant_fullset, dat_plant_dd) %>% select(plate, experiment) %>% distinct() %>%
  group_by(experiment) %>% mutate(platen=paste0("Plate",rank(plate, ties.method = "min")))  %>%
  ungroup() 
gg_platemap_priority <- full_join(dat_plant_fullset, dat_plant_dd) %>%
  left_join(plateconv) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2, pch=factor(transfer_h))) +
  facet_grid(experiment~platen)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkseagreen4", N2C3="darkorange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkseagreen4", N2C3="darkorange"), na.value = "black") +
  scale_shape_manual(values=c(`0`=21, `0.01`=19, `3`=17, `6`=15, `24`=4, `48`=8))+
  labs(shape="Time delay (h)", fill="Strain1", col="Strain2")+
  # theme_bw() +
  ylab("row") + xlab("column")
gg_platemap_priority
ggsave("03_analysis/gg_platemap_priority.png", gg_platemap_priority, width=16, height=6)

## Platemap for sampling?
allCFU %>%
  filter(experiment %in% c("2024_03_Peffect_Trial5", "2024-02_Peffect_Trial4")) %>%
  filter(plant=="col0") %>%
  select(experiment, transfer_h, strain) %>% 
  table()
#### Make figure ####

### Look at raw neon and crimson for N2C3 alone only to see if fluorescence is noticably dimming over time
gg_effectdelay_onfluorloss <- dat_plant_peffect %>%
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
  theme_bw() +
  scale_color_manual(values=c("MOCK"="black", WCS365="darkseagreen4", N2C3="darkorange"))
gg_effectdelay_onfluorloss
ggsave("03_analysis/gg_effectdelay_onfluorloss.png", gg_effectdelay_onfluorloss, height=4, width=8)
  

## only controls (mono)
dat_plant_peffect %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  select(StrainMix) %>% distinct()
gg_priority_controls <- dat_plant_peffect %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400, domStrain = ifelse(WN_FC_fluor_offset>0, "WCS365 more\nabundant", "N2C3 more\nabundant")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first", paste0(transfer_h," hours"))) %>%
  mutate(Delay = factor(Delay, levels=c("Dipped in first", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=domStrain)) +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0.2) +
  facet_grid(.~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-MOCK"="WCS365, then mock", "MOCK-WCS365"="Mock, then WCS365", "N2C3-MOCK"="N2C3, then mock", "MOCK-N2C3"="Mock, then N2C3"))) +
  # scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3")) +
  geom_hline(aes(yintercept=0)) +
  ylab("Which strain dominates\naccording to fluorescence")+xlab("Time between first and second strain") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(col="Plant healthy\n(>400)", shape="Experiment") 
gg_priority_controls
ggsave("03_analysis/gg_priority_controls.png", gg_priority_controls, height=4, width=8)

# 
# 
# gg_priority_controls_v3 <- dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset)) %>%
#   filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
#   mutate(fluor_strain1=ifelse(strain1=="MOCK", NA, fluor_strain1), fluor_strain2=ifelse(strain2=="MOCK", NA, fluor_strain2)) %>%
#   mutate(Alive = Healthiness_hsv>400, domStrain = ifelse(WN_FC_fluor_offset>0, "WCS365 more\nabundant", "N2C3 more\nabundant")) %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first\nonly", paste0(transfer_h," hours"))) %>%
#   mutate(Delay = factor(Delay, levels=c("Dipped in first\nonly", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
#   mutate(fluor_of_strain = ifelse(!is.na(fluor_strain1), fluor_strain1, 
#                                   ifelse(!is.na(fluor_strain2), fluor_strain2, 
#                                          NA))) %>%
#   mutate(fluor_of_oppositestrain = ifelse(fluor_of_strain=="crim", "neon", 
#                                   ifelse(fluor_of_strain=="neon", "crim", 
#                                          NA))) %>%
#   mutate(fluor_val_of_strain = ifelse(fluor_of_strain=="crim", crim_raw_blanked, 
#                                       ifelse(fluor_of_strain=="neon", neon_raw_blanked,
#                                              NA))) %>%
#   mutate(fluor_val_of_oppositestrain = ifelse(fluor_of_strain=="crim", neon_raw_blanked, 
#                                       ifelse(fluor_of_strain=="neon", crim_raw_blanked,
#                                              NA))) %>%
#   ggplot(aes(x=Delay)) +
#   geom_col(aes(y=(fluor_val_of_strain), fill=fluor_of_strain), position=position_dodge2()) +
#   geom_col(aes(y=-(fluor_val_of_oppositestrain), fill=fluor_of_oppositestrain), position=position_dodge2()) +
#   facet_grid(experiment~StrainMix, scales = "free", drop=TRUE, space="free_x"
#              , labeller = labeller(StrainMix=c("WCS365-MOCK"="WCS365, then mock", "MOCK-WCS365"="Mock, then WCS365", "N2C3-MOCK"="N2C3, then mock", "MOCK-N2C3"="Mock, then N2C3"))) +
#   scale_fill_manual(values=c(neon="green", crim = "darkred")) +
#   geom_hline(aes(yintercept=0)) +
#   ylab("Positive bars = Raw fluorescence (blanked) of strain\nNegative bars = Raw fluorescence (blanked) of other ex/em spectra")+xlab("Time between first and second strain") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   labs(fill="Fluorescence\nof strain") 
# gg_priority_controls_v3
# ggsave("03_analysis/gg_priority_controls_v3.png", gg_priority_controls_v3, height=6, width=10)
# 
# gg_priority_mixedonly <- dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
#   filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
#   mutate(Alive = Healthiness_hsv>400) %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
#   facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
#              , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dipped in first", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3")) +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First Strain - Second Strain") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
#   labs(col="Plant healthy\n(>400)", shape="Experiment") +
#   ylim(-10,10)
# gg_priority_mixedonly
# ggsave("03_analysis/gg_priority_mixedonly.png", gg_priority_mixedonly, height=6, width=10)

gg_priority_mixedonly <- dat_plant_peffect  %>%
  filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01, "    Dipped in\nfirst strain", paste0(transfer_h," hours"))) %>%
  mutate(Delay = factor(Delay, levels=c("    Dipped in\nfirst strain", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  # geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  geom_jitter(aes(col=Alive), width=0.1, height=0) +
  facet_grid(.~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-N2C3"="WCS365, then N2C3", "N2C3-WCS365"="N2C3, then WCS365"))) +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3")) +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change fluorescence")+xlab("Delay between first and second strain") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(col="Plant healthy\n(>400)", shape="Experiment") +
  ylim(-10,10)
gg_priority_mixedonly
ggsave("03_analysis/gg_priority_mixedonly.png", gg_priority_mixedonly, height=4, width=6)

#### all dips
gg_priority_mixed_dd <- dat_plant_peffect  %>%
  full_join(dat_plant_dd %>% filter(plant_quality=="good")) %>%
  filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01 & Type=="Priority", "Dipped in first strain;\nput in second strain", 
                        ifelse(transfer_h == 0.01 & Type == "doubledip", "   Dipped in both strains;\nput in MOCK", paste0(transfer_h," hours")))) %>%
  mutate(Delay = factor(Delay, levels=c( "   Dipped in both strains;\nput in MOCK","Dipped in first strain;\nput in second strain", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  # geom_jitter(aes(col=Alive), width=0.1, height=0) +
  facet_grid(.~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-N2C3"="WCS365, then N2C3", "N2C3-WCS365"="N2C3, then WCS365"))) +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3")) +
  scale_shape_discrete(guide="none") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change fluorescence")+xlab("Delay between first and second strain") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(col="Plant healthy\n(>400)", shape="Experiment") +
  ylim(-10,10) 
gg_priority_mixed_dd
ggsave("03_analysis/gg_priority_mixed_dd.png", gg_priority_mixed_dd, height=4, width=7)

### Just the last tow
dat_plant_fullset  %>%
  full_join(dat_plant_dd) %>%
  mutate(strain0 = ifelse(is.na(strain0), "none", strain0)) %>%
  filter(experiment %in% c("2024-02-14","2024-03-20")) %>%
  filter(!is.na(WN_FC_fluor_offset), transfer_h!=0) %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  mutate(Delay = ifelse(transfer_h == 0.01 & Type=="Priority", "Dipped in first strain;\nput in second strain", 
                        ifelse(transfer_h == 0.01 & Type == "doubledip", "   Dipped in both strains;\nput in MOCK", paste0(transfer_h," hours")))) %>%
  mutate(Delay = factor(Delay, levels=c( "   Dipped in both strains;\nput in MOCK","Dipped in first strain;\nput in second strain", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x=Delay, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=strain0), width=0.1, height=0) +
  # geom_jitter(aes(col=Alive), width=0.1, height=0) +
  facet_grid(experiment~StrainMix, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(StrainMix=c("WCS365-N2C3"="WCS365, then N2C3", "N2C3-WCS365"="N2C3, then WCS365"))) +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3")) +
  # scale_shape_discrete(guide="none") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change fluorescence")+xlab("Delay between first and second strain") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(col="Plant healthy\n(>400)", shape="Experiment") +
  ylim(-10,10) 
# 
# dat_plant %>%
#   filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
#   mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
#   ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
#   geom_jitter(aes(col=factor(transfer_h)), width=0, height=0.1) +
#   geom_smooth(aes(col=factor(transfer_h)),method="glm", method.args=c(family="binomial"), se=FALSE) +
#   facet_grid(strain1~.) +
#   scale_color_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))


##### BRMS ####

dat_plant_forbrms <- dat_plant_peffect %>%
  filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("WCS365-N2C3", "N2C3-WCS365"))) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  # mutate(FirstStrain = ifelse(transfer_h==0, "Both", as.character(strain1))) %>%
  select(Alive, StrainMix, WN_FC_fluor_offset, transfer_h, experiment, Type) %>%
  filter(transfer_h!=0, (Type!="MIX")) %>%
  mutate(WN_FC_fluor_offset = ifelse(StrainMix=="MOCK-MOCK",0,WN_FC_fluor_offset)) %>%
  drop_na()
dat_plant_forbrms %>% 
  select(transfer_h, StrainMix) %>% table()
dat_plant_forbrms %>% 
  select(experiment, transfer_h) %>% table()

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
# # 
# ### BRMS (with ALL data version) ####
# dat_plant_fullset_forbrms <- dat_plant_fullset %>%
#   filter(transfer_h!=0) %>%
#   filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("WCS365-N2C3", "N2C3-WCS365"))) %>%
#   mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
#   # mutate(FirstStrain = ifelse(transfer_h==0, "Both", as.character(strain1))) %>%
#   select(Alive, StrainMix, WN_FC_fluor_offset, transfer_h, experiment, Type) %>%
#   # filter(transfer_h!=0, (Type!="MIX")) %>%
#   mutate(WN_FC_fluor_offset = ifelse(StrainMix=="MOCK-MOCK",0,WN_FC_fluor_offset)) %>%
#   drop_na()
# 
# 
# brms_mixedonlyALL <- brm(Alive ~ WN_FC_fluor_offset + StrainMix + StrainMix:transfer_h, data=dat_plant_fullset_forbrms
#                       , seed=52093
#                       , iter=4000
#                       , family="bernoulli"
#                       , file="03_analysis/brms_mixedonlyALL"
# )
# brms_mixedonlyALL
# 
# brms_mixedonly_nostrainmixALL <- brm(Alive ~ WN_FC_fluor_offset, data=dat_plant_fullset_forbrms
#                                   , seed=52093
#                                   , family="bernoulli"
#                                   , file="03_analysis/brms_mixedonly_nostrainmixALL"
# )
# brms_mixedonly_nostrainmixALL
# 
# brms_mixedonly_nofluorALL <- brm(Alive ~ StrainMix + StrainMix:transfer_h, data=dat_plant_fullset_forbrms
#                               , seed=52093
#                               , family="bernoulli"
#                               , file="03_analysis/brms_mixedonly_nofluorALL"
# )
# brms_mixedonly_nofluorALL
# 
# 
# LOO_compare <- LOO(brms_mixedonlyALL, brms_mixedonly_nostrainmixALL,brms_mixedonly_nofluorALL)
# LOO_compare$diffs
# 
# LOO_compare2 <- LOO(brms_mixedonly_nostrainmixALL,brms_mixedonly_nofluorALL)
# LOO_compare2$diffs
# 
# # LOO_compare_exp2 <- LOO(brms_mixedonly_exp2, brms_mixedonly_exp2_nostrainmix,brms_mixedonly_exp2_nofluor)
# # LOO_compare_exp2$diffs
# 
# sink("03_analysis/model_loo_compare.txt")
# print("both experimetnts 1")
# LOO_compare$diffs
# 
# print("only single models")
# LOO_compare2$diffs
# 
# print("both model")
# brms_mixedonly
# print("fluor only")
# brms_mixedonly_nostrainmix
# print("Strainmix only")
# brms_mixedonly_nofluor
# 
# sink()
# 

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
  left_join(dat_plant_fullset)
# left_join(dat_plant)
factor(dat_mixedonly_exp2_withpredictions$transfer_h) %>% unique()
gg_brmfit_all <- dat_mixedonly_exp2_withpredictions %>%
  filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
  # mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
  geom_jitter(aes(col=factor(transfer_h)), width=0, height=0.1) +
  geom_ribbon(aes(fill=factor(transfer_h), ymin=Q2.5, ymax=Q97.5), alpha=0.1, show.legend = FALSE)+
  geom_line(aes(col=factor(transfer_h), y=Q50))+
  # geom_hline(aes(yintercept=0.5), col="grey", lty=2) +
  geom_vline(aes(xintercept=0), col="black", lty=2) +
  facet_grid(.~strain1,
             labeller=labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3="N2C3, then WCS365"))
             ) +
  # scale_color_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))+
  # scale_fill_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))+
  scale_color_manual(values=c(`0`="grey",`0.01` ="lightblue",`3` ="turquoise3",`6` ="skyblue3",`24` ="royalblue4",`48`="black" ))+
  scale_fill_manual(values=c(`0`="grey",`0.01` ="lightblue",`3` ="turquoise3",`6` ="skyblue3",`24` ="royalblue4",`48`="black" ))+
  theme_bw() +
  ylab("Probability of healthy plant") +
  xlab("WCS365:N2C3 (Log2 fold change fluorescence)")+ 
  labs(col="Delay in\nsecond strain (h)", fill="Delay in\nsecond strain (h)")+
  scale_y_continuous(breaks=c(0,0.5,1))
gg_brmfit_all
ggsave("03_analysis/gg_brmfit_all.png", gg_brmfit_all, height=4, width=8)


# 
# dat_mixedonly_exp2_withpredictionsALL <- posterior_linpred(brms_mixedonlyALL, transform=TRUE) %>%
#   apply(MARGIN=2, FUN=function(x) c(mean = mean(x), sd = sd(x), quantile(x,prob=0.5),quantile(x, prob=0.05), quantile(x, prob=0.95)
#                                     , quantile(x, prob=0.025), quantile(x, prob=0.975))) %>%
#   t() %>% as.data.frame() %>%
#   rename(pred_mean = mean
#          , pred_sd = sd
#          , Q50 = `50%`
#          , Q5 = `5%`
#          , Q95 = `95%`
#          , Q2.5 = `2.5%`
#          , Q97.5 = `97.5%`) %>%
#   bind_cols(brms_mixedonlyALL$data) %>%
#   left_join(dat_plant_fullset)
# # left_join(dat_plant)
# 
# gg_brmfit_ALL <- dat_mixedonly_exp2_withpredictionsALL %>%
#   filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
#   # mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
#   ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
#   geom_jitter(aes(col=factor(transfer_h)), width=0, height=0.1) +
#   geom_ribbon(aes(fill=factor(transfer_h), ymin=Q2.5, ymax=Q97.5), alpha=0.1, show.legend = FALSE)+
#   geom_line(aes(col=factor(transfer_h), y=Q50))+
#   geom_hline(aes(yintercept=0.5), col="grey", lty=2) +
#   geom_vline(aes(xintercept=0), col="grey", lty=2) +
#   facet_grid(.~strain1,
#              labeller=labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3="N2C3, then WCS365"))
#   ) +
#   # scale_color_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))+
#   # scale_fill_manual(values=c("black","darkgrey","yellow","gold","goldenrod2","goldenrod4"))+
#   scale_color_manual(values=c(`0`="grey",`0.01` ="lightblue",`3` ="turquoise3",`6` ="skyblue3",`24` ="royalblue4",`48`="black" ))+
#   scale_fill_manual(values=c(`0`="grey",`0.01` ="lightblue",`3` ="turquoise3",`6` ="skyblue3",`24` ="royalblue4",`48`="black" ))+
#   theme_bw() +
#   ylab("Probability of healthy plant") +
#   xlab("WCS365 : N2C3\nLog2 fold change")+ 
#   labs(col="Delay in\nsecond strain (h)", fill="Delay in\nsecond strain (h)")+
#   scale_y_continuous(breaks=c(0,0.5,1))
# gg_brmfit_ALL
# ggsave("03_analysis/gg_brmfit_ALL.png", gg_brmfit_ALL, height=4, width=8)


### Panelled plot for single dips ####

gg_priority_mixedonly_panelled <- plot_grid(gg_priority_mixedonly+labs(title=""), gg_brmfit_all+labs(title=""), nrow=2, align="v", rel_heights = c(1.2,1))
gg_priority_mixedonly_panelled
ggsave(filename="03_analysis/gg_priority_mixedonly_panelled.png", gg_priority_mixedonly_panelled, height=6, width=6)

gg_priority_mixed_dd_panelled <- plot_grid(gg_priority_mixed_dd+labs(title=""), gg_brmfit_all+labs(title=""), nrow=2, align="v", rel_heights = c(1.4,1))
gg_priority_mixed_dd_panelled
ggsave(filename="03_analysis/gg_priority_mixed_dd_panelled.png", gg_priority_mixed_dd_panelled, height=7, width=6)


gg_priority_mixed_dd
#### Control plots from brms ####

dat_plant_controls <- dat_plant_peffect %>%
  filter(StrainMix %in% c("WCS365-N2C3", "N2C3-WCS365")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("WCS365-N2C3", "N2C3-WCS365"))) %>%
  mutate(Alive = as.factor(Healthiness_hsv>400)) %>%
  # mutate(FirstStrain = ifelse(transfer_h==0, "Both", as.character(strain1))) %>%
  select(Alive, Healthiness_hsv, StrainMix, WN_FC_fluor_offset, transfer_h, experiment, Type, strain2_od, ratio) %>%
  filter(transfer_h==0 | Type=="MIX") %>%
  mutate(WN_FC_fluor_offset = ifelse(StrainMix=="MOCK-MOCK",0,WN_FC_fluor_offset)) 

gg_priority_mixedcontrols <- dat_plant_controls %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(!StrainMix %in% c("WCS365-N2C3","N2C3-WCS365")) %>%
  # mutate(Alive = Healthiness_hsv>400, domStrain = ifelse(WN_FC_fluor_offset>0, "WCS365 more\nabundant", "N2C3 more\nabundant")) %>%
  # mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  # mutate(Delay = ifelse(transfer_h == 0.01, "Dipped in first", paste0(transfer_h," hours"))) %>%
  # mutate(Delay = factor(Delay, levels=c("Dipped in first", "3 hours","6 hours", "24 hours", "48 hours"))) %>%
  ggplot(aes(x="Simultaneous\ninoculation", y=WN_FC_fluor_offset)) +
  geom_jitter(aes(col=factor(Alive), pch=experiment), width=0.1, height=0.2) +
  # facet_grid(.~ratio, scales = "free", drop=TRUE, space="free_x") +
  # scale_shape_manual(values=c(21,19))+
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3")) +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365:N2C3\nLog2 fold change fluorescence")+xlab("Time between first and second strain") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(col="Plant healthy\n(>400)", shape="Experiment") 
gg_priority_mixedcontrols
ggsave("03_analysis/gg_priority_mixedcontrols.png", gg_priority_mixedcontrols, height=4, width=4)

# stats
dat_plant_controls %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  summarise(n = n(), perc = sum(Healthiness_hsv>400))
dat_plant_controls %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  summarise(n = n(), perc = sum(WN_FC_fluor_offset>0))


#### Double dip plot ####
# Get invisible point to make plots look symmetric
extradata <- data.frame(StrainMix=c("WCS365-N2C3"), WN_FC_fluor_offset=c(10, -10), opaq=FALSE)
dat_plant_dd %>%
  select(strain0, strain1, strain2, StrainMix) %>% distinct()

dat_forddplot <- dat_plant_dd %>%
  # filter(experiment == "2024-02-14") %>%
  mutate(opaq=TRUE) %>%
  full_join(extradata) %>%
  full_join( dat_plant_controls %>%
               filter(!is.na(WN_FC_fluor_offset)) ) %>%
  mutate(Alive=Healthiness_hsv>400) %>%
  mutate(Control=ifelse(StrainMix %in% c("WCS365-MOCK", "N2C3-MOCK"), "Controls\n(Monoculture)", "Mixed\ninoculant")) %>%
  mutate(StrainMix2 = ifelse(is.na(strain2_od), gsub("-.*$","",StrainMix), StrainMix)) %>%
  mutate(StrainMixsimple = ifelse(StrainMix == "MOCK-MOCK", "MOCK"
                                  , ifelse(transfer_h==0, "Simultaneous\ninoculation"
                                    , ifelse(StrainMix == "WCS365-MOCK", "Dipped in WCS365;\nthen dipped in MOCK;\nplaced in MOCK well",
                                           ifelse(StrainMix == "N2C3-MOCK", "Dipped in N2C3;\nthen dipped in MOCK;\nplaced in MOCK well"
                                                  , ifelse(StrainMix == "WCS365-N2C3", "Dipped in WCS365;\nthen dipped in N2C3;\nplaced in MOCK well"
                                                           , ifelse(StrainMix == "N2C3-WCS365", "Dipped in N2C3;\nthen dipped in WCS365;\nplaced in MOCK well",NA))))))) %>%
  mutate(StrainMixsimple = factor(StrainMixsimple, levels=c(
    "MOCK"
    , "Dipped in WCS365;\nthen dipped in MOCK;\nplaced in MOCK well"
    , "Dipped in N2C3;\nthen dipped in MOCK;\nplaced in MOCK well"
    , "Simultaneous\ninoculation"
    , "Dipped in WCS365;\nthen dipped in N2C3;\nplaced in MOCK well"
    ,"Dipped in N2C3;\nthen dipped in WCS365;\nplaced in MOCK well"
  ))) %>%
  mutate(DippedTwice = !is.na(strain2_od)) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK", "WCS365-MOCK", "N2C3-MOCK", "WCS365-N2C3", "N2C3-WCS365"))) %>%
  filter(StrainMixsimple!="MOCK") %>%
  mutate(WN_FC_fluor_offset_adj = ifelse(Control == "Controls\n(Monoculture)", sign(WN_FC_fluor_offset), WN_FC_fluor_offset)) %>%
  mutate(WN_FC_fluor_offset_adj2 = ifelse(Control == "Controls\n(Monoculture)", NA, WN_FC_fluor_offset)) %>%
  mutate(plant_quality = ifelse(is.na(plant_quality), "good", plant_quality))

gg_doubledip_fluor_good  <- dat_forddplot %>%
  filter(plant_quality=="good") %>%
  ggplot(aes(x=StrainMixsimple)) +
  geom_boxplot(aes(y=WN_FC_fluor_offset_adj2), outlier.colour = NA) +
  geom_jitter(aes(y=WN_FC_fluor_offset_adj, col=Alive), height=0, width=0.2) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(.~Control, drop=TRUE, scales = "free")+
  # scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("WCS365:N2C3\nLog2 fold change") + xlab ("Dipping treatment") +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3"), na.value = rgb(0,0,0,0))  +
  labs(col="Plant healthy\n(>400)") +
  scale_y_continuous(breaks=c(-10,-5,0,5,10))
gg_doubledip_fluor_good
ggsave("03_analysis/gg_doubledip_fluor_good.png", gg_doubledip_fluor_good, height=4, width=6)

gg_doubledip_fluor_bothgoodandbad  <- dat_forddplot %>%
  # filter(plant_quality=="good") %>%
  ggplot(aes(x=StrainMixsimple)) +
  geom_boxplot(aes(y=WN_FC_fluor_offset_adj2), outlier.colour = NA) +
  geom_jitter(aes(y=WN_FC_fluor_offset_adj, col=Alive), height=0, width=0.2) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(.~Control, drop=TRUE, scales = "free")+
  # scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("WCS365:N2C3\nLog2 fold change") + xlab ("Dipping treatment") +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3"), na.value = rgb(0,0,0,0))  +
  labs(col="Plant healthy\n(>400)") +
  scale_y_continuous(breaks=c(-10,-5,0,5,10))
gg_doubledip_fluor_bothgoodandbad
ggsave("03_analysis/gg_doubledip_fluor_bothgoodandbad.png", gg_doubledip_fluor_bothgoodandbad, height=4, width=6)


gg_doubledip_fluor_good_nostrain0  <- dat_forddplot %>%
  filter(plant_quality=="good", !(!is.na(strain0) & experiment=="2024-03-20") ) %>%
  ggplot(aes(x=StrainMixsimple)) +
  geom_boxplot(aes(y=WN_FC_fluor_offset_adj2), outlier.colour = NA) +
  geom_jitter(aes(y=WN_FC_fluor_offset_adj, col=Alive), height=0, width=0.2) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(.~Control, drop=TRUE, scales = "free")+
  # scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("WCS365:N2C3\nLog2 fold change") + xlab ("Dipping treatment") +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3"), na.value = rgb(0,0,0,0))  +
  labs(col="Plant healthy\n(>400)") +
  scale_y_continuous(breaks=c(-10,-5,0,5,10))
gg_doubledip_fluor_good_nostrain0
ggsave("03_analysis/gg_doubledip_fluor_good_nostrain0.png", gg_doubledip_fluor_good_nostrain0, height=4, width=6)

gg_doubledip_fluor_nostrain0  <- dat_forddplot %>%
  filter(!(!is.na(strain0) & experiment=="2024-03-20") ) %>%
  ggplot(aes(x=StrainMixsimple)) +
  geom_boxplot(aes(y=WN_FC_fluor_offset_adj2), outlier.colour = NA) +
  geom_jitter(aes(y=WN_FC_fluor_offset_adj, col=Alive), height=0, width=0.2) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(.~Control, drop=TRUE, scales = "free")+
  # scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("WCS365:N2C3\nLog2 fold change") + xlab ("Dipping treatment") +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3"), na.value = rgb(0,0,0,0))  +
  labs(col="Plant healthy\n(>400)") +
  scale_y_continuous(breaks=c(-10,-5,0,5,10))
gg_doubledip_fluor_nostrain0
ggsave("03_analysis/gg_doubledip_fluor_nostrain0.png", gg_doubledip_fluor_nostrain0, height=4, width=6)


gg_doubledip_fluor_nostrain0_nosiminoc  <- dat_forddplot %>%
  filter(!(!is.na(strain0) & experiment=="2024-03-20") ) %>%
  filter(StrainMixsimple !="Simultaneous\ninoculation") %>%
  mutate(Control = ifelse(Control=="Mixed\ninoculant", "Sequentially\ndipped", Control)) %>%
  ggplot(aes(x=StrainMixsimple)) +
  geom_boxplot(aes(y=WN_FC_fluor_offset_adj2), outlier.colour = NA) +
  geom_jitter(aes(y=WN_FC_fluor_offset_adj, col=Alive), height=0, width=0.2) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(.~Control, drop=TRUE, scales = "free")+
  # scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("WCS365:N2C3\nLog2 Fold change fluorescence") + xlab ("Dipping treatment") +
  scale_color_manual(values=c(`TRUE`="darkgreen", `FALSE` = "yellow3"), na.value = rgb(0,0,0,0))  +
  labs(col="Plant healthy\n(>400)") +
  scale_y_continuous(breaks=c(-10,-5,0,5,10))
gg_doubledip_fluor_nostrain0_nosiminoc
ggsave("03_analysis/gg_doubledip_fluor_nostrain0_nosiminoc.png", gg_doubledip_fluor_nostrain0_nosiminoc, height=4, width=6)

dat_forddplot %>% 
  filter(strain2=="MOCK") %>% View()

#### CFU stuff ####
allCFU <- allCFU %>%
  mutate(experiment = ifelse(experiment == "2023-12_priorityeffects", "2023-11-08",
                             ifelse(experiment == "2023-12_priorityeffects2", "2023-12-06", 
                                    ifelse(experiment == "2024-02_Peffect_Trial3", "2024-01-30", 
                                           ifelse(experiment == "2024-02_Peffect_Trial4", "2024-02-14", 
                                                  ifelse(experiment == "2024_03_Peffect_Trial5","2024-03-20",experiment))))))

allCFU %>% 
  filter(!is.na(abs_600))
allCFU %>% 
  filter(experiment=="2023-12-06") 

allCFU_sum <- allCFU %>%
  # filter(experiment == "2024-02-14") %>% 
  group_by(experiment, transfer_h, strain, plant, abs_600, plant_quality) %>%
  summarise(meanCount = mean(CFU_well), sd = sd(CFU_well), n()) %>%
  mutate(abs_600 = ifelse(is.na(abs_600), "root", as.character(abs_600))) %>%
  mutate(treatment = ifelse(is.na(abs_600), "Plant root", paste0("Abs600=", abs_600,"\nsecond strain inoculant"))) %>% ungroup()

allCFU_sum %>%
  ggplot() +
  geom_pointrange(aes(x=factor(transfer_h), y=log10(meanCount), ymin=log10(meanCount-sd), ymax=log10(meanCount+sd)
                      , col=factor(treatment), pch=experiment))+
  facet_grid(.~strain)

allCFU_sum %>%
  ggplot() +
  geom_pointrange(aes(x=factor(transfer_h), y=log10(meanCount), ymin=log10(meanCount-sd), ymax=log10(meanCount+sd)
                      , col=factor(treatment), pch=plant_quality))+
  facet_grid(.~strain)

### WORKING FROM HERE: Keep the good/bad plant qualifiers in
# JUST overnights
justON_cfu <- allCFU_sum %>%
  filter(abs_600 !="root") %>%
  unite(transfer_h, abs_600, col="overnightType") %>%
  select(overnightType, strain, meanCount, sd ) %>%
  rename(meanCount_strain2 = meanCount, sd_strain2 = sd, strain2 = strain) %>%
  mutate(strain1 = ifelse(strain2=="WCS365", "N2C3", 
                          ifelse(strain2=="N2C3", "WCS365", NA))) 
# JUST plant roots; merge with well data
cfu_log2fc <- allCFU_sum %>%
  filter(abs_600 =="root") %>%
  rename(meanCount_strain1 = meanCount, sd_strain1 = sd) %>%
  rename(strain1=strain) %>%
  select(experiment, transfer_h, strain1, meanCount_strain1, sd_strain1) %>%
  full_join(justON_cfu, relationship = "many-to-many") %>%
  mutate(ratio_mean_21 = meanCount_strain2/meanCount_strain1
         , log2fc_mean_S21 = log2(ratio_mean_21)) %>%
  mutate(log2fc_mean_WN = ifelse(strain2=="WCS365", log2fc_mean_S21,
                                 ifelse(strain2=="N2C3", -log2fc_mean_S21, NA))) %>%
  separate(overnightType, into=c("remove","strain2_od"), sep="_", remove=FALSE) %>%
  select(-remove) %>%mutate(strain2_od = as.numeric(strain2_od)) %>%
  mutate(transfer_h = as.numeric(transfer_h))
cfu_log2fc # THIS IS the log2fc of root vs well cells

# merge with fluor ratios and plant data
dat_plant_fullset %>%
  select(experiment, transfer_h, strain1, strain2, strain2_od, WN_FC_fluor_offset, Healthiness_hsv) %>%
  filter(experiment %in% c("2024-02-14","2023-12-06" ), strain1 %in% c("N2C3", "WCS365"), strain2 %in% c("WCS365", "N2C3")) %>%
  select(-WN_FC_fluor_offset, -Healthiness_hsv) %>% distinct()
cfu_log2fc %>% select(experiment, transfer_h, strain1, strain2, strain2_od,overnightType) %>% distinct()

merged_cfu_plant <- dat_plant_fullset %>%
  select(experiment, transfer_h, strain1, strain2, strain2_od, WN_FC_fluor_offset, Healthiness_hsv) %>%
  filter(experiment %in% c("2024-02-14","2023-12-06" ), strain1 %in% c("N2C3", "WCS365"), strain2 %in% c("WCS365", "N2C3")) %>%
  right_join(cfu_log2fc %>% select(experiment, transfer_h, strain1, strain2, strain2_od, log2fc_mean_S21, log2fc_mean_WN, overnightType), relationship="many-to-many") %>% 
  filter(!is.na(WN_FC_fluor_offset)) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400))
# dat_plant_fullset %>% filter(strain2_od==0.002) %>% View()

gg_inoccfu_vs_cluof_6_0.001 <- merged_cfu_plant %>%
  filter(overnightType!="6_0.001") %>%
  filter(experiment == "2024-02-14") %>%
  ggplot(aes(x=log2fc_mean_WN, y=WN_FC_fluor_offset, col=factor(strain1))) +
  # geom_point(aes(col=factor(Alive)))+
  geom_point(aes(pch=factor(transfer_h)))+
  geom_smooth(method="lm")+
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  # facet_grid(experiment~., labeller = labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3 = "N2C3, then WCS365")))+
  xlab("WCS365:N2C3 at time of inoculation\n(Log2 fold change CFU counts)") +
  ylab("WCS365:N2C3 after 1 week plant growth\n(Log2 fold change fluorescence)")+
  scale_color_manual(values=c(WCS365="darkseagreen4", N2C3="darkorange"))+
  labs(shape="Inoculation\ndelay (h)", col="First strain", title="Used CFU counts from 6h")


gg_inoccfu_vs_cluof_24_0.001 <- merged_cfu_plant %>%
  filter(overnightType!="24_0.001") %>%
  filter(experiment == "2024-02-14") %>%
  ggplot(aes(x=log2fc_mean_WN, y=WN_FC_fluor_offset, col=factor(strain1))) +
  # geom_point(aes(col=factor(Alive)))+
  geom_point(aes(pch=factor(transfer_h)))+
  geom_smooth(method="lm")+
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  # facet_grid(experiment~., labeller = labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3 = "N2C3, then WCS365")))+
  xlab("WCS365:N2C3 at time of inoculation\n(Log2 fold change CFU counts)") +
  ylab("WCS365:N2C3 after 1 week plant growth\n(Log2 fold change fluorescence)")+
  scale_color_manual(values=c(WCS365="darkseagreen4", N2C3="darkorange"))+
  labs(shape="Inoculation\ndelay (h)", col="First strain", title="Used CFU counts from 24h")

gg_cfu_root_vs_well <- plot_grid(gg_inoccfu_vs_cluof_6_0.001, gg_inoccfu_vs_cluof_24_0.001, align = "h")
gg_cfu_root_vs_well
ggsave(filename = "03_analysis/gg_cfu_root_vs_well.png", gg_cfu_root_vs_well, width = 10, height=4)

## Look at simultaneous inoculation ##
# Get cells only

justON_expand <- justON_cfu %>%
  select(-strain1) %>%
  rename(strain1 = strain2, meanCount_strain1 = meanCount_strain2, sd_strain1 = sd_strain2, overnightType1 = overnightType) %>%
  mutate(strain2 = ifelse(strain1=="WCS365","N2C3",
                          ifelse(strain1=="N2C3", "WCS365", NA))) %>%
  full_join(justON_cfu %>% rename(overnightType2 = overnightType), relationship = "many-to-many") %>%
  separate(overnightType1, into=c("type1","strain1_od"), sep="_", remove=FALSE) %>%
  separate(overnightType2, into = c("type2","strain2_od"), sep="_", remove=FALSE) %>%
  mutate(strain1_od = as.numeric(strain1_od), strain2_od = as.numeric(strain2_od)) %>%
  filter(!(strain1_od==0.001 & strain2_od==0.001 & type1!=type2))

merged_mixedinoc_cfu <- dat_plant_fullset %>%
  filter(Type=="MIX") %>%
  # filter(strain2_od !=0.002) %>%
  select(experiment, transfer_h, strain1, strain1_od, strain2, strain2_od, WN_FC_fluor_offset, Healthiness_hsv) %>%
  left_join(justON_expand, relationship = "many-to-many") %>%
  mutate(ratio_inoc_S21 = meanCount_strain2/meanCount_strain1, 
         log2fc_mean_S21 = log2(ratio_inoc_S21)) %>%
  mutate(log2fc_mean_WN = ifelse(strain2=="WCS365", log2fc_mean_S21
                                 , ifelse(strain2=="N2C3", -log2fc_mean_S21, NA))) 

gg_cfu_well_vs_well_6_0.001 <- merged_mixedinoc_cfu%>%
  filter(!(type1=="24" & type2=="24")) %>%
  filter(experiment == "2024-02-14") %>%
  ggplot(aes(x=log2fc_mean_WN, y=WN_FC_fluor_offset, col=factor(strain1))) +
  # geom_point(aes(col=factor(Alive)))+
  geom_point(aes(pch=factor(transfer_h)))+
  geom_smooth(method="lm")+
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  # facet_grid(experiment~., labeller = labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3 = "N2C3, then WCS365")))+
  xlab("WCS365:N2C3 at time of inoculation\n(Log2 fold change CFU counts)") +
  ylab("WCS365:N2C3 after 1 week plant growth\n(Log2 fold change fluorescence)")+
  scale_color_manual(values=c(WCS365="darkseagreen4", N2C3="darkorange"))+
  labs(shape="Inoculation\ndelay (h)", col="First strain", title="Used CFU counts from 6h")
gg_cfu_well_vs_well_6_0.001
ggsave(filename = "03_analysis/gg_cfu_well_vs_well_6_0.001.png", gg_cfu_well_vs_well_6_0.001, height=4, width=5.5)


gg_cfu_well_vs_well_24_0.001 <- merged_mixedinoc_cfu%>%
  filter(!(overnightType1=="6_0.001" & overnightType2=="6_0.001")) %>%
  filter(experiment == "2024-02-14") %>%
  ggplot(aes(x=log2fc_mean_WN, y=WN_FC_fluor_offset, col=factor(strain1))) +
  # geom_point(aes(col=factor(Alive)))+
  geom_point(aes(pch=factor(transfer_h)))+
  geom_smooth(method="lm")+
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  # facet_grid(experiment~., labeller = labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3 = "N2C3, then WCS365")))+
  xlab("WCS365:N2C3 at time of inoculation\n(Log2 fold change CFU counts)") +
  ylab("WCS365:N2C3 after 1 week plant growth\n(Log2 fold change fluorescence)")+
  scale_color_manual(values=c(WCS365="darkseagreen4", N2C3="darkorange"))+
  labs(shape="Inoculation\ndelay (h)", col="First strain", title="Used CFU counts from 24h")
gg_cfu_well_vs_well_24_0.001
ggsave(filename = "03_analysis/gg_cfu_well_vs_well_24_0.001.png", gg_cfu_well_vs_well_24_0.001, height=4, width=5.5)

##### Merged #######

merged_mixedinoc_cfu_filt <- merged_mixedinoc_cfu%>%
  # filter(!(overnightType1=="6_0.001" & overnightType2=="6_0.001")) %>%
  filter(!(overnightType1=="24_0.001" & overnightType2=="24_0.001")) %>%
  filter(experiment == "2024-02-14") %>%
  select(transfer_h, strain1, strain2, WN_FC_fluor_offset, Healthiness_hsv, log2fc_mean_WN, overnightType1, , overnightType2) %>%
  mutate(Type = "Simultaneous\nInoculation")

merged_cfu_plant_filt <- merged_cfu_plant %>%
  # filter(overnightType!="6_0.001") %>%
  filter(overnightType!="24_0.001") %>%
  filter(experiment == "2024-02-14") %>%
  select(transfer_h, strain1, strain2, WN_FC_fluor_offset, Healthiness_hsv, log2fc_mean_WN, overnightType) %>%
  mutate(Type="PlantvsWell")

gg_priority_vs_cellload <- full_join(merged_mixedinoc_cfu_filt,merged_cfu_plant_filt) %>%
  mutate(Type_strain1 = ifelse(Type=="PlantvsWell", paste0(strain1, " on root first"), Type)) %>%
  filter(transfer_h!=0) %>%
  mutate(transfer_h_adj = ifelse(Type_strain1=="Simultaneous\nInoculation", "Simultaneous", paste0(transfer_h,"h"))) %>%
  ggplot(aes(x=log2fc_mean_WN, y=WN_FC_fluor_offset, col=factor(Type_strain1))) +
  # geom_point(aes(col=factor(Alive)))+
  geom_point(aes(pch=factor(transfer_h_adj)))+
  geom_smooth( method="lm")+
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  # facet_grid(strain1~., labeller = labeller(strain1=c(WCS365="WCS365, then N2C3", N2C3 = "N2C3, then WCS365")))+
  xlab("WCS365:N2C3 at time of inoculation\n(Log2 fold change CFU counts)") +
  ylab("WCS365:N2C3 after 1 week plant growth\n(Log2 fold change fluorescence)")+
  # scale_color_manual(values=c(MixedCulture="black", PlantvsWell="forestgreen"))+
  scale_color_manual(values=c("WCS365 on root first"="darkseagreen4", "N2C3 on root first"="darkorange", "Simultaneous\nInoculation"="black"))+
  scale_shape_manual(values=c("Simultaneous"=21, "24h"=15, "6h"=17)) +
  labs(shape="Inoculation\n delay", col="Order of arrival")+
  theme_bw()
gg_priority_vs_cellload
ggsave(filename = "03_analysis/gg_priority_vs_cellload.png", gg_priority_vs_cellload, height=4, width=6)


