#!bin/bash
library(cowplot)
library(gridExtra)
library(lubridate)
library(ggstance)
# library(brms)
# library(lme4)
# library(lmerTest)
library(tidyverse)
library(brms)

dir.create("04_ratio_experiment")
########## Load ########
load("03_EDA_and_adj/allDat_final.RData")
load("03_EDA_and_adj/plant_colours.RData")

allDat_final <- allDat_final %>%
  mutate(lacZ_path = ifelse(is.na(lacZ_path), "WT",lacZ_path)) %>%
  mutate(lacZ_path = factor(lacZ_path, levels=c("WT","L"))) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400))


unique(allDat_final$experiment)

# dat_ac <- allDat_final %>%
#   filter(experiment %in% c("2023-06-14_age_concentration_Pf5_PROPER", "2023-07-26_ac_ccpwn"
#                            ,"2023-08-02_ac_ccpwn")) 
dat_ac <- allDat_final %>%
  filter(experiment %in% c("2023-06-14_age_concentration_Pf5_PROPER"
                           , "2023-06-21_ac_ccpwn"
                           , "2023-06-27_ac_ccpwn"
                           , "2023-07-26_ac_ccpwn"
                           ,"2023-08-02_ac_ccpwn"))%>%
    # filter(experiment %in% c("2023-06-14_age_concentration_Pf5_PROPER", "2023-07-26_ac_ccpwn"
    #                          ,"2023-08-02_ac_ccpwn"))  %>%
  filter(lacZ_path=="WT") %>% # keep just singles from first lacz exp
  filter(is.na(protect2)) %>%
  mutate(protect_v2 = ifelse(is.na(protect_v2), protect, protect_v2))

dat_ac %>%
  select(path_od, protect_od) %>% table()
dat_ac %>% select(lacZ_path, experiment) %>% table()

#### Plant alive/dead stats ####
sink("04_ratio_experiment/plant_summaries.txt")
print("protect path ratio\n")
dat_ac %>%
  group_by(protect, path, ratio) %>%
  summarise(totalP = n(), alive = sum(Alive), propAlive = sum(Alive)/n()) %>% 
  as.data.frame()
print("\nratio\n")
dat_ac %>%
  group_by(ratio) %>%
  summarise(totalP = n(), alive = sum(Alive), propAlive = sum(Alive)/n())
sink()

#### Labellers #####
plant_age_labeller <- c(`5`="Inoculated at\n5 days old", `6`="Inoculated at\n6 days old", `7`="Inoculated at\n7 days old")
N2C3_labeller <- c(`0`="No Pathogen",`1`="N2C3 ~ 10 cells",`2`="N2C3 ~ 10^2 cells", `3`="N2C3 ~ 10^3 cells", `4`="N2C3 ~ 10^4 cells", `5`="N2C3 ~ 10^5 cells",`6`="N2C3 ~ 10^6 cells", `7`="N2C3 ~ 10^7 cells", `8`="N2C3 ~ 10^8 cells")
WCS_labeller <- c(`0`="No Protective",`1`="WCS365 ~ 10 cells",`2`="WCS365 ~ 10^2 cells", `3`="WCS365 ~ 10^3 cells", `4`="WCS365 ~ 10^4 cells", `5`="WCS365 ~ 10^5 cells",`6`="WCS365 ~ 10^6 cells", `7`="WCS365 ~ 10^7 cells", `8`="WCS365 ~ 10^8 cells")
protect_labeller = c(MOCK = "MOCK", WCS365="WCS365", CHAO = "CHA0", CH267 = "CH267", PF5="Pf5")

####### Tile plots #######
# sdr <- data.frame(plant_age=c(5,6,7), lr = c(0.875, 0.817, 0.939))
dat_ccpw_col0_summarised <- dat_ac %>%
  group_by(protect,protect_cells_log, path_cells_log, plant_age) %>%
  summarise(medHealth = median(Healthiness_hsv), propAlive = sum(Alive)/n()) %>% ungroup()

gg_allstrains_competition <- dat_ac %>% 
  # filter(protect=="WCS365") %>%
  ggplot(aes(x=factor(round(protect_cells_log)), factor(round(path_cells_log))))+
  geom_point(data=dat_ccpw_col0_summarised, mapping=aes(x=factor(round(protect_cells_log)), factor(round(path_cells_log))
                                                        , fill=medHealth),
             size=25, pch=22)+
  geom_jitter(aes(col=UniqueID, size=all_plant_pixels), show.legend = FALSE,  position = position_jitterdodge(jitter.width = 0.2, jitter.height=0.2, seed=1)) +
  # scale_fill_gradient(low="black",high="white")+
  scale_fill_gradientn(colours = c("black","white") )+
  # scale_fill_gradientn(colours = c("black","white"), limits=c(0,1))+
  # scale_fill_manual(values=c("black","white"))+
  scale_color_manual(values = col_all) +
  facet_grid(protect~plant_age, labeller = labeller(plant_age = plant_age_labeller)) +
  xlab("Protect cells (~ log10(+1))") + ylab("N2C3 cells (~ log10(+1))") +
  labs(fill="Health Score")
gg_allstrains_competition
ggsave(filename = "04_ratio_experiment/gg_allstrains_competition.png", gg_allstrains_competition,
       height=15, width=14)


#### flipped version of plot ####
gg_allstrains_competition_FLIPPED <- dat_ac %>% 
  # filter(protect=="WCS365") %>%
  ggplot(aes(x=factor(round(protect_cells_log)), factor(round(path_cells_log))))+
  geom_point(data=dat_ccpw_col0_summarised, mapping=aes(x=factor(round(protect_cells_log)), factor(round(path_cells_log))
                                                        , fill=medHealth),
             size=25, pch=22)+
  geom_jitter(aes(col=UniqueID, size=all_plant_pixels), show.legend = FALSE,  position = position_jitterdodge(jitter.width = 0.2, jitter.height=0.2, seed=1)) +
  # scale_fill_gradient(low="black",high="white")+
  scale_fill_gradientn(colours = c("black","white"))+
  # scale_fill_gradientn(colours = c("black","white"), limits=c(0,1))+
  # scale_fill_manual(values=c("black","white"))+
  scale_color_manual(values = col_all) +
  facet_grid(plant_age~protect, labeller = labeller(plant_age = plant_age_labeller), scales="free_x", space = "free_x") +
  xlab("Protect cells (~ log10(+1))") + ylab("N2C3 cells (~ log10(+1))") +
  labs(fill="Health Score")
gg_allstrains_competition_FLIPPED
ggsave(filename = "04_ratio_experiment/gg_allstrains_competition_FLIPPED.png", gg_allstrains_competition_FLIPPED,
       height=9, width=17)
# 
# #### BRMS ratio vs abundance ############
# 
# dat_ac_monocultures_nolacZ <- dat_ac %>%
#   filter(protect =="MOCK" | path == "MOCK") %>%
#   # filter(!(protect=="MOCK"&path=="MOCK")) %>%
#   mutate(Strain = ifelse(protect == "MOCK", as.character(path), as.character(protect))) %>%
#   select(Alive, Healthiness_hsv, Strain, protect, path, plant_age, total_cells_log, plate, experiment, path_cells_log, protect_cells_log) %>%
#   mutate(Strain = factor(Strain, levels=c("MOCK","N2C3","WCS365","CHAO","CH267","PF5"))) %>%
#   mutate(plant_age_adj = plant_age-5)
# # mutate(Strain = factor(Strain, levels=c("N2C3","WCS365","CHAO","CH267","PF5"))) 
# 
# dat_ac_monocultures_nolacZ %>%
#   filter(protect == "MOCK" & path=="MOCK") %>%
#   ggplot()+
#   geom_jitter(aes(x=protect, y=Alive))
# # dat_ac_nopath_nolacz <- dat_ac %>%
# # filter(path =="MOCK", protect !="MOCK") 
# 
# # Get baselines for monocultures
# brm_monocultures <- brm(Alive ~ Strain + Strain:total_cells_log + plant_age + (1 | plate), data=dat_ac_monocultures_nolacZ
#                                        , seed=12498215
#                                        , family="bernoulli"
#                                        , prior = c(set_prior("constant(0)", class = "b", coef = "StrainMOCK:total_cells_log"))
#                                        , file="04_ratio_experiment/brm_monocultures"
#                                        , iter=4000)
# brm_monocultures
# 
# # Plot mock estimates
# draws_monoculture  <- as_draws_df(brm_monocultures) %>%
#   rownames_to_column(var="iter") %>%
#   select(iter, starts_with("b_")) %>%
#   pivot_longer(-iter, names_to="coef", values_to="draw") %>%
#   mutate(coef = gsub("b_","",coef)) %>%
#   group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
#   mutate(mod = "Monoculture") %>%
#   mutate(Coefficient = str_replace_all(coef,c(plant_age="MOCK:Effect of plant age\nat inoculation",Intercept="MOCK:Presence of strain"))) %>%
#   separate(Coefficient, sep=":", into=c("Strain","Coefficient"), fill = "right") %>%
#   mutate(Strain = gsub("Strain","",Strain)) %>%
#   mutate(Coefficient = ifelse(is.na(Coefficient), "Presence of strain", ifelse(
#     Coefficient=="total_cells_log", "Cell load", Coefficient))) %>%
#   mutate(Strain = factor(Strain, levels=c("MOCK","N2C3","WCS365","CHAO","CH267","PF5"))) %>%
#   mutate(Coefficient = factor(Coefficient, levels=c("Effect of plant age\nat inoculation","Presence of strain","Cell load")))
# 
# gg_monoculture_forpriors <- draws_monoculture %>%
#   filter(Strain !="MOCK") %>%
#   # filter(!(Strain=="MOCK" & Coefficient == "Cell density")) %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~Strain) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
#   ylab("Estimate posterior")+ 
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# gg_monoculture_forpriors
# 
# 
# # Look at estimates
# mpriors <- brm_monocultures$fit %>%
#   as.data.frame() %>%
#   select(starts_with("b_")) %>%
#   apply(MARGIN=2, FUN=function(x) c(m = mean(x), std = sd(x))) %>%
#   as_tibble() %>%
#   rename_all(~gsub("b_","",.)) %>% as.data.frame()
# mpriors
# 
# ###### competition only ###
# dat_ac_comp_nolacz <- dat_ac %>%
#   filter(protect !="MOCK" & path !="MOCK") %>%
#   mutate(ratiofc = log10(protect_cells/path_cells)) %>%
#   select(Alive, Healthiness_hsv, protect, path, plant_age, total_cells_log, plate, experiment, path_cells_log, protect_cells_log,ratiofc)
# 
# 
# brm_comp <- brm(Alive ~  plant_age + path_cells_log + protect_cells_log:protect +  ratiofc:protect +(1 | plate) , data=dat_ac_comp_nolacz
#                 , seed=12498215
#                 , family="bernoulli"
#                 , file="04_ratio_experiment/brm_comp"
#                 , prior = c(
#                   set_prior(paste0("normal(",mpriors[,"StrainWCS365:total_cells_log"][1],",",mpriors[,"StrainWCS365:total_cells_log"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                   set_prior(paste0("normal(",mpriors[,"StrainCHAO:total_cells_log"][1],",",mpriors[,"StrainCHAO:total_cells_log"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                   set_prior(paste0("normal(",mpriors[,"StrainCH267:total_cells_log"][1],",",mpriors[,"StrainCH267:total_cells_log"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                   set_prior(paste0("normal(",mpriors[,"StrainPF5:total_cells_log"][1],",",mpriors[,"StrainPF5:total_cells_log"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                   set_prior(paste0("normal(",mpriors[,"StrainN2C3:total_cells_log"][1],",",mpriors[,"StrainN2C3:total_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                   set_prior(paste0("normal(",mpriors[,"plant_age"][1],",",mpriors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                   # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                 ),iter=4000)
# brm_comp
# 
# # Plot mock estimates
# draws_allstrains_bin_wpriors_RATIO  <- as_draws_df(brm_grad_allstrains_bin_wpriors_RATIO) %>%
#   rownames_to_column(var="iter") %>%
#   select(iter, starts_with("b_")) %>%
#   pivot_longer(-iter, names_to="coef", values_to="draw") %>%
#   mutate(coef = gsub("b_","",coef)) %>%
#   group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
#   mutate(coef = gsub("protect_cells_log:protect",":protect_cells_log__",coef)) %>%
#   separate(coef, sep="__", into=c("coef","protect")) %>%
#   mutate(coef = str_replace_all(coef, c(`^path_cells_log$`="N2C3 cells",plant_age="Effect of plant age\nat inoculation",
#                                         `^:protect_cells_log$`="Protective cells", 
#                                         `path_cells_log::protect_cells_log`="    Interaction between\nN2C3 and protective"))) %>%
#   mutate(coef = factor(coef, levels=c("Effect of plant age\nat inoculation","N2C3 cells","Protective cells","    Interaction between\nN2C3 and protective"))) %>%
#   mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   mutate(protect = str_replace_all(protect, c(CHAO='CHA0', PF5='Pf5'))) %>%
#   mutate(protect = factor(protect, levels=c("Shared effects\nbetween all\nprotective strains","WCS365","CHA0","CH267","Pf5")))
# 
# # fixef(brm_grad_allstrains_bin_wpriors, probs = c(0.025, 0.05, 0.95, 0.975)) %>% View()
# # fixef(brm_grad_mockonly_bin_forpriors) %>% View()
# 
# gg_allstrains_wpriors_RATIO <- draws_allstrains_bin_wpriors_RATIO %>%
#   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   filter(coef!="Intercept") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~protect, scales="free_x") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         # , axis.title.x = element_blank()
#         , axis.title.y = element_blank()
#   )+
#   xlab("Coefficients")+
#   labs(col="Probability of\ndirection")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# gg_allstrains_wpriors_RATIO
# 
# 
# ##### simplified brms coef plots ######
# 
# gg_mockonly_forpriors_simplified <- draws_mockonly_bin_forpriors %>%
#   filter(coef!="Intercept") %>%
#   filter(coef!="plant_age") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = FALSE) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~mod) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         , axis.title.x = element_blank())+
#   ylab("Estimate posterior")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
#   ylim(-4, (0.4/0.25)*4)
# gg_mockonly_forpriors_simplified
# 
# ### Ratio version
# gg_allstrains_wpriors_RATIO_simplified <- draws_allstrains_bin_wpriors_RATIO %>%
#   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   filter(coef!="Intercept") %>%
#   filter(protect !="Shared effects\nbetween all\nprotective strains") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~protect, scales="free_x") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         # , axis.title.x = element_blank()
#         , axis.title.y = element_blank()
#   )+
#   xlab("Coefficients")+
#   labs(col="Probability of\ndirection")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
#   ylim(-0.25, 0.4)
# # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# gg_allstrains_wpriors_RATIO_simplified
# 
# gg_bayesmodel_nolacz_bin_RATIO_simp <- plot_grid(gg_mockonly_forpriors_simplified, gg_allstrains_wpriors_RATIO_simplified
#                                                  , nrow=1, align="h", axis = "tb", rel_widths = c(1,4)) 
# gg_bayesmodel_nolacz_bin_RATIO_simp
# ggsave(filename = "04_ratio_experiment/gg_bayesmodel_nolacz_bin_RATIO_simp.png"
#        , gg_bayesmodel_nolacz_bin_RATIO_simp, height=4, width=8)
# 
# ### now for plant age only
# 
# gg_plantage_RATIO_simplified  <- draws_mockonly_bin_forpriors %>%
#   filter(coef == "plant_age") %>%
#   full_join(draws_allstrains_bin_wpriors_RATIO %>% filter(coef == "Effect of plant age\nat inoculation") %>%
#               rename(Coefficient = coef)) %>%
#   mutate(protect = ifelse(is.na(protect), "Pathogen only", "Protective and \npathogen together\n(across all strains)")) %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=protect, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~"Pooled effect\nacross all treatments") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#   )+
#   ylab("Estimate posterior")+xlab("Coefficient")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# gg_plantage_RATIO_simplified
# ggsave(filename = "04_ratio_experiment/gg_plantage_RATIO_simplified.png", gg_plantage_RATIO_simplified, height=3, width=4)
# 
# ## ratio verison
# brm_grad_allstrains_bin_wpriors_NOINTER_RATIO <- brm(Alive ~  plant_age + path_cells_log + protect_cells_log:protect +  (1 | plate) , data=dat_ac_comp_nolacz
#                                                      , seed=12498215
#                                                      , family="bernoulli"
#                                                      , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_NOINTER_RATIO"
#                                                      , prior = c(
#                                                        # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                                        # set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                                        set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                                        set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                                        set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                                        set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                                        set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                                        set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                                        # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                                      ),iter=4000)
# brm_grad_allstrains_bin_wpriors_NOINTER_RATIO
# 
# brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO <- brm(Alive ~  plant_age + protect_cells_log:path_cells_log:protect +  (1 | plate) , data=dat_ac_comp_nolacz
#                                                        , seed=12498215
#                                                        , family="bernoulli"
#                                                        , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO"
#                                                        , prior = c(
#                                                          # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                                          # set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                                          # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                                          # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                                          # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                                          # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                                          # set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                                          set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                                          # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                                        ),iter=4000)
# brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO
# loo_of_cellvsratio_RATIO <- loo(brm_grad_allstrains_bin_wpriors_RATIO, brm_grad_allstrains_bin_wpriors_NOINTER_RATIO, brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO)
# loo_of_cellvsratio_RATIO
# sink("loo_of_cellvsinteraction_RATIO.txt")
# loo_of_cellvsratio_RATIO
# sink()
# 
# 
# #### BRMS all strains, but with prior for mocks, nolacz ############
# 
# dat_ac_onlymock_nolacZ <- dat_ac %>%
#   filter(protect =="MOCK") 
# # dat_ac_nopath_nolacz <- dat_ac %>%
#   # filter(path =="MOCK", protect !="MOCK") 
# dat_ac_nomock_nolacz <- dat_ac %>%
#   filter(protect !="MOCK") 
# # dat_ac_nomock_nolacz2<- dat_ac %>%
#   # filter(protect !="MOCK", path!="MOCK") 
# dat_ac_monoculture_nolacz <- dat_ac %>%
#   filter(path =="MOCK", protect !="MOCK") 
# # Get baselines for pathogens
# brm_grad_mockonly_bin_forpriors <- brm(Alive ~ plant_age + path + path_cells_log + (1 | plate), data=dat_ac_onlymock_nolacZ
#                                        , seed=12498215
#                                        , family="bernoulli"
#                                        , file="04_ratio_experiment/brm_grad_onlymock_bin_forpriors_nolacZ"
#                                        , iter=4000)
# 
# # Plot mock estimates
# draws_mockonly_bin_forpriors  <- as_draws_df(brm_grad_mockonly_bin_forpriors) %>%
#   rownames_to_column(var="iter") %>%
#   select(iter, starts_with("b_")) %>%
#   pivot_longer(-iter, names_to="coef", values_to="draw") %>%
#   mutate(coef = gsub("b_","",coef)) %>%
#   group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
#   mutate(mod = "Model with\nno protectives") %>%
#   mutate(Coefficient = str_replace_all(coef,c(plant_age="Effect of plant age\nat inoculation",pathN2C3="Presence of N2C3", path_cells_log="N2C3 cells", lacZ_pathL="Effect of\nlacZ insertion",Intercept="Plant with no\ninoculant"))) %>%
#   mutate(Coefficient = factor(Coefficient, levels=c("Plant with no\ninoculant", "Effect of plant age\nat inoculation","Effect of\nlacZ insertion", "Presence of N2C3", "N2C3 cells"))) 
# 
# gg_mockonly_forpriors <- draws_mockonly_bin_forpriors %>%
#   filter(coef!="Intercept") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = FALSE) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~mod) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         , axis.title.x = element_blank())+
#   ylab("Estimate posterior")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# gg_mockonly_forpriors
# 
# 
# # Look at estimates
# estimate_for_priors <- brm_grad_mockonly_bin_forpriors$fit %>%
#   as.data.frame() %>%
#   select(starts_with("b_")) %>%
#   apply(MARGIN=2, FUN=function(x) c(m = mean(x), std = sd(x))) %>%
#   as_tibble() %>%
#   rename_all(~gsub("b_","",.)) %>% as.data.frame()
# estimate_for_priors
# estimate_for_priors[,"pathN2C3"][1]
# 
# ## Commensals in monoculture
# brm_grad_monoculture_bin_wpriors <- brm(Alive ~  -1 + plant_age + protect_cells_log:protect + (1 | plate) , data=dat_ac_monoculture_nolacz
#                                        , seed=12498215
#                                        , family="bernoulli"
#                                        , file="04_ratio_experiment/brm_grad_monoculture_bin_wpriors"
#                                        , prior = c(
#                                          # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                          # set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                          # set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                          set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                          # set_prior("normal(0,2.5)", class="b", coef="protectWCS365"),
#                                        ),iter=4000)
# brm_grad_monoculture_bin_wpriors
# # Plot monoculture estimates
# draws_monoculture_bin_forpriors  <- as_draws_df(brm_grad_monoculture_bin_wpriors) %>%
#   rownames_to_column(var="iter") %>%
#   select(iter, starts_with("b_")) %>%
#   pivot_longer(-iter, names_to="coef", values_to="draw") %>%
#   mutate(coef = gsub("b_","",coef)) %>%
#   group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
#   mutate(mod = "Monoculture Effects") %>%
#   mutate(Coefficient = gsub("protect_cells_log:protect","",coef)) %>%
#   mutate(Coefficient = str_replace_all(Coefficient,c(plant_age="Effect of plant age\nat inoculation",WCS365="Effect of WCS365\ncell density\nin monoculture", CHAO="Effect of CHA0\ncell density\nin monoculture", 
#                                               CH267="Effect of CH267\ncell density\nin monoculture", PF5="Effect of Pf5\ncell density\nin monoculture"))) %>%
#   mutate(Coefficient = factor(Coefficient, levels=c("Effect of plant age\nat inoculation"
#                                                     ,"Effect of WCS365\ncell density\nin monoculture"
#                                                     ,"Effect of CHA0\ncell density\nin monoculture"
#                                                     ,"Effect of CH267\ncell density\nin monoculture"
#                                                     ,"Effect of Pf5\ncell density\nin monoculture")))
# 
# gg_monoculture_forpriors <- draws_monoculture_bin_forpriors %>%
#   # filter(coef!="Intercept") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = FALSE) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~mod) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         , axis.title.x = element_blank())+
#   ylab("Estimate posterior")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# gg_monoculture_forpriors
# 
# # Look at estimates
# estimate_for_priors_mono <- brm_grad_monoculture_bin_wpriors$fit %>%
#   as.data.frame() %>%
#   select(starts_with("b_")) %>%
#   apply(MARGIN=2, FUN=function(x) c(m = mean(x), std = sd(x))) %>%
#   as_tibble() %>%
#   rename_all(~gsub("b_","",.)) %>% as.data.frame()
# estimate_for_priors_mono
# # estimate_for_priors[,"pathN2C3"][1]
# 
# 
# ### commensal vs path
# # brm_grad_allstrains_bin_wpriors <- brm(Alive ~  plant_age + path+ path_cells_log + protect_cells_log:protect + protect_cells_log:path_cells_log:protect+ (1 | plate) , data=dat_ac_nomock_nolacz
# brm_grad_allstrains_bin_wpriors <- brm(Alive ~  plant_age + path + path_cells_log + protect_cells_log:protect + protect_cells_log:path_cells_log:protect +  (1 | plate) , data=dat_ac_nomock_nolacz
#                                        , seed=12498215
#                                        , family="bernoulli"
#                                        , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors"
#                                        , prior = c(
#                                          # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                          set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                          set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                          set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                          # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                        ),iter=4000)
# 
# # save(brm_grad_allstrains_bin_wpriors, file = "04_ratio_experiment/WORKINGBRMS.RData")
# 
# # Plot mock estimates
# draws_allstrains_bin_wpriors  <- as_draws_df(brm_grad_allstrains_bin_wpriors) %>%
#   rownames_to_column(var="iter") %>%
#   select(iter, starts_with("b_")) %>%
#   pivot_longer(-iter, names_to="coef", values_to="draw") %>%
#   mutate(coef = gsub("b_","",coef)) %>%
#   group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
#   mutate(coef = gsub("protect_cells_log:protect",":protect_cells_log__",coef)) %>%
#   separate(coef, sep="__", into=c("coef","protect")) %>%
#   mutate(coef = str_replace_all(coef, c(`^path_cells_log$`="N2C3 cells",plant_age="Effect of plant age\nat inoculation",
#                                         `^:protect_cells_log$`="Protective cells", 
#                                         `path_cells_log::protect_cells_log`="    Interaction between\nN2C3 and protective"))) %>%
#   mutate(coef = factor(coef, levels=c("Effect of plant age\nat inoculation","N2C3 cells","Protective cells","    Interaction between\nN2C3 and protective"))) %>%
#   mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   mutate(protect = str_replace_all(protect, c(CHAO='CHA0', PF5='Pf5'))) %>%
#   mutate(protect = factor(protect, levels=c("Shared effects\nbetween all\nprotective strains","WCS365","CHA0","CH267","Pf5")))
# 
# # fixef(brm_grad_allstrains_bin_wpriors, probs = c(0.025, 0.05, 0.95, 0.975)) %>% View()
# # fixef(brm_grad_mockonly_bin_forpriors) %>% View()
# 
# gg_allstrains_wpriors <- draws_allstrains_bin_wpriors %>%
#   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   filter(coef!="Intercept") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~protect, scales="free_x") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         # , axis.title.x = element_blank()
#         , axis.title.y = element_blank()
#   )+
#   xlab("Coefficients")+
#   labs(col="Probability of\ndirection")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# gg_allstrains_wpriors
# 
# gg_bayesmodel_nolacz_bin <- plot_grid(gg_mockonly_forpriors, gg_allstrains_wpriors
#                                       , nrow=1, align="h", axis = "tb", rel_widths = c(1,4)) 
# gg_bayesmodel_nolacz_bin
# ggsave(filename = "04_ratio_experiment/gg_bayesmodel_nolacz_bin.png"
#        , gg_bayesmodel_nolacz_bin, height=5, width=10)
# 
# ###### VERSION TWO 
# dat_ac_comp_nolacz <- dat_ac_nomock_nolacz %>%
#   filter(path!="MOCK", protect!="MOCK") %>%
#   mutate(ratiofc = log10(protect_cells/path_cells))
# brm_grad_allstrains_bin_wpriors_RATIO <- brm(Alive ~  plant_age + path_cells_log +protect_cells_log:protect +  protect_cells_log:path_cells_log:protect +(1 | plate) , data=dat_ac_comp_nolacz
#                                        , seed=12498215
#                                        , family="bernoulli"
#                                        , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_RATIO"
#                                        , prior = c(
#                                          # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                          # set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                          set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                          set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                          set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                          # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                        ),iter=4000)
# 
# # Plot mock estimates
# draws_allstrains_bin_wpriors_RATIO  <- as_draws_df(brm_grad_allstrains_bin_wpriors_RATIO) %>%
#   rownames_to_column(var="iter") %>%
#   select(iter, starts_with("b_")) %>%
#   pivot_longer(-iter, names_to="coef", values_to="draw") %>%
#   mutate(coef = gsub("b_","",coef)) %>%
#   group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
#   mutate(coef = gsub("protect_cells_log:protect",":protect_cells_log__",coef)) %>%
#   separate(coef, sep="__", into=c("coef","protect")) %>%
#   mutate(coef = str_replace_all(coef, c(`^path_cells_log$`="N2C3 cells",plant_age="Effect of plant age\nat inoculation",
#                                         `^:protect_cells_log$`="Protective cells", 
#                                         `path_cells_log::protect_cells_log`="    Interaction between\nN2C3 and protective"))) %>%
#   mutate(coef = factor(coef, levels=c("Effect of plant age\nat inoculation","N2C3 cells","Protective cells","    Interaction between\nN2C3 and protective"))) %>%
#   mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   mutate(protect = str_replace_all(protect, c(CHAO='CHA0', PF5='Pf5'))) %>%
#   mutate(protect = factor(protect, levels=c("Shared effects\nbetween all\nprotective strains","WCS365","CHA0","CH267","Pf5")))
# 
# # fixef(brm_grad_allstrains_bin_wpriors, probs = c(0.025, 0.05, 0.95, 0.975)) %>% View()
# # fixef(brm_grad_mockonly_bin_forpriors) %>% View()
# 
# gg_allstrains_wpriors_RATIO <- draws_allstrains_bin_wpriors_RATIO %>%
#   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   filter(coef!="Intercept") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~protect, scales="free_x") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         # , axis.title.x = element_blank()
#         , axis.title.y = element_blank()
#   )+
#   xlab("Coefficients")+
#   labs(col="Probability of\ndirection")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# gg_allstrains_wpriors_RATIO
# 
# 
# ##### simplified brms coef plots ######
# 
# gg_mockonly_forpriors_simplified <- draws_mockonly_bin_forpriors %>%
#   filter(coef!="Intercept") %>%
#   filter(coef!="plant_age") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = FALSE) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~mod) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         , axis.title.x = element_blank())+
#   ylab("Estimate posterior")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
#   ylim(-4, (0.4/0.25)*4)
# gg_mockonly_forpriors_simplified
# 
# gg_allstrains_wpriors_simplified <- draws_allstrains_bin_wpriors %>%
#   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   filter(coef!="Intercept") %>%
#   filter(protect !="Shared effects\nbetween all\nprotective strains") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~protect, scales="free_x") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         # , axis.title.x = element_blank()
#         , axis.title.y = element_blank()
#   )+
#   xlab("Coefficients")+
#   labs(col="Probability of\ndirection")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
#   ylim(-0.25, 0.4)
# # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# gg_allstrains_wpriors_simplified
# 
# gg_bayesmodel_nolacz_bin_simp <- plot_grid(gg_mockonly_forpriors_simplified, gg_allstrains_wpriors_simplified
#                                       , nrow=1, align="h", axis = "tb", rel_widths = c(1,4)) 
# gg_bayesmodel_nolacz_bin_simp
# ggsave(filename = "04_ratio_experiment/gg_bayesmodel_nolacz_bin_simp.png"
#        , gg_bayesmodel_nolacz_bin_simp, height=4, width=8)
# 
# ### Ratio version
# gg_allstrains_wpriors_RATIO_simplified <- draws_allstrains_bin_wpriors_RATIO %>%
#   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
#   filter(coef!="Intercept") %>%
#   filter(protect !="Shared effects\nbetween all\nprotective strains") %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~protect, scales="free_x") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         # , axis.title.x = element_blank()
#         , axis.title.y = element_blank()
#   )+
#   xlab("Coefficients")+
#   labs(col="Probability of\ndirection")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
#   ylim(-0.25, 0.4)
# # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# gg_allstrains_wpriors_RATIO_simplified
# 
# gg_bayesmodel_nolacz_bin_RATIO_simp <- plot_grid(gg_mockonly_forpriors_simplified, gg_allstrains_wpriors_RATIO_simplified
#                                            , nrow=1, align="h", axis = "tb", rel_widths = c(1,4)) 
# gg_bayesmodel_nolacz_bin_RATIO_simp
# ggsave(filename = "04_ratio_experiment/gg_bayesmodel_nolacz_bin_RATIO_simp.png"
#        , gg_bayesmodel_nolacz_bin_RATIO_simp, height=4, width=8)
# 
# ### now for plant age only
# 
# gg_plantage_simplified  <- draws_mockonly_bin_forpriors %>%
#   filter(coef == "plant_age") %>%
#   full_join(draws_allstrains_bin_wpriors %>% filter(coef == "Effect of plant age\nat inoculation") %>%
#               rename(Coefficient = coef)) %>%
#   mutate(protect = ifelse(is.na(protect), "Pathogen only", "Protective and \npathogen together\n(across all strains)")) %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=protect, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~"Pooled effect\nacross all treatments") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         )+
#   ylab("Estimate posterior")+xlab("Coefficient")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# gg_plantage_simplified
# ggsave(filename = "04_ratio_experiment/gg_plantage_simplified.png", gg_plantage_simplified, height=3, width=4)
# 
# 
# 
# gg_plantage_RATIO_simplified  <- draws_mockonly_bin_forpriors %>%
#   filter(coef == "plant_age") %>%
#   full_join(draws_allstrains_bin_wpriors_RATIO %>% filter(coef == "Effect of plant age\nat inoculation") %>%
#               rename(Coefficient = coef)) %>%
#   mutate(protect = ifelse(is.na(protect), "Pathogen only", "Protective and \npathogen together\n(across all strains)")) %>%
#   mutate(SigEffect = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
#   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
#                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
#   ggplot() +
#   geom_pointrange(aes(x=protect, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_grid(.~"Pooled effect\nacross all treatments") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#   )+
#   ylab("Estimate posterior")+xlab("Coefficient")+
#   labs(col="Two-sided\nPD < 0.05")+
#   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
# gg_plantage_RATIO_simplified
# ggsave(filename = "04_ratio_experiment/gg_plantage_RATIO_simplified.png", gg_plantage_RATIO_simplified, height=3, width=4)
# 
# # 
# # gg_allstrains_wpriors_simplified <- draws_allstrains_bin_wpriors %>%
# #   # mutate(protect=ifelse(is.na(protect), "Shared effects\nbetween all\nprotective strains",protect)) %>%
# #   filter(coef!="Intercept") %>%
# #   filter(protect !="Shared effects\nbetween all\nprotective strains") %>%
# #   # filter(protect!="Pathogen only") %>%
# #   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
# #   mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
# #   mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
# #                                ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
# #   ggplot() +
# #   geom_pointrange(aes(x=coef, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll)) +
# #   geom_hline(aes(yintercept=0)) +
# #   facet_grid(.~protect, scales="free_x") +
# #   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
# #         # , axis.title.x = element_blank()
# #         , axis.title.y = element_blank()
# #   )+
# #   xlab("Coefficients")+
# #   labs(col="Probability of\ndirection")+
# #   scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
# #   ylim(-0.25, 0.4)
# # # scale_color_manual(values=c(`FALSE`="black",`TRUE`="red"))
# # gg_allstrains_wpriors_simplified
# # 
# # gg_bayesmodel_nolacz_bin_simp <- plot_grid(gg_mockonly_forpriors_simplified, gg_allstrains_wpriors_simplified
# #                                            , nrow=1, align="h", axis = "tb", rel_widths = c(1,4)) 
# # gg_bayesmodel_nolacz_bin_simp
# # ggsave(filename = "04_ratio_experiment/gg_bayesmodel_nolacz_bin_simp.png"
# #        , gg_bayesmodel_nolacz_bin_simp, height=5, width=10)
# 
# ####### LOO compare with model without interaction ####
# brm_grad_allstrains_bin_wpriors_NOINTER <- brm(Alive ~  plant_age + path + path_cells_log + protect_cells_log:protect +  (1 | plate) , data=dat_ac_nomock_nolacz
#                                                , seed=12498215
#                                                , family="bernoulli"
#                                                , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_NOINTER"
#                                                , prior = c(
#                                                  # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                                  set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                                  set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                                  set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                                  # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                                ),iter=4000)
# brm_grad_allstrains_bin_wpriors_NOINTER
# 
# brm_grad_allstrains_bin_wpriors_ONLYINTER <- brm(Alive ~  plant_age + path + protect_cells_log:path_cells_log:protect +  (1 | plate) , data=dat_ac_nomock_nolacz
#                                                  , seed=12498215
#                                                  , family="bernoulli"
#                                                  , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_ONLYINTER"
#                                                  , prior = c(
#                                                    # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                                    set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                                    set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                                    # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                                  ),iter=4000)
# brm_grad_allstrains_bin_wpriors_ONLYINTER
# loo_of_cellvsratio <- loo(brm_grad_allstrains_bin_wpriors, brm_grad_allstrains_bin_wpriors_NOINTER, brm_grad_allstrains_bin_wpriors_ONLYINTER)
# loo_of_cellvsratio
# sink("loo_of_cellvsinteraction.txt")
# loo_of_cellvsratio
# sink()
# 
# ## ratio verison
# brm_grad_allstrains_bin_wpriors_NOINTER_RATIO <- brm(Alive ~  plant_age + path_cells_log + protect_cells_log:protect +  (1 | plate) , data=dat_ac_comp_nolacz
#                                                , seed=12498215
#                                                , family="bernoulli"
#                                                , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_NOINTER_RATIO"
#                                                , prior = c(
#                                                  # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                                  # set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                                  set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                                  set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                                  set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                                  # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                                ),iter=4000)
# brm_grad_allstrains_bin_wpriors_NOINTER_RATIO
# 
# brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO <- brm(Alive ~  plant_age + protect_cells_log:path_cells_log:protect +  (1 | plate) , data=dat_ac_comp_nolacz
#                                                  , seed=12498215
#                                                  , family="bernoulli"
#                                                  , file="04_ratio_experiment/brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO"
#                                                  , prior = c(
#                                                    # set_prior("normal(0,2.5)", class="b", coef="protect_cells_log:protectWCS365"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors[,"pathN2C3"][1],",",estimate_for_priors[,"pathN2C3"][2],")"), class="b", coef="pathN2C3"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectWCS365"][2],")"), class="b", coef="protect_cells_log:protectWCS365"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCHAO"][2],")"), class="b", coef="protect_cells_log:protectCHAO"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectCH267"][2],")"), class="b", coef="protect_cells_log:protectCH267"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][1],",",estimate_for_priors_mono[,"protect_cells_log:protectPF5"][2],")"), class="b", coef="protect_cells_log:protectPF5"),
#                                                    # set_prior(paste0("normal(",estimate_for_priors[,"path_cells_log"][1],",",estimate_for_priors[,"path_cells_log"][2],")"), class="b", coef="path_cells_log"),
#                                                    set_prior(paste0("normal(",estimate_for_priors[,"plant_age"][1],",",estimate_for_priors[,"plant_age"][2],")"), class="b", coef="plant_age")
#                                                    # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
#                                                  ),iter=4000)
# brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO
# loo_of_cellvsratio_RATIO <- loo(brm_grad_allstrains_bin_wpriors_RATIO, brm_grad_allstrains_bin_wpriors_NOINTER_RATIO, brm_grad_allstrains_bin_wpriors_ONLYINTER_RATIO)
# loo_of_cellvsratio_RATIO
# sink("loo_of_cellvsinteraction_RATIO.txt")
# loo_of_cellvsratio_RATIO
# sink()
# 
# ##### Putting bayes model onto plots ######
# 
# ## Binomial 
# # Get predictions
# newdat_x_mock <- dat_ac_onlymock_nolacZ %>% 
#   select(protect, plant_age, path_cells_log,path,  protect_cells_log, path, Healthiness_hsv, Alive) %>%
#   mutate(path=as.character(path)) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate")
# newdat_x_nomock <- dat_ac_nomock_nolacz %>% 
#   select(protect, plant_age, path_cells_log,  protect_cells_log, path, Healthiness_hsv, Alive) %>%
#   mutate(path=as.character(path)) %>%
#   mutate(protect = as.factor(as.vector(protect))) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate")
# newdat_x_monoculture <- dat_ac_nomock_nolacz %>% 
#   select(protect, plant_age, protect_cells_log, path, Healthiness_hsv, Alive) %>%
#   mutate(protect = as.factor(as.vector(protect))) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate")
# newdat_x_comp <- dat_ac_nomock_nolacz %>% 
#   select(protect, plant_age, path_cells_log, protect_cells_log, Healthiness_hsv, Alive) %>%
#   mutate(protect = as.factor(as.vector(protect))) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate")
# 
# mean_predict_onlymock_nolacz <- apply(posterior_linpred(brm_grad_mockonly_bin_forpriors, newdata = newdat_x_mock, allow_new_levels=TRUE, transform=TRUE), MARGIN = 2, mean)
# mean_predict_nomock_nolacz <- apply(posterior_linpred(brm_grad_allstrains_bin_wpriors, newdata = newdat_x_nomock, allow_new_levels=TRUE, transform=TRUE), MARGIN = 2, mean)
# mean_predict_monoculture_nolacz <- apply(posterior_linpred(brm_grad_monoculture_bin_wpriors, newdata = newdat_x_monoculture, allow_new_levels=TRUE, transform=TRUE), MARGIN = 2, mean)
# mean_predict_comp_nolacz <- apply(posterior_linpred(brm_grad_allstrains_bin_wpriors_RATIO, newdata = newdat_x_comp, allow_new_levels=TRUE, transform=TRUE), MARGIN = 2, mean)
# newdat_x_mock$prediction <- mean_predict_onlymock_nolacz
# newdat_x_nomock$prediction <- mean_predict_nomock_nolacz
# newdat_x_monoculture$prediction <- mean_predict_monoculture_nolacz
# newdat_x_comp$prediction <- mean_predict_comp_nolacz
# predictions_nolacZ <- full_join(newdat_x_mock, newdat_x_nomock)
# predictions_nolacZ_RATIO <- full_join(newdat_x_mock, newdat_x_monoculture ) %>%
#   full_join(newdat_x_comp)
# 
# allPredictions_justmocks <- predictions_nolacZ %>%
#   filter(path=="MOCK" | protect=="MOCK") %>%
#   group_by(protect, path, plant_age) %>%
#   mutate(prediction = median(prediction)) %>%
#   ungroup()
# allPredictions_all <- predictions_nolacZ %>%
#   filter(!(path=="MOCK" | protect=="MOCK")) %>%
#   group_by(protect, path, protect_cells_log, path_cells_log, plant_age) %>%
#   summarise(prediction = mean(prediction)) %>%
#   bind_rows(allPredictions_justmocks) %>%
#   ungroup() %>%
#   filter(protect!="MOCK")
# 
# allPredictions_mockmock_nolacZ <- predictions_nolacZ %>%
#   filter(protect=="MOCK") %>%
#   mutate(protect="No protective") %>%
#   select(protect, path, path_cells_log, protect_cells_log, prediction, plant_age) %>%
#   distinct()
# 
# gg_allprotect_bin_nolacZ_fit <- dat_ac %>%
#   filter(protect!="MOCK") %>%
#   ggplot() +
#   geom_jitter(aes(x=round(protect_cells_log), y=Alive, col=as.factor(plant_age)), height=0.1, width=0.1,show.legend = FALSE) +
#   geom_line(data=allPredictions_all, aes(x=round(protect_cells_log), y=prediction, col=as.factor(plant_age),group=as.factor(plant_age))) +
#   scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         , axis.title.y = element_blank()
#         , axis.ticks.y = element_line(colour="white")
#         , axis.text.y = element_text(colour="white"))+
#   facet_grid(round(path_cells_log)~protect, scales="free", drop = TRUE, space="free"
#              , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
#                                    , protect = protect_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1))+
#   xlab("Protective cells (log10(+1))") + labs(col="Plant age at\ninoculation")
# gg_allprotect_bin_nolacZ_fit
# 
# gg_mocksonly_bin_nolacZ_fit <- dat_ac %>%
#   filter(protect=="MOCK") %>%
#   mutate(protect="No protective") %>%
#   ggplot() +
#   geom_jitter(aes(x=round(protect_cells_log), y=Alive, col=as.factor(plant_age)), height=0.1, width=0.1,show.legend = FALSE) +
#   geom_hline(data=allPredictions_mockmock_nolacZ, aes(yintercept=prediction, col=as.factor(plant_age)), show.legend = FALSE) +
#   scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_blank()
#         , strip.text.y.right = element_blank()
#         , axis.ticks.x = element_blank())+
#   facet_grid(round(path_cells_log)~protect, scales="free"
#              , drop = TRUE, space="free"
#              , labeller = labeller(plant_age = plant_age_labeller
#                                    , `round(path_cells_log)`=N2C3_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1)) +
#   scale_y_continuous(breaks=c(0,0.5, 1)) +
#   ylab("Probability of healthy plant") + xlab("") 
# gg_mocksonly_bin_nolacZ_fit
# 
# gg_col0_collapsedage_bayes_nolacZ_bin <- plot_grid(gg_mocksonly_bin_nolacZ_fit, gg_allprotect_bin_nolacZ_fit, nrow=1, rel_widths = c(1,4), align = "h")
# gg_col0_collapsedage_bayes_nolacZ_bin
# ggsave("04_ratio_experiment/gg_col0_collapsedage_bayes_nolacZ_bin.png", gg_col0_collapsedage_bayes_nolacZ_bin, height=6, width=9)
# 
# ##### Simplified brms plots ####
# range_pathcellslog <- seq(range(dat_ac_onlymock_nolacZ %>% pull(path_cells_log))[1], range(dat_ac_onlymock_nolacZ$path_cells_log)[2], length.out=100)
# newdat_expanded_mock <- dat_ac_onlymock_nolacZ %>%
#   # filter(path=="N2C3") %>%
#   select(protect, plant_age,path,  protect_cells_log, path) %>%
#   mutate(path=as.character(path)) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate") %>%
#   full_join(data.frame(path="N2C3", path_cells_log=range_pathcellslog), relationship="many-to-many") %>%
#   mutate(path=ifelse(path_cells_log == 0, "MOCK", path)) %>%
#   drop_na()
# 
# newdat_expanded_nomock <- dat_ac_nomock_nolacz%>%
#   filter(protect%in% c("WCS365","CHAO","CH267","PF5")) %>%
#   mutate(protect = factor(protect, levels=c("WCS365","CHAO","CH267","PF5"))) %>%
#   select(protect, plant_age,path,  protect_cells_log, path) %>%
#   mutate(path=as.character(path)) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate") %>%
#   full_join(data.frame(path="N2C3", path_cells_log=range_pathcellslog), relationship="many-to-many") %>%
#   mutate(path=ifelse(path_cells_log == 0, "MOCK", path)) %>%
#   drop_na()
# 
# # Version where we have 95% prediction intervals
# mean_pred95_predict_onlymock_cont <- t(apply(posterior_linpred(brm_grad_mockonly_bin_forpriors, newdata = newdat_expanded_mock, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
# mean_pred95_predict_nomock_cont <- t(apply(posterior_linpred(brm_grad_allstrains_bin_wpriors, newdata = newdat_expanded_nomock, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
# newdat_expanded_mock <- bind_cols(newdat_expanded_mock,mean_pred95_predict_onlymock_cont)
# newdat_expanded_nomock <- bind_cols(newdat_expanded_nomock,mean_pred95_predict_nomock_cont)
# predictions_cont <- full_join(newdat_expanded_mock, newdat_expanded_nomock)
# 
# gg_allprotect_bin_nolacZ_fit_simp <- dat_ac %>%
#   # filter(protect!="MOCK") %>%
#   filter(plant_age==5) %>%
#   ggplot() +
#   geom_point(aes(x=round(path_cells_log), y=Alive, group=(round(protect_cells_log)) , col=(round(protect_cells_log))), show.legend = FALSE, position=position_jitterdodgev(jitter.height=0, jitter.width=0.5, dodge.height=0.3)) +
#   geom_smooth(data=predictions_cont %>% filter(plant_age==5), aes(x=round(path_cells_log), y=mean, col=(round(protect_cells_log)), group=factor(round(protect_cells_log))), se=FALSE) +
#   # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#         )+
#   facet_grid(.~protect, scales="free", drop = TRUE, space="free"
#              , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
#                                    , protect = protect_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1)) +
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
#   ylab("Probability of healthy plant") +
#   xlab("N2C3 (pathogen) cells (log10(+1))") + labs(col="Protective cells\n(log10(+1))") +
#   scale_color_gradient(low="gold", high="darkgreen")
#   # scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
# gg_allprotect_bin_nolacZ_fit_simp
# ggsave("04_ratio_experiment/gg_allprotect_bin_nolacZ_fit_simp.png", gg_allprotect_bin_nolacZ_fit_simp, height=3, width=8)
# 
# 
# ### Ratio version of simplified prediction plot ####
# range_pathcellslog <- seq(range(dat_ac_onlymock_nolacZ %>% pull(path_cells_log))[1], range(dat_ac_onlymock_nolacZ$path_cells_log)[2], length.out=100)
# newdat_expanded_mock_ratio <- dat_ac_onlymock_nolacZ %>%
#   # filter(path=="N2C3") %>%
#   select(protect, plant_age,path,  protect_cells_log, path) %>%
#   mutate(path=as.character(path)) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate") %>%
#   full_join(data.frame(path="N2C3", path_cells_log=range_pathcellslog), relationship="many-to-many") %>%
#   mutate(path=ifelse(path_cells_log == 0, "MOCK", path)) %>%
#   drop_na()
# # commensal
# newdat_expanded_monoculture<-  newdat_x_monoculture %>%
#   select(plant_age, protect, protect_cells_log, plate) %>% distinct()
# range_pathcellslog2 <- seq(1, range(dat_ac_onlymock_nolacZ$path_cells_log)[2], length.out=100)
# newdat_expanded_comp <- dat_ac_comp_nolacz%>%
#   filter(protect%in% c("WCS365","CHAO","CH267","PF5")) %>%
#   mutate(protect = factor(protect, levels=c("WCS365","CHAO","CH267","PF5"))) %>%
#   select(protect, plant_age,path,  protect_cells_log, path) %>%
#   mutate(path=as.character(path)) %>%
#   mutate(plant_age = as.numeric(plant_age)) %>%
#   distinct() %>%
#   mutate(plate="newplate") %>%
#   full_join(data.frame(path="N2C3", path_cells_log=range_pathcellslog2), relationship="many-to-many") %>%
#   mutate(path=ifelse(path_cells_log == 0, "MOCK", path)) %>%
#   filter(path!="MOCK") %>%
#   drop_na() 
# 
# # Version where we have 95% prediction intervals
# mean_pred95_predict_onlymock_cont <- t(apply(posterior_linpred(brm_grad_mockonly_bin_forpriors, newdata = newdat_expanded_mock_ratio, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
# mean_pred95_predict_monoculture_cont <- t(apply(posterior_linpred(brm_grad_monoculture_bin_wpriors, newdata = newdat_expanded_monoculture, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
# mean_pred95_predict_comp_cont <- t(apply(posterior_linpred(brm_grad_allstrains_bin_wpriors_RATIO, newdata = newdat_expanded_comp, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
# newdat_expanded_mock_ratio <- bind_cols(newdat_expanded_mock_ratio,mean_pred95_predict_onlymock_cont)
# newdat_expanded_monoculture <- bind_cols(newdat_expanded_monoculture,mean_pred95_predict_monoculture_cont)
# newdat_expanded_comp <- bind_cols(newdat_expanded_comp,mean_pred95_predict_comp_cont)
# predictions_cont_RATIO <- full_join(newdat_expanded_mock_ratio, newdat_expanded_comp) 
# predictions_mocks <- newdat_expanded_monoculture 
# 
# gg_allprotect_bin_nolacZ_fit_RATIO_simp <- dat_ac %>%
#   # filter(protect!="MOCK") %>%
#   filter(plant_age==5) %>%
#   ggplot() +
#   geom_point(aes(x=(path_cells_log), y=Alive, group=(round(protect_cells_log)) , col=(round(protect_cells_log))), show.legend = FALSE, position=position_jitterdodgev(jitter.height=0, jitter.width=0.5, dodge.height=0.3)) +
#   geom_line(data=predictions_cont_RATIO %>% filter(plant_age==5, path=="N2C3"), aes(x=(path_cells_log), y=mean, col=(round(protect_cells_log)), group=factor(round(protect_cells_log)))) +
#   geom_hline(data=predictions_mocks %>% filter(plant_age==5), aes(yintercept=mean, col=(round(protect_cells_log)), group=factor(round(protect_cells_log)) )) +
#   # geom_smooth(data=predictions_cont_RATIO %>% filter(plant_age==5), aes(x=round(path_cells_log), y=mean, col=(round(protect_cells_log)), group=factor(round(protect_cells_log))), se=FALSE) +
#   # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#   )+
#   facet_grid(.~protect, scales="free", drop = TRUE, space="free"
#              , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
#                                    , protect = protect_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1)) +
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
#   ylab("Probability of healthy plant") +
#   xlab("N2C3 (pathogen) cells (log10(+1))") + labs(col="Protective cells\n(log10(+1))") +
#   scale_color_gradient(low="gold", high="darkgreen")
# # scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
# gg_allprotect_bin_nolacZ_fit_RATIO_simp
# ggsave("04_ratio_experiment/gg_allprotect_bin_nolacZ_fit_RATIO_simp.png", gg_allprotect_bin_nolacZ_fit_RATIO_simp, height=3, width=8)
# 
# ## Get age plot
# gg_allprotect_plantage_simp <- dat_ac %>%
#   # filter(protect!="MOCK") %>%
#   filter((protect_cells_log>4& protect_cells_log<5) | protect_cells_log==0) %>%
#   ggplot() +
#   geom_point(aes(x=round(path_cells_log), y=Alive, group=as.character(plant_age) , col=as.character(plant_age)), show.legend = FALSE, position=position_jitterdodgev(jitter.height=0, jitter.width=0.5, dodge.height=0.3)) +
#   geom_smooth(data=predictions_cont %>% filter((protect_cells_log>4& protect_cells_log<5) | protect_cells_log==0), aes(x=round(path_cells_log), y=mean, col=as.character(plant_age), group=as.character(plant_age)), se=FALSE) +
#   # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#   )+
#   facet_grid(.~protect, scales="free", drop = TRUE, space="free"
#              , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
#                                    , protect = protect_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1)) +
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
#   ylab("Probability of healthy plant") +
#   xlab("N2C3 (pathogen) cells (log10(+1))") + labs(col="Plant age (days)\nat inoculation") +
#   scale_color_manual(values=c("lightblue","blue","darkblue"))
# # scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
# gg_allprotect_plantage_simp
# ggsave("04_ratio_experiment/gg_allprotect_plantage_simp.png", gg_allprotect_plantage_simp, height=3, width=8)
# 
# gg_allprotect_plantage_simp_RATIO <- dat_ac %>%
#   # filter(protect!="MOCK") %>%
#   filter((protect_cells_log>4& protect_cells_log<5) | protect_cells_log==0) %>%
#   ggplot() +
#   geom_point(aes(x=round(path_cells_log), y=Alive, group=as.character(plant_age) , col=as.character(plant_age)), show.legend = FALSE, position=position_jitterdodgev(jitter.height=0, jitter.width=0.5, dodge.height=0.3)) +
#   geom_smooth(data=predictions_cont_RATIO %>% filter((protect_cells_log>4& protect_cells_log<5) | protect_cells_log==0), aes(x=round(path_cells_log), y=mean, col=as.character(plant_age), group=as.character(plant_age)), se=FALSE) +
#   # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#   )+
#   facet_grid(.~protect, scales="free", drop = TRUE, space="free"
#              , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
#                                    , protect = protect_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1)) +
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
#   ylab("Probability of healthy plant") +
#   xlab("N2C3 (pathogen) cells (log10(+1))") + labs(col="Plant age (days)\nat inoculation") +
#   scale_color_manual(values=c("lightblue","blue","darkblue"))
# # scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
# gg_allprotect_plantage_simp_RATIO
# ggsave("04_ratio_experiment/gg_allprotect_plantage_simp_RATIO.png", gg_allprotect_plantage_simp_RATIO, height=3, width=8)
# 

# #### Flipped
# gg_allprotect_bin_nolacZ_fit_flipped <- dat_ac %>%
#   filter(lacZ_path=="WT") %>%
#   filter(protect!="MOCK") %>%
#   ggplot() +
#   geom_jitter(aes(x=round(protect_cells_log), y=Alive, col=as.factor(plant_age)), height=0.1, width=0.1,show.legend = FALSE) +
#   # geom_line(data=allPredictions_avelacz, aes(x=round(protect_cells_log), y=prediction_bin2, col=as.factor(plant_age),group=as.factor(plant_age))) +
#   geom_smooth(data=allPredictions_all, aes(x=round(protect_cells_log), y=prediction, col=as.factor(plant_age),group=as.factor(plant_age)), method="glm", method.args=c(family="binomial"), se=FALSE) +
#   scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
#   )+
#   facet_grid(protect~round(path_cells_log), scales="free", drop = TRUE, space="free", labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller))+
#   scale_x_continuous(breaks=seq(0,8,1))+
#   scale_y_continuous(breaks=c(0,0.5, 1)) +
#   xlab("Protective cells (log10(+1))") + labs(col="Plant age at\ninoculation")+ylab("Probability of healthy plant\n(Health score > 400)")
# gg_allprotect_bin_nolacZ_fit_flipped
# ggsave("04_ratio_experiment/gg_allprotect_bin_nolacZ_fit_flipped.png", gg_allprotect_bin_nolacZ_fit_flipped, height=6, width=8)

#### Brms with FULL model, priors ####
###### Monoculture for priors ######

dat_ac_monocultures_nolacZ <- dat_ac %>%
  filter(protect =="MOCK" | path == "MOCK") %>%
  # filter(!(protect=="MOCK"&path=="MOCK")) %>%
  mutate(Strain = ifelse(protect == "MOCK", as.character(path), as.character(protect))) %>%
  select(Alive, Healthiness_hsv, Strain, protect, path, plant_age, total_cells_log, plate, experiment, path_cells_log, protect_cells_log) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","N2C3","WCS365","CHAO","CH267","PF5"))) %>%
  mutate(plant_age_adj = plant_age-5)
# mutate(Strain = factor(Strain, levels=c("N2C3","WCS365","CHAO","CH267","PF5"))) 

# Get baselines for monocultures
brm_monocultures <- brm(Alive ~ Strain + Strain:total_cells_log + plant_age + (1 | plate), data=dat_ac_monocultures_nolacZ
                        , seed=12498215
                        , family="bernoulli"
                        , prior = c(set_prior("constant(0)", class = "b", coef = "StrainMOCK:total_cells_log"))
                        , file="04_ratio_experiment/brm_monocultures"
                        , iter=4000)
brm_monocultures

# Plot mock estimates
draws_monoculture  <- as_draws_df(brm_monocultures) %>%
  rownames_to_column(var="iter") %>%
  select(iter, starts_with("b_")) %>%
  pivot_longer(-iter, names_to="coef", values_to="draw") %>%
  mutate(coef = gsub("b_","",coef)) %>%
  group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
  mutate(mod = "Monoculture") %>%
  mutate(Coefficient = str_replace_all(coef,c(plant_age="MOCK:Effect of plant age\nat inoculation",Intercept="MOCK:Presence of strain"))) %>%
  separate(Coefficient, sep=":", into=c("Strain","Coefficient"), fill = "right") %>%
  mutate(Strain = gsub("Strain","",Strain)) %>%
  mutate(Coefficient = ifelse(is.na(Coefficient), "Presence of strain", ifelse(
    Coefficient=="total_cells_log", "Cell load", Coefficient))) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","N2C3","WCS365","CHAO","CH267","PF5"))) %>%
  mutate(Coefficient = factor(Coefficient, levels=c("Effect of plant age\nat inoculation","Presence of strain","Cell load")))

gg_monoculture_forpriors <- draws_monoculture %>%
  filter(Strain !="MOCK") %>%
  # filter(!(Strain=="MOCK" & Coefficient == "Cell density")) %>%
  # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~Strain) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Estimate posterior")+ 
  labs(col="Two-sided\nPD < 0.05")+
  scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
gg_monoculture_forpriors
ggsave("04_ratio_experiment/gg_monoculture_forpriors.png", gg_monoculture_forpriors, width=8, height=4)

# Look at estimates
mpriors <- brm_monocultures$fit %>%
  as.data.frame() %>%
  select(starts_with("b_")) %>%
  apply(MARGIN=2, FUN=function(x) c(m = mean(x), std = sd(x))) %>%
  as_tibble() %>%
  rename_all(~gsub("b_","",.)) %>% as.data.frame()
mpriors


###### Full model, with priors ######
brm_full_sdpriors <- brm(Alive ~  plant_age + protect:protect_cells_log + path:path_cells_log + protect_cells_log:path_cells_log:protect  +  (1 | plate) , data=dat_ac
               , seed=12498215
               , family="bernoulli"
               , file="04_ratio_experiment/brm_full_sdpriors"
               , prior = c(
                 set_prior(paste0("normal(",mpriors[,"Intercept"][1],",",mpriors[,"Intercept"][2],")"), class="Intercept"),
                 # set_prior(paste0("normal(",mpriors[,"StrainN2C3"][1],",",mpriors[,"StrainN2C3"][2],")"), class="b", coef="pathN2C3"),
                 # set_prior(paste0("normal(",mpriors[,"StrainWCS365"][1],",",mpriors[,"StrainWCS365"][2],")"), class="b", coef="protectWCS365"),
                 # set_prior(paste0("normal(",mpriors[,"StrainCHAO"][1],",",mpriors[,"StrainCHAO"][2],")"), class="b", coef="protectCHAO"),
                 # set_prior(paste0("normal(",mpriors[,"StrainCH267"][1],",",mpriors[,"StrainCH267"][2],")"), class="b", coef="protectCH267"),
                 # set_prior(paste0("normal(",mpriors[,"StrainPF5"][1],",",mpriors[,"StrainPF5"][2],")"), class="b", coef="protectPF5"),
                 set_prior("constant(0)", class="b", coef="pathMOCK:path_cells_log"),
                 set_prior("constant(0)", class="b", coef="protectMOCK:protect_cells_log"),
                 set_prior("constant(0)", class="b", coef="protectMOCK:protect_cells_log:path_cells_log"),
                 set_prior(paste0("normal(",mpriors[,"StrainWCS365:total_cells_log"][1],",",mpriors[,"StrainWCS365:total_cells_log"][2],")"), class="b", coef="protectWCS365:protect_cells_log"),
                 set_prior(paste0("normal(",mpriors[,"StrainCHAO:total_cells_log"][1],",",mpriors[,"StrainCHAO:total_cells_log"][2],")"), class="b", coef="protectCHAO:protect_cells_log"),
                 set_prior(paste0("normal(",mpriors[,"StrainCH267:total_cells_log"][1],",",mpriors[,"StrainCH267:total_cells_log"][2],")"), class="b", coef="protectCH267:protect_cells_log"),
                 set_prior(paste0("normal(",mpriors[,"StrainPF5:total_cells_log"][1],",",mpriors[,"StrainPF5:total_cells_log"][2],")"), class="b", coef="protectPF5:protect_cells_log"),
                 set_prior(paste0("normal(",mpriors[,"StrainN2C3:total_cells_log"][1],",",mpriors[,"StrainN2C3:total_cells_log"][2],")"), class="b", coef="pathN2C3:path_cells_log"),
                 set_prior(paste0("normal(",mpriors[,"plant_age"][1],",",mpriors[,"plant_age"][2],")"), class="b", coef="plant_age")
                 # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
               ),iter=4000)
brm_full_sdpriors

###### LOO with full model ######

brm_nointer_sdpriors <- brm(Alive ~  plant_age + protect:protect_cells_log + path:path_cells_log  + (1 | plate) , data=dat_ac
                         , seed=12498215
                         , family="bernoulli"
                         , file="04_ratio_experiment/brm_nointer_sdpriors"
                         , prior = c(
                           set_prior(paste0("normal(",mpriors[,"Intercept"][1],",",mpriors[,"Intercept"][2],")"), class="Intercept"),
                           # set_prior(paste0("normal(",mpriors[,"StrainN2C3"][1],",",mpriors[,"StrainN2C3"][2],")"), class="b", coef="pathN2C3"),
                           # set_prior(paste0("normal(",mpriors[,"StrainWCS365"][1],",",mpriors[,"StrainWCS365"][2],")"), class="b", coef="protectWCS365"),
                           # set_prior(paste0("normal(",mpriors[,"StrainCHAO"][1],",",mpriors[,"StrainCHAO"][2],")"), class="b", coef="protectCHAO"),
                           # set_prior(paste0("normal(",mpriors[,"StrainCH267"][1],",",mpriors[,"StrainCH267"][2],")"), class="b", coef="protectCH267"),
                           # set_prior(paste0("normal(",mpriors[,"StrainPF5"][1],",",mpriors[,"StrainPF5"][2],")"), class="b", coef="protectPF5"),
                           set_prior("constant(0)", class="b", coef="pathMOCK:path_cells_log"),
                           set_prior("constant(0)", class="b", coef="protectMOCK:protect_cells_log"),
                           # set_prior("constant(0)", class="b", coef="protectMOCK:path_cells_log:protect_cells_log"),
                           set_prior(paste0("normal(",mpriors[,"StrainWCS365:total_cells_log"][1],",",mpriors[,"StrainWCS365:total_cells_log"][2],")"), class="b", coef="protectWCS365:protect_cells_log"),
                           set_prior(paste0("normal(",mpriors[,"StrainCHAO:total_cells_log"][1],",",mpriors[,"StrainCHAO:total_cells_log"][2],")"), class="b", coef="protectCHAO:protect_cells_log"),
                           set_prior(paste0("normal(",mpriors[,"StrainCH267:total_cells_log"][1],",",mpriors[,"StrainCH267:total_cells_log"][2],")"), class="b", coef="protectCH267:protect_cells_log"),
                           set_prior(paste0("normal(",mpriors[,"StrainPF5:total_cells_log"][1],",",mpriors[,"StrainPF5:total_cells_log"][2],")"), class="b", coef="protectPF5:protect_cells_log"),
                           set_prior(paste0("normal(",mpriors[,"StrainN2C3:total_cells_log"][1],",",mpriors[,"StrainN2C3:total_cells_log"][2],")"), class="b", coef="pathN2C3:path_cells_log"),
                           set_prior(paste0("normal(",mpriors[,"plant_age"][1],",",mpriors[,"plant_age"][2],")"), class="b", coef="plant_age")
                           # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
                         ),iter=4000)

brm_onlyinter_sdpriors <- brm(Alive ~  plant_age + path:path_cells_log + protect:protect_cells_log:path_cells_log +(1 | plate) , data=dat_ac
                            , seed=12498215
                            , family="bernoulli"
                            , file="04_ratio_experiment/brm_onlyinter_sdpriors"
                            , prior = c(
                              set_prior(paste0("normal(",mpriors[,"Intercept"][1],",",mpriors[,"Intercept"][2],")"), class="Intercept"),
                              # set_prior(paste0("normal(",mpriors[,"StrainN2C3"][1],",",mpriors[,"StrainN2C3"][2],")"), class="b", coef="pathN2C3"),
                              # set_prior(paste0("normal(",mpriors[,"StrainWCS365"][1],",",mpriors[,"StrainWCS365"][2],")"), class="b", coef="protectWCS365"),
                              # set_prior(paste0("normal(",mpriors[,"StrainCHAO"][1],",",mpriors[,"StrainCHAO"][2],")"), class="b", coef="protectCHAO"),
                              # set_prior(paste0("normal(",mpriors[,"StrainCH267"][1],",",mpriors[,"StrainCH267"][2],")"), class="b", coef="protectCH267"),
                              # set_prior(paste0("normal(",mpriors[,"StrainPF5"][1],",",mpriors[,"StrainPF5"][2],")"), class="b", coef="protectPF5"),
                              set_prior("constant(0)", class="b", coef="pathMOCK:path_cells_log"),
                              # set_prior("constant(0)", class="b", coef="protectMOCK:protect_cells_log"),
                              # set_prior("constant(0)", class="b", coef="protectMOCK:path_cells_log:protect_cells_log"),
                              set_prior("constant(0)", class="b", coef="path_cells_log:protectMOCK:protect_cells_log"),
                              # set_prior(paste0("normal(",mpriors[,"StrainWCS365:total_cells_log"][1],",",mpriors[,"StrainWCS365:total_cells_log"][2],")"), class="b", coef="protectWCS365:protect_cells_log"),
                              # set_prior(paste0("normal(",mpriors[,"StrainCHAO:total_cells_log"][1],",",mpriors[,"StrainCHAO:total_cells_log"][2],")"), class="b", coef="protectCHAO:protect_cells_log"),
                              # set_prior(paste0("normal(",mpriors[,"StrainCH267:total_cells_log"][1],",",mpriors[,"StrainCH267:total_cells_log"][2],")"), class="b", coef="protectCH267:protect_cells_log"),
                              # set_prior(paste0("normal(",mpriors[,"StrainPF5:total_cells_log"][1],",",mpriors[,"StrainPF5:total_cells_log"][2],")"), class="b", coef="protectPF5:protect_cells_log"),
                              set_prior(paste0("normal(",mpriors[,"StrainN2C3:total_cells_log"][1],",",mpriors[,"StrainN2C3:total_cells_log"][2],")"), class="b", coef="pathN2C3:path_cells_log"),
                              set_prior(paste0("normal(",mpriors[,"plant_age"][1],",",mpriors[,"plant_age"][2],")"), class="b", coef="plant_age")
                              # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
                            ),iter=4000)
brm_onlyinter_sdpriors


brm_onlyinter2_sdpriors <- brm(Alive ~  plant_age + protect:protect_cells_log:path_cells_log +(1 | plate) , data=dat_ac
                              , seed=12498215
                              , family="bernoulli"
                              , file="04_ratio_experiment/brm_onlyinter2_sdpriors"
                              , prior = c(
                                set_prior(paste0("normal(",mpriors[,"Intercept"][1],",",mpriors[,"Intercept"][2],")"), class="Intercept"),
                                # set_prior(paste0("normal(",mpriors[,"StrainN2C3"][1],",",mpriors[,"StrainN2C3"][2],")"), class="b", coef="pathN2C3"),
                                # set_prior(paste0("normal(",mpriors[,"StrainWCS365"][1],",",mpriors[,"StrainWCS365"][2],")"), class="b", coef="protectWCS365"),
                                # set_prior(paste0("normal(",mpriors[,"StrainCHAO"][1],",",mpriors[,"StrainCHAO"][2],")"), class="b", coef="protectCHAO"),
                                # set_prior(paste0("normal(",mpriors[,"StrainCH267"][1],",",mpriors[,"StrainCH267"][2],")"), class="b", coef="protectCH267"),
                                # set_prior(paste0("normal(",mpriors[,"StrainPF5"][1],",",mpriors[,"StrainPF5"][2],")"), class="b", coef="protectPF5"),
                                # set_prior("constant(0)", class="b", coef="pathMOCK:path_cells_log"),
                                # set_prior("constant(0)", class="b", coef="protectMOCK:protect_cells_log"),
                                set_prior("constant(0)", class="b", coef="protectMOCK:protect_cells_log:path_cells_log"),
                                # set_prior(paste0("normal(",mpriors[,"StrainWCS365:total_cells_log"][1],",",mpriors[,"StrainWCS365:total_cells_log"][2],")"), class="b", coef="protectWCS365:protect_cells_log"),
                                # set_prior(paste0("normal(",mpriors[,"StrainCHAO:total_cells_log"][1],",",mpriors[,"StrainCHAO:total_cells_log"][2],")"), class="b", coef="protectCHAO:protect_cells_log"),
                                # set_prior(paste0("normal(",mpriors[,"StrainCH267:total_cells_log"][1],",",mpriors[,"StrainCH267:total_cells_log"][2],")"), class="b", coef="protectCH267:protect_cells_log"),
                                # set_prior(paste0("normal(",mpriors[,"StrainPF5:total_cells_log"][1],",",mpriors[,"StrainPF5:total_cells_log"][2],")"), class="b", coef="protectPF5:protect_cells_log"),
                                # set_prior(paste0("normal(",mpriors[,"StrainN2C3:total_cells_log"][1],",",mpriors[,"StrainN2C3:total_cells_log"][2],")"), class="b", coef="pathN2C3:path_cells_log"),
                                set_prior(paste0("normal(",mpriors[,"plant_age"][1],",",mpriors[,"plant_age"][2],")"), class="b", coef="plant_age")
                                # set_prior("normal(0.74,0.4)", class="b", coef="lacZ_pathL")
                              ),iter=4000)
brm_onlyinter2_sdpriors


loo_of_sdpriors <- loo(brm_full_sdpriors, brm_onlyinter_sdpriors, brm_nointer_sdpriors,brm_onlyinter2_sdpriors)
loo_of_sdpriors
sink("04_ratio_experiment/loo_of_sdpriors.txt")
loo_of_sdpriors
sink()

###### plotting bayes results ######

# Plot mock estimates
draws_fullsdpriors  <- as_draws_df(brm_full_sdpriors) %>%
  rownames_to_column(var="iter") %>%
  select(iter, starts_with("b_")) %>%
  pivot_longer(-iter, names_to="coef", values_to="draw") %>%
  mutate(coef = gsub("b_","",coef)) %>%
  group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
  mutate(mod = "Full Model") %>%
  filter(med!=0) %>%
  mutate(coef_adj = str_replace_all(coef
                                       ,c(plant_age="Pooled effects:Effect of plant age\nat inoculation"
                                          ,`^pathN2C3$`="N2C3:Presence of strain"
                                          ,Intercept="Pooled effects:Plant with no\ninoculant"
                                          ,`path_cells_log:protect_cells_log`="Interaction between \npathogen and protective\ncell densities"
                                          ,`protect_cells_log:path_cells_log`="Interaction between \npathogen and protective\ncell densities"
                                          ,`protect_cells_log` = "Inoculation concentration"
                                          , `path_cells_log` = "Inoculation concentration"))) %>%
  separate(coef_adj, sep=":", into=c("Strain","Coefficient"), fill="right") %>%
  mutate(Coefficient = ifelse(is.na(Coefficient), "Presence of strain", Coefficient)) %>%
  mutate(Strain = gsub("path|protect", "",Strain)) %>%
  mutate(Coefficient = factor(Coefficient
                              , levels=c("Plant with no\ninoculant"
                                         , "Effect of plant age\nat inoculation"
                                         , "Presence of strain", "Inoculation concentration", "Interaction between \npathogen and protective\ncell densities"))) %>%
  mutate(Strain = factor(Strain, levels=c("Pooled effects", "N2C3","WCS365","CHAO","CH267","PF5")))

gg_fullsdpriors <- draws_fullsdpriors %>%
  # filter(coef!="Intercept") %>%
  # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q5*Q95)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~Strain, drop=TRUE, scales="free_x") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
        , axis.title.x = element_blank())+
  ylab("Estimate posterior")+
  labs(col="Two-sided\nPD < 0.05")+
  scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
gg_fullsdpriors
ggsave("04_ratio_experiment/gg_fullsdpriors.png", gg_fullsdpriors, width=6, height=3)

###### simplified stats plots ######
gg_fullsdpriors_simp <- draws_fullsdpriors %>%
  filter(!Strain%in%c("N2C3","Pooled effects")) %>%
  # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~Strain, drop=TRUE, scales="free_x") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
        , axis.title.x = element_blank())+
  ylab("Estimate posterior")+
  labs(col="Two-sided\nPD < 0.05")+
  scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
gg_fullsdpriors_simp
ggsave("04_ratio_experiment/gg_fullsdpriors_simp.png", gg_fullsdpriors_simp, width=6, height=3)

gg_monoculture_forpriors_simp <- draws_monoculture %>%
  filter(Strain !="MOCK") %>%
  # filter(!(Strain=="MOCK" & Coefficient == "Cell density")) %>%
  # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~Strain) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Coefficient Estimate\nfrom Bayesian model using \nhealthy/unhealthy response metric")+ 
  labs(col="Two-sided\nPD < 0.05")+
  scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
gg_monoculture_forpriors_simp
ggsave("04_ratio_experiment/gg_monoculture_forpriors_simp.png", gg_monoculture_forpriors_simp, height=4, width=8)
###### simplified bayes predictions ######
range_pathcellslog <- seq(range(dat_ac %>% pull(path_cells_log))[1], range(dat_ac$path_cells_log)[2], length.out=100)

newdat_expanded_full <- dat_ac%>%
  # filter(protect%in% c("WCS365","CHAO","CH267","PF5")) %>%
  mutate(protect = factor(protect, levels=c("MOCK","WCS365","CHAO","CH267","PF5"))) %>%
  select(protect, plant_age, path,  protect_cells_log, path) %>%
  mutate(path=as.character(path)) %>%
  mutate(plant_age = as.numeric(plant_age)) %>%
  distinct() %>%
  mutate(plate="newplate") %>%
  full_join(data.frame(path="N2C3", path_cells_log=range_pathcellslog), relationship="many-to-many") %>%
  mutate(path=ifelse(path_cells_log == 0, "MOCK", path)) %>%
  drop_na()

# Version where we have 95% prediction intervals
mean_pred95_predict_full <- t(apply(posterior_linpred(brm_full_sdpriors, newdata = newdat_expanded_full, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
newdat_expanded_full <- bind_cols(newdat_expanded_full,mean_pred95_predict_full)

gg_fullmodelpredictions_simp <- dat_ac %>%
  # filter(protect!="MOCK") %>%
  filter(plant_age==5) %>%
  ggplot() +
  geom_point(aes(x=round(path_cells_log), y=Alive, group=(round(protect_cells_log)) , col=(round(protect_cells_log)))
             , show.legend = FALSE, position=position_jitterdodgev(jitter.height=0.02, jitter.width=0.5, dodge.height=0.3)
             , cex=0.1) +
  geom_smooth(data=newdat_expanded_full %>% filter(plant_age==5), aes(x=round(path_cells_log), y=mean, col=(round(protect_cells_log)), group=factor(round(protect_cells_log))), se=FALSE) +
  # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
  # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  # )+
  facet_grid(.~protect, scales="free", drop = TRUE, space="free"
             , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
                                   , protect = protect_labeller))+
  scale_x_continuous(breaks=seq(0,8,1)) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  ylab("Probability of healthy plant") +
  xlab("N2C3 (pathogen) cells (log10(+1))") + labs(col="Protective cells\n(log10(+1))") +
  scale_color_gradient(low="gold", high="darkgreen")
# scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
gg_fullmodelpredictions_simp
ggsave("04_ratio_experiment/gg_fullmodelpredictions_simp.png", gg_fullmodelpredictions_simp, height=3, width=7)


## Stats save
sink("04_ratio_experiment/brms_summaries.txt")
print("Monocultures\n")
draws_monoculture
print("\nFull model\n")
draws_fullsdpriors
sink()


#### Plot experimental design for sanity check ####
#re-arranging plate IDs
exp_plate <- dat_ac %>%
  select(experiment, plate, plant_age) %>% distinct() %>%
  group_by(experiment, plant_age) %>% mutate(platenum=rank(plate)) %>% ungroup() %>%
  unite(plant_age, platenum, col="age_plate")

## Protective setup
gg_platesmap_ratio <- dat_ac %>%
  mutate(protect_cells_log = ifelse(protect=="MOCK", 6, protect_cells_log),
         path_cells_log = ifelse(path=="MOCK", 6, path_cells_log)) %>%
  left_join(exp_plate) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(experiment_date = gsub("_.*$","",experiment)) %>%
  ggplot(aes(x=col, y=row)) +
  geom_tile(aes(fill=protect, alpha=protect_cells_log)) +
  geom_point(aes(col=path, alpha=path_cells_log),pch=19) +
  # geom_text(aes(label=ratio))+
  # facet_grid(experiment~platenum) +
  facet_grid(age_plate~experiment_date, labeller = labeller(age_plate=c(`5_1`="5-day-old\nplant"
                                                                   ,`5_2`="5-day-old\nplant"
                                                                   ,`5_3`="5-day-old\nplant"
                                                                   ,`5_4`="5-day-old\nplant"
                                                                   ,`6_1`="6-day-old\nplant"
                                                                   ,`6_2`="6-day-old\nplant"                                                                   ,`6_1`="6-day-old\nplant"
                                                                   ,`7_1`="7-day-old\nplant"                                                                   ,`6_1`="6-day-old\nplant"
                                                                   ,`7_2`="7-day-old\nplant"
                                                                   
  ))) +
  scale_color_manual(values=c(MOCK="grey", N2C3="red")) +
  scale_fill_manual(values=c(MOCK="white", WCS365="darkgreen", CHAO="green", CH267="goldenrod", PF5="purple")) +
  labs(col="Pathogen treatment", fill="Commensal/protective\n strain", alpha = "Approx cell\nabundances (log10)")+
  xlab("Columns") + ylab("Rows")
  # theme_bw()
gg_platesmap_ratio
ggsave("04_ratio_experiment/gg_platesmap_ratio.png",gg_platesmap_ratio, height=10, width=16 )

#### brms for plant age ####

brm_plantage <- brm(Alive ~ plant_age + plant_age:path + plant_age:protect + plant_age:path:protect, data=dat_ac
                    , file = "04_ratio_experiment/brm_plantage"
                    , family = "bernoulli"
                    , iter=4000)

###### make stats plots ######
# Plot mock estimates
draws_plantage  <- as_draws_df(brm_plantage) %>%
  rownames_to_column(var="iter") %>%
  select(iter, starts_with("b_")) %>%
  pivot_longer(-iter, names_to="coef", values_to="draw") %>%
  mutate(coef = gsub("b_","",coef)) %>%
  group_by(coef) %>% summarise(Q2.5=quantile(draw, 0.025), Q97.5=quantile(draw, 0.975),  Q5=quantile(draw, 0.05), Q95=quantile(draw, 0.95), med = quantile(draw, 0.5)) %>%
  mutate(mod = "Plant Age Model") %>%
  filter(med!=0) %>%
  mutate(coef_adj = str_replace_all(coef
                                    ,c(Intercept = "Intercept::",
                                      `plant_age$`="plant_age::",
                                       `pathN2C3$`="pathN2C3:",
                                      `plant_age:protect` = "plant_age::protect"))) %>%
  separate(coef_adj, sep=":", into=c("plantage","path","protect"), fill="left") %>%
  mutate(path = gsub("path","",path), protect = gsub("protect","",protect)) %>%
  mutate(path = ifelse(is.na(path)|path=="","MOCK", path)
         , protect = ifelse(is.na(protect)|protect=="","MOCK", protect)) %>%
  mutate(path = factor(path, levels=c("MOCK","N2C3")), protect = factor(protect, levels=c("MOCK","WCS365","CHAO","CH267","PF5"))) %>%
  unite(protect, path, remove=FALSE, sep="+", col="StrainMix") %>%
  mutate(plantage = ifelse(is.na(plantage)|plantage=="Intercept", "No interaction","Interaction with\nplantage")) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK+MOCK"
                                                ,"MOCK+N2C3"
                                                ,"WCS365+MOCK"
                                                , "CHAO+MOCK"
                                                , "CH267+MOCK"
                                                , "PF5+MOCK"
                                                ,"WCS365+N2C3"
                                                , "CHAO+N2C3"
                                                , "CH267+N2C3"
                                                , "PF5+N2C3")))

gg_plantage <- draws_plantage %>%
  filter(coef!="Intercept") %>%
  mutate(comp = ifelse(path == "MOCK" & protect == "MOCK", "Plant\nalone", ifelse(path=="MOCK" | protect == "MOCK", "Monoculture", "Competition between\npathagen and protective"))) %>%
  mutate(Coefficient = ifelse(path == "MOCK" & protect == "MOCK", "Effect of\nplant age",ifelse(path=="MOCK", paste0("Interaction with\n",protect), ifelse(protect == "MOCK", paste0("Interaction with\n",path), paste0("Interaction with\n",StrainMix))))) %>%
  mutate(comp = factor(comp, levels=c("Plant\nalone", "Monoculture","Competition between\npathagen and protective"))) %>%
  arrange(path, protect) %>%
  mutate(Coefficient = factor(Coefficient, levels=unique(Coefficient))) %>%
   # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~comp, drop=TRUE, scales="free_x", space="free_x") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Estimate posterior")+
  labs(col="Two-sided\nPD < 0.05")+
  scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))
gg_plantage
ggsave("04_ratio_experiment/gg_plantage.png", gg_plantage, height=4, width=7)

###### prediction plot ######
newx_plantage = seq(5,7, length.out=13)
newdat_plantage <- distinct(brm_plantage$data) %>% mutate(plant_age=as.numeric(plant_age)) %>%
  select(-plant_age) %>% mutate(plate="newplate") %>%
  full_join(data.frame(plant_age = newx_plantage, plate="newplate"), relationship="many-to-many")
mean_pred95_predict_plantage <- t(apply(posterior_linpred(brm_plantage, newdata=newdat_plantage, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
predictions_plantage <- bind_cols(newdat_plantage, mean_pred95_predict_plantage)


gg_plantage_predictions <- predictions_plantage %>%
  ggplot() +
  geom_ribbon(aes(x=plant_age, ymin=Q2.5, ymax=Q97.5, fill=path), alpha=0.2)+
  geom_line(aes(x=plant_age, y=mean, col=path))+
  facet_grid(.~ protect) +
  scale_x_continuous(breaks=c(5,6,7))+
  scale_color_manual(values=c(N2C3="orange", MOCK="black")) +
  scale_fill_manual(values=c(N2C3="orange", MOCK="black")) +
  xlab("Plant age (days)") +ylab("Probability of healthy plant\n(95% PI)") +
  labs(col="Pathogen added", fill="Pathogen added")
gg_plantage_predictions
ggsave("04_ratio_experiment/gg_plantage_predictions.png",gg_plantage_predictions, height=3, width=7)


predictions_plantage %>%
  filter(plant_age %in% c(5,6,7)) %>%
  ggplot() +
  geom_pointrange(aes(x=plant_age, y=mean, ymin=Q2.5, ymax=Q97.5, col=path), alpha=0.2)+
  geom_line(aes(x=plant_age, y=mean,col=path))+
  facet_grid(.~ protect) +
  scale_x_continuous(breaks=c(5,6,7))+
  scale_color_manual(values=c(N2C3="darkorange", MOCK="black")) +
  xlab("Plant age (days)") +ylab("Posterior prediction intervals\n(95% PI)")+
  labs(col="Pathogen added")
  
### STATISTICAL SUMMARIES
sink("04_ratio_experiment/brms_results_text.txt")
print("\n Monoculture age\n")
fixef(brm_monocultures)
print("\n Full model age\n")
fixef(brm_full_sdpriors)
print("\nplant age\n")
fixef(brm_plantage)
sink()



dat_ac %>%
  group_by(protect, path, plant_age) %>%
  reframe(propAlive = mean(Alive)) 

dat_ac %>%
  group_by(protect, path) %>%
  reframe(propAlive = mean(Alive)) 

dat_ac %>%
  group_by(protect, ratio, plant_age) %>%
  reframe(propAlive = mean(Alive)) 
