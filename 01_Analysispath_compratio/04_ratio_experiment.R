#!bin/bash
library(cowplot)
library(gridExtra)
library(lubridate)
library(ggstance)
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

#### Brms with FULL model, priors ####
###### Monoculture for priors ######

dat_ac_monocultures_nolacZ <- dat_ac %>%
  filter(protect =="MOCK" | path == "MOCK") %>%
  # filter(!(protect=="MOCK"&path=="MOCK")) %>%
  mutate(Strain = ifelse(protect == "MOCK", as.character(path), as.character(protect))) %>%
  select(UniqueID, Alive, Healthiness_hsv, Strain, protect, path, plant_age, total_cells_log, plate, experiment, path_cells_log, protect_cells_log, all_plant_pixels) %>%
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
    Coefficient=="total_cells_log", "Strain dose", Coefficient))) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","N2C3","WCS365","CHAO","CH267","PF5"))) %>%
  mutate(Coefficient = factor(Coefficient, levels=c("Effect of plant age\nat inoculation","Presence of strain","Strain dose")))
# 
# gg_monoculture_forpriors <- draws_monoculture %>%
#   filter(Strain !="MOCK") %>%
#   # filter(!(Strain=="MOCK" & Coefficient == "Cell density")) %>%
#   # filter(protect!="Pathogen only") %>%
#   mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
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
# ggsave("04_ratio_experiment/gg_monoculture_forpriors.png", gg_monoculture_forpriors, width=8, height=4)

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
                                          ,`path_cells_log:protect_cells_log`="Interaction with\nN2C3 dose"
                                          ,`protect_cells_log:path_cells_log`="Interaction with\nN2C3 dose"
                                          ,`protect_cells_log` = "Inoculation dose\nof protective"
                                          , `path_cells_log` = "Inoculation dose\nof protective"))) %>%
  separate(coef_adj, sep=":", into=c("Strain","Coefficient"), fill="right") %>%
  mutate(Coefficient = ifelse(is.na(Coefficient), "Presence of strain", Coefficient)) %>%
  mutate(Strain = gsub("path|protect", "",Strain)) %>%
  mutate(Coefficient = factor(Coefficient
                              , levels=c("Plant with no\ninoculant"
                                         , "Effect of plant age\nat inoculation"
                                         , "Presence of strain", "Inoculation dose\nof protective", "Interaction with\nN2C3 dose"))) %>%
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
  mutate(Strain = ifelse(Strain=="CHAO","CHA0", ifelse(Strain=="PF5","Pf5", as.character(Strain)))) %>%
  mutate(Strain2 = paste0(Strain, "\neffects")) %>%
  mutate(Strain2 = factor(Strain2, levels=c("WCS365\neffects", "CHA0\neffects", "CH267\neffects", "Pf5\neffects"))) %>%
  # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  # geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5), col="black") +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q5, ymax=Q95), col=rgb(0.5, 0.5, 0.5, 0.5), lwd=2)+
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~Strain2, drop=TRUE, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Coefficient Estimate\n(>0 means positive effect\non plant health)")+
  labs(col="Two-sided\nPD < 0.05")+
  # scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="lightpink2")) +
  xlab("Predictor in model")
gg_fullsdpriors_simp
ggsave("04_ratio_experiment/gg_fullsdpriors_simp.png", gg_fullsdpriors_simp, width=6, height=4)

gg_monoculture_forpriors_simp <- draws_monoculture %>%
  filter(Strain !="MOCK") %>%
  mutate(Strain = ifelse(Strain=="CHAO","CHA0", ifelse(Strain=="PF5","Pf5", as.character(Strain)))) %>%
  mutate(Strain2 = paste0(Strain, "\neffects")) %>%
  mutate(Strain2 = factor(Strain2, levels=c("N2C3\neffects","WCS365\neffects", "CHA0\neffects", "CH267\neffects", "Pf5\neffects"))) %>%
  # filter(!(Strain=="MOCK" & Coefficient == "Cell density")) %>%
  # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  # geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5), col="black") +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q5, ymax=Q95), col=rgb(0.5, 0.5, 0.5, 0.5), lwd=2)+geom_hline(aes(yintercept=0)) +
  facet_grid(.~Strain2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Coefficient Estimate\n(>0 means positive effect\non plant health)")+ 
  labs(col="Two-sided\nPD < 0.05")+
  # scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink")) +
  xlab("Predictor in model")
gg_monoculture_forpriors_simp
ggsave("04_ratio_experiment/gg_monoculture_forpriors_simp.png", gg_monoculture_forpriors_simp, height=4, width=6)


gg_monoculture_rawdat <- dat_ac_monocultures_nolacZ %>%
  # filter(Strain !="MOCK") %>%
  filter(plant_age==5) %>%
  mutate(Strain = ifelse(Strain=="CHAO","CHA0", ifelse(Strain=="PF5","Pf5", as.character(Strain)))) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","WCS365","CHA0","CH267","Pf5","N2C3"))) %>%
  # mutate(Strain2 = paste0(Strain, "\neffects")) %>%
  # mutate(Strain2 = factor(Strain2, levels=c("N2C3\neffects","WCS365\neffects", "CHA0\neffects", "CH267\neffects", "Pf5\neffects"))) %>%
  # filter(!(Strain=="MOCK" & Coefficient == "Cell density")) %>%
  # filter(protect!="Pathogen only") %>%
  # mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  # mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  # mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               # ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot(aes(x=total_cells_log, y=Healthiness_hsv)) +
  # geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_jitter(aes(fill=UniqueID, size=all_plant_pixels),alpha = 0.8, col='black', pch=21, height=0, width=0.2) +
  # geom_boxplot() +
  geom_hline(aes(yintercept=400)) +
  facet_grid(.~Strain, scales="free_x",drop = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5) )+
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_fill_manual(values=col_all) +
  scale_radius(range=c(0.5,5)) +
  ylab("Health Score") +
  xlab("Inoculation dose (log10(+1))")+
  guides(fill='none')+
  labs(cex="Plant size\n(pixels)")
gg_monoculture_rawdat
ggsave("04_ratio_experiment/gg_monoculture_rawdat.png", gg_monoculture_rawdat, height=3, width=8)


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

healthColRamp <- colorRampPalette(c("gold", "darkseagreen4"))

gg_fullmodelpredictions_simp <- dat_ac %>%
  # filter(protect!="MOCK") %>%
  filter(plant_age==5) %>%
  ggplot() +
  geom_point(aes(x=round(path_cells_log), y=Alive, group=(round(protect_cells_log)) , col=factor(round(protect_cells_log)))
             , show.legend = FALSE, position=position_jitterdodgev(jitter.height=0.02, jitter.width=0.5, dodge.height=0.3)
             , cex=0.1) +
  geom_smooth(data=newdat_expanded_full %>% filter(plant_age==5), aes(x=round(path_cells_log), y=mean, col=factor(round(protect_cells_log)), group=factor(round(protect_cells_log))), se=FALSE) +
  # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
  # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  # )+
  facet_grid(.~protect, scales="free", drop = TRUE, space="free"
             , labeller = labeller(plant_age = plant_age_labeller, `round(path_cells_log)`=N2C3_labeller
                                   , protect = protect_labeller))+
  scale_x_continuous(breaks=seq(0,8,1)) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  ylab("Probability of healthy plant") +
  xlab("N2C3 (pathogen) inoculation\ndose (log10(+1))") + labs(col="Protective inoculation\n dose (log10(+1))"
                                                               # , title="Effect of inoculation dose on\nplant health outcome (against N2C3)"
                                                               ) +
  # scale_color_gradient(low="gold", high="darkseagreen4") +
  scale_color_manual(values=c(healthColRamp(6)))+
  theme_bw()
# scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
gg_fullmodelpredictions_simp
ggsave("04_ratio_experiment/gg_fullmodelpredictions_simp.png", gg_fullmodelpredictions_simp, height=3, width=6)

## For monoculture
range_totalcellslog <- seq(range(dat_ac %>% pull(protect_cells_log))[1], range(dat_ac$protect_cells_log)[2], length.out=100)

brm_monocultures

newdat_expanded_mono <- dat_ac_monocultures_nolacZ%>%
  # filter(protect%in% c("WCS365","CHAO","CH267","PF5")) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","N2C3","WCS365","CHAO","CH267","PF5"))) %>%
  select(Strain, plant_age) %>%
  distinct() %>%
  mutate(plant_age = as.numeric(plant_age)) %>%
  mutate(plate="newplate") %>%
  full_join(data.frame(plate="newplate", total_cells_log=range_pathcellslog), relationship="many-to-many") %>%
  mutate(total_cells_log=ifelse(Strain=="MOCK",0, total_cells_log)) %>%
  # mutate(Strain=ifelse(total_cells_log == 0, "MOCK", as.character(Strain))) %>%
  drop_na()

# Version where we have 95% prediction intervals
mean_pred95_predict_mono <- t(apply(posterior_linpred(brm_monocultures, newdata = newdat_expanded_mono, allow_new_levels=TRUE, transform = TRUE), MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
newdat_expanded_mono <- bind_cols(newdat_expanded_mono,mean_pred95_predict_mono)
newdat_expanded_mono_forplot <- newdat_expanded_mono %>% filter(plant_age==5) %>%
  rowwise() %>%
  mutate(total_cells_log = ifelse(Strain=="MOCK",rnorm(1, 0, 0.7), total_cells_log ) ) %>% ungroup() %>%
  filter((total_cells_log>2 & Strain !="MOCK") | Strain == "MOCK")%>%
  mutate(Strain = ifelse(Strain=="CHAO","CHA0", ifelse(Strain=="PF5","Pf5", as.character(Strain)))) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","WCS365","CHA0","CH267","Pf5","N2C3")))

gg_monomodelpredictions_simp <- dat_ac_monocultures_nolacZ %>%
  # filter(Strain!="MOCK") %>%
  mutate(Strain = ifelse(Strain=="CHAO","CHA0", ifelse(Strain=="PF5","Pf5", as.character(Strain)))) %>%
  mutate(Strain = factor(Strain, levels=c("MOCK","WCS365","CHA0","CH267","Pf5","N2C3"))) %>%
  filter(plant_age==5) %>%
  ggplot() +
  geom_jitter(aes(x=total_cells_log, y=Alive)
             , show.legend = FALSE, height=0.1, width=0.1
             , cex=0.1) +
  geom_ribbon(data=newdat_expanded_mono_forplot, aes(x=(total_cells_log), ymin=Q2.5, ymax=Q97.5), alpha=0.25) +
  geom_line(data=newdat_expanded_mono_forplot, aes(x=(total_cells_log), y=mean)) +
  # scale_color_manual(values = c('5'="darkolivegreen1",`6` = "olivedrab3", `7`="darkolivegreen")) +
  # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  # )+
  facet_grid(.~Strain, scales="free", drop = TRUE, space="free"
             , labeller = labeller(plant_age = plant_age_labeller, `round(total_cells_log)`=N2C3_labeller
                                   ))+
  scale_x_continuous(breaks=c(0,3,4,5,6,7,8)) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  ylab("Probability of healthy plant") +
  xlab("Inoculation dose (log10(+1))") + 
  # scale_color_gradient(low="gold", high="darkseagreen4") +
  scale_color_manual(values=c(healthColRamp(6)))+
  theme_bw()
# scale_color_manual(values=c(`3`="goldenrod", `4`="goldenrod4", `5`="brown", `6`="darkred", `7`="red"))
gg_monomodelpredictions_simp
ggsave("04_ratio_experiment/gg_monomodelpredictions_simp.png", gg_monomodelpredictions_simp, height=3, width=7)



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
brm_plantage
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
  mutate(comp = ifelse(path == "MOCK" & protect == "MOCK", "Plant\nalone", ifelse(path=="MOCK" | protect == "MOCK", "Plant age interaction\nwith strains in monoculture", "Plant age interaction with\npathogen and protective"))) %>%
  mutate(Coefficient = ifelse(path == "MOCK" & protect == "MOCK", "Effect of\nplant age",ifelse(path=="MOCK", paste0("Interaction with\n",protect), ifelse(protect == "MOCK", paste0("Interaction with\n",path), paste0("Interaction with\n",StrainMix))))) %>%
  mutate(comp = factor(comp, levels=c("Plant\nalone", "Plant age interaction\nwith strains in monoculture","Plant age interaction with\npathogen and protective"))) %>%
  arrange(path, protect) %>%
  mutate(Coefficient = factor(Coefficient, levels=unique(Coefficient))) %>%
   # filter(protect!="Pathogen only") %>%
  mutate(SigEffect = 0<sign(Q2.5*Q97.5)) %>%
  mutate(SigEffect2 = 0<sign(Q5*Q95)) %>%
  mutate(SigEffectAll = ifelse(SigEffect, "PD<0.025 (95% CI)", 
                               ifelse(SigEffect2, "PD<0.05 (90% CI)", NA))) %>%
  ggplot() +
  # geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5, col=SigEffectAll), show.legend = TRUE) +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q2.5, ymax=Q97.5), col="black") +
  geom_pointrange(aes(x=Coefficient, y=med, ymin=Q5, ymax=Q95), col=rgb(0.5, 0.5, 0.5, 0.5), lwd=2) +
  geom_hline(aes(yintercept=0)) +
  facet_grid(.~comp, drop=TRUE, scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  # scale_color_manual(values=c(`PD<0.025 (95% CI)`="red", `PD<0.05 (90% CI)`="pink"))+
  # labs(col="Two-sided\nPD < 0.05")+
  ylab("Coefficient Estimate\n(>0 means positive effect\non plant health)")
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
  geom_point(data=dat_ac, aes(x=plant_age, y=Alive, group=path , col=path)
             , position=position_jitterdodgev(jitter.height=0.1, jitter.width=0.5, dodge.height=-0.3)
             , cex=0.1)+
  # geom_jitter(data=dat_ac, aes(x=plant_age, y=as.numeric(Alive), col=path), width=0.25, height=0.1, cex=0.1) +
  geom_ribbon(aes(x=plant_age, ymin=Q2.5, ymax=Q97.5, fill=path), alpha=0.2)+
  geom_line(aes(x=plant_age, y=mean, col=path))+
  facet_grid(.~ protect, labeller = labeller(protect=protect_labeller)) +
  scale_x_continuous(breaks=c(5,6,7))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(values=c(N2C3="darkorange", MOCK="black")) +
  scale_fill_manual(values=c(N2C3="darkorange", MOCK="black")) +
  xlab("Plant age (days)") +ylab("Probability of healthy plant") +
  labs(col="Pathogen added", fill="Pathogen added")+
  theme_bw() +
  theme(legend.position = "top")
gg_plantage_predictions
ggsave("04_ratio_experiment/gg_plantage_predictions.png",gg_plantage_predictions, height=3.5, width=7)


## Idea: what if we subtracted average MOCK effect from each individual strain
## Then we subtracted the interaction between WCSC and N2C3 for each
posterior_linpred(brm_plantage, newdata=newdat_plantage, allow_new_levels=TRUE, transform = TRUE)
mean_pred95_predict_plantage <- t(apply(, MARGIN = 2, function(x) c(mean=mean(x), Q2.5=as.numeric(quantile(x, 0.025)), Q97.5=as.numeric(quantile(x, 0.975)))))
predictions_plantage <- bind_cols(newdat_plantage, mean_pred95_predict_plantage)

mockMeans <- predictions_plantage %>% filter(path=="MOCK", protect=="MOCK") %>%
  select(plant_age, mean) %>%
  rename(mean_mock = mean)
monocultureMeans <- predictions_plantage %>% filter(path=="MOCK", protect!="MOCK") %>%
  select(protect, plant_age, mean) %>%
  rename(mean_monoprotect = mean)
predictions_plantage %>%
  filter(path=="MOCK" | protect=="MOCK") %>%
  filter(! (path=="MOCK" & protect == "MOCK")) %>%
  left_join(mockMeans, relationship = "many-to-many") %>%
  mutate(adjMean = mean-mean_mock, adjQ2.5 = Q2.5-mean_mock, adjQ97.5 = Q97.5-mean_mock) %>%
  mutate(Strain = ifelse(!is.na(path) & path!="MOCK", as.character(path), ifelse(!is.na(protect) & protect !="MOCK", as.character(protect), NA))) %>%
  mutate(Strain = factor(Strain, levels=c("WCS365","CHAO","CH267","PF5","N2C3"))) %>%
  ggplot() +
  geom_line(aes(x=plant_age, y=adjMean)) +
  geom_ribbon(aes(x=plant_age, ymin=adjQ2.5, ymax=adjQ97.5), alpha=0.25, fill="yellow2") +
  facet_grid(.~Strain)+
  ylab("Decrease in probability of healthy plant\nrelative to MOCK treatment")+
  xlab("Plant age at inoculation") +
  scale_x_continuous(breaks=c(5,6,7)) +
  scale_y_continuous(breaks=c(0,-0.2, 0.4, -0.6,-0.8))+
  theme_bw()

predictions_plantage %>%
  filter(protect!="MOCK") %>%
  # filter(! (path=="MOCK" & protect == "MOCK")) %>%
  left_join(monocultureMeans, relationship = "many-to-many") %>%
  mutate(adjMean = mean-mean_monoprotect, adjQ2.5 = Q2.5-mean_monoprotect, adjQ97.5 = Q97.5-mean_monoprotect) %>%
  # mutate(Strain = ifelse(!is.na(path) & path!="MOCK", as.character(path), ifelse(!is.na(protect) & protect !="MOCK", as.character(protect), NA))) %>%
  # mutate(Strain = factor(Strain, levels=c("WCS365","CHAO","CH267","PF5","N2C3"))) %>%
  ggplot() +
  geom_line(aes(x=plant_age, y=adjMean, group=path, col=path)) +
  geom_ribbon(aes(x=plant_age, ymin=adjQ2.5, ymax=adjQ97.5, group=path, fill=path), alpha=0.25) +
  facet_grid(.~protect)+
  ylab("Decrease in probability of healthy plant\nrelative to monoculture treatment")+
  xlab("Plant age at inoculation") +
  scale_x_continuous(breaks=c(5,6,7)) +
  scale_y_continuous(breaks=c(0,-0.2, 0.4, -0.6,-0.8))+
  scale_color_manual(values=c(N2C3="darkorange", MOCK="yellow2")) +
  scale_fill_manual(values=c(N2C3="darkorange", MOCK="yellow2")) +
  theme_bw()

predictions_plantage %>%
  # filter(path=="MOCK" | protect=="MOCK") %>%
  # filter(! (path=="MOCK" & protect == "MOCK")) %>%
  left_join(mockMeans, relationship = "many-to-many") %>%
  mutate(adjMean = mean-mean_mock, adjQ2.5 = Q2.5-mean_mock, adjQ97.5 = Q97.5-mean_mock) %>%
  mutate(Strain = ifelse(!is.na(path) & path!="MOCK", as.character(path), ifelse(!is.na(protect) & protect !="MOCK", as.character(protect), NA))) %>%
  mutate(Strain = factor(Strain, levels=c("WCS365","CHAO","CH267","PF5","N2C3"))) %>%
  ggplot() +
  geom_line(aes(x=plant_age, y=adjMean, group=path, col=path)) +
  geom_ribbon(aes(x=plant_age, ymin=adjQ2.5, ymax=adjQ97.5, group=path, fill=path), alpha=0.25) +
  facet_grid(.~protect)+
  ylab("Decrease in probability of healthy plant\nrelative to MOCK treatment")+
  xlab("Plant age at inoculation") +
  scale_x_continuous(breaks=c(5,6,7)) +
  scale_y_continuous(breaks=c(0,-0.2, 0.4, -0.6,-0.8))+
  scale_color_manual(values=c(N2C3="darkorange", MOCK="yellow2")) +
  scale_fill_manual(values=c(N2C3="darkorange", MOCK="yellow2")) +
  theme_bw()


### verison where it's monoculture first
predictions_plantage %>%
  mutate(Strain = ifelse(!is.na(path) & path!="MOCK" & protect == "MOCK", as.character(path), ifelse(!is.na(protect) & protect !="MOCK", as.character(protect), ifelse(protect=="MOCK" & path=="MOCK","MOCK",NA)))) %>%
  filter(!is.na(Strain), !(protect!="MOCK" & path!="MOCK")) %>%
  ggplot() +
  geom_point(data=dat_ac, aes(x=plant_age, y=Alive)
             # , position=position_jitterdodgev(jitter.height=0.1, jitter.width=0.5, dodge.height=-0.3)
             , cex=0.1)+
  # geom_jitter(data=dat_ac, aes(x=plant_age, y=as.numeric(Alive), col=path), width=0.25, height=0.1, cex=0.1) +
  geom_ribbon(aes(x=plant_age, ymin=Q2.5, ymax=Q97.5), alpha=0.2)+
  geom_line(aes(x=plant_age, y=mean))+
  facet_grid(.~ Strain, labeller = labeller(protect=protect_labeller)) +
  scale_x_continuous(breaks=c(5,6,7))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  # scale_color_manual(values=c(N2C3="darkorange", MOCK="black")) +
  # scale_fill_manual(values=c(N2C3="darkorange", MOCK="black")) +
  xlab("Plant age (days)") +ylab("Probability of healthy plant") +
  labs(col="Pathogen added", fill="Pathogen added")+
  theme_bw() +
  theme(legend.position = "top")



# 
# predictions_plantage %>%
#   filter(plant_age %in% c(5,6,7)) %>%
#   ggplot() +
#   geom_pointrange(aes(x=plant_age, y=mean, ymin=Q2.5, ymax=Q97.5, col=path), alpha=0.2)+
#   geom_line(aes(x=plant_age, y=mean,col=path))+
#   facet_grid(.~ protect) +
#   scale_x_continuous(breaks=c(5,6,7))+
#   scale_color_manual(values=c(N2C3="darkorange", MOCK="black")) +
#   xlab("Plant age (days)") +ylab("Posterior prediction intervals\n(95% PI)")+
#   labs(col="Pathogen added")
#   
### STATISTICAL SUMMARIES
sink("04_ratio_experiment/brms_results_text.txt")
print("\n Monoculture age\n")
fixef(brm_monocultures)
print("\n Full model age\n")
fixef(brm_full_sdpriors)
print("\nplant age\n")
fixef(brm_plantage)
sink()

#### Age and ratio panelled plot ####

## Join with other
gg_ratio_and_age_panelled <- plot_grid(gg_fullmodelpredictions_simp, gg_plantage_predictions+
            theme(legend.position = "right"), nrow=2, align="v", axis="lr")
gg_ratio_and_age_panelled
ggsave(filename = "04_ratio_experiment/gg_ratio_and_age_panelled.png", gg_ratio_and_age_panelled, height=7, width=7)
