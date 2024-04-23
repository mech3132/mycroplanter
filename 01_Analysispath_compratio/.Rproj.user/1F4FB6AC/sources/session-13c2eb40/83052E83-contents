#!bin/bash
library(cowplot)
library(gridExtra)
library(lubridate)
library(ggstance)
library(tidyverse)
library(brms)

dir.create("05_plant_mutant")
########## Load ########
load("03_EDA_and_adj/allDat_final.RData")
load("03_EDA_and_adj/plant_colours.RData")

allDat_final <- allDat_final %>%
  mutate(lacZ_path = ifelse(is.na(lacZ_path), "WT",lacZ_path)) %>%
  mutate(lacZ_path = factor(lacZ_path, levels=c("WT","L"))) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400))


dat_plantmutant <- allDat_final %>% filter(experiment=="2023-07-05_plantmutants") %>%
  filter(plant!='nej')
# What we want is a ratio where col0 is protected, but the other mutants may not be
# Therefore, we are going to search for the lowest protect:path ratio that still allows statistically
# protection against N2C3

#### EDA plots ####

dat_plantmutant %>% 
  ggplot(aes(x=plant, y=Healthiness_hsv)) +
  geom_boxplot()+
  geom_jitter(aes(col=factor(path_od))) +
  facet_grid(path~protect)+
  # scale_colour_manual(values=col_all)+
  ylab("Health Score") +
  xlab("Plant genotype")

dat_plantmutant %>% 
  # filter(path!="MOCK") %>%
  ggplot(aes(x=plant, y=Healthiness_hsv)) +
  geom_boxplot()+
  geom_jitter(aes(col=factor(path_od))) +
  facet_grid(ratio~protect)+
  # scale_colour_manual(values=col_all)+
  ylab("Health Score") +
  xlab("Plant genotype")

dat_plantmutant%>% 
  filter(path_od %in% c(0,0.0001)) %>%
  ggplot(aes(x=plant, y=Healthiness_hsv)) +
  geom_boxplot()+
  geom_jitter(aes(col=factor(path_od))) +
  facet_grid(path~protect)+
  # scale_colour_manual(values=col_all)+
  ylab("Health Score") +
  xlab("Plant genotype")

### Plot all
dat_plantmutant %>%
  # filter(!col%in% c(1,12)) %>%
  mutate(ratio = factor(ratio, levels=c("1-0","10-1","1-1","1-10"))) %>%
  ggplot() +
  geom_jitter(aes(x=plant, y=Healthiness_hsv, col=plate), width=0.2, height=0) +
  facet_grid(protect ~ path_od)

dat_plantmutant %>%
  mutate(ratio = factor(ratio, levels=c("1-0","10-1","1-1","1-10"))) %>%
  group_by(plant, protect, path_od) %>%
  summarise(nAlive = sum(Alive),total=n()) %>%
  ungroup() %>%
  mutate(propAlive = nAlive/total) %>%
  # select(plant, protect, path_od, nAlive, total)
  # mutate(Alive = Healthiness_hsv>400) %>%
  # filter(!col%in% c(1,12)) %>%
  # filter(!protect %in% c("MOCK")) %>%
  ggplot() +
  geom_bar(aes(x=plant, y=propAlive), stat="identity") +
  facet_grid(protect ~ path_od)

#### Singel BRMS and plots ####
### Do separate brms for all
dat_plantmutant_forbrms <- dat_plantmutant %>%
  filter(protect !="CH267") %>%
  # filter(ratio !="1-10", path_od != 0.01) %>%
  select(Healthiness_hsv, Alive, path_od, path, protect, plant, path_cells_log, ratio) %>%
  mutate(ratio = factor(ratio, levels=c("1-0","10-1","1-1","1-10"))) %>%
  mutate(plant = factor(plant, levels=c("col0","bbc","bik1"))) %>%
  mutate(protect = factor(protect, levels=c("MOCK","WCS365","CHAO","PF5"))) %>%
  mutate(path = factor(path, levels=c("MOCK","N2C3")))


dat_plantmutant_summarise <- dat_plantmutant_forbrms %>%
  group_by(path, ratio, path_cells_log, plant,path_od,  protect) %>%
  summarise(propAlive = sum(Alive)/n()) %>%
  ungroup() 
dat_plantmutant_summarise%>%
  ggplot() +
  geom_bar(aes(x=path_cells_log, y=propAlive, group=plant, fill=plant), stat='identity', position = position_dodge())+
  facet_grid(protect~path, scales = "free_x", space = c("free_x"))+
  scale_fill_manual(values=c(col0="darkblue", bbc="pink", bik1="magenta"))

dat_plantmutant_forbrms_onlyn2c3 <- dat_plantmutant_forbrms %>%
  filter(path=="N2C3")
brm_alldat_bin <- brm(Alive ~ plant*protect*path_cells_log , data = dat_plantmutant_forbrms_onlyn2c3
                  , family="bernoulli"
                  , seed = 1245
                  , control = list(max_treedepth=15)
                  , iter=4000
                  , file="05_plant_mutant/brm_alldatnomock_bin"
)
fixef(brm_alldat_bin)

dat_plantmutant_forbrms_onlymock <- dat_plantmutant_forbrms %>%
  filter(path=="MOCK")
brm_alldat_mockonly_bin <- brm(Alive ~ plant*protect , data = dat_plantmutant_forbrms_onlymock
                      , family="bernoulli"
                      , seed = 1245
                      , control = list(max_treedepth=15)
                      , iter=4000
                      , file="05_plant_mutant/brm_alldatonlymock_bin"
)
brm_alldat_mockonly_bin

sink("05_plant_mutant/stats_and_numbers.txt")
brm_alldat_bin

print("")
print("")
brm_alldat_mockonly_bin
sink()

dat_forpred<- dat_plantmutant_forbrms_onlyn2c3 %>%
  select( plant, protect, path_cells_log, path, path_od, ratio) %>%
  distinct() 
alldatbin_predict <- posterior_linpred(brm_alldat_bin, newdata = dat_forpred, transform=TRUE)
dat_withpred <- alldatbin_predict %>% as.data.frame() %>%
  apply(MARGIN=2, function(x) c(med = quantile(x,0.5), lwr = quantile(x, 0.025), upr = quantile(x, 0.975), lwr5 = quantile(x, 0.05), upr95 = quantile(x,0.95))) %>%
  t() %>%
  cbind(dat_forpred) 

dat_forpred_mock <- dat_plantmutant_forbrms_onlymock %>%
  select( plant, protect, path_cells_log, path, path_od,ratio) %>%
  distinct() 
alldatbin_predict_mock <- posterior_linpred(brm_alldat_mockonly_bin, newdata = dat_forpred_mock, transform=TRUE)
dat_withpred_mock <- alldatbin_predict_mock %>% as.data.frame() %>%
  apply(MARGIN=2, function(x) c(med = quantile(x,0.5), lwr = quantile(x, 0.025), upr = quantile(x, 0.975), lwr5 = quantile(x, 0.05), upr95 = quantile(x,0.95))) %>%
  t() %>%
  cbind(dat_forpred_mock) 

# Combine
alldat_withpred <- full_join(dat_withpred, dat_withpred_mock)

dat_plantmutant_summarise <- dat_plantmutant_summarise %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10 (commensal:N2C3)", 
                         ifelse(factor(round(path_cells_log))==6, "1-1 (commensal:N2C3)",
                                ifelse(factor(round(path_cells_log))==5, "10-1 (commensal:N2C3)",
                                       ifelse(factor(round(path_cells_log))==0, "No pathogen", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No pathogen", "10-1 (commensal:N2C3)","1-1 (commensal:N2C3)","1-10 (commensal:N2C3)"))) 
# 
# gg_brms_plantmutant_alldat <- alldat_withpred %>%
#   mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10", 
#                          ifelse(factor(round(path_cells_log))==6, "1-1",
#                                 ifelse(factor(round(path_cells_log))==5, "10-1",
#                                        ifelse(factor(round(path_cells_log))==0, "No pathogen", NA))))) %>%
#   mutate(ratio2 = factor(ratio2, levels=c("No pathogen", "10-1","1-1","1-10"))) %>%
#   ggplot()+
#   geom_bar(dat_plantmutant_summarise, mapping=aes(x=ratio2, y=propAlive, group=plant, fill=plant),alpha=0.25 , stat='identity', position = position_dodge())+
#   # geom_pointrange(aes(x=factor(round(path_cells_log)), y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=plant), position = position_dodge2(width=0.9)) +
#   # geom_pointrange(aes(x=factor(round(path_cells_log)), y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=plant), position = position_dodge2(width=0.9), lwd=2) +
#   geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=plant), position = position_dodge2(width=0.9)) +
#   # geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=plant), position = position_dodge2(width=0.9), lwd=2) +
#   facet_grid(protect ~ path, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
#   # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
#   scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
#   scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
# labs(col="Plant\ngenotype", fill="Plant\ngenotype") +
#   ylab("Probability of healthy plant\n(Health score > 400)")+ 
#   xlab("Protective : Pathogen ratio")
# gg_brms_plantmutant_alldat
# ggsave(filename="05_plant_mutant/gg_brms_plantmutant_alldat.png", gg_brms_plantmutant_alldat, height=8, width=10)
# 
# ## reduced version
# gg_brms_plantmutant_simp <- alldat_withpred %>%
#   mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10", 
#                          ifelse(factor(round(path_cells_log))==6, "1-1",
#                                 ifelse(factor(round(path_cells_log))==5, "10-1",
#                                        ifelse(factor(round(path_cells_log))==0, "No pathogen", NA))))) %>%
#   mutate(ratio2 = factor(ratio2, levels=c("No pathogen", "10-1","1-1","1-10"))) %>%
#   filter(ratio2 !="1-10") %>%
#   ggplot()+
#   geom_bar(dat_plantmutant_summarise %>% filter(ratio2 !="1-10"), mapping=aes(x=ratio2, y=propAlive, group=plant, fill=plant),alpha=0.25 , stat='identity', position = position_dodge())+
#   # geom_pointrange(aes(x=factor(round(path_cells_log)), y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=plant), position = position_dodge2(width=0.9)) +
#   # geom_pointrange(aes(x=factor(round(path_cells_log)), y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=plant), position = position_dodge2(width=0.9), lwd=2) +
#   geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=plant, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.9), size=0) +
#   geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=plant, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.9), size=0) +
#   facet_grid(protect ~ path, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
#   # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
#   scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
#   scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
#   scale_linewidth_manual(values=c(2, 0.5)) +
#   labs(col="Plant genotype\n(Model estimates)", fill="Plant genotype\n(Bars represent\nobserved proportions)", linewidth="Model estimates") +
#   ylab("Probability of healthy plant\n(Health score > 400)")+ 
#   xlab("Protective : Pathogen ratio")+
#   guides(color="none")
# gg_brms_plantmutant_simp
# ggsave(filename="05_plant_mutant/gg_brms_plantmutant_simp.png", gg_brms_plantmutant_simp, height=4, width=6)


#### Plot predictions on actualy data ####
# get data
dat_plantmutant_forplotting <- dat_plantmutant_forbrms %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10 (commensal:N2C3)", 
                         ifelse(factor(round(path_cells_log))==6, "1-1 (commensal:N2C3)",
                                ifelse(factor(round(path_cells_log))==5, "10-1 (commensal:N2C3)",
                                       ifelse(factor(round(path_cells_log))==0, "No\nN2C3", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No\nN2C3", "1-10 (commensal:N2C3)","1-1 (commensal:N2C3)","10-1 (commensal:N2C3)"))) %>%
  filter(ratio2 !="1-10") %>%
  mutate(pathogen=factor(round(path_cells_log)))
  
gg_brms_plantmutant_simp_wdat <- alldat_withpred %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10 (commensal:N2C3)", 
                         ifelse(factor(round(path_cells_log))==6, "1-1 (commensal:N2C3)",
                                ifelse(factor(round(path_cells_log))==5, "10-1 (commensal:N2C3)",
                                       ifelse(factor(round(path_cells_log))==0, "No\nN2C3", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No\nN2C3", "1-10 (commensal:N2C3)","1-1 (commensal:N2C3)","10-1 (commensal:N2C3)"))) %>%
  filter(ratio2 !="1-10") %>%
  ggplot()+
  # geom_jitter(data=dat_plantmutant_forplotting,
             # aes(x=ratio2, y=Alive), cex=0.1, height=0.2, width=0.2)+
  geom_point(data=dat_plantmutant_forplotting
             , aes(x=ratio2, y=Alive, col=plant, group=plant)
             , show.legend = TRUE, position=position_jitterdodge(jitter.height=0.1, jitter.width=0.5, dodge=0.5)
             , cex=0.1) +
  geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=plant, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.5), size=0) +
  geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=plant, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.5), size=0) +
  facet_grid(protect ~ path, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
  # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
  # scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
  scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
  scale_linewidth_manual(values=c(2, 0.5)) +
  labs(col="Plant genotype",  linewidth="Model estimates") +
  ylab("Probability of healthy plant")+ 
  xlab("Protective : Pathogen ratio")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_bw() 
gg_brms_plantmutant_simp_wdat
ggsave(filename="05_plant_mutant/gg_brms_plantmutant_simp_wdat.png", gg_brms_plantmutant_simp_wdat, height=4, width=5)


gg_brms_plantmutant_simp_wdat_h <- alldat_withpred %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10", 
                         ifelse(factor(round(path_cells_log))==6, "1-1",
                                ifelse(factor(round(path_cells_log))==5, "10-1",
                                       ifelse(factor(round(path_cells_log))==0, "No\nN2C3", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No\nN2C3", "1-10","1-1","10-1"))) %>%
  filter(ratio2 !="1-10") %>%
  ggplot()+
  # geom_jitter(data=dat_plantmutant_forplotting,
  # aes(x=ratio2, y=Alive), cex=0.1, height=0.2, width=0.2)+
  geom_point(data=dat_plantmutant_forplotting
             , aes(x=ratio2, y=Alive, col=plant, group=plant)
             , show.legend = TRUE, position=position_jitterdodge(jitter.height=0.1, jitter.width=0.5, dodge=0.6)
             , cex=0.1) +
  geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=plant, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  geom_pointrange(aes(x=ratio2, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=plant, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  facet_grid(.~protect, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
  # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
  # scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
  scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
  scale_linewidth_manual(values=c(2, 0.5)) +
  labs(col="Plant genotype",  linewidth="Model estimates") +
  ylab("Probability of healthy plant")+ 
  xlab("Protective : Pathogen ratio")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_bw()+
  theme(legend.position = "bottom") +
  guides(lwd=guide_legend(nrow=2,byrow=TRUE))
gg_brms_plantmutant_simp_wdat_h
ggsave(filename="05_plant_mutant/gg_brms_plantmutant_simp_wdat_h.png", gg_brms_plantmutant_simp_wdat_h, height=4, width=7)

gg_plantmutant_v2 <- alldat_withpred %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1:10", 
                         ifelse(factor(round(path_cells_log))==6, "1:1",
                                ifelse(factor(round(path_cells_log))==5, "1:0.1",
                                       ifelse(factor(round(path_cells_log))==0, "1:0 (no N2C3)", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("1:0 (no N2C3)", "1:10","1:1","1:0.1"))) %>%
  filter(ratio2 !="1:10") %>%
  ggplot()+
  # geom_jitter(data=dat_plantmutant_forplotting,
  # aes(x=ratio2, y=Alive), cex=0.1, height=0.2, width=0.2)+
  geom_point(data=dat_plantmutant_forplotting %>%
               mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1:10", 
                                      ifelse(factor(round(path_cells_log))==6, "1:1",
                                             ifelse(factor(round(path_cells_log))==5, "1:0.1",
                                                    ifelse(factor(round(path_cells_log))==0, "1:0 (no N2C3)", NA))))) %>%
               mutate(ratio2 = factor(ratio2, levels=c("1:0 (no N2C3)", "1:10","1:1","1:0.1"))) %>% filter(ratio2 !="1:10") 
             , aes(x=plant, y=Alive, col=ratio2, group=ratio2)
             , show.legend = TRUE, position=position_jitterdodge(jitter.height=0.1, jitter.width=0.5, dodge=0.6)
             , cex=0.1) +
  geom_point(aes(x=plant, y=`med.50%`, col=ratio2), size=3, position = position_dodge2(width=0.6), size=0) +
  geom_pointrange(aes(x=plant, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=ratio2, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  geom_pointrange(aes(x=plant, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=ratio2, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  facet_grid(.~protect, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="No microbiota", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
  # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
  # scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
  # scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
  scale_color_manual(values=c("1:0 (no N2C3)"="darkgrey", "1:1"="darkorange", "1:0.1"="gold2")) +
  scale_linewidth_manual(values=c(2, 0.5)) +
  labs(col="Microbiota:Pathogen\ninoculation ratio",  linewidth="Model estimates") +
  ylab("Probability of healthy plant")+ 
  xlab("Plant genotype")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_bw()
gg_plantmutant_v2
ggsave("05_plant_mutant/gg_plantmutant_v2.png",gg_plantmutant_v2, height=4, width=8)

###### SHORTCUT TO STATS #####
alldat_withpred %>%
  filter(protect %in% c("WCS365","CHAO"), ratio %in% c("1-1"), plant%in% c("col0","bbc")) %>%
  arrange(protect)
alldat_withpred %>%
  filter(protect %in% c("WCS365","CHAO"), ratio %in% c("10-1"), plant=="bbc") %>%
  arrange(protect)
alldat_withpred %>%
  filter(protect %in% c("WCS365","CHAO"), path=="MOCK",  plant=="bbc") %>%
  arrange(protect)

alldat_withpred %>%
  filter(protect %in% c("WCS365","CHAO"), path=="N2C3",  ratio %in% c("1-1"), plant%in% c("col0","bik1")) %>%
  arrange(protect)

gg_plantmutant_v3 <- alldat_withpred %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1:10 (commensal:N2C3)", 
                         ifelse(factor(round(path_cells_log))==6, "1:1 (commensal:N2C3)",
                                ifelse(factor(round(path_cells_log))==5, "1:0.1 (commensal:N2C3)",
                                       ifelse(factor(round(path_cells_log))==0, "No\nN2C3", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No\nN2C3", "1:10 (commensal:N2C3)","1:1 (commensal:N2C3)","1:0.1 (commensal:N2C3)"))) %>%
  filter(ratio2 !="1:10 (commensal:N2C3)") %>%
  filter(protect!="MOCK") %>%
  ggplot()+
  # geom_jitter(data=dat_plantmutant_forplotting,
  # aes(x=ratio2, y=Alive), cex=0.1, height=0.2, width=0.2)+
  geom_point(data=dat_plantmutant_forplotting%>%filter(protect!="MOCK") 
             , aes(x=plant, y=Alive, col=ratio2, group=ratio2)
             , show.legend = TRUE, position=position_jitterdodge(jitter.height=0.1, jitter.width=0.5, dodge=0.6)
             , cex=0.1) +
  geom_point(aes(x=plant, y=`med.50%`, col=ratio2), size=3, position = position_dodge2(width=0.6)) +
  geom_pointrange(aes(x=plant, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=ratio2, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  geom_pointrange(aes(x=plant, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=ratio2, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  facet_grid(.~protect, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
  # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
  # scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
  # scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
  scale_color_manual(values=c("No\nN2C3"="black", "1:1 (commensal:N2C3)"="darkorange", "1:0.1 (commensal:N2C3)"="gold2")) +
  scale_linewidth_manual(values=c(2, 0.5)) +
  labs(col="Protective:Pathogen\ninoculation ratio",  linewidth="Model estimates") +
  ylab("Probability of healthy plant")+ 
  xlab("Plant genotype")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_bw()
gg_plantmutant_v3
ggsave("05_plant_mutant/gg_plantmutant_v3.png",gg_plantmutant_v3, height=4, width=7)


gg_plantmutant_wcs <- alldat_withpred %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10 (commensal:N2C3)", 
                         ifelse(factor(round(path_cells_log))==6, "1-1 (commensal:N2C3)",
                                ifelse(factor(round(path_cells_log))==5, "10-1 (commensal:N2C3)",
                                       ifelse(factor(round(path_cells_log))==0, "No\nN2C3", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No\nN2C3", "1-10 (commensal:N2C3)","1-1 (commensal:N2C3)","10-1 (commensal:N2C3)"))) %>%
  filter(ratio2 !="1-10 (commensal:N2C3)") %>%
  filter(protect %in% c("MOCK", "WCS365") )%>%
  ggplot()+
  # geom_jitter(data=dat_plantmutant_forplotting,
  # aes(x=ratio2, y=Alive), cex=0.1, height=0.2, width=0.2)+
  geom_point(data=dat_plantmutant_forplotting %>% filter(protect %in% c("MOCK", "WCS365") )
             , aes(x=plant, y=Alive, col=ratio2, group=ratio2)
             , show.legend = TRUE, position=position_jitterdodge(jitter.height=0.1, jitter.width=0.5, dodge=0.6)
             , cex=0.1) +
  geom_point(aes(x=plant, y=`med.50%`, col=ratio2), size=3, position = position_dodge2(width=0.6)) +
  geom_pointrange(aes(x=plant, y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=ratio2, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  geom_pointrange(aes(x=plant, y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=ratio2, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  facet_grid(.~protect, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
  # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
  # scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
  # scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
  scale_color_manual(values=c("No\nN2C3"="black", "1-1 (commensal:N2C3)"="darkorange", "10-1 (commensal:N2C3)"="gold2")) +
  scale_linewidth_manual(values=c(2, 0.5)) +
  labs(col="Protective:Pathogen\ninoculation ratio",  linewidth="Model estimates") +
  ylab("Probability of healthy plant")+ 
  xlab("Plant genotype")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_bw()
gg_plantmutant_wcs
ggsave("05_plant_mutant/gg_plantmutant_wcs.png",gg_plantmutant_wcs, height=4, width=6)


alldat_withpred %>%
  mutate(ratio2 = ifelse(factor(round(path_cells_log)) == 7, "1-10", 
                         ifelse(factor(round(path_cells_log))==6, "1-1",
                                ifelse(factor(round(path_cells_log))==5, "10-1",
                                       ifelse(factor(round(path_cells_log))==0, "No\nN2C3", NA))))) %>%
  mutate(ratio2 = factor(ratio2, levels=c("No\nN2C3", "1-10","1-1","10-1"))) %>%
  filter(ratio2 !="1-10") %>%
  ggplot()+
  # geom_jitter(data=dat_plantmutant_forplotting,
  # aes(x=ratio2, y=Alive), cex=0.1, height=0.2, width=0.2)+
  geom_point(data=dat_plantmutant_forplotting
             , aes(x=factor(path_od), y=Alive, col=protect, group=factor(path_od))
             , show.legend = TRUE, position=position_jitterdodge(jitter.height=0.1, jitter.width=0.5, dodge=0.6)
             , cex=0.1) +
  geom_point(aes(x=factor(path_od), y=`med.50%`, col=protect), size=3, position = position_dodge2(width=0.6)) +
  geom_pointrange(aes(x=factor(path_od), y=`med.50%`, ymin=`lwr.2.5%`, ymax = `upr.97.5%` , col=protect, lwd="95% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  geom_pointrange(aes(x=factor(path_od), y=`med.50%`, ymin=`lwr5.5%`, ymax = `upr95.95%` , col=protect, lwd="90% prediction\ninterval"), position = position_dodge2(width=0.6), size=0) +
  facet_grid(.~plant, scales="free_x", space = "free_x", labeller = labeller(protect=c(MOCK="MOCK", WCS365="WCS365",CHAO="CHA0", PF5 = "Pf5")))+
  # scale_color_manual(values=c(MOCK="black", N2C3="orange"))+
  # scale_fill_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) + 
  # scale_color_manual(values=c(col0="black", bbc="magenta", bik1="salmon")) +
  # scale_color_manual(values=c("No\nN2C3"="black", "1-1"="darkorange", "10-1"="gold2")) +
  scale_linewidth_manual(values=c(2, 0.5)) +
  labs(col="Inoculation ratio",  linewidth="Model estimates") +
  ylab("Probability of healthy plant")+ 
  xlab("Protective : Pathogen ratio")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_bw()+
  theme(legend.position = "bottom") +
  guides(lwd=guide_legend(nrow=2,byrow=TRUE))


sink("05_plant_mutant/stats_and_numbers.txt") 
print("actualy raw data")
dat_plantmutant_summarise %>% as.data.frame()

alldat_withpred %>%
  select(plant, protect, path, ratio, path_cells_log, `lwr.2.5%`, `upr.97.5%`)

brm_alldat_bin

brm_alldat_mockonly_bin

sink()

##### Plot experimental design for sanity check ########
#re-arranging plate IDs
exp_plate <- dat_plantmutant %>%
  select(experiment, plate, plant) %>% distinct() %>%
  group_by(experiment, plant) %>% mutate(platenum=rank(plate)) %>% ungroup() %>%
  unite(plant, platenum, col="plant_plate", remove=FALSE)

## Protective setup
# gg_plates_plantpath <- dat_plantmutant %>%
#   left_join(exp_plate) %>%
#   mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
#   ggplot(aes(x=col, y=row)) +
#   geom_tile(aes(fill=path, alpha=path_cells_log)) +
#   # geom_tile(aes(fill=path, alpha=path_cells_log)) +
#   # geom_text(aes(label=ratio))+
#   # facet_grid(experiment~platenum) +
#   facet_grid(platenum~plant) +
#   # scale_fill_manual(values=c(col0="darkgreen", bbc="gold", bik1="goldenrod"))
#   scale_fill_manual(values=c(MOCK="grey", N2C3="darkorange"))
# 
#   # scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", CHAO="green", CH267="goldenrod", PF5="purple"))
# gg_plates_plantpath
# ggsave("05_plant_mutant/gg_plates_plantpath.png",gg_plates_plantpath, height=5, width=10 )

gg_plates_plantprotect <- dat_plantmutant %>%
  mutate(protect_cells_log = ifelse(protect=="MOCK", 7, protect_cells_log),
         path_cells_log = ifelse(path=="MOCK", 7, path_cells_log)) %>%
  left_join(exp_plate) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot(aes(x=col, y=row)) +
  geom_tile(aes(fill=protect)) +
  # geom_point(aes(col=path, alpha=factor(round(path_cells_log)), cex=factor(round(path_cells_log)))) +
  geom_point(aes(col=path, alpha=factor(round(path_cells_log)))) +
  # geom_text(aes(label=ratio))+
  # facet_grid(experiment~platenum) +
  facet_grid(platenum~plant, labeller=labeller(platenum=c(`1`="5-day-old\nplant",`2`="5-day-old\nplant"))) +
  scale_fill_manual(values=c(MOCK="white", WCS365="darkgreen", CHAO="green", CH267="goldenrod", PF5="purple")) +
  scale_color_manual(values=c(MOCK="grey", N2C3="red")) +
  # scale_size_manual(values=c(0.1, 0.5, 1, 2))+
  # theme_bw() +
  labs(col="Pathogen treatment", fill="Protective strain", alpha = "Pathogen cells\n(log10)", size = "Pathogen cells\n(log10)") +
  ylab("Row") + xlab("Column")
gg_plates_plantprotect
ggsave("05_plant_mutant/gg_plates_plantprotect.png",gg_plates_plantprotect, height=5, width=10 )




