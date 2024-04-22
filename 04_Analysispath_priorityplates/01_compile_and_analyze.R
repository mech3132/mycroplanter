#!bin/bash
library(tidyverse)

dir.create("01_compile_and_analyze")
dat1 <- read.csv("00_raw_data/process_athal_2022-04-19.csv")
dat1a <- read.csv("00_raw_data/2022-08-07_platecount_invivo_comp_lacz.csv")
dat1b <- read.csv("00_raw_data/2022-08-07_platecount_noplantcom_lacz.csv")

meta1a <- read.csv("00_meta/inoc_plant_plates.csv")
meta1b <- read.csv("00_meta/inoc_dot_plates.csv")

dat_plant <- full_join(dat1, meta1a) %>% 
  full_join(dat1a) %>%
  unite(strain1, strain2, col="strains", remove=FALSE, sep="+") %>%
  rename(dilution_plant = dilution, white_plant=white, blue_plant =blue)

dat_dots <- full_join(dat1b, meta1b) %>% 
  filter(!(is.na(white)&is.na(blue))) %>%
  unite(ratio1, ratio2, col = "ratio", remove=FALSE, sep=":") %>%
  unite(strain1, strain2, col="strains", remove=FALSE, sep="+") %>%
  unite(strains,ratio,col=dotplate,remove=FALSE) %>%
  select(dotplate, row, col, ratio1, ratio2, ratio, dilution, white, blue, strain1, strain2, strains, date_inoc, date_process) %>%
  rename(date_process_dot = date_process, dilution_dot = dilution, white_dot=white, blue_dot =blue)

alldat <- full_join(dat_plant, dat_dots)

# For lacz counts only
invivo_dat <- alldat %>% filter(!is.na(white_plant)) %>%
  mutate(dotplate=NA) %>%
  select(plate, dotplate, strains, ratio, white_plant, blue_plant,delay)  %>%
  rename(white=white_plant, blue=blue_plant) %>%
  mutate(type="invivo")
invitro_dat <- alldat %>% filter(!is.na(white_dot)) %>%
  select(plate, dotplate, strains, ratio, white_dot, blue_dot,delay) %>%
  rename(white=white_dot, blue=blue_dot) %>%
  mutate(type="invitro")

alldat_lacz_long <- rbind(invivo_dat, invitro_dat) %>%
  mutate(propwhite = white/(white+blue), propblue = blue/(white+blue)) %>%
  filter(!is.na(propwhite)) %>%
  unite(strains, ratio, delay, col="treatment", remove=FALSE) %>% 
  mutate(delay=ifelse(is.na(delay),"-",delay)) %>%# Manual change
  mutate(ratio=ifelse(is.na(ratio),"1:1",ratio)) # Manual change

alldat_lacz_summarize <- alldat_lacz_long %>%
  group_by(type, treatment, strains, ratio, delay) %>%
  summarise(white = mean(propwhite), sdpropwhite = sd(propwhite), basewhite = mean(propwhite)) %>%
  ungroup() %>%
  mutate(blue = 1-white) %>%
  pivot_longer(c(white, blue), names_to="colony", values_to = "prop") 

##### Plant growth ####

alldat %>% 
  filter(!is.na(treatment), treatment!="Mock") %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  rowwise() %>%
  mutate(Treatment = treatment) %>%
  # mutate(Treatment = ifelse(ratio=="1:1","W+N (1:1)", Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  mutate(Treatment = gsub("L","", Treatment)) %>%
  mutate(Treatment = gsub("N[+]N","N", Treatment)) %>%
  mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
                                                ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=weight_mg)) +
  geom_jitter(aes(x=Treatment, y=weight_mg)) +
  facet_grid(.~Timing, drop=TRUE, scales="free") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Plant weight (mg)") + xlab("Treatments\n6uL @ 0.001OD, 1:1 for concurrent\n6uL @ 0.0005OD each for delayed")

gg_plantgrowth <- alldat %>% 
  filter(!is.na(treatment), treatment!="Mock") %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  rowwise() %>%
  mutate(Treatment = treatment) %>%
  # mutate(Treatment = ifelse(ratio=="1:1","W+N (1:1)", Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  mutate(Treatment = gsub("L","", Treatment)) %>%
  mutate(Treatment = gsub("N[+]N","N", Treatment)) %>%
  mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
                                                ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=weight_mg)) +
  geom_jitter(aes(x=Treatment, y=weight_mg)) +
  facet_grid(.~Timing, drop=TRUE, scales="free") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Plant weight (mg)") + xlab("Treatments\n6uL @ 0.001OD, 1:1 for concurrent\n6uL @ 0.0005OD each for delayed")
gg_plantgrowth
ggsave(filename="01_compile_and_analyze/gg_plantgrowth.png",gg_plantgrowth,height=5, width=10)

##### lacZ ####
alldat_lacz_summarize %>% 
  rowwise() %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  mutate(type=ifelse(type=="invitro","in vitro", "in vivo")) %>%
  mutate(type=factor(type, levels=c("in vivo","in vitro"))) %>%
  mutate(Treatment=strains) %>%
  mutate(Treatment = ifelse(delay=="++",gsub("[+]","++",Treatment), Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  # mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
  #                                               ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  ungroup() %>%
  ggplot() + 
  geom_bar(aes(x=Treatment, y=prop, fill=colony, group=type), stat="identity", position="stack") +
  geom_errorbar(aes(x=Treatment, ymin=basewhite-sdpropwhite, max = basewhite+sdpropwhite)) +
  facet_grid(type~Timing, drop=TRUE,scales = "free")+
  scale_fill_manual(values=c(blue="blue", white="grey"))

gg_lacz <- alldat_lacz_summarize %>% 
  rowwise() %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  mutate(type=ifelse(type=="invitro","in vitro", "in vivo")) %>%
  mutate(type=factor(type, levels=c("in vivo","in vitro"))) %>%
  mutate(Treatment=strains) %>%
  mutate(Treatment = ifelse(delay=="++",gsub("[+]","++",Treatment), Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  # mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
  #                                               ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  ungroup() %>%
  filter(Treatment%in%c("W+NL","WL+N","NSSPP+NL","W++NL","WL++N","N++WL","NL++W")) %>%
  mutate(Treatment=factor(Treatment, levels=c("W+NL","WL+N","NSSPP+NL","W++NL","WL++N","N++WL","NL++W"))) %>%
  ggplot() + 
  geom_bar(aes(x=Treatment, y=prop, fill=colony, group=type), stat="identity", position="stack") +
  geom_errorbar(aes(x=Treatment, ymin=basewhite-sdpropwhite, max = basewhite+sdpropwhite)) +
  facet_grid(type~Timing, drop=TRUE,scales = "free")+
  scale_fill_manual(values=c(blue="blue", white="grey"))+
  xlab("Treatment") + ylab("Proportion of blue/white colonies")
gg_lacz
ggsave(filename="01_compile_and_analyze/gg_lacz_barplots.png",gg_lacz,height=5, width=8)




##### Save ####


write.table(alldat, file = "01_compile_and_analyze/merged_2022-04-04.txt", quote=FALSE, row.names = FALSE, sep="\t")

###### Only WN #####


##### Plant growth ####

gg_plantgrowth_justWN <- alldat %>% 
  filter(!is.na(treatment), treatment!="Mock") %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  rowwise() %>%
  mutate(Treatment = treatment) %>%
  # mutate(Treatment = ifelse(ratio=="1:1","W+N (1:1)", Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  mutate(Treatment = gsub("L","", Treatment)) %>%
  mutate(Treatment = gsub("N[+]N","N", Treatment)) %>%
  mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
                                                ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  filter(Timing == "Delayed inoculation") %>%
  mutate(Treatment = ifelse(Treatment == "W++0", "WCS365\nthen mock"
                            , ifelse(Treatment == "N++0", "N2C3\nthen mock"
                                     , ifelse(Treatment == "0++N", "Mock\nthen N2C3"
                                              , ifelse(Treatment == "N++W", "N2C3\nthen WCS365"
                                                       , ifelse(Treatment == "W++N", "WCS365\nthen N2C3", NA)))))) %>%
  mutate(Treatment = factor(Treatment
                            , levels=c("WCS365\nthen mock","N2C3\nthen mock","Mock\nthen N2C3", "WCS365\nthen N2C3", "N2C3\nthen WCS365"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=weight_mg)) +
  geom_jitter(aes(x=Treatment, y=weight_mg)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Plant weight (mg)") + xlab("Treatments\n6uL @ 0.001OD, 1:1 for concurrent\n6uL @ 0.0005OD each for delayed")
gg_plantgrowth_justWN
ggsave(filename="01_compile_and_analyze/gg_plantgrowth+WNONLY.png",gg_plantgrowth_justWN,height=3.5, width=5)


unique(alldat$treatment)
gg_plantgrowth_justWN_timing <- alldat %>% 
  # filter(!is.na(treatment), treatment!="Mock") %>%
  filter(!is.na(treatment)) %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  mutate(Timing = ifelse(treatment=="Mock", "Concurrent inoculation", Timing)) %>%
  rowwise() %>%
  mutate(Treatment = treatment) %>%
  # mutate(Treatment = ifelse(ratio=="1:1","W+N (1:1)", Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  mutate(Treatment = gsub("L","", Treatment)) %>%
  mutate(Treatment = gsub("N[+]N","N", Treatment)) %>%
  mutate(Treatment = factor(Treatment, levels=c("Mock","W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
                                                ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  # filter(Timing == "Delayed inoculation") %>%
  filter(Treatment %in% c("W","N","W+N","W+NL","W++0","W++N","N++W","N++0","0++N","Mock")) %>%
  mutate(Treatment = ifelse(Treatment == "W++0", "WCS365\nthen mock"
                            , ifelse(Treatment == "N++0", "N2C3\nthen mock"
                                     , ifelse(Treatment == "0++N", "Mock\nthen N2C3"
                                              , ifelse(Treatment == "N++W", "N2C3\nthen WCS365"
                                                       , ifelse(Treatment == "W++N", "WCS365\nthen N2C3", 
                                                                ifelse(Treatment == "W","WCS365",
                                                                       ifelse(Treatment == "N", "N2C3",
                                                                              ifelse(Treatment=="Mock", "Mock\n(MgSO4)",
                                                                                     ifelse(Treatment == "W+N", "WCS365+N2C3")))))))))) %>%
  mutate(Treatment = factor(Treatment
                            , levels=c("WCS365\nthen mock","N2C3\nthen mock"
                                       ,"Mock\nthen N2C3", "WCS365\nthen N2C3"
                                       , "N2C3\nthen WCS365", 
                                       "Mock\n(MgSO4)","WCS365","N2C3","WCS365+N2C3"))) %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=weight_mg)) +
  geom_jitter(aes(x=Treatment, y=weight_mg)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Plant weight (mg)") + xlab("Treatments\n6uL @ 0.001OD, 1:1 for concurrent\n6uL @ 0.0005OD each for delayed")+
  facet_grid(.~Timing, scales="free_x", drop=TRUE)
gg_plantgrowth_justWN_timing
ggsave(filename="01_compile_and_analyze/gg_plantgrowth_WNONLY_wtiming.png",gg_plantgrowth_justWN_timing,height=4, width=6)

##### lacZ ####
gg_priority_WN <- alldat_lacz_summarize %>% 
  rowwise() %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  mutate(type=ifelse(type=="invitro","in vitro", "in vivo")) %>%
  mutate(type=factor(type, levels=c("in vivo","in vitro"))) %>%
  mutate(Treatment=strains) %>%
  mutate(Treatment = ifelse(delay=="++",gsub("[+]","++",Treatment), Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  # mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
  #                                               ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  ungroup() %>%
  filter(Timing == "Delayed inoculation", type == "in vivo") %>%
  select(Treatment, prop, colony) %>%
  separate(Treatment, into = c("First","Second"), sep="[+][+]", remove=FALSE) %>% 
  mutate(Isolate = ifelse(colony=="white"& (First =="N" | Second =="N"), "N2C3",
                          ifelse(colony=="blue" & (First == "NL"|Second == "NL"), "N2C3",
                                 ifelse(colony=="white" & (First=="W" | Second == "W"), "WCS365",
                                        ifelse(colony=="blue" & (First =="WL" | Second =="WL"), "WCS365",NA))))) %>%
  mutate(Treatment = gsub("L","",Treatment)) %>%
  group_by(Isolate, Treatment) %>%
  summarise(meanprop = mean(prop), minprop = min(prop), maxprop = max(prop)) %>%
  ungroup() %>%
  mutate(minprop = ifelse(Isolate=="WCS365", NA, minprop),
         maxprop = ifelse(Isolate=="WCS365", NA, maxprop)) %>%
  mutate(Treatment = ifelse(Treatment == "W++N", "WCS365 first", "N2C3 first")) %>%
  ggplot() + 
  geom_bar(aes(x=Treatment, y=meanprop, fill=Isolate), stat="identity", position="stack") +
  geom_errorbar(aes(x=Treatment, ymin=minprop, max = maxprop)) +
  # facet_grid(type~Timing, drop=TRUE,scales = "free")+
  scale_fill_manual(values=c(WCS365="darkgreen", N2C3="orange")) +
  ylab("Proportion of\nprotective/pathogenic colonies")
gg_priority_WN
ggsave(filename="01_compile_and_analyze/gg_lacz_barplots_WNONLY.png",gg_priority_WN,height=3, width=4)

### Barplots with both concurrent and delayed inoculations

gg_lacz_WN <- alldat_lacz_summarize %>% 
  rowwise() %>%
  mutate(Timing=ifelse(delay=="-","Concurrent inoculation","Delayed inoculation")) %>%
  mutate(type=ifelse(type=="invitro","in vitro", "in vivo")) %>%
  mutate(type=factor(type, levels=c("in vivo","in vitro"))) %>%
  mutate(Treatment=strains) %>%
  mutate(Treatment = ifelse(delay=="++",gsub("[+]","++",Treatment), Treatment)) %>%
  mutate(Treatment = gsub("Q","SSPP",gsub("S","SS",gsub("P","PP",Treatment)))) %>%
  # mutate(Treatment = factor(Treatment, levels=c("W","NSSPP","NSS","NPP","N","W+N","NSSPP+N","NSS+N","NPP+N"
  #                                               ,"W++0","N++0","0++N","N++W","W++N"))) %>%
  ungroup() %>%
  filter(Treatment%in%c("W+NL","WL+N","W++NL","WL++N","N++WL","NL++W")) %>%
  mutate(Treatment = gsub("+","\nwith\n",gsub("++","\nthen\n",gsub("L","-LacZ",gsub("N","N2C3",gsub("W","WCS365", Treatment))), fixed=TRUE), fixed=TRUE)) %>%
  # mutate(Treatment=factor(Treatment, levels=c("W+NL","WL+N","W++NL","WL++N","N++WL","NL++W"))) %>%
  ggplot() + 
  geom_bar(aes(x=Treatment, y=prop, fill=colony, group=type), stat="identity", position="stack") +
  geom_errorbar(aes(x=Treatment, ymin=basewhite-sdpropwhite, max = basewhite+sdpropwhite)) +
  facet_grid(type~Timing, drop=TRUE,scales = "free")+
  scale_fill_manual(values=c(blue="blue", white="grey"))+
  xlab("Treatment") + ylab("Proportion of blue/white colonies")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))+
  labs(fill="Colony\ncolor")
gg_lacz_WN
ggsave(filename="01_compile_and_analyze/gg_lacz_barplots_WNonly.png",gg_lacz_WN,height=5, width=8)



















