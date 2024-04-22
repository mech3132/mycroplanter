#!bin/bash
library(brms)
library(tidyverse)
library(cowplot)

load("02_compile_data/allCFU_data_processed.RData")

dir.create("04_cfu_compare")

### Brms ###
allCFU_filt <- allCFU_data_processed %>%
  filter(!is.na(Treatment), !is.na(WN_FC_fluor_offset), !is.na(WN_FC_inocCFU) ) %>%
  filter(experiment == "2024-02-14") %>%
  filter(transfer_h>0) %>%
  mutate(Treatment = factor(Treatment, levels=c("Simultaneous inoculation", "WCS365 first", "N2C3 first")))


brm_prior_vs_inoc <- brm(WN_FC_fluor_offset ~ Treatment*WN_FC_inocCFU, data=allCFU_filt
                          , seed=52093
                         , iter=4000
                         , file="04_cfu_compare/brm_prior_vs_inoc")
brm_prior_vs_inoc

dat_wn_pred_plant <- allCFU_data_processed %>%
  filter(!is.na(Treatment), !is.na(WN_FC_fluor_offset), !is.na(WN_FC_inocCFU), !is.na(Alive)) %>%
  filter(experiment == "2024-02-14") %>%
  filter(transfer_h>0) %>%
  mutate(Treatment = factor(Treatment, levels=c("Simultaneous inoculation", "WCS365 first", "N2C3 first")))


brm_wn_plant <- brm(Alive ~ Treatment*WN_FC_fluor_offset, data=dat_wn_pred_plant
                         , seed=52093
                    , family="bernoulli"
                         , iter=4000
                         , file="04_cfu_compare/brm_wn_plant")
brm_wn_plant

sink("04_cfu_compare/brms_stats.txt")
brm_prior_vs_inoc
print("")
brm_wn_plant
sink()

#### Plot with brms ####
newX <- seq(min(allCFU_filt$WN_FC_inocCFU), max(allCFU_filt$WN_FC_inocCFU), length.out=20)
newX_pred <- allCFU_filt %>% select(Treatment) %>%
  mutate(newDat=TRUE) %>%
  full_join(data.frame(newDat=TRUE, WN_FC_inocCFU=newX), relationship = "many-to-many")

dat_prior_vs_inoc_withpred <- posterior_predict(brm_prior_vs_inoc, newdata = newX_pred, transform=TRUE) %>%
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
  bind_cols(newX_pred) 
  # distinct() %>%
  # left_join(allCFU_filt)

newX_plant <- seq(min(dat_wn_pred_plant$WN_FC_fluor_offset), max(dat_wn_pred_plant$WN_FC_fluor_offset), length.out=20)
newX_plantpred <- dat_wn_pred_plant %>% select(Treatment) %>%
  mutate(newDat=TRUE) %>%
  full_join(data.frame(newDat=TRUE, WN_FC_fluor_offset=newX_plant), relationship = "many-to-many")
dat_wn_plant_wpred <- posterior_linpred(brm_wn_plant, newdata = newX_plantpred, transform=TRUE) %>%
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
  bind_cols(newX_plantpred) 

##### Plot with brms pred #####
gg_cfu_priorityvsratio_wbrm <- allCFU_filt %>%
  ggplot() +
  geom_jitter(aes(x=WN_FC_inocCFU, y=WN_FC_fluor_offset, group=Treatment, col=Treatment,pch=factor(transfer_h)), height=0.1, width=0.1)  +
  geom_ribbon(data=dat_prior_vs_inoc_withpred, mapping = aes(x=WN_FC_inocCFU, ymin=Q2.5, ymax =Q97.5, group=Treatment, fill=Treatment), col=NA, alpha=0.25) +
  geom_line(data=dat_prior_vs_inoc_withpred, mapping = aes(x=WN_FC_inocCFU,y=pred_mean , group=Treatment, col=Treatment)) +
  scale_color_manual(values=c(`N2C3 first` = "darkorange", `WCS365 first` = "darkseagreen4", `Simultaneous inoculation` = "black")) +
  scale_fill_manual(values=c(`N2C3 first` = "darkorange", `WCS365 first` = "darkseagreen4", `Simultaneous inoculation` = "black")) +
  scale_shape_manual(values=c(`0` = 4, `6` = 19, `24` = 21)) +
  labs(shape = "Inoculation delay (h)", col= "Order of inoculation", fill="Order of inoculation") +
  xlab("WCS365:N2C3 fold-change at inoculation\n(CFU counts)") +
  ylab("WCS365:N2C3 fold-change in final community\n(Fluoresence-based estimates)") +
  theme_bw()
# geom_smooth(method="lm")
# facet_grid(.~experiment)
gg_cfu_priorityvsratio_wbrm
ggsave(filename = "04_cfu_compare/gg_cfu_priorityvsratio_wbrm.png", gg_cfu_priorityvsratio_wbrm, height=4.5, width=6.5)

# 
# gg_cfu_priorityvsratio_6only <- allCFU_data_processed %>%
#   filter(experiment == "2024-02-14") %>%
#   filter(transfer_h==6) %>%
#   ggplot(aes(x=WN_FC_inocCFU, y=WN_FC_fluor_offset, group=Treatment, col=Treatment)) +
#   geom_jitter(aes(), height=0.1, width=0.1)  +
#   # geom_smooth(span=50) +
#   geom_smooth(method="lm") +
#   scale_color_manual(values=c(`N2C3 first` = "darkorange", `WCS365 first` = "darkseagreen4", `Simultaneous inoculation` = "black")) +
#   # scale_shape_manual(values=c(`0` = 4, `6` = 19, `24` = 21)) +
#   labs(shape = "Inoculation delay (h)", col= "Order of inoculation") +
#   xlab("WCS365:N2C3 fold-change at inoculation\n(CFU counts)") +
#   ylab("WCS365:N2C3 fold-change in final community\n(Fluoresence-based estimates)") +
#   theme_bw() +
#   ylim(-5, 8)
# # geom_smooth(method="lm")
# # facet_grid(.~experiment)
# gg_cfu_priorityvsratio_6only
# ggsave(filename = "02_compile_data/gg_cfu_priorityvsratio_6only.png", gg_cfu_priorityvsratio_6only, height=4.5, width=6.5)

## Try plant health now
gg_planthealth_fluor <- dat_wn_pred_plant %>%
  ggplot() +
  geom_jitter(aes(x=WN_FC_fluor_offset, y=Alive, col=Treatment, group=Treatment), width=0, height=0.1, show.legend = FALSE) +
  geom_ribbon(data=dat_wn_plant_wpred, mapping = aes(x=WN_FC_fluor_offset, ymin=Q2.5, ymax =Q97.5, group=Treatment, fill=Treatment), col=NA, alpha=0.25, show.legend = FALSE) +
  geom_line(data=dat_wn_plant_wpred, mapping = aes(x=WN_FC_fluor_offset,y=pred_mean , group=Treatment, col=Treatment), show.legend = FALSE) +
  # geom_smooth(method="glm", method.args = c(family="binomial")) +
  scale_color_manual(values=c(`N2C3 first` = "darkorange", `WCS365 first` = "darkseagreen4", `Simultaneous inoculation` = "black")) +
  scale_fill_manual(values=c(`N2C3 first` = "darkorange", `WCS365 first` = "darkseagreen4", `Simultaneous inoculation` = "black")) +
  xlab("WCS365:N2C3\n(log2 fold change fluorescence)") + ylab("Probability of healthy plant") +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  theme_bw()
gg_planthealth_fluor
ggsave("04_cfu_compare/gg_planthealth_fluor.png", gg_planthealth_fluor, height=4, width=6)

### Multipanel

gg_panelled <- plot_grid(gg_cfu_priorityvsratio,gg_planthealth_fluor, align="h",rel_widths = c(2,1.2))
gg_panelled
ggsave("04_cfu_compare/gg_panelled.png", gg_panelled, width=10, height=4.25)
