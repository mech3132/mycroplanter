#!bin/bash

dir.create("03a-1_HealthScore_Assessment")
library("plotwidgets")
library('plyr') # need to load plyr first
library(tidyverse)
library(cowplot)


load("03_EDA_and_adj/allDat_final.RData") 
load("03_EDA_and_adj/plant_colours.RData")
reference <- allDat_final %>%
  # mutate(Healthiness_rgb = ((diff_gr_median*2 + diff_gb_median*1)*(diff_gr_median/10 + log10(all_plant_pixels+1)))) %>%
  # mutate(Healthiness_hsv = (hue)*(2*saturation+log10(all_plant_pixels+1)/100)) %>%
  select(UniqueID, Healthiness_rgb, Healthiness_hsv)
datsurv <- read.csv("03a_HealthScore_assessment/plants_for_survey.txt")
dat <- read.csv("03a_HealthScore_assessment/assessment_data.csv", na.strings = c("NA","","na")) %>%
  filter(!is.na(rank))
planttypes <- read.csv("03a_HealthScore_assessment/plant_types.csv")

###### First, let's check if all sets are "complete"

dat %>%
  select(person, id) %>%table()
# Make person a code, rather than an initial
person_legend <- dat %>%
  select(person) %>% distinct() %>% rename(personname=person) %>%
  mutate(person = paste0("HUMAN", rank(personname)))

# Get actual rankings
reference <- datsurv %>% select(UniqueID, row) %>%
  rename(id=row) %>% left_join(reference)
colnames(dat)
alg_temp_wide <- dat %>% select(set, id) %>% distinct() %>% arrange(set) %>%
  left_join(reference) %>%
  mutate(healthy_hsv = Healthiness_hsv>400, healthy_rgb = Healthiness_rgb>400) %>%
  group_by(set) %>%
  mutate(HealthScore_HSV = rank(Healthiness_hsv), HealthScore_RGB = rank(Healthiness_rgb)) %>%
  ungroup() 
alg_ranks <-   alg_temp_wide %>% select(set, HealthScore_HSV, HealthScore_RGB, id) %>%
  pivot_longer(-c(set,id), names_to="person", values_to="rank")
alg_healthy <- alg_temp_wide %>% select(set, healthy_hsv, healthy_rgb, id) %>%
  rename(HealthScore_HSV = healthy_hsv, HealthScore_RGB = healthy_rgb) %>%
  pivot_longer(-c(set,id), names_to="person", values_to="healthy")
dat_walg <- dat %>% rename(personname = person) %>% left_join(person_legend) %>%
  full_join(full_join(alg_ranks,alg_healthy))

dat_all  <- allDat_final %>%
  select(UniqueID, Healthiness_hsv, Healthiness_rgb) %>%
  right_join(datsurv %>% select(-Healthiness)) %>%
  rename(id=row) %>%
  left_join(dat_walg) %>%
  select(person, set, rank, id, healthy) %>%
  filter(!is.na(person)) %>%
  left_join(planttypes)


######## Look at all rankings #########
factorPlants <- alg_temp_wide %>%
  arrange(Healthiness_hsv) %>%
  filter(!is.na(id)) %>%
  pull(id) %>% unique()
gg_rawranks <- dat_all %>%
  mutate(person=factor(person, levels=unique(c("HealthScore_HSV","HealthScore_RGB", person)))) %>%
  mutate(id=factor(id, levels=factorPlants)) %>%
  # filter(person!="HealthScore_RGB") %>%
  ggplot() +
  geom_tile(aes(x=person, y=factor(id), fill=rank), col="black") +
  facet_grid(.~factor(set), drop=TRUE, scales="free")+
  scale_fill_gradient(low="orange", high="green") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position = "top")+
  xlab("Person") +
    ylab("Random plant ID")+
  labs(fill="Rank") +
  coord_flip() 
gg_rawranks

ggsave("03a-1_HealthScore_Assessment/gg_rawranks.png", gg_rawranks, height=3, width=12)

# dat_all %>%
#   mutate(person=factor(person, levels=c("HealthScore_HSV","HealthScore_RGB",  "AB", "LF", "MYC", "SH","SS", "ZMM"))) %>%
#   mutate(id=factor(id, levels=factorPlants)) %>%
#   filter(group!=2) %>%
#   # filter(person!="HealthScore_RGB") %>%
#   ggplot() +
#   geom_tile(aes(x=person, y=factor(id), fill=rank)) +
#   facet_grid(factor(set)~., drop=TRUE, scales="free")+
#   scale_fill_gradient(low="orange", high="green") +
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
#   xlab("Person") +
#   ylab("Random plant ID")+
#   labs(fill="Rank")

######## Get correlations ######
#TESTING
# dat_all <- dat_all %>%
#   # filter(id!=7713) %>%
#   # filter(person!="HealthScore_RGB") %>%
#   filter(group!=2)

allCorMat <- list()
allCorDF <- data.frame(p1=NA)
for ( s in 1:10) {
  allCorMat[[s]] <- filter(dat_all, set==s) %>%
    select(id, rank, person) %>%
    pivot_wider(names_from=person, values_from=rank) %>%
    select(-id) %>% cor(method="kendall")
  allCorDF <- filter(dat_all, set==s) %>%
    select(id, rank, person) %>%
    pivot_wider(names_from=person, values_from=rank) %>%
    select(-id) %>% cor(method="kendall") %>%
    as.data.frame() %>% rownames_to_column(var="p1") %>%
    pivot_longer(-p1, names_to="p2", values_to="tau") %>%
    mutate(set=s) %>%
    full_join(allCorDF)
}
allCorMat
allCorDF <- filter(allCorDF, !is.na(p1)) %>% 
  mutate(p1 = factor(p1, levels=unique(c(c("HealthScore_RGB", "HealthScore_HSV", p1))))) %>%
  mutate(p2 = factor(p2, levels=unique(c(c("HealthScore_RGB", "HealthScore_HSV", p2))))) 

# just alg
dat_alg <- dat_all %>% filter(person %in% c("HealthScore_RGB", "HealthScore_HSV")) %>%
  pivot_wider(names_from="person", values_from="rank")
dat_alg_long <- dat_all %>%  filter(person %in% c("HealthScore_RGB", "HealthScore_HSV")) %>%
  rename(rank_alg = rank, algorithm=person)

allCorDF %>%
  filter(p1!=p2) %>%
   ggplot(aes(x=p1, y=tau)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(.~p2, drop=TRUE) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

forHeatmap <- allCorDF %>%
  group_by(p1, p2) %>%
  summarise(meanTau = mean(tau), sdTau = sd(tau)) %>%
  ungroup() 

forHeatmap %>% ggplot() +
  geom_tile(aes(x=p1, y=p2, fill=meanTau))
forHeatmap_df <- forHeatmap %>%
  select(-sdTau) %>%
  pivot_wider(names_from=p2, values_from=meanTau)
forHeatmap_mat <- as.matrix(forHeatmap_df[,-1])
rownames(forHeatmap_mat) <- pull(forHeatmap_df[,1])

png(filename="03a-1_HealthScore_Assessment/heatmap_tau_pairedkendall.png"
    , height=8, width=8,res=600,  units = "in"
    )
heatmap(1-forHeatmap_mat, scale = "none", margins = c(12,12))
dev.off()
# 
# dat_all %>%
#   filter(!person %in% c("HealthScore_HSV","HealthScore_RGB")) %>%
#   group_by(set, id) %>%
#   summarise(meanRank = mean(rank))  %>%
#   mutate(newRank = rank(meanRank)) %>% ungroup() %>%
#   left_join(dat_alg_long) %>%
#   ggplot(aes(x=rank_alg, y=newRank)) +
#   geom_jitter(aes(col=factor(set)), width=0.25, height=0.25) +
#   geom_smooth(method="lm") +
#   facet_grid(.~algorithm)

### What if we do drop-out calculations to see how each one predicts others' average?
allPersons <- dat_all$person %>% unique()

droponeRanks <- data.frame()
for ( p in allPersons ) {
  newMeanRanks <- dat_all %>% filter(person!=p) %>%
    group_by(id, set) %>%
    summarise(droppedMeanRank = mean(rank)) %>% 
    ungroup() %>% 
    left_join(dat_all %>% filter(person==p))
  droponeRanks <- rbind(droponeRanks, newMeanRanks)
}


gg_dropone_fits <- droponeRanks %>%
  mutate(person = factor(person, levels = unique(c("HealthScore_RGB", "HealthScore_HSV", person)))) %>%
  ggplot() + 
  geom_point(aes(x=rank, y=droppedMeanRank)) +
  geom_abline(aes(intercept=0, slope=1), col="red") +
  facet_wrap(.~person) +
  ylab("Mean Rank of all others") + xlab("Rank of single person/metric")
gg_dropone_fits

ggsave(filename = "03a-1_HealthScore_Assessment/gg_dropone_ranks.png", gg_dropone_fits, height=8, width=8)

droponeRanks %>%
  group_by(person) %>%
  summarise(tau = cor.test(rank, droppedMeanRank, method="spearman")$estimate
            ,z = cor.test(rank, droppedMeanRank, method="spearman")$statistic
            ,p = cor.test(rank, droppedMeanRank, method="spearman")$p.value) %>%
  arrange(tau)

droponeRanks %>%
  group_by(person) %>%
  summarise(tau = cor.test(rank, droppedMeanRank, method="kendall")$estimate
            ,z = cor.test(rank, droppedMeanRank, method="kendall")$statistic
            ,p = cor.test(rank, droppedMeanRank, method="kendall")$p.value) %>%
  arrange(tau)

### Are people consistent with themselves? ###

# 7713, 5158, 2603, 60 are consistent ones

healthy_not_agree <- dat_all %>%
  filter(id %in% c(7713, 5158, 2603, 60)) %>%
  rowwise() %>%
  mutate(sign10=ifelse(healthy, 1, -1)) %>% ungroup() %>%
  group_by(person, id) %>%
  summarise(healthy_percent = sum(healthy)/n(), signbigger = sum(sign10)) %>% ungroup() %>%
  mutate(healthy_percent_posneg = ifelse(signbigger <0, -1*healthy_percent, healthy_percent)
         , fullAgreement = ifelse(healthy_percent %in% c(0,1), TRUE, FALSE))

######### Fix the math behind this. Something not right.
healthy_not_agree <- dat_all %>%
  filter(id %in% c(7713, 5158, 2603, 60)) %>%
  rowwise() %>%
  mutate(sign10=ifelse(healthy, 1, -1)) %>% ungroup() %>%
  group_by(person, id) %>%
  summarise(healthy_percent = sum(healthy)/n(), signbigger = sum(sign10), sumTotal=n(), Healthy=sum(healthy), NotHealthy=sum(!healthy)) %>% ungroup() %>%
  mutate( percent_agree = ifelse(healthy_percent< 0.5, (1-healthy_percent)*(-1), 
                                 ifelse(healthy_percent==0, -1, healthy_percent))
         , fullAgreement = ifelse(percent_agree %in% c(-1,1), TRUE, FALSE))

healthy_not_agree %>%
  rename(plant_id = id) %>%
  ggplot() +
  geom_bar(aes(x=factor(plant_id), y=percent_agree, fill=fullAgreement), stat="identity") +
  facet_wrap(.~person) 

gg_agreement <- healthy_not_agree %>%
  select(person, id, Healthy, NotHealthy) %>%
  pivot_longer(c(Healthy, NotHealthy), names_to="Healthy", values_to="count") %>%
  rename(plant_id = id) %>%
  ggplot() +
  geom_bar(aes(x=factor(plant_id), y=count, fill=Healthy), stat="identity") +
  facet_wrap(.~person) +
  xlab("Random Plant ID") + ylab("Count") + labs(fill="Healthy?") +
  scale_fill_manual(values=c(Healthy="darkgreen", NotHealthy="orange")) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_agreement

ggsave("03a-1_HealthScore_Assessment/gg_agreement.png", gg_agreement, 
       width=6, height=5)

######## Remove half coloured ones ##############

######## Get correlations ######
#TESTING
dat_all2 <- dat_all %>%
  # filter(id!=7713) %>%
  # filter(person!="HealthScore_RGB") %>%
  filter(group!=2)

allCorMat2 <- list()
allCorDF2 <- data.frame(p1=NA)
for ( s in 1:10) {
  allCorMat2[[s]] <- filter(dat_all2, set==s) %>%
    select(id, rank, person) %>%
    pivot_wider(names_from=person, values_from=rank) %>%
    select(-id) %>% cor(method="kendall")
  allCorDF2 <- filter(dat_all2, set==s) %>%
    select(id, rank, person) %>%
    pivot_wider(names_from=person, values_from=rank) %>%
    select(-id) %>% cor(method="kendall") %>%
    as.data.frame() %>% rownames_to_column(var="p1") %>%
    pivot_longer(-p1, names_to="p2", values_to="tau") %>%
    mutate(set=s) %>%
    full_join(allCorDF2)
}
allCorMat2
allCorDF2 <- filter(allCorDF2, !is.na(p1)) %>% 
  mutate(p1 = factor(p1, levels=unique(c("HealthScore_RGB", "HealthScore_HSV", p1)))) %>%
  mutate(p2 = factor(p2, levels=unique(c("HealthScore_RGB", "HealthScore_HSV", p2))))

# just alg
# dat_alg2 <- dat_all2 %>% filter(person %in% c("HealthScore_RGB", "HealthScore_HSV")) %>%
  # pivot_wider(names_from="person", values_from="rank")
# dat_alg_long2 <- dat_all2 %>%  filter(person %in% c("HealthScore_RGB", "HealthScore_HSV")) %>%
  # rename(rank_alg = rank, algorithm=person)

allCorDF2 %>%
  filter(p1!=p2) %>%
  ggplot(aes(x=p1, y=tau)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(.~p2, drop=TRUE) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

forHeatmap2 <- allCorDF2 %>%
  group_by(p1, p2) %>%
  summarise(meanTau = mean(tau), sdTau = sd(tau)) %>%
  ungroup() 

forHeatmap2 %>% ggplot() +
  geom_tile(aes(x=p1, y=p2, fill=meanTau))
forHeatmap_df2 <- forHeatmap2 %>%
  select(-sdTau) %>%
  pivot_wider(names_from=p2, values_from=meanTau)
forHeatmap_mat2 <- as.matrix(forHeatmap_df2[,-1])
rownames(forHeatmap_mat2) <- pull(forHeatmap_df2[,1])

png(filename="03a-1_HealthScore_Assessment/heatmap_tau_pairedkendall_nohalfcol.png"
    , height=8, width=8,res=600,  units = "in"
)
heatmap(1-forHeatmap_mat2, scale = "none", margins = c(12,12))
dev.off()
# 
# dat_all2 %>%
#   filter(!person %in% c("HealthScore_HSV","HealthScore_RGB")) %>%
#   group_by(set, id) %>%
#   summarise(meanRank = mean(rank))  %>%
#   mutate(newRank = rank(meanRank)) %>% ungroup() %>%
#   left_join(dat_alg_long) %>%
#   ggplot(aes(x=rank_alg, y=newRank)) +
#   geom_jitter(aes(col=factor(set)), width=0.25, height=0.25) +
#   geom_smooth(method="lm") +
#   facet_grid(.~algorithm)

### What if we do drop-out calculations to see how each one predicts others' average?
allPersons <- dat_all2$person %>% unique()

droponeRanks2 <- data.frame()
for ( p in allPersons ) {
  newMeanRanks <- dat_all %>% filter(person!=p) %>%
    group_by(id, set) %>%
    summarise(droppedMeanRank = mean(rank)) %>% 
    ungroup() %>% 
    left_join(dat_all %>% filter(person==p))
  droponeRanks2 <- rbind(droponeRanks2, newMeanRanks)
}


gg_dropone_fits2 <- droponeRanks2 %>%
  mutate(person = factor(person, levels = unique(c("HealthScore_RGB", "HealthScore_HSV", person)))) %>%
  ggplot() + 
  geom_point(aes(x=rank, y=droppedMeanRank)) +
  geom_abline(aes(intercept=0, slope=1), col="red") +
  facet_wrap(.~person) +
  ylab("Mean Rank of all others") + xlab("Rank of single person/metric")
gg_dropone_fits2

ggsave(filename = "03a-1_HealthScore_Assessment/gg_dropone_ranks_nohalfcol.png", gg_dropone_fits2, height=8, width=8)

droponeRanks2 %>%
  group_by(person) %>%
  summarise(tau = cor.test(rank, droppedMeanRank, method="spearman")$estimate
            ,z = cor.test(rank, droppedMeanRank, method="spearman")$statistic
            ,p = cor.test(rank, droppedMeanRank, method="spearman")$p.value) %>%
  arrange(tau)

droponeRanks2 %>%
  group_by(person) %>%
  summarise(tau = cor.test(rank, droppedMeanRank, method="kendall")$estimate
            ,z = cor.test(rank, droppedMeanRank, method="kendall")$statistic
            ,p = cor.test(rank, droppedMeanRank, method="kendall")$p.value) %>%
  arrange(tau)

######## Find ideal healthy/not healthy division ######
dat_all
plant_hsv_only <- reference %>% select(id, Healthiness_hsv) %>%
  arrange(Healthiness_hsv) %>%
  rownames_to_column(var="rn") %>% filter(!is.na(Healthiness_hsv))
onlyPeople <- dat_all  %>%
  filter(!person %in% c("HealthScore_HSV","HealthScore_RGB"))

contingencyHealthiness <- data.frame(person=NA)
for ( b in 1:nrow(plant_hsv_only)) {
  # b=4
  currentDat <- plant_hsv_only %>%
    mutate(B = as.numeric(rn)>b) %>%
    right_join(onlyPeople) %>%
    select(id, person, Healthiness_hsv,healthy, B)
  curB <- plant_hsv_only %>%
    filter(rn %in% c(b, b+1)) %>%
    select(Healthiness_hsv) %>% pull() %>%mean()
  for ( p in unique(onlyPeople$person)) {
    curdat1 <- currentDat %>% filter(person==p) %>%
      select(healthy,B) %>%
      drop_na()
    dimTemp <- dim(table(curdat1))
    if ( any(!dimTemp>1 ) ) {
      newDat <- data.frame(person=p, b=b, HealthScore_HSV = curB,  stat=NA, pval =NA)
    } else {
      cor(curdat1$healthy, curdat1$B)
      table(curdat1$healthy, curdat1$B)
      curTest <- chisq.test(curdat1$healthy, curdat1$B)
      newDat <- data.frame(person=p, b=b, HealthScore_HSV = curB,  stat=curTest$statistic, pval = curTest$p.value, MSE=mean((curTest$residuals)^2))
      contingencyHealthiness = full_join(contingencyHealthiness, newDat)
    }
    
  }
    
}
contingencyHealthiness <- contingencyHealthiness %>%
  filter(!is.na(person))


contingencyHealthiness %>%
  ggplot(aes(x=HealthScore_HSV, y=stat, group=person, col=person)) +
  geom_point() +
  geom_line()

contingencyHealthiness %>%
  ggplot(aes(x=HealthScore_HSV, y=log10(pval), group=person, col=person)) +
  geom_point() +
  geom_line()


gg_MSE_deadalive <- contingencyHealthiness %>%
  ggplot(aes(x=HealthScore_HSV, y=stat, group=person, col=person)) +
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=FALSE)+
  ylab("Test statistic\n(Chi squared test against human classifications\n at each threshold)") +
  xlab("Health Score threshold")
gg_MSE_deadalive
ggsave(filename = "03a-1_HealthScore_Assessment/gg_MSE_deadalive.png", gg_MSE_deadalive, height=4, width=6)

hist(plant_hsv_only$Healthiness_hsv)
# 
# ######## GG BIPLOT PLANTS ###############
# #### Make custom biplot
# ggbiplot_PLANTS <- function (pcobj, colors="black",outline_col= "black", text_col="black", arrow_col="yellow", choices = 1:2, scale = 1, pc.biplot = TRUE, 
#                              obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
#                              ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
#                              alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
#                              varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
#                              ...) 
# {
#   library(ggplot2)
#   library(plyr)
#   library(scales)
#   library(grid)
#   stopifnot(length(choices) == 2)
#   if (inherits(pcobj, "prcomp")) {
#     nobs.factor <- sqrt(nrow(pcobj$x) - 1)
#     d <- pcobj$sdev
#     u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
#     v <- pcobj$rotation
#   }
#   else if (inherits(pcobj, "princomp")) {
#     nobs.factor <- sqrt(pcobj$n.obs)
#     d <- pcobj$sdev
#     u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
#     v <- pcobj$loadings
#   }
#   else if (inherits(pcobj, "PCA")) {
#     nobs.factor <- sqrt(nrow(pcobj$call$X))
#     d <- unlist(sqrt(pcobj$eig)[1])
#     u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
#     v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
#                                                   1]), FUN = "/")
#   }
#   else if (inherits(pcobj, "lda")) {
#     nobs.factor <- sqrt(pcobj$N)
#     d <- pcobj$svd
#     u <- predict(pcobj)$x/nobs.factor
#     v <- pcobj$scaling
#     d.total <- sum(d^2)
#   }
#   else {
#     stop("Expected a object of class prcomp, princomp, PCA, or lda")
#   }
#   choices <- pmin(choices, ncol(u))
#   df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
#                               FUN = "*"))
#   v <- sweep(v, 2, d^var.scale, FUN = "*")
#   df.v <- as.data.frame(v[, choices])
#   names(df.u) <- c("xvar", "yvar")
#   names(df.v) <- names(df.u)
#   if (pc.biplot) {
#     df.u <- df.u * nobs.factor
#   }
#   r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
#   v.scale <- rowSums(v^2)
#   df.v <- r * df.v/sqrt(max(v.scale))
#   if (obs.scale == 0) {
#     u.axis.labs <- paste("standardized PC", choices, sep = "")
#   }
#   else {
#     u.axis.labs <- paste("PC", choices, sep = "")
#   }
#   u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
#                                             100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
#   if (!is.null(labels)) {
#     df.u$labels <- labels
#   }
#   if (!is.null(groups)) {
#     df.u$groups <- groups
#   }
#   if (varname.abbrev) {
#     df.v$varname <- abbreviate(rownames(v))
#   }
#   else {
#     df.v$varname <- rownames(v)
#   }
#   df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
#   df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
#   # MAKE VARIABLE NAME
#   df.u <- df.u %>% as.data.frame() %>%
#     rownames_to_column(var="UniqueIDs")
#   # MAKE COLOURS
#   if (length(colors) != nrow(df.u)) {
#     colors <- rep("black", nrow(df.u))
#     names(colors) <- df.u$UniqueIDs
#   }
#   g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
#     ylab(u.axis.labs[2]) + coord_equal() 
#   if (!is.null(df.u$labels)) {
#     if (!is.null(df.u$groups)) {
#       g <- g + geom_text(aes(color=groups, label = labels), 
#                          size = labels.size)
#     }
#     else {
#       g <- g + geom_text(aes(label = labels), size = labels.size)
#     }
#   }
#   else {
#     if (!is.null(df.u$groups)) {
#       # g <- g + geom_point(aes(color = groups), alpha = alpha, show.legend = FALSE)
#       g <- g + geom_point(aes(fill=UniqueIDs, col=groups), alpha = alpha, show.legend = FALSE, pch=21, size=2) +# Don't colour by groups
#         scale_fill_manual(values=colors)
#     }
#     else {
#       g <- g + geom_point(aes(fill=UniqueIDs), col=outline_col, alpha = alpha, show.legend = FALSE, pch=21, size=2)+
#         scale_fill_manual(values=colors)
#     }
#   }
#   if (!is.null(df.u$groups) && ellipse) {
#     theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
#     circ <- cbind(cos(theta), sin(theta))
#     ell <- ddply(df.u, "groups", function(x) {
#       if (nrow(x) <= 2) {
#         return(NULL)
#       }
#       sigma <- var(cbind(x$xvar, x$yvar))
#       mu <- c(mean(x$xvar), mean(x$yvar))
#       ed <- sqrt(qchisq(ellipse.prob, df = 2))
#       data.frame(sweep(circ %*% chol(sigma) * ed, 2, 
#                        mu, FUN = "+"), groups = x$groups[1])
#     })
#     names(ell)[1:2] <- c("xvar", "yvar")
#     # g <- g + geom_path(data = ell, aes(color = groups, group = groups)) # eventually change so colours are unique here
#     g <- g + geom_path(data = ell, aes(color=groups, group = groups))
#   }
#   if (var.axes) {
#     if (circle) {
#       theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
#                                                 length = 50))
#       circle <- data.frame(xvar = r * cos(theta), yvar = r * 
#                              sin(theta))
#       g <- g + geom_path(data = circle, color = muted("white"), 
#                          size = 1/2, alpha = 1/3)
#     }
#     g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
#                                            xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
#                                                                                                   "picas")), color = muted(arrow_col))
#   }
#   if (var.axes) {
#     g <- g + geom_text(data = df.v, aes(label = varname, 
#                                         x = xvar, y = yvar, angle = angle, hjust = hjust), 
#                        color = text_col, size = varname.size)
#   }
#   return(g)
# }
# 



##### Get final version of health score assessment plots ########
# 
# ## Filter out testing groups
# allDat_final$experiment %>% unique()
# allDat_prelim <- allDat_final %>%
#   filter(experiment %in% c("2023-05-30_ccppwn","2023-06-07_ccpwn", "2023-05-03_age_concentration","2023-05-10_age_concentration","2023-05-16_basiceffects_n2c3"
#                             ,"2023-05-24_age_concentration")) %>%
#   filter(path!="PAO1", protect!="PO6") %>%
#   mutate(path = factor(path, levels=c("MOCK","N2C3"))
#          , protect = factor(protect, levels=c("MOCK","WCS365","CHAO","CH267","PF5")))
# allDat_prelim %>%
#   select(protect, path) %>% table()
# # "2023-05-03_age_concentration","2023-05-10_age_concentration","2023-05-16_basiceffects_n2c3"
# # ,"2023-05-24_age_concentration", ###### Other experiments we could use
# 
# allRGB <- allDat_prelim %>%
#   mutate(logpix = log10(all_plant_pixels)) %>%
#   select(logpix, pix_r_median, pix_g_median, pix_b_median) %>%
#   rename(Number_pixels_log=logpix, red_intensity = pix_r_median, green_intensity = pix_g_median, blue_intensity = pix_b_median) %>%
#   as.matrix()
# rownames(allRGB) <- allDat_prelim$UniqueID
# allHSV <- cbind(allRGB[,1], t(apply(allRGB, 1,function(x) rgb2hsv(r=x[2],g=x[3],b=x[4]))))
# colnames(allHSV) <- c("Number_pixels_log","hue","saturation","value")
# allHSL <- cbind(allRGB[,1], t(apply(allRGB, 1,function(x) rgb2hsl(as.matrix(x[2:4])))))
# colnames(allHSL) <- c("Number_pixels_log","hue","saturation","lightness")
# 
# # Get all colours
# # Get colors
# col_all <- apply(allRGB, 1, FUN = function(x) rgb(x[2]/255,x[3]/255,x[4]/255))
# names(col_all) <- rownames(allRGB)
# 
# 
# 
# ###### make PCA #########
# 
# pca_rgb <- prcomp(allRGB)
# pca_hsv <- prcomp(allHSV)
# pca_hsl <- prcomp(allHSL)
# 
# 
# # PCA
# # pca_all <- prcomp(alldat_mat)
# gg_biplot_rgb <- ggbiplot_PLANTS(pca_rgb, colors = col_all, outline_col=col_all,arrow_col="black", text_col="black",choices = c(1,2))
# gg_biplot_rgb
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_biplot_rgb.png", gg_biplot_rgb
#        , width=5, height=4)
# 
# pca_rgb$x %>% as.data.frame() %>%
#   arrange(-PC2) %>% head()
# 
# gg_biplot_hsv <- ggbiplot_PLANTS(pca_hsv, colors = col_all, outline_col=col_all,arrow_col="black", text_col="black",choices = c(1,2)) +
#   xlim(-5,5)
# gg_biplot_hsv
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_biplot_hsv.png", gg_biplot_hsv
#        , width=5, height=4)
# 
# gg_biplot_hsl <- ggbiplot_PLANTS(pca_hsl, colors = col_all, outline_col=col_all,arrow_col="black", text_col="black",choices = c(1,2))
# gg_biplot_hsl
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_biplot_hsl.png", gg_biplot_hsl
#        , width=5, height=4)
# 
# # Calculate x contrib to each axis; can we omit number of pixels?
# load_coord <- t(apply(pca_rgb$rotation, 1, FUN = function(x)x*pca_rgb$sdev))^2
# x_contrib <- t(apply(load_coord, 1, FUN=function(x) x/colSums(load_coord)))
# # Multiply out by the variance explained by each axis
# var_expl_axis <- pca_rgb$sdev^2/sum(pca_rgb$sdev^2)
# total_var_expl <- apply(x_contrib, 1, function(x) x*var_expl_axis)
# colSums(total_var_expl)
# ## 0.000087 is number pixels; or 0.0087% of variation explained by pixel count
# # The first 2 PCs explain this much variation in dataset:
# sum(var_expl_axis[c(1,2)])
# 
# # ###### Try without size to compare amount of variance explained by size
# 
# pca_rgb_nosize <- prcomp(allRGB[,-1])
# pca_hsv_nosize <- prcomp(allHSV[,-1])
# pca_hsl_nosize <- prcomp(allHSL[,-1])
# 
# 
# # PCA
# # pca_all <- prcomp(alldat_mat)
# gg_biplot_rgb_nosize <- ggbiplot_PLANTS(pca_rgb_nosize, colors = col_all, outline_col=col_all,arrow_col="black", text_col="black",choices = c(1,2))
# gg_biplot_rgb_nosize
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_biplot_rgb_nosize.png", gg_biplot_rgb_nosize
#        , width=5, height=4)
# 
# gg_biplot_hsv_nosize <- ggbiplot_PLANTS(pca_hsv_nosize, colors = col_all, outline_col=col_all,arrow_col="black", text_col="black",choices = c(1,2)) +
#   xlim(-2.5,4)
# gg_biplot_hsv_nosize
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_biplot_hsv_nosize.png", gg_biplot_hsv_nosize
#        , width=5, height=4)
# 
# gg_biplot_hsl_nosize <- ggbiplot_PLANTS(pca_hsl_nosize, colors = col_all, outline_col=col_all,arrow_col="black", text_col="black",choices = c(1,2))
# gg_biplot_hsl_nosize
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_biplot_hsl_nosize.png", gg_biplot_hsl_nosize
#        , width=5, height=4)
# 
# 
# # Note; I tried an NMDS too for rgb, but the shape looks exactly like a PCA and it's so much slowe rto run
# ###
# # 
# # 
# # ### Can a separate this out into MOCK, WCS, and N2C3 treatments only?
# # # unique(allDat_filt$experiment)
# # toKeepIDs <- allDat_filt %>%
# #   # filter(experiment == "2023-06-14_age_concentration_Pf5_PROPER") %>%
# #   # select(protect_od, path_od) %>% table()
# #   filter((protect=="WCS365" & path=="MOCK") | (protect == "MOCK" & path=="N2C3") | (protect=="MOCK" & path=="MOCK")) %>%
# #   filter(protect_od %in% c(0,0.01), path_od %in% c(0,0.01)) %>%
# #   filter(MS==0.5, MES==0.5, pH==5.8, plant=="col0", plant_age==5) %>%
# #   pull(UniqueID)
# # col_all_filt <- col_all[toKeepIDs]
# # 
# # pca_singlestrain <- prcomp(alldat_mat[toKeepIDs,])
# # groups_treat <- allDat_filt %>% filter(UniqueID %in% toKeepIDs) %>% pull(StrainMix)
# # # gg_biplot_singleStrains <- ggbiplot_PLANTS(pca_singlestrain, colors = col_all_filt, choices = c(1,2), groups = groups_treat, ellipse = TRUE)
# # # gg_biplot_singleStrains +
# # #   scale_color_manual(values=c("MOCK-MOCK"="black", "MOCK-N2C3"="darkred", "WCS365-MOCK"="blue"))
# # 
# # gg_pca_wcm <- ggbiplot(pca_singlestrain, choices = c(1,2), groups = groups_treat, ellipse = TRUE)+
# #   scale_color_manual(values=c("MOCK-MOCK"="black", "MOCK-N2C3"="darkred", "WCS365-MOCK"="blue"))+
# #   scale_x_reverse()+
# #   scale_y_reverse()
# # gg_pca_wcm
# # ggsave("03a-1_HealthScore_Assessment/gg_pca_wcm.png",gg_pca_wcm, width=5, height=4)
# 
# ######## Loadings for healthy/stressed, dead/alive #########
# pca_rgb_x <- data.frame(pca_rgb$x) %>%
#   mutate(PC1_adj = (PC1-min(PC1))/max(PC1), PC2_adj = ((-PC2)-min(-PC2))/max(-PC2)) %>%
#   mutate(combPC = PC1_adj * PC2_adj) %>%
#   rownames_to_column(var="UniqueID")
# 
# pca_hsv_x_nosize <- data.frame(pca_hsv_nosize$x) %>%
#   mutate(PC1_adj = (PC1-min(PC1))/max(PC1), PC2_adj = ((-PC2)-min(-PC2))/max(-PC2)) %>%
#   mutate(combPC = PC1_adj * PC2_adj) %>%
#   rownames_to_column(var="UniqueID")
# round(var_expl_axis[c(1)],3)*100
# pca_hsv_x <- data.frame(pca_hsv$x) %>%
#   mutate(PC1_adj = (PC1-min(PC1))/max(PC1), PC2_adj = ((-PC2)-min(-PC2))/max(-PC2)) %>%
#   mutate(combPC = PC1_adj * PC2_adj) %>%
#   rownames_to_column(var="UniqueID")
# round(var_expl_axis[c(1)],3)*100
# 
# gg_allplants_colsize <- pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_point(aes(x=PC1, y=-PC2, fill=UniqueID, size = all_plant_pixels), pch=21,  show.legend = FALSE)+
#   scale_fill_manual(values=col_all) +
#   ylab(paste0("PC2\n",round(var_expl_axis[c(2)],3)*100,"% of variation"))+
#   xlab(paste0("PC1\n",round(var_expl_axis[c(1)],3)*100,"% of variation"))+
#   scale_y_reverse()
# gg_allplants_colsize
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_allplants_rgb_colsize.png",gg_allplants_colsize
#        ,height=4, width=5)
# # Try colouring by ratios?
# 
# gg_allplants_hsv_nosize <- pca_hsv_x_nosize %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_point(aes(x=PC1, y=-PC2, fill=UniqueID, size = all_plant_pixels), pch=21,  show.legend = FALSE)+
#   scale_fill_manual(values=col_all) +
#   ylab(paste0("PC2\n",round(var_expl_axis[c(2)],3)*100,"% of variation"))+
#   xlab(paste0("PC1\n",round(var_expl_axis[c(1)],3)*100,"% of variation"))+
#   scale_y_reverse()
# gg_allplants_hsv_nosize
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_allplants_hsv_nosize.png",gg_allplants_hsv_nosize
#        ,height=4, width=5)
# 
# 
# gg_allplants_hsv <- pca_hsv_x %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_point(aes(x=PC1, y=-PC2, fill=UniqueID, size = all_plant_pixels), pch=21,  show.legend = FALSE)+
#   scale_fill_manual(values=col_all) +
#   ylab(paste0("PC2\n",round(var_expl_axis[c(2)],3)*100,"% of variation"))+
#   xlab(paste0("PC1\n",round(var_expl_axis[c(1)],3)*100,"% of variation"))+
#   scale_y_reverse() +
#   xlim(-1,1.5)
# gg_allplants_hsv
# ggsave(filename = "03a-1_HealthScore_Assessment/gg_allplants_hsv.png",gg_allplants_hsv
#        ,height=4, width=5)
# 
# gg_diff_gb <- pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_point(aes(x=PC1_adj, y=PC2_adj, col=diff_gb_median), show.legend = TRUE)+
#   scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Difference in\nGreen/Blue intensity")
# gg_diff_gb
# ggsave("03a-1_HealthScore_Assessment/gg_diff_gb.png", gg_diff_gb, width=5, height=4)
# 
# gg_diff_gr <- pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_point(aes(x=PC1_adj, y=PC2_adj, col=diff_gr_median), show.legend = TRUE)+
#   scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Difference in\nGreen/Red intensity")
# gg_diff_gr
# ggsave("03a-1_HealthScore_Assessment/gg_diff_gr.png", gg_diff_gr, width=5, height=4)
# 
# gg_pixels <- pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_point(aes(x=PC1_adj, y=PC2_adj, col=log10(all_plant_pixels)), show.legend = TRUE)+
#   scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Plant size (log pixels)")
# gg_pixels
# ggsave("03a-1_HealthScore_Assessment/gg_pixels.png", gg_pixels, width=5, height=4)
# 
# gg_comb <- pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   mutate(comb = (diff_gb_median + diff_gr_median*2)*log10(all_plant_pixels) )%>%
#   ggplot() +
#   geom_point(aes(x=PC1_adj, y=PC2_adj, col=comb), show.legend = TRUE)+
#   scale_y_reverse()
# gg_comb
# ggsave("03a-1_HealthScore_Assessment/gg_comb.png", gg_comb, width=5, height=4)
# 
# ###### Try hsv
# 
# gg_hsv_grad_h <- pca_hsv_x_nosize %>%
#   left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID")) %>%
#   ggplot() +
#   geom_point(aes(x=PC1_adj, y=PC2_adj, col=hue), show.legend = TRUE)+
#   scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Hue")
# gg_hsv_grad_h
# ggsave("03a-1_HealthScore_Assessment/gg_hsv_grad_h.png", gg_hsv_grad_h, width=5, height=4)
# 
# 
# gg_hsv_grad_s <- pca_hsv_x_nosize %>%
#   left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID")) %>%
#   ggplot() +
#   geom_point(aes(x=PC1_adj, y=PC2_adj, col=saturation), show.legend = TRUE)+
#   scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Saturation")
# gg_hsv_grad_s
# ggsave("03a-1_HealthScore_Assessment/gg_hsv_grad_s.png", gg_hsv_grad_s, width=5, height=4)
# 
# 
# # Can this be used as a single variable?
# pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   ggplot() +
#   geom_jitter(aes(x=diff_gb_median, y=1, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
#   scale_color_manual(values=col_all)
# 
# pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   mutate(comb=(diff_gr_median+diff_gb_median)
#          ,comb2 = ((diff_gr_median*2 + diff_gb_median*1)*log10(all_plant_pixels+1)) )%>%
#   ggplot() +
#   geom_jitter(aes(x=comb2, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
#   scale_color_manual(values=col_all)
# 
# pca_rgb_x %>%
#   left_join(allDat_prelim) %>%
#   mutate(comb=(diff_gr_median+diff_gb_median)
#          ,comb2 = ((diff_gr_median*2 + diff_gb_median*1)*(diff_gr_median/10 + log10(all_plant_pixels+1))) )%>%
#   ggplot() +
#   geom_jitter(aes(x=comb2, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
#   scale_color_manual(values=col_all)
# 
# 
# pca_hsv_x_nosize %>%
#   left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
#   left_join(allDat_prelim) %>%
#   mutate(comb=(1*hue)*(3*saturation))%>%
#   ggplot() +
#   geom_jitter(aes(x=comb, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
#   scale_color_manual(values=col_all)
# 
# 
# pca_hsv_x_nosize %>%
#   left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
#   left_join(allDat_prelim) %>%
#   mutate(comb=(hue)*(2*saturation+log10(all_plant_pixels+1)/100))%>%
#   ggplot() +
#   geom_jitter(aes(x=comb, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
#   scale_color_manual(values=col_all)
# 
# 
# pca_hsv_x_nosize %>%
#   left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
#   left_join(allDat_prelim) %>%
#   mutate(comb=(hue+log10(all_plant_pixels+1)/100)*(2*saturation))%>%
#   ggplot() +
#   geom_jitter(aes(x=comb, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
#   scale_color_manual(values=col_all)
# 
# #### UPDATE allDAT WITH NEW METRIC #####
# allDat_prelim <- allDat_prelim %>%
#   left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
#   mutate(
#     # Healthiness_nosize = diff_gb_median + diff_gr_median*2
#     Healthiness_rgb = (diff_gb_median + diff_gr_median*2)*log10(all_plant_pixels+1)
#     , Healthiness_hsv = hue*saturation*log10(all_plant_pixels+1)*1000
#     # , Healthiness_hsv = (hue)*(2*saturation+log10(all_plant_pixels+1)/100)
#     # , protect_cells_log = log10(protect_cells+1)
#     # , path_cells_log = log10(path_cells+1)
#     # , path_cells_adj_log = log10(path_cells_adj+1)
#     , ratio_protect_log = log10(ratio_protect+1)
#     , ratio_path_log = log10(ratio_path + 1)
#     , ratio_protpath_log = (ratio_protect_log+1)/(ratio_path_log+1)-1
#     # , total_cells_log = log10(total_cells+1)
#     # , total_cells_adj_log = log10(total_cells_adj + 1)
#   )
# 
# 
# 
# #### Look at histogram of "healthiness" ####
# ggHist <- allDat_prelim %>%
#   # filter(path=="N2C3", protect=="MOCK") %>%
#   ggplot() +
#   geom_histogram(aes(x=Healthiness_rgb), bins=100) +
#   ylab("Number of plants") +
#   xlab("") +
#   theme(axis.text.x = element_blank())
# ggCol <- allDat_prelim %>%
#   ggplot() +
#   geom_jitter(aes(x=Healthiness_rgb, y=1, col=UniqueID, size=all_plant_pixels), height=0.1, width=0, show.legend=FALSE)+
#   scale_color_manual(values=col_all)+
#   ylab("")+
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
#   xlab("Health score\n((Green:Blue + 2*(Green:Red))*log10(Number of pixels))")
# 
# gg_Healthiness <- plot_grid(ggHist, ggCol, ncol=1, align='v')
# gg_Healthiness
# ggsave(filename = "03_EDA_and_adj/gg_Healthiness_rgb.png",
#        gg_Healthiness, width=12, height=6)
# # 
# # ggHist_nosize <- allDat_filt %>%
# #   # filter(path=="N2C3", protect=="MOCK") %>%
# #   ggplot() +
# #   geom_histogram(aes(x=Healthiness_rgb_now), bins=100) +
# #   ylab("Number of plants") +
# #   xlab("") +
# #   theme(axis.text.x = element_blank())
# # ggCol_nosize <- allDat_filt %>%
# #   ggplot() +
# #   geom_jitter(aes(x=Healthiness_nosize, y=1, col=UniqueID, size=all_plant_pixels), height=0.1, width=0, show.legend=FALSE)+
# #   scale_color_manual(values=col_all)+
# #   ylab("")+
# #   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
# #   xlab("Healthiness metric\n((Green:Blue + 2*(Green:Red)))")
# # 
# # gg_Healthiness_nosize <- plot_grid(ggHist_nosize, ggCol_nosize, ncol=1, align='v')
# # gg_Healthiness_nosize
# # ggsave(filename = "03_EDA_and_adj/gg_Healthiness_nosize.png",
# #        gg_Healthiness_nosize, width=12, height=6)
# 
# 
# ggHist <- allDat_prelim %>%
#   # filter(path=="N2C3", protect=="MOCK") %>%
#   ggplot() +
#   geom_histogram(aes(x=Healthiness_hsv), bins=100) +
#   ylab("Number of plants") +
#   xlab("") +
#   theme(axis.text.x = element_blank())
# ggCol <- allDat_prelim %>%
#   ggplot() +
#   geom_jitter(aes(x=Healthiness_hsv, y=1, col=UniqueID, size=all_plant_pixels), height=0.1, width=0, show.legend=FALSE)+
#   scale_color_manual(values=col_all)+
#   ylab("")+
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
#   xlab("Health score\n(Hue*Saturation*log10(Plant pixels + 1))")
# 
# gg_Healthiness <- plot_grid(ggHist, ggCol, ncol=1, align='v')
# gg_Healthiness
# ggsave(filename = "03_EDA_and_adj/gg_Healthiness_hsv.png",
#        gg_Healthiness, width=12, height=6)
