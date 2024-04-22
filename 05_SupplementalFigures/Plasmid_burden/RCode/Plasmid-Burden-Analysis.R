library(dplyr)
library(ggplot2)
library(lubridate)
library(gcplyr)
library(tidyr)
# ---- Experimental Notes -------
## Volume for all samples was 200 ul
## Inoculation amount for all samples was 0.01 (abs 600nm path length of 100ul on plate reader)
# All trials were run for ~24h at 28C
# ---------------------------------


All2345 <- read.csv("dat/001_Growth_curve_data_compiled_24-03-29.csv") #Compiled raw data from growth curve trials 1-5

# Convert Time to numeric
All2345$Time <- as.numeric(as.character(All2345$Time))

All2345$Well <- factor(All2345$Well,
                       levels = paste(rep(LETTERS[1:8], each = 12), 1:12, sep = ""))


####### Calculating blanks for each trial
blank_values <-  All2345  %>%
  filter(is.na(Strain) | Strain == "Blank") %>%
  filter(Media == "LB" | Media == "0.5MS0.5MES30mMSuccinate") %>%
  #filter(Time=="0") %>%
  group_by(Experiment) %>%
  summarise(blank=mean(Abs600))

All2345 <- All2345 %>%  
  full_join(blank_values) %>% 
  mutate(Abs600_blanked = Abs600 - blank)


### Filtering outliers and bacteria inoculated outside of 0.01 or in MS (MSMES + Sucrose)

All2345_filt <- All2345 %>%
  filter(Media != "MS") %>%
  filter(Media != "NA") %>%
  filter(!(Well == "F7" & Experiment == "Trial_2023_12_14")) %>%
  filter(!(Well == "B5" & Experiment == "Trial5_2024-03-20")) %>%
  filter(Innoculation == "0.01" | is.na(Innoculation) ) 

# Replace NA values with a placeholder label
All2345_filt$Plasmid[is.na(All2345_filt$Plasmid)] <- "NA"

All2345_filt$Plasmid <- factor(All2345_filt$Plasmid, 
                               levels = c("Crim", "Neon", "NA"))

#### GROWTH CURVE PLOT #####
gg_growth_curve_1 <- All2345_filt %>%
  filter(Strain != "Blank" & !is.na(Strain)) %>%
  group_by(Strain,Plasmid, Media, Time) %>%
  summarise(mean_Abs600_blanked = mean(Abs600_blanked),
            sd_Abs600_blanked = sd(Abs600_blanked)) %>%
  ggplot(aes(x = Time, y = mean_Abs600_blanked, colour = Strain)) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_Abs600_blanked - sd_Abs600_blanked,
                    ymax = mean_Abs600_blanked + sd_Abs600_blanked),
                width = 0.1) +
  facet_wrap(Media~ Plasmid, labeller = labeller(Plasmid = c("Crim" = "Crimson", "Neon"="Neon", "NA" = "No Plasmid"), Media = c("LB" = "LB", "0.5MS0.5MES30mMSuccinate" = "MS + Succinate"))) +
  labs(x = "Time", y = "Abs600_blanked")+
  scale_color_manual(values = c(
    "WCS365" = "darkseagreen4", 
    "N2C3" = "darkorange",
    "WCS365-LacZ" = "black"
  ))+
  labs(
    x = "Time (hours)", 
    y = "Abs 600nm"
  ) +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
  ) 

gg_growth_curve_1
ggsave("gg_growth_curve_1.png", gg_growth_curve_1, height=6, width=10)


### -------  CALCULATING SUMMARY STATS --------- #####

All234_filt <- All2345_filt
######Calculating derivatives 
ex_dat_mrg <- All234_filt

ex_dat_mrg <- ex_dat_mrg %>% 
  rename(Measurements = Abs600_blanked)

#Calculate the plain derivative (the slope of the origonal density data *how quickly does the population grow at each time point)
ex_dat_mrg <- mutate(group_by(ex_dat_mrg, Well, Experiment, Name ), deriv = calc_deriv(x = Time, y = Measurements)) #Can also group by more catagories -> (group_by(ex_dat_mrg, Well, Bacteria_strain, Phage)

head(ex_dat_mrg)
#Wells we want to visualize below
sample_wells <- c("B2", "B8", "G2", "G8")


# plot the derivative for some sample_wells
ggplot(data = dplyr::filter(ex_dat_mrg, Well %in% sample_wells), aes(x = Time, y = deriv)) +
  geom_line() +
  facet_wrap(Experiment~Well, scales = "free")

#calcualte the per-capita derivative (the growth rate of the cells)
ex_dat_mrg <- mutate(group_by(ex_dat_mrg, Well), deriv_percap = calc_deriv(x = Time, y = Measurements, percapita = TRUE, blank = 0.2996))

# plot the per-capita derivative
ggplot(data = dplyr::filter(ex_dat_mrg, Well %in% sample_wells), aes(x = Time, y = deriv_percap)) +
  geom_line() +
  facet_wrap(Experiment~Well, scales = "free")

#fitting a linear regression to multiple points, reducing this jumpiness in the derivative
ex_dat_mrg <- mutate(group_by(ex_dat_mrg,Experiment, Well), deriv_percap5 = calc_deriv(x = Time, y = Measurements,
                                                                                       percapita = TRUE, blank = 0, 
                                                                                       window_width_n = 5, trans_y = "log"))

# Now let's plot the smoothed per-capita derivative
ggplot(data = dplyr::filter(ex_dat_mrg, Well %in% sample_wells), aes(x = Time, y = deriv_percap5)) +
  geom_line() +
  facet_wrap(Experiment~Well, scales = "free")

######Converting per-capita growth rates into doubling times
ex_dat_mrg <- mutate(group_by(ex_dat_mrg,Experiment, Well), deriv_percap5 = calc_deriv(x = Time, y = Measurements,
                                                                                       percapita = TRUE, blank = 0,window_width_n = 5, trans_y = "log"), doub_time = doubling_time(y = deriv_percap5))

#View(ex_dat_mrg)

dat_derv <- ex_dat_mrg

#Maximum growth rate and minimum doubling time

ex_dat_mrg_sum <-
  summarize(group_by(ex_dat_mrg, Well, Experiment, Name),
            max_percap = max_gc(deriv_percap5, na.rm = TRUE), max_percap_time = extr_val(Time, which_max_gc(deriv_percap5)), doub_time = doubling_time(y = max_percap))

head(ex_dat_mrg_sum)

dat_MG <- ex_dat_mrg_sum

ggplot(data = dplyr::filter(ex_dat_mrg, Well %in% sample_wells), aes(x = Time, y = deriv_percap5)) +
  geom_line() +
  facet_wrap(Experiment~Well) +
  geom_point(data = dplyr::filter(ex_dat_mrg_sum, Well %in% sample_wells),
             aes(x = max_percap_time, y = max_percap),
             size = 2, color = "red") + coord_cartesian(ylim = c(-1, NA))

#Maximum density
#The maximum bacterial density can be a measure of bacterial growth yield/efficiency. If your bacteria plateau in density, the maximum density can also be a measure of bacterial carrying capacity. 
ex_dat_mrg_sum <-
  summarize(group_by(ex_dat_mrg, Well, Experiment, Name),
            max_dens = max_gc(Measurements, na.rm = TRUE),
            max_time = extr_val(Time, which_max_gc(Measurements)))

dat_MD <- ex_dat_mrg_sum
head(ex_dat_mrg_sum)

ggplot(data = dplyr::filter(ex_dat_mrg, Well %in% sample_wells), aes(x = Time, y = Measurements)) +
  geom_line() +
  facet_wrap(Experiment~Well) +
  geom_point(data = dplyr::filter(ex_dat_mrg_sum, Well %in% sample_wells),
             aes(x = max_time, y = max_dens), size = 2, color = "red")

#Area under the curve - The area under the curve is a common metric of total bacterial growth, for instance in the presence of antagonists like antibiotics or phages

ex_dat_mrg_sum <-
  summarize(group_by(ex_dat_mrg, Well,Experiment, Name),
            auc = auc(x = Time, y = Measurements))

dat_AUC <- ex_dat_mrg_sum

head(ex_dat_mrg_sum)


### Combine all the data
library(purrr)

# List of data frames we want to combine 
data_frames_list <- list(dat_derv, dat_MG, dat_MD, dat_AUC)

# Merging all data frames
all_dat_summary <- reduce(data_frames_list, full_join, by = c("Well", "Experiment", "Name"))

write.table(all_dat_summary, file="002_Growth_curve_data_and_summary_stats_compiled_2024-03-29.txt", row.names = FALSE, quote = FALSE, sep="\t")


### -------  DOUBLING TIME FIGURE  --------- #####
all_dat_summary_filt <-all_dat_summary %>%
  filter(Media != "MS") %>%
  filter(Media != "NA") %>%
  filter(Strain != "NA") %>%
  filter(Strain != "Blank")

gg_doubling_Time <- ggplot(data = all_dat_summary_filt, aes(x=Strain, y=max_percap_time)) +
  geom_boxplot() +
  geom_point(aes(color=Strain, shape = Experiment)) +
  facet_grid(Media~Plasmid, scales = "free_x", drop=TRUE, space="free_x", labeller = labeller(Plasmid = c("Crim" = "Crimson", "Neon"="Neon", "NA" = "No Plasmid"), Media = c("LB" = "LB", "0.5MS0.5MES30mMSuccinate" = "MS + Succinate")))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = " ", y = "Doubling Time")+
  scale_color_manual(values = c(
    "WCS365" = "darkseagreen4", 
    "N2C3" = "darkorange",
    "WCS365-LacZ" = "black"
  ))+
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
  ) 

gg_doubling_Time
ggsave("gg_doubling_Time_1.png", gg_doubling_Time, height=6, width=10)
