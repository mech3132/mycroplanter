#!bin/bash

###### Step one: compile data ###########

library(lubridate) # Better to manage dates with
library(tidyverse) # data wrangling with dplyr package + ggplot2 package for plotting
# NOTE: always load tidyverse last because the command "select" is a common function name and can be masked by other packages

# I name the output directory the same name as the R script so I can match code to outputs easier
dir.create("01_compile_data_forleah")

#### Load in all metadata ####
# Get metadata file names
meta_filenames <- list.files("00_raw_data/meta")
# Iterate through vector of filenames and load them in, joining with each other.
allMeta <- data.frame(plate=NA) # Empty dataframe must have at least one common column name as files you're going to load in
for ( m in meta_filenames) {
  tempmeta <- read.csv(paste0("00_raw_data/meta/",m)) # Adding the folder "meta" to filepath
  allMeta <- full_join(allMeta, tempmeta)
}
# remove all lines where plate is NA
allMeta  <- filter(allMeta, !is.na(plate))
# allMeta %>% filter(experiment=="2024-02-14_doubledip")

#### Load in all pixel data ####
allPixelDir <- list.dirs("00_raw_data/dat/pixel_data")[list.dirs("00_raw_data/dat/pixel_data")!="00_raw_data/dat/pixel_data"]
# Above, the list.dirs returns the original directory ("00_raw_data/dat/pixel_data") so I just remove that by removing the index where the value equals that path

# Loop through each directory and get the file pixel_data.txt
allPixel <- data.frame(subimage=NA)
for (p in allPixelDir) {
  # The directory name follows the formula: "processed_",EXPERIMENTNAME, __, PLATENAME", so we need to extract the platename
  exp__plate <- gsub(".*processed_","",p) # gsub subs first string with second string. .* is "unlimited number of wildcard variables".
  experiment <- strsplit(exp__plate, split="__")[[1]][1]
  plate <- strsplit(exp__plate, split="__")[[1]][2]
  tempdat <- read.delim(paste0(p,"/pixel_data.txt")) %>%
    mutate(experiment = experiment, plate=plate) # mutate creates new columns in your data frame. Adding experiment and plate to each
  allPixel <- full_join(allPixel, tempdat)
}
# remove subimage == NA
allPixel <- filter(allPixel, !is.na(subimage)) 

#### Load in all fluorescence data ####
# Get fluor file names
fluor_filenames <- list.files("00_raw_data/dat/fluor_data/")
# Iterate through vector of filenames and load them in, joining with each other.
allFluor <- data.frame(plate=NA) # Empty dataframe must have at least one common column name as files you're going to load in
for ( f in fluor_filenames) {
  # f = fluor_filenames[2]
  tempfluor <- read.csv(paste0("00_raw_data/dat/fluor_data/",f)) # Adding the folder to filepath with paste0
  allFluor <- full_join(allFluor, tempfluor)
}
# remove all lines where plate is NA
allFluor  <- filter(allFluor, !is.na(plate))

######## Get fluorescence standard curves ############
# lm_offset <- read.delim("../Project_HT_methodspaper_fluoranalysis/02_compile_data/offest_average.txt")
# I usually use the path above ^^^ but I've copied it here for you since you won't have my other directories
lm_offset <- read.delim("offest_average.txt")
# This value is the average offset across exp 2023-09-13 and 2023-09-14
# Slope was froced to 1
# outliers were removed.

####### Adjust files and process/add additional info ############
# Conversion for OD/abs to cells #
# OD = 3.166*abs - 0.268 (but only When 0.3>Abs>0.1)--> formual derived from experiment
# Cells per ml for OD1 is 5e8
# cells_per_Abs1 = (3.166*1-0.268)*5E8 ## OLD ONE, including all Abs600
cells_per_Abs1 = (3.19*1-0.275)*5E8 ## NEW ONE, excluding Abs600<0.125


allMeta_adj  <- allMeta %>%
  # Make all dates lubridate dates, which are easier to handle than basic "R" dates
  mutate(date_germ = as.Date(date_germ)
         , date_ster = as.Date(date_ster)
         , date_inoc = as.Date(date_inoc)
         , date_process = as.Date(date_process)) %>%
  # Separate the strain1 and strain2 names into lacZ and fluor info
  # First, get rid of "mixed" ones
  separate(strain1, into= c("strain1","other"), sep="-") %>%
  separate(other, into= c("lacZ_strain1","fluor_strain1"), sep="_") %>% 
  mutate(lacZ_strain1=ifelse(lacZ_strain1=="",NA,lacZ_strain1),fluor_strain1=ifelse(fluor_strain1=="",NA,fluor_strain1)) %>%
  separate(strain2, into= c("strain2","other"), sep="-") %>%
  separate(other, into= c("lacZ_strain2","fluor_strain2"), sep="_") %>% 
  mutate(lacZ_strain2=ifelse(lacZ_strain2=="",NA,lacZ_strain2),fluor_strain2=ifelse(fluor_strain2=="",NA,fluor_strain2))%>%
  mutate(col=factor(col, levels=c(seq(1,12)))) %>%
  # select(experiment, plate, row, col, strain1, strain2) %>%
  # View()
  unite(strain1, strain2, col="StrainMix", sep = "-", remove=FALSE) %>% # make a strainmix column by uniting the strain2 and strain1 columns
  rowwise() %>% # this tells it to calculate by row, rather that summing across all rows
  mutate(totalOD = strain2_od + strain1_od # Calculate new column that is TOTAL OD
         , ratio_strain2 = strain2_od/min(c(strain2_od, strain1_od)) # ratio of strain2:strain1
         , ratio_strain1 = strain1_od/min(c(strain2_od, strain1_od)) # ratio of strain1:strain2
         # , ratio_strain1_adj = strain1_od_adj/min(c(strain2_od, strain1_od_adj))
  ) %>% 
  # When dividing by zero, remove these.
  mutate(ratio_strain2 = ifelse(is.infinite(ratio_strain2), 1, 
                                ifelse(is.na(ratio_strain2), 0, ratio_strain2))
         , ratio_strain1 = ifelse(is.infinite(ratio_strain1), 1, 
                                  ifelse(is.na(ratio_strain1), 0, ratio_strain1))) %>%
  mutate(ratio = paste0(c(ratio_strain2, ratio_strain1), collapse = "-")) %>% # category of ratios used
  # mutate(ratio_adj = paste0(c(ratio_strain2, ratio_strain1_adj), collapse = "-")) %>%
  ungroup() %>% # undoes the "rowwise" command
  # Below, we calculate the number of cells (and log of number of cells)
  mutate(strain2_cells = strain2_od*cells_per_Abs1,
         strain1_cells = strain1_od*cells_per_Abs1, 
         total_cells = totalOD*cells_per_Abs1,
         strain2_cells_log = log10(strain2_cells+1),
         strain1_cells_log = log10(strain1_cells+1),
         total_cells_log = log10(total_cells+1)) %>%
  mutate(plant_age = date_inoc - date_germ) # Plant age at inoculation

# Optionally, we can manually view the new data frame like this:
# View(allMeta_adj)
allMeta_adj %>%
  select(experiment) %>% table()

# Get hsv
allRGB <- allPixel %>% select(pix_r_median, pix_g_median, pix_b_median)
allHSV <- t(apply(allRGB, 1, function(x) rgb2hsv(r=x[1],g=x[2],b=x[3])))
colnames(allHSV) <- c("hue","saturation","value")
allPixel <- cbind(allPixel, allHSV)
# Adjust the pixel data
allPixel_adj <- allPixel %>% 
  # The subimage name is the experiment__plate__well, so we're going to separate that out
  separate(subimage, into = c("remove","well"), remove=TRUE, sep="_") %>% select(-remove) %>%
  separate(well, into=c("row","col"), remove=FALSE, sep=1) %>% 
  mutate(col=factor(col, levels=c(seq(1,12)))) %>% # Make this a factor, which is a sequence from 1-12 (otherwise the alphabetical order will be 1, 10, 11, 12, 2, 3, 4...) 
  mutate(Healthiness_rgb = (diff_gb_median + diff_gr_median*2)*log10(all_plant_pixels+1),
         Healthiness_hsv = hue*saturation*log10(all_plant_pixels+1)*1000) # This is from my own analysis


##### Adjust fluorescence data
# Going to rename columns to remain consistent with the protect/path wording
# Avoid '.' because those are typically "wildcard" symbols

allFluor_adj <- allFluor %>%
  # rename(neon_raw = Neon.raw,  crim_raw=Crim.raw) %>%
  # rename(neon_CFU = Neon.cfu, crim_CFU = Crim.cfu) %>%
  mutate(col=factor(col, levels=c(seq(1,12)))) # Make this a factor, which is a sequence from 1-12 (otherwise the alphabetical order will be 1, 10, 11, 12, 2, 3, 4...)


##### JOIN and check #####
allDat <- full_join(allMeta_adj,allPixel_adj) %>%
  full_join(allFluor_adj) %>% 
  unite(experiment, plate, well, col="UniqueID", sep="__", remove=FALSE)

### Platecheck 
allDat %>%
  filter(experiment=="2023-12_priorityeffects") %>%
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2)) +
  facet_grid(readtime~plate)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", N2C3="orange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkgreen", N2C3="orange"), na.value = "black")

gg_platemap_1 <- allDat %>%
  filter(experiment=="2023-12_priorityeffects") %>%
  filter(readtime==7) %>%
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2)) +
  facet_grid(readtime~plate)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", N2C3="orange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkgreen", N2C3="orange"), na.value = "black")
gg_platemap_1

gg_platemap_2 <- allDat %>%
  filter(experiment=="2023-12_priorityeffects2") %>%
  # select(plate, row, col, strain1, strain2,Healthiness_hsv,crim_btm, neon_btm ) %>% View()
  # select(plate) %>% table()
  # select(plate, row, col, strain1, strain2, ) %>%
  # View()
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2)) +
  facet_grid(.~plate)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", N2C3="orange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkgreen", N2C3="orange"), na.value = "black")
gg_platemap_2

gg_platemap_dd <- allDat %>%
  filter(experiment=="2024-02-14_doubledip") %>%
  # select(plate, row, col, strain1, strain2,Healthiness_hsv,crim_btm, neon_btm ) %>% View()
  # select(plate) %>% table()
  # select(plate, row, col, strain1, strain2, ) %>%
  # View()
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=strain1)) +
  geom_point(aes(x=col, y=(row), col=strain2)) +
  facet_grid(.~plate)+
  scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", N2C3="orange"), na.value = "white")+
  scale_color_manual(values=c(MOCK="white", WCS365="darkgreen", N2C3="orange"), na.value = "black")
gg_platemap_dd

## Visual check of all fluorescence
### ALL first strains were neon

allDat %>%
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  # select(StrainMix) %>% View()
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK",
                       ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis neon"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=log10(neon_btm+1))) +
  geom_point(aes(x=col, y=(row), col=mock)) +
  facet_grid(readtime~plate)+
  scale_fill_gradient(low="white", high="green")+
  scale_color_manual(values=c("grey","darkgreen","blue"))

allDat %>%
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK",
                       ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis neon"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=log10(crim_btm+1))) +
  geom_point(aes(x=col, y=(row), col=mock)) +
  facet_grid(readtime~plate)+
  scale_fill_gradient(low="white", high="red")+
  scale_color_manual(values=c("grey","darkgreen","black"))

## Let's confirm  that all the controls look okay---
allDat %>%
  # filter(plant =="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(ismock = ifelse(StrainMix=="MOCK-MOCK","TRUE MOCK",
                       ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis one of neon/crim"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=row, fill=log10(neon_top+1))) +
  geom_point(aes(x=col, y=row, col=ismock)) +
  facet_grid(readtime~plate)+
  scale_fill_gradient(low="white", high="green")+
  scale_color_manual(values=c("grey","darkgreen","white"), na.value = "black")

allDat %>%
  # filter(plant =="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(ismock = ifelse(StrainMix=="MOCK-MOCK","TRUE MOCK",
                         ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis one of neon/crim"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=row, fill=log10(crim_top+1))) +
  geom_point(aes(x=col, y=row, col=ismock)) +
  facet_grid(readtime~plate)+
  scale_fill_gradient(low="white", high="darkred")+
  scale_color_manual(values=c("grey","darkgreen","white"), na.value = "black")

allDat %>%
  filter(plant =="col0") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(ismock = ifelse(StrainMix=="MOCK-MOCK","TRUE MOCK",
                         ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis one of neon/crim"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=row, fill=log10(crim_top+1))) +
  geom_point(aes(x=col, y=row, col=StrainMix)) +
  facet_grid(readtime~plate)+
  scale_fill_gradient(low="white", high="darkred")+
  scale_color_manual(values=c("N2C3-WCS365"="grey","WCS365-N2C3"="grey"
                              ,"MOCK-WCS365"="darkgreen","WCS365-MOCK"="darkgreen"
                              ,"MOCK-N2C3"="darkorange","N2C3-MOCK"="darkorange"
                              ,"MOCK-MOCK" = "white"
                              ), na.value = "black")

## These wells looked like they were contaminated
# MANUALLY REMOVE:
man_remove <- c("2023-12_priorityeffects_Nov-10th-transfer-PlateD_cropped_G_11"
                ,"2023-12_priorityeffects_Nov-10th-transfer-PlateD_cropped_H_10"
                ,"2023-12_priorityeffects_Nov-9th-transfer-PlateC_cropped_G_8"
                ,"2023-12_priorityeffects_Nov-9th-transfer-PlateC_cropped_F_8" ## E and F look weird in follow up sanity checks
                ,"2023-12_priorityeffects_Nov-9th-transfer-PlateC_cropped_E_8")


# Are there any mis-alignments (ie: StrainMix is NA or all_plant_pixels is NA, indicating that some plate/experiment/well combo doesn't line up)
allDat %>%
  filter(is.na(StrainMix) | is.na(all_plant_pixels))
# Should have no rows!
# allDat %>% filter(is.na(plate)) %>% View()

# Are there any NAs in StrainMix?
allDat %>% filter((is.na(StrainMix))) 
# Are there any NAs in all_plant_pixels?
allDat %>% filter(is.na(all_plant_pixels))
# Make a table of data to see if number of samples makes sense
allDat %>%
  select(StrainMix, strain1) %>%
  table()
allDat %>%
  select(StrainMix, strain2) %>%
  table()

######### Make a histgram of pixel values to see whether there are any strange outliers ###########

allDat %>%
  ggplot() +
  geom_density(aes(x=pix_g_median), fill=NA, col="green")+
  geom_density(aes(x=pix_r_median), fill=NA, col="red")+
  geom_density(aes(x=pix_b_median), fill=NA, col="blue")
# The values where pixels are close to zero (<50) are actually wells that don't have plants in them (cross-checked)
# Filter so that all green pixels > 100, and red pixels must be >50
allDat_raw <- allDat %>%
  filter((pix_g_median>100| pix_r_median>50))

##### Split dataset into halves #######

#### Top vs Bottom ####
allDat_raw %>%
  ggplot() +
  geom_point(aes(x=log10(crim_top), y=log10(crim_btm)), col="darkred") +
  geom_point(aes(x=log10(neon_top), y=log10(neon_btm)), col="green") +
  geom_abline(aes(intercept=0, slope=1))
# Top and bottom correlate very well; top values are actually a bit higher than bottom values?

## Look at starting values, see which one is more variable
allDat_raw %>% 
  filter(readtime==0)%>%
  ggplot() +
  geom_point(aes(x=log10(crim_top), y=log10(crim_btm)), col="darkred") +
  geom_point(aes(x=log10(neon_top), y=log10(neon_btm)), col="green") +
  geom_abline(aes(intercept=0, slope=1))
# Crimson is not as nice as neon

allDat_raw %>% 
  filter(readtime==0) %>%
  select(row, col, plate, crim_top, neon_top, crim_btm, neon_btm) %>%
  pivot_longer(-c(row, col, plate), names_to="type", values_to="fluor") %>%
  separate(type, into=c("fluortype","tb"), remove=FALSE) %>%
  ggplot() +
  geom_tile(aes(x=col, y=row, fill=log10(fluor))) +
  facet_grid(plate ~ type)

### Top reads are actually better, and since we use top for other exp let's just use top for this one too

######### PROCESS ################
# Check out blanks 
# Check out blanks to subtract

allDat_raw %>%
  # filter(experiment %in% c("2023-09-13_fluor","2023-09-14_fluor")) %>%
  # filter(all_plant_pixels>100) %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  select(experiment, crim_top, neon_top) %>%
  pivot_longer(-c(experiment), names_to="fluor", values_to="raw") %>%
  ggplot() +
  geom_point(aes(x=fluor, y=log10(raw))) +
  facet_wrap(.~experiment)

### Remove crimson values over E4; and remove neon values over 4.5

allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  filter(crim_top<10^4, neon_top< 10^4.5) %>%
  select(experiment, crim_top, neon_top) %>%
  pivot_longer(-c(experiment), names_to="fluor", values_to="raw") %>%
  ggplot() +
  geom_point(aes(x=fluor, y=log10(raw))) +
  facet_wrap(.~experiment)
# Looks better; this filter works

### Now calculate blanks, without the contaminated samples
blank_values <-  allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  filter(crim_top<10^4, neon_top< 10^4.5) %>%
  group_by(experiment) %>%
  summarise(crim_top_blanks=mean(crim_top), neon_top_blanks=mean(neon_top)
  )

contaminated_blanks <-  allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  filter(crim_top>=10^4 |  neon_top >= 10^4.5) %>%
  unite(experiment, plate, row, col, col=UniqueID) %>%
  select(UniqueID) %>% pull() %>% c(man_remove) %>% unique()

### Check mono-inoculations to see if they look right ###
allDat_raw %>%
  filter(plant=="col0") %>%
  filter(readtime==7) %>%
  filter(strain2=="MOCK") %>%
  filter(plate_filled_with !="MIX") %>%
  select(experiment, plate, strain1, fluor_strain1, crim_top, neon_top) %>%
  mutate(strain1_fluor_value = ifelse(fluor_strain1=="crim", crim_top, neon_top)) %>%
  ggplot() +
  geom_point(aes(x=strain1, y=log10(strain1_fluor_value), col=fluor_strain1)) +
  facet_grid(experiment ~ plate)
# SHould be clustered into distinct groups (presence/absence); not bi-model. Confirmed they look decent.

allDat_raw_adj <- allDat_raw %>%
  full_join(blank_values) %>%
  mutate(fluor_offset = lm_offset[1,1]) %>%
  unite(experiment, plate, row, col, col=UniqueID, remove=FALSE) %>%
  filter(!UniqueID %in% contaminated_blanks) %>%
  mutate(crim_raw_blanked = crim_top -crim_top_blanks
         , neon_raw_blanked = neon_top-neon_top_blanks) %>%
  mutate(
    ratio_cn_fluor = crim_raw_blanked/neon_raw_blanked ## TRY THIS-- they are deleting a lot of stuff accidentally if you adjust for bg fluor
  ) %>%
  # Log 2 fold change for crim vs neon
  mutate(CN_FC_fluor = log2(ratio_cn_fluor)
         , CN_FC_fluor_offset = CN_FC_fluor -fluor_offset) %>%
  ### Calculate CFU for those that have it
  rowwise() %>%
  #S12 means Strain 1- Strain 2
  mutate(S12_FC_fluor = ifelse(!is.na(fluor_strain1), 
                               ifelse(fluor_strain1=="crim", CN_FC_fluor, -CN_FC_fluor)
                               ,ifelse(!is.na(fluor_strain2), ifelse(fluor_strain2=="crim", -CN_FC_fluor, CN_FC_fluor), NA))
  ) %>%
  mutate(S12_FC_fluor_offset = ifelse(!is.na(fluor_strain1), 
                                      ifelse(fluor_strain1=="crim", CN_FC_fluor_offset, -CN_FC_fluor_offset)
                                      ,ifelse(!is.na(fluor_strain2), ifelse(fluor_strain2=="crim", -CN_FC_fluor_offset, CN_FC_fluor_offset), NA))
  ) 

dat_plant <- allDat_raw_adj %>%
  filter( plant == "col0") %>% # remove empty wells
  filter(readtime==7) %>% # remove time point 0 reads
  # Here, we remove ones that w and w; or n and n. These were controls, we're not going to plot them below
  filter(!(strain1=="WCS365" & strain2=="WCS365" ), !(strain1=="N2C3" & strain2=="N2C3") ) %>%
  # Mutate the strains so we know which is which
  rowwise() %>%
  mutate(
    WN_FC_fluor = ifelse(strain1=="WCS365", S12_FC_fluor, 
                         ifelse(strain2=="WCS365", -S12_FC_fluor, 
                                ifelse(strain1=="N2C3", -S12_FC_fluor, 
                                       ifelse(strain2=="N2C3", S12_FC_fluor, NA)))),
    WN_FC_fluor_offset = ifelse(strain1=="WCS365", S12_FC_fluor_offset, 
                                ifelse(strain2=="WCS365", -S12_FC_fluor_offset, 
                                       ifelse(strain1=="N2C3", -S12_FC_fluor_offset, 
                                              ifelse(strain2=="N2C3", S12_FC_fluor_offset, NA))))
  ) %>%
  ungroup() %>%
  mutate(strain2 = ifelse(is.na(strain2), "No second strain", strain2)) %>%
  mutate(strain2 = factor(strain2, levels=c("MOCK","WCS365","N2C3", "No second strain"))
         , strain1 = factor(strain1, levels=c("MOCK","WCS365","N2C3")))

### Confirm data looks okay ###
## First, look at raw values when strain 1 is MOCK
dat_plant %>%
  select(experiment, strain1, strain2, fluor_strain1, fluor_strain2, crim_top, neon_top) %>%
  filter(strain1=="MOCK") %>%
  ggplot() +
  geom_point(aes(x=fluor_strain2, y=log10(crim_top)) ) +
  facet_grid(.~experiment)
## Secpmd, look at raw values when strain 2 is MOCK
dat_plant %>%
  select(experiment, strain1, strain2, fluor_strain1, fluor_strain2, crim_top, neon_top) %>%
  filter( strain2=="MOCK" | is.na(strain2)) %>%
  ggplot() +
  geom_point(aes(x=fluor_strain1, y=log10(neon_top))) +
  facet_grid(.~experiment)

### Now, let's look at WCS:N2C3 fold change between all treatments
dat_plant %>%
  ggplot() +
  geom_jitter(aes(x=strain1, y=WN_FC_fluor_offset, col=strain2), width=0.2, height=0)+
  facet_grid(experiment~transfer_h, labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") +
  labs(col="Second inoculated strain") +
  scale_color_manual(values=c(MOCK="black", N2C3="orange", WCS365="green", `No second strain`="grey"))+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

# Version with plant health; collapse all time transfers so it's not an insane sized plot
dat_plant %>%
  ggplot(aes(x=strain1, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Healthiness_hsv, pch=factor(transfer_h)), width=0.1, height=0) +
  facet_grid(.~strain2, labeller = labeller(strain2=c(N2C3="Second strain N2C3", WCS365="Second strain WCS365", MOCK="No second strain"))) +
  scale_color_gradient(low="yellow", high="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") 



gg_ratio_planthealth_exp2 <- dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
  facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
             , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
  scale_shape_manual(values=c(15,19))+
  scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
  geom_hline(aes(yintercept=0)) +
  ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_ratio_planthealth_exp2
ggsave("01_compile_data_forleah/gg_ratio_planthealth_exp2.png", gg_ratio_planthealth_exp2, height=5, width=8)


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
ggsave("01_compile_data_forleah/gg_ratio_planthealth_exp5.png", gg_ratio_planthealth_exp5, height=6, width=10)


###### Save data #########
# Saving as RData file means you keep things like "factor" levels-- useful for plotting. 
# BUT you cannot us RData files between computers (bc sometimes diff versions of R)!! So can't "share" RData files very easily.
# Saving as a txt means you can share with anyone and the data will be the same. But need to re-format dates etc
allDat_final <- allDat_raw_adj

save(allDat_final, file="01_compile_data_forleah/allDat_final.RData")
write.table(allDat_final, file="01_compile_data_forleah/allDat_final.RData", sep="\t", row.names = FALSE, quote=FALSE)


save(dat_plant, file="01_compile_data_forleah/dat_plant.RData")
write.table(dat_plant, file="01_compile_data_forleah/dat_plant.RData", sep="\t", row.names = FALSE, quote=FALSE)

###
gg_doubledip <- allDat_final %>% filter(experiment == "2024-02-14_doubledip") %>%
  filter(plant=="col0") %>%
  mutate(StrainMix2 = ifelse(is.na(strain2_od), gsub("-.*$","",StrainMix), StrainMix)) %>%
  mutate(StrainMixsimple = ifelse(StrainMix == "MOCK-MOCK", "MOCK"
                                  , ifelse(StrainMix == "WCS365-MOCK", "Dipped in WCS365;\nplaced in MOCK well",
                                           ifelse(StrainMix == "N2C3-MOCK", "Dipped in N2C3;\nplaced in MOCK well"
                                                  , ifelse(StrainMix == "WCS365-N2C3", "Dipped in WCS365;\nthen dipped in N2C3;\nplaced in MOCK well"
                                                           , ifelse(StrainMix == "N2C3-WCS365", "Dipped in N2C3;\nthen dipped in WCS365;\nplaced in MOCK well")))))) %>%
 mutate(StrainMixsimple = factor(StrainMixsimple, levels=c(
   "MOCK"
   , "Dipped in WCS365;\nplaced in MOCK well"
   , "Dipped in N2C3;\nplaced in MOCK well"
   , "Dipped in WCS365;\nthen dipped in N2C3;\nplaced in MOCK well"
   ,"Dipped in N2C3;\nthen dipped in WCS365;\nplaced in MOCK well"
 ))) %>%
   mutate(DippedTwice = !is.na(strain2_od)) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK", "WCS365-MOCK", "N2C3-MOCK", "WCS365-N2C3", "N2C3-WCS365"))) %>%
  ggplot(aes(x=StrainMixsimple, y=Healthiness_hsv)) +
  geom_boxplot() +
  geom_jitter( height=0) +
  geom_hline(aes(yintercept=400)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  # scale_color_manual(values=c("lightblue","darkblue")) +
  ylab("Health score") + xlab ("Dipping treatment")
gg_doubledip

ggsave("02_compile_data/gg_doubledip.png", gg_doubledip, height=4, width=6)
# allDat_final %>% filter(experiment == "2024-02-14_doubledip") %>%
#   filter(plant=="col0") %>%
#   mutate(StrainMix2 = ifelse(is.na(strain2_od), gsub("-.*$","",StrainMix), StrainMix)) %>%
#   mutate(DippedTwice = !is.na(strain2_od)) %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK", "WCS365-MOCK", "N2C3-MOCK", "WCS365-N2C3", "N2C3-WCS365"))) %>%
#   mutate(Alive = Healthiness_hsv>400) %>%
#   group_by(StrainMix) %>% summarise(propAlive = sum(Alive)/n()) %>% ungroup() %>%
#   ggplot(aes(x=StrainMix, y=propAlive)) +
#   geom_bar(stat="identity", position = position_stack()) 
