#!bin/bash

###### Step one: compile data ###########

library(lubridate) # Better to manage dates with
library(tidyverse) # data wrangling with dplyr package + ggplot2 package for plotting
# NOTE: always load tidyverse last because the command "select" is a common function name and can be masked by other packages

# I name the output directory the same name as the R script so I can match code to outputs easier
dir.create("02_compile_data")

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
allMeta  <- filter(allMeta, !is.na(plate))  %>% 
  select(-readtime) %>% # remove readtime bc of merging issues
  unite(row, col, col="well", remove=FALSE, sep="") # make a manual well column

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

#### Get CFU data ####

cfu_filenames <- list.files("00_raw_data/dat/cfu_data/")
# Iterate through vector of filenames and load them in, joining with each other.
allCFU <- data.frame(experiment=NA) # Empty dataframe must have at least one common column name as files you're going to load in
for ( f in cfu_filenames) {
  # f = fluor_filenames[2]
  tempCFU <- read.csv(paste0("00_raw_data/dat/cfu_data/",f)) # Adding the folder to filepath with paste0
  allCFU <- full_join(allCFU, tempCFU)
}
# remove all lines where plate is NA
allCFU  <- filter(allCFU, !is.na(experiment)) 

######## Get fluorescence offset ############
lm_offset <- read.delim("../Analysispath_fluoranalysis/02_compile_data/offest_average.txt")

####### Adjust files and process/add additional info ############
# Conversion for OD/abs to cells #
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
  unite(experiment, plate, well, col="UniqueID", sep="__", remove=FALSE) %>%
  filter(readtime !=3)
### REMOVE THINGS THAT I KNOW ARE WRONG ###
allDat <- allDat %>%
  filter(strain2_od !=0.002) %>%
  rowwise() %>%
  mutate(remove = ifelse(length(grep("Strain 2 was innoculated incorectly", notes))>0, FALSE, TRUE)) %>% ungroup() %>%
  filter(remove) %>% select(-remove)

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

gg_platemap_3 <- allDat %>%
  filter(experiment=="2024-02_Peffect_Trial3") %>%
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
gg_platemap_3


gg_platemap_4 <- allDat %>%
  filter(experiment=="2024-02_Peffect_Trial4") %>%
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
gg_platemap_4


# gg_platemap_dd <- allDat %>%
#   filter(experiment=="2024-02-14_doubledip") %>%
#   # select(plate, row, col, strain1, strain2,Healthiness_hsv,crim_btm, neon_btm ) %>% View()
#   # select(plate) %>% table()
#   # select(plate, row, col, strain1, strain2, ) %>%
#   # View()
#   # filter(plant !="NOPLANT") %>%
#   mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
#   ggplot() +
#   geom_tile(aes(x=col, y=(row), fill=strain1)) +
#   geom_point(aes(x=col, y=(row), col=strain2)) +
#   facet_grid(.~plate)+
#   scale_fill_manual(values=c(MOCK="grey", WCS365="darkgreen", N2C3="orange"), na.value = "white")+
#   scale_color_manual(values=c(MOCK="white", WCS365="darkgreen", N2C3="orange"), na.value = "black")
# gg_platemap_dd

## Visual check of all fluorescence
### ALL first strains were neon


allDat %>%
  filter(readtime==7) %>%
  # filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  # select(StrainMix) %>% View()
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK",
                       ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis neon"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=log10(neon_btm+1))) +
  geom_point(aes(x=col, y=(row), col=mock)) +
  facet_wrap(.~plate)+
  scale_fill_gradient(low="white", high="green")+
  scale_color_manual(values=c("grey","darkgreen","blue"))

allDat %>%
  filter(readtime==7) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK",
                       ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis neon"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=log10(crim_btm+1))) +
  geom_point(aes(x=col, y=(row), col=mock)) +
  facet_wrap(.~plate)+
  scale_fill_gradient(low="white", high="red")+
  scale_color_manual(values=c("grey","darkgreen","black"))

allDat %>%
  filter(readtime==7) %>%
  filter(plant !="NOPLANT") %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK",
                       ifelse(strain1=="MOCK", "First strain\nis mock", "First strain\nis neon"))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=Healthiness_hsv)) +
  geom_point(aes(x=col, y=(row), col=strain1)) +
  facet_wrap(.~plate)+
  scale_fill_gradient(low="white", high="darkgreen")
  # scale_color_manual(values=c("grey","darkgreen","black"))

## Look at top 
# LOOKING FOR CONTAMINATION: basically any white dots that are on dark neon colours
allDat %>%
  filter(readtime==7) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK","Not mock")) %>%
  ggplot() +
  geom_tile(aes(x=col, y=row, fill=log10(neon_top+1))) +
  geom_point(aes(x=col, y=row, col=mock)) +
  facet_wrap(.~plate)+
  scale_fill_gradient(low="white", high="green")+
  scale_color_manual(values=c("white","black"), na.value = "black")
# LOOKING FOR CONTAMINATION: basically any white dots that are on dark red colours
allDat %>%
  filter(readtime==7) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  mutate(mock = ifelse(StrainMix=="MOCK-MOCK","MOCK","Not mock")) %>%
  ggplot() +
  geom_tile(aes(x=col, y=row, fill=log10(crim_top+1))) +
  geom_point(aes(x=col, y=row, col=mock)) +
  facet_wrap(.~plate)+
  scale_fill_gradient(low="white", high="darkred")+
  scale_color_manual(values=c("white","black"), na.value = "black")

# MANUALLY REMOVE:
man_remove <- c("2023-12_priorityeffects_Nov-10th-transfer-PlateD_cropped_G_11"
                ,"2023-12_priorityeffects_Nov-10th-transfer-PlateD_cropped_H_10"
                ,"2023-12_priorityeffects_Nov-9th-transfer-PlateC_cropped_G_8"
                ,"2023-12_priorityeffects_Nov-9th-transfer-PlateC_cropped_F_8" ## E and F look weird in follow up sanity checks
                ,"2023-12_priorityeffects_Nov-9th-transfer-PlateC_cropped_E_8"
                , "2024-02_Peffects_Trial3_PlateC_T3h_Trial4_G_2")


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
# Let's look at green values that are >50, <100
allDat %>%
  filter(pix_g_median<50)

# Filter so that all green pixels > 100, and red pixels must be >50
allDat_raw <- allDat %>%
  # The weird transferh 6 stuff is to keep in data where there was a dropped plate
  filter(((pix_g_median>100| pix_r_median>50)|(experiment=="2024-02_Peffect_Trial4" & transfer_h==6 & strain2_od == 0.001)) | plant=="NOPLANT")


##### Split dataset into halves #######
# remove non plants
# allDat_plants <- allDat_raw %>%
#   filter(plant=="col0")
# allDat_noplants <- allDat_raw %>%
#   filter(plant=="NOPLANT")

#### Top vs Bottom ####
gg_topvsbottom <- allDat_raw %>%
  ggplot() +
  geom_point(aes(x=log10(crim_top), y=log10(crim_btm)), col="darkred") +
  geom_point(aes(x=log10(neon_top), y=log10(neon_btm)), col="green") +
  geom_abline(aes(intercept=0, slope=1))
gg_topvsbottom
ggsave("02_compile_data/gg_topvsbottom.png", gg_topvsbottom, height=4, width=6)

## Look at starting values, see which one is more variable
strain1only <- allDat_raw %>% 
  filter(readtime==7, plant=="col0")%>%
  filter(!is.na(strain1), !is.na(strain2), !is.na(fluor_strain1), !is.na(fluor_strain2)) %>%
  mutate(fluor_top_new = ifelse(fluor_strain1 == "crim", crim_top, 
                               ifelse(fluor_strain1 == "neon", neon_top, NA)),
         fluor_btm_new = ifelse(fluor_strain1 == "crim", crim_btm, 
                                ifelse(fluor_strain1 == "neon", neon_btm, NA)),
         fluor = fluor_strain1, 
         strain = strain1,
         oneortwo="strain1") %>%
  select(fluor_top_new, fluor_btm_new, fluor, strain, oneortwo,StrainMix)
strain2only <- allDat_raw %>% 
  filter(readtime==7, plant=="col0")%>%
  filter(!is.na(strain1), !is.na(strain2), !is.na(fluor_strain1), !is.na(fluor_strain2)) %>%
  mutate(fluor_top_new = ifelse(fluor_strain2 == "crim", crim_top, 
                                ifelse(fluor_strain2 == "neon", neon_top, NA)),
         fluor_btm_new = ifelse(fluor_strain2 == "crim", crim_btm, 
                                   ifelse(fluor_strain2 == "neon", neon_btm, NA)),
         fluor = fluor_strain2,
         strain = strain2,
         oneortwo = "strain2") %>%
  select(fluor_top_new, fluor_btm_new, fluor, strain, oneortwo, StrainMix)

gg_fluortopvsbottom_detailed <- bind_rows(strain1only, strain2only) %>%
  ggplot() +
  geom_point(aes(x=log10(fluor_top_new), y=log10(fluor_btm_new), col=fluor)) +
  # geom_point(aes(x=log10(neon_top), y=log10(neon_btm)), col="green") +
  geom_abline(aes(intercept=0, slope=1))+
  facet_grid(strain ~ StrainMix)+
  scale_color_manual(values=c(crim="darkred", neon="green"))+
  xlab("Fluor from top (log10)") + ylab("Fluor from bottom (log10)")
gg_fluortopvsbottom_detailed
ggsave("02_compile_data/gg_fluortopvsbottom_detailed.png", gg_fluortopvsbottom_detailed, height=6, width=10)

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
# Looks better; need to filter out certain ones

### Now calculate blanks
# filter(experiment %in% c("2023-09-13_fluor","2023-09-14_fluor")) %>%
# filter(all_plant_pixels>100) %>%
allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  select(experiment, crim_top, neon_top) %>%
  pivot_longer(-c(experiment), names_to="fluor", values_to="raw") %>%
  ggplot() +
  geom_point(aes(x=fluor, y=log10(raw))) +
  facet_wrap(.~experiment)

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

## Get st curves
# 
# allConv_wide <- allConv %>%
#   pivot_wider(names_from=fluor, values_from = c(od_conv_slope, od_conv_intercept) ) %>%
#   filter(experiment != "allLog10_winter") %>%
#   rename(stcurve= experiment) %>%
#   mutate(fc_cn = log2(od_conv_slope_crim/od_conv_slope_neon)) %>%
#   mutate(MS=0.5, MES=0.5)
# 
# allConv_bg_wide <- allConv_bg %>%
#   filter(version=="winter") %>%
#   pivot_wider(names_from = fluor, values_from=c("bg_conv_intercept","bg_conv_slope"))%>%
#   mutate(MS=0.5, MES=0.5)

### Check mono-inoculations ###
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
         # ,CN_FC_estAbs = log2(ratio_cn_estAbs)
         , CN_FC_fluor_offset = CN_FC_fluor -fluor_offset) %>%
  ### Calculate CFU for those that have it
  rowwise() %>%
  mutate(S12_FC_fluor = ifelse(!is.na(fluor_strain1), 
                               ifelse(fluor_strain1=="crim", CN_FC_fluor, -CN_FC_fluor)
                               ,ifelse(!is.na(fluor_strain2), ifelse(fluor_strain2=="crim", -CN_FC_fluor, CN_FC_fluor), NA))
                               ) %>%
  mutate(S12_FC_fluor_offset = ifelse(!is.na(fluor_strain1), 
                                ifelse(fluor_strain1=="crim", CN_FC_fluor_offset, -CN_FC_fluor_offset)
                                ,ifelse(!is.na(fluor_strain2), ifelse(fluor_strain2=="crim", -CN_FC_fluor_offset, CN_FC_fluor_offset), NA))
  ) %>% ungroup() %>%
  mutate(experiment = ifelse(experiment == "2023-12_priorityeffects", "2023-11-08",
                             ifelse(experiment == "2023-12_priorityeffects2", "2023-12-06", 
                                    ifelse(experiment == "2024-02_Peffect_Trial3", "2024-01-30", 
                                           ifelse(experiment == "2024-02_Peffect_Trial4", "2024-02-14", experiment)))))

dat_plant_all <- allDat_raw_adj %>%
  filter( plant == "col0") %>%
  filter(readtime==7) %>%
  filter(!(strain1=="WCS365" & strain2=="WCS365" ), !(strain1=="N2C3" & strain2=="N2C3") ) %>%
  # Mutate the strains so we know which is which
  rowwise() %>%
  mutate(
    ###### FIX THIS-- line up wcs cs n2c3
    WN_FC_fluor = ifelse(strain1=="WCS365", S12_FC_fluor, 
                         ifelse(strain2=="WCS365", -S12_FC_fluor, 
                                ifelse(strain1=="N2C3", -S12_FC_fluor, 
                                       ifelse(strain2=="N2C3", S12_FC_fluor, NA)))),
    WN_FC_fluor_offset = ifelse(strain1=="WCS365", S12_FC_fluor_offset, 
                         ifelse(strain2=="WCS365", -S12_FC_fluor_offset, 
                                ifelse(strain1=="N2C3", -S12_FC_fluor_offset, 
                                       ifelse(strain2=="N2C3", S12_FC_fluor_offset, NA)))),
  ) %>%
  ungroup() %>%
  mutate(strain2 = ifelse(is.na(strain2), "No second strain", strain2)) %>%
  mutate(strain2 = factor(strain2, levels=c("MOCK","WCS365","N2C3", "No second strain"))
         , strain1 = factor(strain1, levels=c("MOCK","WCS365","N2C3")))
# Separate into double dip and not
dat_plant_dd <- dat_plant_all %>%
  filter(Type=="doubledip")
dat_plant <- dat_plant_all %>%
  filter(Type!="doubledip"| is.na(Type))

### SANITY CHECK; not adjusted yet, but will be in future

## First, look at raw values
dat_plant %>%
  select(experiment, strain1, strain2, fluor_strain1, fluor_strain2, crim_top, neon_top) %>%
  filter(strain1=="MOCK") %>%
  ggplot() +
  geom_point(aes(x=fluor_strain2, y=log10(crim_top)) ) +
  facet_grid(strain2~experiment)
# Two contaminants? Removed them above (plate C, E8, F8)

dat_plant %>%
  select(experiment, strain1, strain2, fluor_strain1, fluor_strain2, crim_top, neon_top) %>%
  filter( strain2=="MOCK" | is.na(strain2)) %>%
  ggplot() +
  geom_point(aes(x=fluor_strain1, y=log10(neon_top))) +
  facet_grid(.~experiment)

### just at the one plate
# S12
# dat_plant %>%
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2, transfer_h, WN_FC_fluor, WN_FC_CFU) %>%
#   # distinct() %>%
#   filter(experiment=="2023-12_priorityeffects2") %>%
#   ggplot() +
#   geom_jitter(aes(x=strain1, y=S12_FC_fluor_offset, pch=strain2, col=Alive), width=0.2, height=0)+
#   facet_grid(experiment~transfer_h, labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   # ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") +
#   # labs(col="Second inoculated strain") +
#   # scale_color_manual(values=c(MOCK="black", N2C3="orange", WCS365="green", `No second strain`="grey"))+
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
# 
# 
# dat_plant %>%
#   # select(UniqueID, strain1, strain2, transfer_h, WN_FC_fluor, WN_FC_CFU) %>%
#   # distinct() %>%
#   # filter(experiment=="2023-12_priorityeffects") %>%
#   ggplot() +
#   geom_jitter(aes(x=strain1, y=WN_FC_fluor_offset, col=strain2), width=0.2, height=0)+
#   facet_grid(experiment~transfer_h, labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") +
#   labs(col="Second inoculated strain") +
#   scale_color_manual(values=c(MOCK="black", N2C3="orange", WCS365="green", `No second strain`="grey"))+
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
# 
# dat_plant %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   ggplot(aes(x=strain1, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Healthiness_hsv, pch=factor(transfer_h)), width=0.1, height=0) +
#   facet_grid(.~strain2) +
#   scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") 
# 
# dat_plant %>%
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   filter(experiment=="2023-12_priorityeffects2") %>%
#   ggplot(aes(x=strain1, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Healthiness_hsv), width=0.1, height=0) +
#   facet_grid(transfer_h~strain2) +
#   # scale_shape_manual(values=c(`TRUE` = 19, `FALSE` =21))+
#   # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
#   scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") 

# 
# dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset)) %>%
#   # filter(StrainMix=="MOCK-MOCK")
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   filter(experiment=="2023-12_priorityeffects2") %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Healthiness_hsv), width=0.1, height=0) +
#   facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x") +
#   # scale_shape_manual(values=c(`TRUE` = 19, `FALSE` =21))+
#   # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
#   scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("First inoculated strain") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


# 
# gg_ratio_planthealth_exp2 <- dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset)) %>%
#   # filter(StrainMix=="MOCK-MOCK")
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   # filter(experiment=="2023-12_priorityeffects2") %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
#   ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Alive, pch=experiment), width=0.1, height=0) +
#   facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
#              , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   # scale_shape_manual(values=c(15,19))+
#   scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
#   # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
#   # scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# gg_ratio_planthealth_exp2
# ggsave("02_compile_data/gg_ratio_planthealth_exp2.png", gg_ratio_planthealth_exp2, height=5, width=8)
# 
# 
# gg_ratio_planthealth_exp3 <- dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset)) %>%
#   # filter(StrainMix=="MOCK-MOCK")
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   # filter(experiment=="2023-12_priorityeffects2") %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
#   ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Healthiness_hsv, pch=experiment), width=0.1, height=0) +
#   facet_grid(.~transfer_h, scales = "free", drop=TRUE, space="free_x"
#              , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   # scale_shape_manual(values=c(15,19))+
#   # scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
#   # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
#   scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# gg_ratio_planthealth_exp3
# ggsave("02_compile_data/gg_ratio_planthealth_exp3.png", gg_ratio_planthealth_exp3, height=5, width=8)
# # 
# gg_ratio_planthealth_exp4 <- dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset)) %>%
#   # filter(StrainMix=="MOCK-MOCK")
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   # filter(experiment=="2023-12_priorityeffects2") %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
#   ggplot(aes(x=StrainMix, y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Healthiness_hsv, pch=experiment), width=0.1, height=0) +
#   facet_grid(experiment~transfer_h, scales = "free", drop=TRUE, space="free_x"
#              , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   # scale_shape_manual(values=c(15,19))+
#   # scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
#   # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
#   scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# gg_ratio_planthealth_exp4
# ggsave("02_compile_data/gg_ratio_planthealth_exp4.png", gg_ratio_planthealth_exp4, height=6, width=10)
# 
# # 
# gg_ratio_planthealth_exp5 <- dat_plant %>%
#   filter(!is.na(WN_FC_fluor_offset)) %>%
#   # filter(StrainMix=="MOCK-MOCK")
#   mutate(Alive = Healthiness_hsv>400) %>%
#   # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
#   # distinct() %>%
#   # filter(experiment=="2023-12_priorityeffects2") %>%
#   mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
#   filter(WN_FC_fluor_offset<15, WN_FC_fluor_offset>-15) %>%
#   ggplot(aes(x=factor(transfer_h), y=WN_FC_fluor_offset)) +
#   geom_boxplot() +
#   geom_jitter(aes(col=Healthiness_hsv, pch=experiment), width=0.1, height=0) +
#   facet_grid(experiment~StrainMix, scales = "free", drop=TRUE, space="free_x"
#              , labeller = labeller(transfer_h=c(`0`="Simultaneous Inoc",`0.01`="Dip", `3`="3h apart", `6`="6h apart",  `24`="24h apart", `48`="48h apart"))) +
#   # scale_shape_manual(values=c(15,19))+
#   # scale_color_manual(values=c(`TRUE`="green", `FALSE` = "orange")) +
#   # scale_color_manual(values=c("grey","yellow","green","darkgreen","blue")) +
#   scale_color_gradient(low="yellow", high="darkgreen") +
#   geom_hline(aes(yintercept=0)) +
#   ylab("WCS365 : N2C3\n Log2 Fold change")+xlab("Firsttrain-SecondStrain") +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# gg_ratio_planthealth_exp5
# ggsave("02_compile_data/gg_ratio_planthealth_exp5.png", gg_ratio_planthealth_exp5, height=6, width=10)
# 
# 

dat_plant %>%
  filter(!is.na(WN_FC_fluor_offset)) %>%
  # filter(StrainMix=="MOCK-MOCK")
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  # select(UniqueID, strain1, strain2,fluor_strain1, fluor_strain2, transfer_h,S12_FC_fluor, WN_FC_fluor, WN_FC_CFU, Healthiness_hsv, ratio) %>%
  # distinct() %>%
  filter(experiment=="2023-12-06") %>%
  filter(StrainMix %in% c("WCS365-N2C3","N2C3-WCS365"), transfer_h>0) %>%
  mutate(StrainMix = factor(StrainMix, levels=c("MOCK-MOCK","WCS365-MOCK","MOCK-WCS365","N2C3-MOCK","MOCK-N2C3","WCS365-N2C3","N2C3-WCS365"))) %>%
  ggplot(aes(x=WN_FC_fluor_offset, y=Alive)) +
  geom_jitter(aes(col=StrainMix), width=0, height=0.2)+
  geom_smooth(method="glm",
              method.args=list(family="binomial"), col="black")+
  # facet_grid(transfer_h~StrainMix)
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

###### CFU processing #######

allCFU_adj <- allCFU %>%
  mutate(experiment = ifelse(experiment == "2023-12_priorityeffects", "2023-11-08",
                             ifelse(experiment == "2023-12_priorityeffects2", "2023-12-06", 
                                    ifelse(experiment == "2024-02_Peffect_Trial3", "2024-01-30", 
                                           ifelse(experiment == "2024-02_Peffect_Trial4", "2024-02-14", experiment))))) 
# REMOVE ALL EXPERIMENTS TO OMIT
onCFU_summary <- allCFU_adj  %>% #filter(experiment %in% c("2024-02-14")) %>%
  filter(plant=="overnight") %>%
  unite(experiment, transfer_h, col="exp_trans", remove=TRUE) %>%
  mutate(batchtest = ifelse(exp_trans=="2024-02-14_6" & abs_600==0.001, "B","A")) %>%
  group_by(batchtest, strain, abs_600) %>%
  summarise(meanCount = mean(CFU_well), sd = sd(CFU_well), rep_on = n()) %>%
  ungroup() %>%
  mutate(treatment = ifelse(is.na(abs_600), "Plant root", paste0("Abs600=", abs_600,"\nsecond strain inoculant"))) %>%
  rename(strain2 = strain, strain2_od = abs_600, strain2_meanCount = meanCount, strain2_sd = sd) 

# Root summary; also add in 
rootCFU_summary <- allCFU_adj %>% filter(plant=="col0")  %>%
  group_by(experiment, strain, transfer_h, abs_600) %>%
  summarise(meanCount = mean(CFU_well), sd = sd(CFU_well), rep_root= n()) %>%
  mutate(abs_600 = ifelse(is.na(abs_600), "root", as.character(abs_600))) %>%
  mutate(treatment = ifelse(is.na(abs_600), "Plant root", paste0("Abs600=", abs_600,"\nsecond strain inoculant"))) %>% ungroup() %>%
  # mutate(batchtest = ifelse(experiment == "2024-02-14" & transfer_h == 6, "B","A")) %>%
  rename(strain1 = strain, strain1_meanCount = meanCount, strain1_sd = sd) %>%
  left_join(dat_plant_all) %>%
  select(experiment, transfer_h, strain1, strain2, strain2_od, WN_FC_fluor_offset, Healthiness_hsv, strain1_meanCount, strain1_sd) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  mutate(batchtest = ifelse(experiment == "2024-02-14" & transfer_h == 6 & (strain2_od==0.001), "B","A")) 

# Check data
onCFU_summary %>%
  ggplot() +
  geom_pointrange(aes(x=(strain2_od), y=log10(strain2_meanCount)
                      , ymin=log10(strain2_meanCount-strain2_sd), ymax=log10(strain2_meanCount+strain2_sd)
                      , col=strain2))+
  # facet_grid(.~strain2) +
  scale_x_log10() +
  xlab("ABS600 of strain 2") + ylab("CFU count of strain2 (log10)") +
  labs(col="Strain 2") + 
  scale_color_manual(values=c(N2C3="darkorange", WCS365 = "darkseagreen4"))

gg_rootCFU <- rootCFU_summary %>% 
  unite(experiment, strain1, col="grp", remove=FALSE ) %>%
  ggplot(aes(x=(transfer_h), y=log10(strain1_meanCount), col=strain1, group=grp)) +
  geom_pointrange(aes(ymin=log10(strain1_meanCount-strain1_sd), ymax=log10(strain1_meanCount+strain1_sd)
                      , col=strain1, pch=experiment)) +
  geom_line(aes()) +
  scale_color_manual(values=c(N2C3="darkorange", WCS365="darkseagreen4")) +
  labs(col="Strain", shape="Experiment") +xlab("Colonization time \n(Strain 1)")+
  ylab("CFUs on root (log10)")
gg_rootCFU
ggsave("02_compile_data/gg_rootCFU.png", gg_rootCFU, height=3, width=5)

### viewing cfu with other data
allCFU_plants <- left_join(rootCFU_summary, onCFU_summary, relationship = "many-to-many") %>%
  # filter(strain1!="MOCK", strain2!="MOCK") %>%
  mutate(ratioCFU_S12 = strain1_meanCount/strain2_meanCount
         , S12_FC_inocCFU = log2(ratioCFU_S12)
         , WN_FC_inocCFU = ifelse(strain1=="WCS365" & strain2 == "N2C3", S12_FC_inocCFU
                                  , ifelse(strain1=="N2C3" & strain2 == "WCS365", -S12_FC_inocCFU, NA))) %>%
  mutate(Treatment= ifelse(strain1=="N2C3" & strain2 == "WCS365", "N2C3 first",
                               ifelse(strain1=="WCS365" & strain2 == "N2C3", "WCS365 first", NA)))
# Check
allCFU_plants %>%
  filter(experiment == "2024-02-14") %>%
  # filter(transfer_h>0) %>%
  ggplot(aes(x=WN_FC_inocCFU, y=WN_FC_fluor_offset, group=Treatment, col=Treatment)) +
  geom_jitter(aes(pch=factor(transfer_h)), height=0.1, width=0.1)  +
  geom_smooth(span=100)

## Get cfu data for both mixes
onCFU_summary1 <- onCFU_summary %>%
  filter(batchtest == "A" & strain2_od!=0.001) %>%
  mutate(batchtest = "B") %>%
  full_join(onCFU_summary) %>%
  select(-c(rep_on, treatment))
onCFU_summary2 <-onCFU_summary1 %>%
    rename(strain1_meanCount = strain2_meanCount
           , strain1_sd = strain2_sd
           , strain1=strain2
           , strain1_od= strain2_od) 

######### WORKING, NOT JOINING PROPERLY, MISSING MIXED #########

allCFU_mixedtreat <- dat_plant_all %>%
  filter(Type=="MIX") %>% 
  mutate(batchtest = ifelse(experiment == "2024-02-14" & transfer_h == 6 & (strain1_od==0.001 | strain2_od==0.001), "B","A")) %>%
  select(experiment, transfer_h, strain1, strain1_od
         , strain2,  strain2_od, WN_FC_fluor_offset, Healthiness_hsv, batchtest, Type) %>%
  left_join(onCFU_summary1) %>%
  left_join(onCFU_summary2) %>%
  mutate(ratioCFU_S12 = strain1_meanCount/strain2_meanCount
                  , S12_FC_inocCFU = log2(ratioCFU_S12)
                  , WN_FC_inocCFU = ifelse(strain1=="WCS365" & strain2 == "N2C3", S12_FC_inocCFU
                                           , ifelse(strain1=="N2C3" & strain2 == "WCS365", -S12_FC_inocCFU, NA))) %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  mutate(Treatment = "Simultaneous inoculation") %>%
  select(one_of(colnames(allCFU_plants))) 

## 

allCFU_data_processed <- full_join(allCFU_mixedtreat, allCFU_plants) %>%
  filter(strain1!=strain2)

###### Save data #########
# Saving as RData file means you keep things like "factor" levels-- useful for plotting. 
# BUT you cannot us RData files between computers (bc sometimes diff versions of R)!! So can't "share" RData files very easily.
# Saving as a txt means you can share with anyone and the data will be the same. But need to re-format dates etc
allDat_final <- allDat_raw_adj 

save(allDat_final, file="02_compile_data/allDat_final.RData")
write.table(allDat_final, file="02_compile_data/allDat_final.RData", sep="\t", row.names = FALSE, quote=FALSE)


save(dat_plant, file="02_compile_data/dat_plant.RData")
write.table(dat_plant, file="02_compile_data/dat_plant.txt", sep="\t", row.names = FALSE, quote=FALSE)


save(dat_plant_dd, file="02_compile_data/dat_plant_dd.RData")
write.table(dat_plant_dd, file="02_compile_data/dat_plant_dd.txt", sep="\t", row.names = FALSE, quote=FALSE)

save(allCFU, file = "02_compile_data/allCFU.RData")
write.table(allCFU, file="02_compile_data/allCFU.txt", sep="\t", row.names = FALSE, quote=FALSE)


save(allCFU_data_processed, file="02_compile_data/allCFU_data_processed.RData")
write.table(allCFU_data_processed, file="02_compile_data/allCFU_data_processed.txt", sep="\t", row.names = FALSE, quote=FALSE)

#### Double dip peak ####
