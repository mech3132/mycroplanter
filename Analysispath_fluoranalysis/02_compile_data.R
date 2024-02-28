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
allMeta  <- filter(allMeta, !is.na(plate))





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
allFluor <- data.frame(plate_alias=NA) # Empty dataframe must have at least one common column name as files you're going to load in
for ( f in fluor_filenames) {
  tempfluor <- read.csv(paste0("00_raw_data/dat/fluor_data/",f)) # Adding the folder to filepath with paste0
  allFluor <- full_join(allFluor, tempfluor)
}
# remove all lines where plate is NA
allFluor  <- filter(allFluor, !is.na(plate_alias))

######## Get fluorescence standard curves ############
# Before, I was converting by standard curves, but now we're using the constant method so I'm not doing that anymore
# allConv <- read.delim("01_standard_curve/conversion_models.txt")
# allConv_bg <- read.delim("01_standard_curve/conversion_models_bg.txt")


####### Adjust files and process/add additional info ############
# Conversion for OD/abs to cells #
# OD = 3.166*abs - 0.268 (but only When 0.3>Abs>0.1)--> formual derived from experiment
# Cells per ml for OD1 is 5e8
# cells_per_Abs1 = (3.166*1-0.268)*5E8 ## OLD ONE, including all Abs600
cells_per_Abs1 = (3.19*1-0.275)*5E8 ## NEW ONE, excluding Abs600<0.125. This equation is from a different script output

allMeta_adj <- allMeta %>%
  # Make all dates lubridate dates, which are easier to handle than basic "R" dates
  mutate(date_germ = as.Date(date_germ)
         , date_ster = as.Date(date_ster)
         , date_inoc = as.Date(date_inoc)
         , date_process = as.Date(date_process)) %>%
  # Separate the path and protect names into lacZ and fluor info
  separate(path, into= c("path","other"), sep="-") %>%
  separate(other, into= c("lacZ_path","fluor_path"), sep="_") %>% 
  mutate(lacZ_path=ifelse(lacZ_path=="",NA,lacZ_path),fluor_path=ifelse(fluor_path=="",NA,fluor_path)) %>%
  separate(protect, into= c("protect","other"), sep="-") %>%
  separate(other, into= c("lacZ_protect","fluor_protect"), sep="_") %>% 
  mutate(lacZ_protect=ifelse(lacZ_protect=="",NA,lacZ_protect),fluor_protect=ifelse(fluor_protect=="",NA,fluor_protect))%>%
  mutate(col=factor(col, levels=c(seq(1,12)))
         , protect = factor(protect, levels=c("MOCK","WCS365")) # Factors make it so MOCK is automatically first when running models or making plots (default is alphabetical)
         , path = factor(path, levels=c("MOCK","N2C3"))) %>%
  unite(protect, path, col="StrainMix", sep = "-", remove=FALSE) %>% # make a strainmix column by uniting the protect and path columns
  rowwise() %>% # this tells it to calculate by row, rather that summing across all rows
  mutate(totalOD = protect_od + path_od # Calculate new column that is TOTAL OD
         , ratio_protect = protect_od/min(c(protect_od, path_od)) # ratio of protect:path
         , ratio_path = path_od/min(c(protect_od, path_od)) # ratio of path:protect
         # , ratio_path_adj = path_od_adj/min(c(protect_od, path_od_adj))
  ) %>% 
  # When dividing by zero, remove these.
  mutate(ratio_protect = ifelse(is.infinite(ratio_protect), 1, 
                                ifelse(is.na(ratio_protect), 0, ratio_protect))
         , ratio_path = ifelse(is.infinite(ratio_path), 1, 
                               ifelse(is.na(ratio_path), 0, ratio_path))) %>%
  mutate(ratio = paste0(c(ratio_protect, ratio_path), collapse = "-")) %>% # category of ratios used
  # mutate(ratio_adj = paste0(c(ratio_protect, ratio_path_adj), collapse = "-")) %>%
  ungroup() %>% # undoes the "rowwise" command
  # Below, we calculate the number of cells (and log of number of cells)
  mutate(protect_cells = protect_od*cells_per_Abs1,
         path_cells = path_od*cells_per_Abs1, 
         total_cells = totalOD*cells_per_Abs1,
         protect_cells_log = log10(protect_cells+1),
         path_cells_log = log10(path_cells+1),
         total_cells_log = log10(total_cells+1)) %>%
  mutate(plant_age = date_inoc - date_germ) # Plant age at inoculation


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
  rename(neon_raw = Neon.raw,  crim_raw=Crim.raw) %>%
  rename(neon_CFU = Neon.cfu, crim_CFU = Crim.cfu) %>%
  mutate(col=factor(col, levels=c(seq(1,12)))) # Make this a factor, which is a sequence from 1-12 (otherwise the alphabetical order will be 1, 10, 11, 12, 2, 3, 4...)


##### JOIN and check #####
allDat <- full_join(allMeta_adj,allPixel_adj) %>%
  full_join(allFluor_adj) # Note: with pipling, the first item is assumed to be what you pipe in. 

### Sanity check ###
## Visual check of everything
allDat %>%
  filter(experiment %in% c("2023-08-08_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=crim_raw)) +
  facet_wrap(.~plate)
# THROW OUT ALL MOCKS HERE-- they look contaminated
allDat %>%
  filter(experiment %in% c("2023-08-08_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=neon_raw)) +
  facet_wrap(.~plate)

allDat %>%
  filter(experiment %in% c("2023-08-16_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=crim_raw)) +
  facet_wrap(.~plate)
allDat %>%
  filter(experiment %in% c("2023-08-16_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=neon_raw)) +
  facet_wrap(.~plate)
# "2023-08-16_fluor_scan_NneonWcrim_halflacZ_004_cropped_G_12" looks contaminated

allDat %>%
  filter(experiment %in% c("2023-09-13_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=crim_raw)) +
  facet_wrap(.~plate)
allDat %>%
  filter(experiment %in% c("2023-09-13_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=neon_raw)) +
  facet_wrap(.~plate)


allDat %>%
  filter(experiment %in% c("2023-09-14_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=crim_raw)) +
  facet_wrap(.~plate)
# "2023-09-14_fluor_Ncrim_WneonLacZ_A_001_crop_D_3"
# "2023-09-14_fluor_Ncrim_WneonLacZ_B_002_crop_D_4"
# These two look contaminated


allDat %>%
  filter(experiment %in% c("2023-11-28_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=crim_raw)) +
  facet_wrap(.~plate)
allDat %>%
  filter(experiment %in% c("2023-11-28_fluor")) %>%
  mutate(row=factor(row, levels=rev(sort(unique(row))))) %>%
  ggplot() +
  geom_tile(aes(x=col, y=(row), fill=neon_raw)) +
  facet_wrap(.~plate)



# MANUALLY REMOVE:
man_remove <- c("2023-09-14_fluor_Ncrim_WneonLacZ_A_001_crop_D_3", "2023-09-14_fluor_Ncrim_WneonLacZ_B_002_crop_D_4"
  ,"2023-08-16_fluor_scan_NneonWcrim_halflacZ_004_cropped_G_12")


# Are there any mis-alignments (ie: StrainMix is NA or all_plant_pixels is NA, indicating that some plate/experiment/well combo doesn't line up)
allDat %>%
  filter(is.na(StrainMix) | is.na(all_plant_pixels))
# Should have no rows!
allDat %>% filter(is.na(plate))

# Are there any NAs in StrainMix?
allDat %>% filter((is.na(StrainMix)))
# Are there any NAs in all_plant_pixels?
allDat %>% filter(is.na(all_plant_pixels))
# Make a table of data to see if number of samples makes sense
allDat %>%
  select(StrainMix, protect_od) %>%
  table()

allDat %>%
  select(StrainMix, path_od) %>%
  table()

######### Make a histgram of pixel values to see whether there are any strange outliers ###########

allDat %>%
  ggplot() +
  geom_density(aes(x=pix_g_median), fill=NA, col="green")+
  geom_density(aes(x=pix_r_median), fill=NA, col="red")+
  geom_density(aes(x=pix_b_median), fill=NA, col="blue")
# The values where pixels are close to zero (<50) are actually wells that don't have plants in them (cross-checked)

# Filter so that all green pixels > 100, and red pixels must be >50
# Keeping "no plant" ones for now so we can compare communities when there is no host
allDat_raw <- allDat %>%
  filter((pix_g_median>100| pix_r_median>50)) 

######### PROCESS ################


allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  select(experiment, crim_raw, neon_raw, fluor_protect) %>%
  pivot_longer(-c(experiment,fluor_protect), names_to="fluor", values_to="raw") %>%
  ggplot() +
  geom_point(aes(x=fluor, y=log10(raw), col=fluor_protect)) +
  facet_wrap(.~experiment)

### Remove crimson values over E4; and remove neon values over 4.5
# Looks like there is a lot of contamination  in 8-8 and 8-16
# Let's just use blanks from the the second set to replace the first set blanks


### Look at the raw values of blank wells; they should look clustered together. Remove the ones that are outliers.
allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  mutate(not_contam = (crim_raw<10^4 & neon_raw < 10^4.5)) %>%
  select(experiment, crim_raw, neon_raw, fluor_protect, not_contam) %>%
  pivot_longer(-c(experiment,fluor_protect, not_contam), names_to="fluor", values_to="raw") %>%
  ggplot() +
  geom_jitter(aes(x=fluor, y=log10(raw), col=fluor_protect, pch=not_contam), width=0.1, height=0) +
  facet_wrap(.~experiment) +
  scale_shape_manual(values=c(21, 19))

blank_values <-  allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  filter(crim_raw<10^4, neon_raw < 10^4.5) %>%
  group_by(experiment) %>%
  summarise(crim_raw_blanks=mean(crim_raw), neon_raw_blanks=mean(neon_raw)
  )
blank_values <- rbind(blank_values, data.frame(experiment="2023-08-08_fluor"
                                               , blank_values%>%filter(experiment=="2023-08-16_fluor")%>%select(-experiment)))
contaminated_blanks <-  allDat_raw %>%
  filter(StrainMix=="MOCK-MOCK") %>%
  filter(crim_raw>=10^4 |  neon_raw >= 10^4.5) %>%
  unite(experiment, plate, row, col, col=UniqueID) %>%
  select(UniqueID) %>% pull() %>% c(man_remove) %>% unique()

# ## Don't need this anymore; no longer converting with standard curves
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

allDat_raw_adj <- allDat_raw %>%
  full_join(blank_values) %>%
  # full_join(allConv_wide, relationship = "many-to-many") %>%
  # full_join(allConv_bg_wide) %>%
  unite(experiment, plate, row, col, col=UniqueID, remove=FALSE) %>%
  filter(!UniqueID %in% contaminated_blanks) %>%
  mutate(crim_raw_blanked = crim_raw -crim_raw_blanks
         , neon_raw_blanked = neon_raw-neon_raw_blanks) %>%
  # Since not using standard curves, don't need to adjust with blanked neon values. Commenting out all below
  # mutate(neon_raw_blanked_adj = neon_raw_blanked-(crim_raw_blanked*bg_conv_slope_neon+bg_conv_intercept_neon)) %>%
  # mutate(crim_estAbs = (crim_raw_blanked-od_conv_intercept_crim)/od_conv_slope_crim
         # , neon_estAbs = (neon_raw_blanked_adj-od_conv_intercept_neon)/od_conv_slope_neon) %>%
  # mutate(ratio_cn_fluor = crim_raw_blanked/neon_raw_blanked_adj
         mutate(ratio_cn_fluor = crim_raw_blanked/neon_raw_blanked
                # , ratio_cn_estAbs = crim_estAbs/neon_estAbs
                ) %>%
  # Log 2 fold change for crim vs neon
  mutate(CN_FC_fluor = log2(ratio_cn_fluor)
         # ,CN_FC_estAbs = log2(ratio_cn_estAbs)
         ) %>%
  ### Calculate CFU for those that have it
  rowwise() %>%
  # Below: calculated actualy CFU ratios
  mutate(ratio_cn_CFU = crim_CFU/neon_CFU) %>%
  mutate(CN_FC_CFU = log2(ratio_cn_CFU)) 

### Try fitting the previous curves ###
dat_sanitycheck_raw <- allDat_raw_adj %>%
  filter( experiment %in% c("2023-09-13_fluor")
          ,StrainMix %in% c("MOCK-MOCK","WCS365-MOCK")) %>%
  select(StrainMix, fluor_protect, crim_raw, neon_raw) %>% distinct() %>%
  pivot_longer(-c(StrainMix,fluor_protect), names_to = "fluor_cells",  values_to = "raw")
dat_sanitycheck_blanked <- allDat_raw_adj %>%
  filter( experiment %in% c("2023-09-13_fluor","2023-09-14_fluor","2023-11-28_fluor")
          ,StrainMix %in% c("MOCK-MOCK","WCS365-MOCK")) %>%
  select(StrainMix, fluor_protect, crim_raw_blanked, neon_raw_blanked) %>%distinct() %>%
  pivot_longer(-c(StrainMix,fluor_protect), names_to = "fluor_cells",  values_to = "raw_blanked")
gg_wcs_vs_mock_fluor_raw <- dat_sanitycheck_raw %>%
  ggplot() +
  geom_jitter(aes(x=StrainMix, y=log10(raw), col=fluor_protect), height=0, width=0.1) +
  facet_grid(.~fluor_cells, labeller = labeller(fluor_cells=c(crim_raw = "Crimson ex/em", neon_raw = "Neon ex/em"))) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Raw (blanked) fluorescence\nlog10") + labs(col="WCS365 \nstrain")+
  scale_color_manual(values=c(crim="orange", neon="green"))
gg_wcs_vs_mock_fluor_raw
ggsave("02_compile_data/gg_wcs_vs_mock_fluor_raw.png", gg_wcs_vs_mock_fluor_raw, height=4, width=8)

gg_wcs_vs_mock_fluor_blanked <- dat_sanitycheck_blanked %>%
  ggplot() +
  geom_jitter(aes(x=StrainMix, y=log10(raw_blanked), col=fluor_protect), height=0, width=0.1) +
  facet_grid(.~fluor_cells, labeller = labeller(fluor_cells=c(crim_raw_blanked = "Crimson ex/em", neon_raw_blanked = "Neon ex/em"))) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Raw (blanked) fluorescence\nlog10") + labs(col="WCS365 \nstrain")+
  scale_color_manual(values=c(crim="orange", neon="green"))
gg_wcs_vs_mock_fluor_blanked
ggsave("02_compile_data/gg_wcs_vs_mock_fluor_blanked.png", gg_wcs_vs_mock_fluor_blanked, height=4, width=8)


gg_wcs_vs_mock_fluor_blanked_onlymatching <-  dat_sanitycheck_blanked %>%
  filter(StrainMix %in% c("MOCK-MOCK","WCS365-MOCK")) %>%
  separate(fluor_cells, into=c("fluor","other","other2"), sep="_", remove=FALSE) %>%
  filter(fluor_protect == fluor) %>% select(-c("other","other2")) %>%
  ggplot() +
  geom_jitter(aes(x=StrainMix, y=log10(raw_blanked), col=fluor_protect), height=0, width=0.1) +
  # facet_grid(.~fluor_cells, labeller = labeller(fluor_cells=c(crim_raw_blanked = "Crimson ex/em", neon_raw_blanked = "Neon ex/em"))) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Raw (blanked) fluorescence\nlog10") + labs(col="WCS365 \nstrain")+
  scale_color_manual(values=c(crim="orange", neon="green"))
gg_wcs_vs_mock_fluor_blanked_onlymatching
ggsave("02_compile_data/gg_wcs_vs_mock_fluor_blanked_onlymatching.png", gg_wcs_vs_mock_fluor_blanked_onlymatching, height=4, width=4)

# 
# gg_curveadjusted_WCSvsmock <- allDat_raw_adj %>%
#   filter( experiment %in% c("2023-09-13_fluor")
#           ,StrainMix %in% c("MOCK-MOCK","WCS365-MOCK")) %>%
#   select(stcurve, StrainMix, fluor_protect, crim_estAbs, neon_estAbs) %>%
#   pivot_longer(-c(stcurve, StrainMix,fluor_protect), names_to = "fluor_cells",  values_to = "estAbs") %>%
#   separate(fluor_cells, into=c("fluor","other","other2"), sep="_", remove=FALSE) %>%
#   filter(fluor_protect == fluor) %>% select(-c("other","other2")) %>%
#   ggplot() +
#   geom_jitter(aes(x=StrainMix, y=log10(estAbs), col=fluor_protect), height=0, width=0.1) +
#   facet_grid(.~stcurve, labeller = labeller(fluor_cells=c(crim_raw_blanked = "Crimson ex/em", neon_raw_blanked = "Neon ex/em"))) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   ylab("Estimated Abs600\n(based on different st. curves)") + labs(col="WCS365 \nstrain")+
#   scale_color_manual(values=c(crim="orange", neon="green"))
# gg_curveadjusted_WCSvsmock
# ggsave("02_compile_data/gg_curveadjusted_WCSvsmock.png", gg_curveadjusted_WCSvsmock, height=4, width=8)

## Calculate conversion factor
# Look at CFU vs fluor data
allDat_raw_adj %>%
  # Remove 2023-08-08 because those are both n2c3 
  filter(experiment %in% c("2023-09-13_fluor","2023-09-14_fluor", "2023-11-28_fluor")) %>%
  filter(!is.na(crim_CFU), !is.na(neon_CFU)) %>%
  select(experiment, fluor_protect, fluor_path,CN_FC_fluor, CN_FC_CFU ) %>%distinct() %>%
  filter(is.finite(CN_FC_fluor), is.finite(CN_FC_CFU)) %>%
  ggplot() +
  geom_point(aes(x=CN_FC_CFU, y=CN_FC_fluor, col=experiment))+
  geom_abline(aes(slope=1, intercept=0))

## Get "general" conversion factor for other experiments
# Here, filtering our data to include only those that are within a "reasonable" range of CFU counts 
# Going to remove CFU counts that are outliers here
dat_for_conv <- allDat_raw_adj %>%
  filter(experiment %in% c("2023-09-13_fluor","2023-09-14_fluor")) %>%
  filter(!is.na(CN_FC_fluor), !is.na(CN_FC_CFU), is.finite(CN_FC_fluor), is.finite(CN_FC_CFU)) %>%
  select(experiment, fluor_protect, plate, CN_FC_fluor,CN_FC_CFU ) %>%
  filter(!(CN_FC_CFU>4 & CN_FC_fluor<(-2)), CN_FC_fluor>-5, CN_FC_CFU>-5)
dat_for_conv %>% ggplot() +
  geom_point(aes(x=CN_FC_CFU, y=CN_FC_fluor, col=experiment, pch=fluor_protect)) +
  geom_abline(aes(intercept=0, slope=1))
# View(dat_for_conv)
lm_forconv <- lm(CN_FC_fluor ~ offset(CN_FC_CFU),data=dat_for_conv)
lm_forconv
##  Save this model
write.table(lm_forconv$coefficients, file = "02_compile_data/offest_average.txt", row.names = FALSE, quote = FALSE, sep="\t")

### This is to calculate offsets by exp, but we're not going to do that here; we're going to use same offset for all exp
# yoffset <- allDat_raw_adj %>%
#   filter(!is.na(crim_CFU), !is.na(neon_CFU)) %>%
#   select(experiment, fluor_protect, fluor_path,CN_FC_fluor, CN_FC_CFU ) %>%distinct() %>%
#   filter(is.finite(CN_FC_fluor), is.finite(CN_FC_CFU)) %>%
#   group_by(experiment,fluor_protect) %>%
#   summarise(yoffset_FC = lm(CN_FC_fluor~CN_FC_CFU)$coefficients[1]) %>%
#   select(experiment, fluor_protect, yoffset_FC) %>%
#   ungroup()

allDat_cleaned <- allDat_raw_adj %>%
  # full_join(yoffset) %>%
  mutate(yoffset_ave=coef(lm_forconv)[1]) %>%
  # mutate(CN_FC_fluor_offset = CN_FC_fluor-yoffset_FC)
mutate(CN_FC_fluor_offset = CN_FC_fluor-yoffset_ave)

### Below: some lines are commented out because these were calculated using the standard curves, which we aren't using anymore
# Experiments 1 and 2, where both isolates were probably N2C3
dat_12 <- allDat_cleaned %>%
  filter( experiment %in% c("2023-08-08_fluor","2023-08-16_fluor")) %>%
  filter( plant == "col0") %>%
  mutate(NN_FC_fluor = ifelse(fluor_protect == "crim",CN_FC_fluor, -CN_FC_fluor )
         ,NN_FC_fluor_offset = ifelse(fluor_protect == "crim",CN_FC_fluor_offset, -CN_FC_fluor_offset )
         ,NN_FC_CFU= ifelse(fluor_protect == "crim",CN_FC_CFU, -CN_FC_CFU )
         # ,NN_FC_estAbs = ifelse(fluor_protect == "crim", CN_FC_estAbs, -CN_FC_estAbs)
  ) %>%
  ungroup() 

dat_34 <- allDat_cleaned %>%
  filter( experiment %in% c("2023-09-13_fluor","2023-09-14_fluor")) %>%
  filter( plant == "col0") %>%
  mutate(WN_FC_fluor = ifelse(fluor_protect == "crim",CN_FC_fluor, -CN_FC_fluor )
         ,WN_FC_fluor_offset = ifelse(fluor_protect == "crim",CN_FC_fluor_offset, -CN_FC_fluor_offset )
         ,WN_FC_CFU= ifelse(fluor_protect == "crim",CN_FC_CFU, -CN_FC_CFU )
         # ,WN_FC_estAbs = ifelse(fluor_protect == "crim", CN_FC_estAbs, -CN_FC_estAbs)
  ) %>%
  ungroup() 

dat_5 <- allDat_cleaned %>%
  filter( experiment %in% c("2023-11-28_fluor")) %>%
  filter( plant == "col0") %>%
  mutate(WN_FC_fluor = ifelse(fluor_protect == "crim",CN_FC_fluor, -CN_FC_fluor )
         ,WN_FC_fluor_offset = ifelse(fluor_protect == "crim",CN_FC_fluor_offset, -CN_FC_fluor_offset )
         ,WN_FC_CFU= ifelse(fluor_protect == "crim",CN_FC_CFU, -CN_FC_CFU )
         # ,WN_FC_estAbs = ifelse(fluor_protect == "crim", CN_FC_estAbs, -CN_FC_estAbs)
  ) %>%
  ungroup() 
## Looking at dat 5 temperoarily
dat_5 %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  ggplot(aes(x=StrainMix, y=WN_FC_fluor)) +
  geom_boxplot() +
  geom_point(aes(col=Alive)) +
  scale_color_manual(values=c("orange","green"))
dat_5 %>%
  filter(ratio=="1-1") %>%
  mutate(Alive = as.numeric(Healthiness_hsv>400)) %>%
  filter(StrainMix == "WCS365-N2C3") %>%
  ggplot(aes(x=WN_FC_fluor, y=Alive)) +
  geom_point() +
  geom_smooth(method="glm", method.args=list(family="binomial"))
dat_5
# Wow! Looks great. Let's keep exp 3, 4, and 5

# ## This is old stuff from comparing standard curves
# ## Check if this makes sense
# gg_comparing_estAbsMethods_NN <- dat_12 %>%
#   select(NN_FC_CFU, NN_FC_fluor_offset, fluor_protect) %>%
#   rename(NN_FC_estAbs = NN_FC_fluor_offset) %>%
#   mutate(stcurve = "Constant Offset Method") %>% distinct() %>%
#   full_join(dat_12 %>% select(NN_FC_CFU, NN_FC_estAbs, stcurve, fluor_protect)) %>%
#   ggplot() +
#   geom_point(aes(x=NN_FC_CFU, y=NN_FC_estAbs, col=fluor_protect)) +
#   geom_abline(aes(slope=1, intercept=0), lty=2, col="grey") +
#   facet_wrap(.~stcurve, labeller=labeller(stcurve=c(exponential_200="Standard curve made with\nExponential culture, 200uL\nIntercept forced through 0"
#                                                     ,exponential_200_winter="Standard curve made with\nExponential culture, 200uL\nFree intercept"
#                                                     ,exponential_275="Standard curve made with\nExponential culture, 275uL\nIntercept forced through 0"
#                                                     ,exponential_275_winter="Standard curve made with\nExponential culture, 275uL\nFree intercept"
#                                                     ,overnight_200="Standard curve made with\nOvernight culture, 200uL\nIntercept forced through 0"
#                                                     ,overnight_200_winter="Standard curve made with\nOvernight culture, 200uL\nFree intercept"
#                                                     ,overnight_275="Standard curve made with\nOvernight culture, 275uL\nIntercept forced through 0"
#                                                     ,overnight_275_winter="Standard curve made with\nOvernight culture, 275uL\nFree intercept"
#                                                     ,`Constant Offset Method`="Constant Offset Method"
#   ))) +
#   xlab("Log2 fold-change of \nknown CFU count ratios (N1:N2)") +
#   ylab("Log2 fold-change of \nfluoresence-based ratio estimates (N1:N2)") +
#   labs(col="Fluoresence of N1") +
#   scale_color_manual(values=c("orange","green"))
# gg_comparing_estAbsMethods_NN
# ggsave("02_compile_data/gg_comparing_estAbsMethods_NN.png", gg_comparing_estAbsMethods_NN, height=9, width=10)
# 
# 
# gg_comparing_estAbsMethods_WN <- dat_34 %>%
#   select(WN_FC_CFU, WN_FC_fluor_offset, fluor_protect) %>%
#   rename(WN_FC_estAbs = WN_FC_fluor_offset) %>%
#   mutate(stcurve = "Constant Offset Method") %>% distinct() %>%
#   full_join(dat_34 %>% select(WN_FC_CFU, WN_FC_estAbs, stcurve, fluor_protect)) %>%
#   ggplot() +
#   geom_point(aes(x=WN_FC_CFU, y=WN_FC_estAbs, col=fluor_protect)) +
#   geom_abline(aes(slope=1, intercept=0), lty=2, col="grey") +
#   facet_wrap(.~stcurve, labeller=labeller(stcurve=c(exponential_200="Standard curve made with\nExponential culture, 200uL\nIntercept forced through 0"
#                                                     ,exponential_200_winter="Standard curve made with\nExponential culture, 200uL\nFree intercept"
#                                                     ,exponential_275="Standard curve made with\nExponential culture, 275uL\nIntercept forced through 0"
#                                                     ,exponential_275_winter="Standard curve made with\nExponential culture, 275uL\nFree intercept"
#                                                     ,overnight_200="Standard curve made with\nOvernight culture, 200uL\nIntercept forced through 0"
#                                                     ,overnight_200_winter="Standard curve made with\nOvernight culture, 200uL\nFree intercept"
#                                                     ,overnight_275="Standard curve made with\nOvernight culture, 275uL\nIntercept forced through 0"
#                                                     ,overnight_275_winter="Standard curve made with\nOvernight culture, 275uL\nFree intercept"
#                                                     ,`Constant Offset Method` ="Constant Offset Method"))) +
#   xlab("Log2 fold-change of \nknown CFU count ratios (WCS365:N2C3)") +
#   ylab("Log2 fold-change of \nfluoresence-based ratio estimates (WCS365:N2C3)") +
#   labs(col="Fluoresence of\nWCS365") +
#   scale_color_manual(values=c("orange","green"))
# gg_comparing_estAbsMethods_WN
# ggsave("02_compile_data/gg_comparing_estAbsMethods_WN.png", gg_comparing_estAbsMethods_WN, height=9, width=10)
# 

# Tried looking at what those 4 are: 3 are from a single plate in exp 2023-09-14 in column 6 (A,B,F). The other is from a
# different experiment. I think there might have been something wonky with those samples.
# 
# dat_wn_fluor %>%
#   filter(as.numeric(col)>3, as.numeric(col)<9) %>%
#   ggplot() +
#   geom_point(aes(x=WN_FC_CFU, y=WN_FC_fluor_offset, col=factor(col)))
# dat_wn_fluor %>%
#   ggplot() +
#   geom_point(aes(x=WN_FC_CFU, y=WN_FC_fluor_offset, col=factor(row)))
# dat_wn_fluor %>%
#   filter(experiment=="2023-09-14_fluor", plate=="Nneon_WcrimLacZ_A_003_crop") %>%
#   ggplot() +
#   geom_point(aes(x=WN_FC_CFU, y=WN_FC_fluor_offset, col=factor(col)))
# 
# dat_34 %>%
#   filter(is.finite(WN_FC_fluor_offset), is.finite(WN_FC_CFU)) %>%
#   ggplot() +
#   geom_point(aes(x=WN_FC_CFU, y=WN_FC_estAbs, pch=fluor_protect, col=stcurve)) 


#### Let's take a peak at final data
dat_wn_fluor %>%
  filter(is.finite(WN_FC_fluor_offset)) %>%
  select(ratio,WN_FC_fluor_offset,Healthiness_hsv, StrainMix  ) %>% distinct() %>%
  mutate(Alive = Healthiness_hsv>400) %>%
  ggplot(aes(x=ratio, y=WN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=Alive), width=0.1, height=0) +
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x", space = "free")+
  scale_color_manual(values=c("orange","darkgreen")) +
  geom_hline(aes(yintercept=0)) +
  xlab("Inoculation Ratio (WCS:N2C3)") +
  ylab("LogFold change (WCS:N2C3)")

## Okay, so offset average doesn't seem to work too well for different strains...?
dat_12 %>%
  filter(is.finite(NN_FC_fluor_offset)) %>%
  select(ratio,NN_FC_fluor_offset,Healthiness_hsv, StrainMix,fluor_protect  ) %>% distinct() %>%
  ggplot(aes(x=ratio, y=NN_FC_fluor_offset)) +
  geom_boxplot() +
  geom_jitter(aes(col=fluor_protect), width=0.1, height=0, alpha=0.5) +
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x", space = "free")+
  xlab("Inoculation Ratio\n(N2C3 (Strain 1):N2C3 (Strain2))") +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values=c("orange","darkgreen")) +
  labs(col="Fluor of N2C3\nStrain1") +
  ylab("LogFold change (N1 : N2)")

dat_n2c3_vs_n2c3_col0 <- dat_12 %>%
  filter(plant=="col0") 
# %>%
  # Remove all the testing curves we did-- this was already removed/ not included above so commenting out here
  # select(-NN_FC_estAbs, -starts_with("od_"), -ratio_cn_estAbs, -stcurve,
         # -CN_FC_estAbs, -crim_estAbs, -neon_estAbs, -fc_cn) %>% distinct()

dat_wcs_vs_n2c3_col0 <- allDat_cleaned %>%
  filter( experiment %in% c("2023-09-13_fluor","2023-09-14_fluor","2023-11-28_fluor")) %>%
  filter(plant=="col0") %>%
  mutate(WN_FC_fluor = ifelse(fluor_protect == "crim",CN_FC_fluor, -CN_FC_fluor )
         ,WN_FC_fluor_offset = ifelse(fluor_protect == "crim",CN_FC_fluor_offset, -CN_FC_fluor_offset )
         ,WN_FC_CFU= ifelse(fluor_protect == "crim",CN_FC_CFU, -CN_FC_CFU )
         # ,WN_FC_estAbs = ifelse(fluor_protect == "crim", CN_FC_estAbs, -CN_FC_estAbs)
  ) %>%
  ungroup() 
# %>%
  # Remove all the testing curves we did
  # select(-WN_FC_estAbs, -starts_with("od_"), -ratio_cn_estAbs, -stcurve,
         # -CN_FC_estAbs, -crim_estAbs, -neon_estAbs, -fc_cn) %>% distinct()


write.table(dat_n2c3_vs_n2c3_col0, file="02_compile_data/dat_n2c3_vs_n2c3.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(dat_wcs_vs_n2c3_col0, file="02_compile_data/dat_wcs_vs_n2c3.txt", row.names = FALSE, quote = FALSE, sep="\t")
