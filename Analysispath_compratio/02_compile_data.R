library(plotwidgets)
library(lubridate)
library(ggbiplot)
library(tidyverse)


dir.create("02_compile_data")

#### GET METADATA ####
# alldirs <- list.dirs("00_Experiments",recursive = TRUE)[grep("/meta", list.dirs("00_Experiments",recursive = TRUE), fixed=TRUE)]
# allMetaFP <- list.files(path=alldirs, full.names = TRUE)[grep('_meta.csv',list.files(path=alldirs, full.names = TRUE), fixed=TRUE)]

# Load and join
# Get metadata file names
meta_filenames <- list.files("00_raw_data/meta")
# Iterate through vector of filenames and load them in, joining with each other.
allMeta <- data.frame(plate=NA) # Empty dataframe must have at least one common column name as files you're going to load in
for ( m in meta_filenames) {
  tempmeta <- read.csv(paste0("00_raw_data/meta/",m)) %>% mutate(notes = as.character(notes))
  allMeta <- full_join(allMeta, tempmeta)
}
# remove all lines where plate is NA
allMeta  <- filter(allMeta, !is.na(plate))

#### Conversion for OD/abs to cells #########
# OD = 3.166*abs - 0.268 (but only When Abs>0.1)
# Cells per ml for OD1 is 5e8
# Cells for Abs 0.1 = 3.166*0.1-0.268*5E8
# cells_per_Abs1 = (3.166*1-0.268)*5E8 ## OLD ONE, including all Abs600
cells_per_Abs1 = (3.19*1-0.275)*5E8 ## NEW ONE, excluding Abs600<0.125


allMeta_adj <- allMeta %>%
  # filter(!(experiment%in% c("2023-07-18_rhizoscreen", "2023-07-26_rhizoscreen"))) %>%
  mutate(date_germ = as.Date(date_germ)
         , date_ster = as.Date(date_ster)
         , date_inoc = as.Date(date_inoc)
         , date_process = as.Date(date_process)) %>%
  separate(path, into= c("path","other"), sep="-") %>%
  separate(other, into= c("lacZ_path","fluor_path"), sep="_") %>% 
  mutate(lacZ_path=ifelse(lacZ_path=="",NA,lacZ_path),fluor_path=ifelse(fluor_path=="",NA,fluor_path)) %>%
  separate(protect, into= c("protect","other"), sep="-") %>%
  separate(other, into= c("lacZ_protect","fluor_protect"), sep="_") %>% 
  mutate(lacZ_protect=ifelse(lacZ_protect=="",NA,lacZ_protect),fluor_protect=ifelse(fluor_protect=="",NA,fluor_protect))%>%
  mutate(protect = str_to_upper(protect), protect_v2 = str_to_upper(protect_v2)) %>%
  mutate(col=factor(col, levels=c(seq(1,12)))
         , protect = factor(protect, levels=c("MOCK","WCS365","CHAO","CH267","PF5","PO6"))
         , path = factor(path, levels=c("MOCK","N2C3","PAO1"))) %>%
  filter(!experiment%in% c("48pilot_fluor","")) %>%
  unite(protect, path, col="StrainMix", sep = "-", remove=FALSE) %>%
  rowwise() %>%
  mutate(path_od_adj = ifelse(path_od == 9e-3, 1e-2,
                              ifelse(path_od == 9e-5, 1e-4
                                     , ifelse(path_od == 9e-7, 1e-6,path_od)))) %>%
  ungroup() %>%
  mutate(path_od_actual = path_od, path_od = path_od_adj, protect_od = protect_od) %>%
  rowwise() %>%
  mutate(totalOD = protect_od + path_od
         # , totalOD_adj = protect_od + path_od_adj
         , ratio_protect = protect_od/min(c(protect_od, path_od))
         # , ratio_protect_adj = protect_od/min(c(protect_od, path_od_adj))
         , ratio_path = path_od/min(c(protect_od, path_od))
         # , ratio_path_adj = path_od_adj/min(c(protect_od, path_od_adj))
  ) %>% 
  mutate(ratio_protect = ifelse(is.infinite(ratio_protect), 1, 
                            ifelse(is.na(ratio_protect), 0, ratio_protect))
         , ratio_path = ifelse(is.infinite(ratio_path), 1, 
                                  ifelse(is.na(ratio_path), 0, ratio_path))) %>%
  mutate(ratio = paste0(c(ratio_protect, ratio_path), collapse = "-")) %>%
  # mutate(ratio_adj = paste0(c(ratio_protect, ratio_path_adj), collapse = "-")) %>%
  ungroup() %>%
  mutate(protect_cells = protect_od*cells_per_Abs1*0.275,
         path_cells = path_od*cells_per_Abs1*0.275, 
         total_cells = totalOD*cells_per_Abs1*0.275,
         protect_cells_log = log10(protect_cells+1),
         path_cells_log = log10(path_cells+1),
         total_cells_log = log10(total_cells+1)) %>%
  mutate(plant_age = date_inoc - date_germ) %>%
  select(where(function(x) any(!is.na(x))))
# allMeta_adj %>%
#   filter(fluor_protect=="neon") %>% View()
# allMeta_rhizo <- allMeta %>%
#   filter((experiment%in% c("2023-07-18_rhizoscreen", "2023-07-26_rhizoscreen"))) %>%
#   mutate(date_germ = as.Date(date_germ)
#          , date_ster = as.Date(date_ster)
#          , date_inoc = as.Date(date_inoc)
#          , date_process = as.Date(date_process)) %>%
#   mutate(col=factor(col, levels=c(seq(1,12)))
#          , protect = factor(protect, levels=unique(c("MOCK",unique(protect))))) %>%
#   filter(!experiment%in% c("48pilot_fluor","")) %>%
#   mutate(protect_od = protect_od) 
# NOTE: There will still be NAs in blank wells but this will be filtered out when joined with pixel data
# View(allMeta_adj)

# allMeta_adj %>% View()

#### GET PIXEL DATA ####
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
    mutate(experiment = experiment, plate=plate)  # mutate creates new columns in your data frame. Adding experiment and plate to each
  # old code had an error in it where pix_b_median was accidentally named pix_b_ so I'm fixing this here
  if ("pix_b_" %in% colnames(tempdat)) {
    tempdat <- tempdat %>% rename("pix_b_median" = "pix_b_")
  }
  
  allPixel <- full_join(allPixel, tempdat)
}
# remove subimage == NA
dat_adj <- allPixel %>%
  filter(!is.na(subimage)) %>%
  separate(subimage, into = c("remove","well"), remove=TRUE, sep="_") %>% select(-remove) %>%
  separate(well, into=c("row","col"), remove=FALSE, sep=1) %>% 
  # separate(experiment__plate, into=c("experiment","plate"), sep="__", remove=TRUE) %>%
  mutate(col=factor(col, levels=c(seq(1,12)))) 


# 
# dat <- read.delim("01_Processed_Scans/combined_pixel_data.txt")
# dat_adj <- dat %>% 
#   separate(subimage, into = c("remove","well"), remove=TRUE, sep="_") %>% select(-remove) %>%
#   separate(well, into=c("row","col"), remove=FALSE, sep=1) %>% 
#   separate(experiment__plate, into=c("experiment","plate"), sep="__", remove=TRUE) %>%
#   mutate(col=factor(col, levels=c(seq(1,12)))) 


##### JOIN and check #####
allDat <- full_join(allMeta_adj,dat_adj) 

### Sanity check ; make sure overlap is adequate ###
# allDat %>% View()
allDat %>% filter((is.na(StrainMix)& all_plant_pixels>50)) 
allDat %>% filter(is.na(all_plant_pixels)) 
# Double checking it looks okay
allDat %>%
  select(StrainMix, protect_od) %>%
  table()

allDat %>%
  select(StrainMix, path_od) %>%
  table()

#### Histograms of pixel values ####
## plant pixels
ggplot(allDat) +
  geom_histogram(aes(x=all_plant_pixels))

allDat %>%
  ggplot() +
  geom_density(aes(x=pix_g_median), fill=NA, col="green")+
  geom_density(aes(x=pix_r_median), fill=NA, col="red")+
  geom_density(aes(x=pix_b_median), fill=NA, col="blue")

allDat %>%
  filter(pix_g_median<50)
# Filter out all cases where green pixels is less than 100. means there is no plant.
allDat %>%
  filter(pix_r_median<50) 
# Cases where red pixels are less than 50 are the same problem.
allDat %>%
  filter(pix_b_median < 60, pix_b_median>40)
# Keep the blue ones with pix<50; those are real
allDat %>%
  filter( pix_b_median>100) 
# Also checked blue > 100; those are real too.

## Check out those that have no protective or pathogen
allDat %>%
  filter(is.na(protect)| is.na(path))
# They all have >50 pixels

allDat <- allDat %>%
  filter(!is.na(protect), !is.na(path))%>%
  filter(pix_g_median>100, pix_r_median>50)

allDat %>%
  ggplot() +
  geom_density(aes(x=pix_g_median), fill=NA, col="green")+
  geom_density(aes(x=pix_r_median), fill=NA, col="red")+
  geom_density(aes(x=pix_b_median), fill=NA, col="blue")

allDat %>%
  ggplot() +
  geom_density(aes(x=diff_gb_median), fill=NA, col="green")+
  geom_density(aes(x=diff_rb_median), fill=NA, col="red")+
  geom_density(aes(x=diff_gr_median), fill=NA, col="brown")

###### Save data #########
save(allDat, file="02_compile_data/allDat.RData")
write.table(allDat, file="02_compile_data/allDat.txt", sep="\t", row.names = FALSE, quote=FALSE)
