
library(MASS)
library(plotwidgets) # for hsl
library(lubridate) # dates
library(lme4) # lmm models
library(lmerTest)
library(ggbiplot) # pca plotting
library(tidyverse)
library(cowplot) # panels


#### Load #######
dir.create("03_EDA_and_adj")
load("02_compile_data/allDat.RData")

########## CUSTOM BIPLOT ##############
#### Make custom biplot
ggbiplot_PLANTS <- function (pcobj, colors="black",outline_col= "black", text_col="black", arrow_col="yellow", choices = 1:2, scale = 1, pc.biplot = TRUE, 
                             obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                             ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                             alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                             varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                             ...) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  # MAKE VARIABLE NAME
  df.u <- df.u %>% as.data.frame() %>%
    rownames_to_column(var="UniqueIDs")
  # MAKE COLOURS
  if (length(colors) != nrow(df.u)) {
    colors <- rep("black", nrow(df.u))
    names(colors) <- df.u$UniqueIDs
  }
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal() 
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(color=groups, label = labels), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      # g <- g + geom_point(aes(color = groups), alpha = alpha, show.legend = FALSE)
      g <- g + geom_point(aes(fill=UniqueIDs, col=groups), alpha = alpha, show.legend = FALSE, pch=21, size=2) +# Don't colour by groups
        scale_fill_manual(values=colors)
    }
    else {
      g <- g + geom_point(aes(fill=UniqueIDs), col=outline_col, alpha = alpha, show.legend = FALSE, pch=21, size=2)+
        scale_fill_manual(values=colors)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circ <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circ %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    # g <- g + geom_path(data = ell, aes(color = groups, group = groups)) # eventually change so colours are unique here
    g <- g + geom_path(data = ell, aes(color=groups, group = groups))
  }
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = muted(arrow_col))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = text_col, size = varname.size)
  }
  return(g)
}


# #### Look at mock ####
# mock_col_5 <- allDat %>%
#   filter(protect=="MOCK" & path == "MOCK"
#          , plant=="col0"
#          , plant_age==5
#          , MS==0.5
#          , MES == 0.5
#          , pH ==5.8) 
#   
# ### PCA of colours ####
# UniqueIDs_mock <- mock_col_5 %>%
#   filter(!is.na(pix_b_median), !is.na(pix_r_median), !is.na(pix_g_median), pix_g_median>100, pix_r_median>50) %>%
#   mutate(ID="ID") %>%
#   unite(ID,experiment, plate, row, col, col="UniqueID", sep="__") %>%
#   select(UniqueID) %>% pull()
# 
# rgb_mat_mock <- mock_col_5 %>%
#   filter(!is.na(pix_b_median), !is.na(pix_r_median), !is.na(pix_g_median), pix_g_median>100, pix_r_median>50) %>%
#   mutate(logpix = log(all_plant_pixels)) %>%
#   select(logpix, pix_r_median, pix_g_median,pix_b_median) %>%
#   as.matrix() 
# rownames(rgb_mat_mock) <- UniqueIDs_mock
# 
# ratio_mat_mock <-  mock_col_5 %>%
#   filter(!is.na(pix_b_median), !is.na(pix_r_median), !is.na(pix_g_median), pix_g_median>100, pix_r_median>50) %>%
#   mutate(logpix = log(all_plant_pixels)) %>%
#   select(logpix, ratio_r_median, ratio_g_median, ratio_b_median) %>%
#   as.matrix() 
# rownames(ratio_mat_mock) <- UniqueIDs_mock
# 
# ### Get colours for each unique ID
# 
# # Get colors
# colours <- apply(rgb_mat_mock, 1, FUN = function(x) rgb(x[2]/255,x[3]/255,x[4]/255))
# names(colours) <- rownames(rgb_mat_mock)
# 
# pca_rgb <- prcomp(rgb_mat_mock)
# ggbiplot_PLANTS(pca_rgb, colors = colours)
# 
# ggbiplot(pca_rgb, labels = rownames(pca_rgb$x))+
#   xlim(-10,2)
# # Check out three outliers:
# "2023-05-30_ccppwn_CH267_CHAO_001_b__D__5"
# "2023-05-10_age_concentration__N2C3_IH001_cropped__A__6"
# "2023-05-03_age_concentration__age_concentration_n2c3001_cropped__E__12"
# # They are all normal plants, but they are split between yellow and green.
# 
# 
# 
# ### Now, try ratios
# pca_ratio <- prcomp(ratio_mat_mock)
# ggbiplot_PLANTS(pca_ratio, colors = colours)

### Try with ALL plants??

allDat_filt <- allDat %>%
  mutate(ID="ID") %>%
  unite(ID,experiment, plate, row, col, col="UniqueID", sep="__", remove=FALSE) %>%
  mutate(UniqueID = gsub("-","_", UniqueID)) %>%
  # filter(!is.na(pix_b_median), !is.na(pix_r_median), !is.na(pix_g_median), pix_g_median>100, pix_r_median>50) %>%
  # mutate(totalOD_adj = path_od_adj + protect_od) %>%
  # mutate(protect_cells = (protect_od)*5e8*0.275, path_cells = (path_od)*5e8*0.275, total_cells = (totalOD)*5e8*0.275) %>%
  # mutate(path_cells_adj = (path_od_adj/0.3)*2E9, total_cells_adj = (totalOD_adj/0.3)*2E9) %>%
  # mutate(protect_cells_log = log10(protect_cells+1)
         # , path_cells_log = log10(path_cells+1)
         # , path_cells_adj_log = log10(path_cells_adj+1)
         # ) %>%
  mutate(plant = factor(plant, levels=c("col0","bbc","bik1","nej"))
         ) %>%
  filter(protect !="PO6", path!="PAO1") 

UniqueIDs_all <- allDat_filt %>%
  select(UniqueID) %>% pull()
  

######### Different colour theories ###########
## Get a subset of experiments to include for the PCA
allDat_filt %>% select(experiment) %>% unique()
## Choose ONLY those that are not being used in actual experiments, so we have independence
toKeepExp <- c("2023-05-03_age_concentration", "2023-05-10_age_concentration"
               ,"2023-05-16_basiceffects_n2c3"
               # , "2023-05-24_age_concentration" ### Used this for OD stuff
  , "2023-05-30_ccppwn"
  ,"2023-06-07_ccpwn"
  # ,"2023-08-08_fluor","2023-08-16_fluor"
  # , "2023-09-06_ccpwncombos"
  )
allDat_forpca <- allDat_filt %>%
  filter(experiment %in% toKeepExp)
UniqueIDs_forpca <- allDat_forpca$UniqueID
length(UniqueIDs_forpca)
allDat_forpca %>%
  select(StrainMix, plant_age) %>% table()

allDat_filt[which(is.na(allDat_filt$pix_b_median)),] %>% View()

allRGB <- allDat_filt %>%
  mutate(logpix = log10(all_plant_pixels)) %>%
  select(logpix, pix_r_median, pix_g_median, pix_b_median) %>%
  rename(Number_pixels_log=logpix, red_intensity = pix_r_median, green_intensity = pix_g_median, blue_intensity = pix_b_median) %>%
  as.matrix()
rownames(allRGB) <- UniqueIDs_all
allHSV <- cbind(allRGB[,1], t(apply(allRGB, 1,function(x) rgb2hsv(r=x[2],g=x[3],b=x[4]))))
colnames(allHSV) <- c("Number_pixels_log","hue","saturation","value")
allHSL <- cbind(allRGB[,1], t(apply(allRGB, 1,function(x) rgb2hsl(as.matrix(x[2:4])))))
colnames(allHSL) <- c("Number_pixels_log","hue","saturation","lightness")

# Get all colours
# Get colors
col_all <- apply(allRGB, 1, FUN = function(x) rgb(x[2]/255,x[3]/255,x[4]/255))
names(col_all) <- rownames(allRGB)

###### make PCA #########
## Get tabel of all treatments for PCA
sink("03_EDA_and_adj/healthiness_pca_info.txt")
allDat_forpca %>%
  select(path, protect) %>% table()
allDat_forpca %>%
  select(experiment) %>% unique()
allDat_forpca %>%
  select(protect, protect_od) %>% table()
allDat_forpca %>%
  select(path, path_od) %>% table()
sink()
# Make some platemaps
plateNumMap <- allDat_forpca %>%
  select(experiment, plate) %>% unique() %>%
  group_by(experiment) %>%
  mutate(platenum = rank(plate)) %>%
  ungroup()
gg_platemap_for_pca <- allDat_forpca %>%
  mutate(path_cells_log = ifelse(path=="MOCK",6, path_cells_log), protect_cells_log = ifelse(protect=="MOCK", 6, protect_cells_log)) %>%
  select(UniqueID, experiment, plate, row, col, StrainMix, protect, path, protect_od, path_od,
         protect_cells_log, path_cells_log, total_cells_log, plant, plant_age) %>%
  mutate(row = factor(row, levels=rev(LETTERS[1:8]))) %>%
  left_join(plateNumMap) %>%
  ggplot() +
  # geom_tile(aes(x=col, y=row), fill="white") +
  geom_tile(aes(x=col, y=row, fill=protect, alpha=protect_cells_log)) +
  # geom_point(aes(x=col, y=row), col="grey") +
  geom_point(aes(x=col, y=row, col=path, alpha = path_cells_log, pch=factor(plant_age))) +
  facet_grid(platenum~experiment) +
  scale_color_manual(values=c(MOCK="grey",N2C3="darkred")) +
  scale_fill_manual(values=c(MOCK="white", WCS365= "darkgreen", CHAO="goldenrod",CH267 ="blue", PF5="purple")) +
  scale_shape_manual(values=c(19, 17, 15 )) +
  labs(col="Pathogen added", fill="Commensal added", alpha="Approx cells (log)", shape="Plant age\nat inoculation") +
  xlab("Columns in 96-well plate") + ylab("Rows in 96-well plate")
gg_platemap_for_pca
ggsave("03_EDA_and_adj/gg_platemap_for_pca.png", gg_platemap_for_pca, height=10, width=15)

pca_rgb <- prcomp(allRGB[allDat_forpca$UniqueID,])
pca_hsv <- prcomp(allHSV[allDat_forpca$UniqueID,])
pca_hsl <- prcomp(allHSL[allDat_forpca$UniqueID,])

col_all


# PCA
# pca_all <- prcomp(alldat_mat)
gg_biplot_rgb <- ggbiplot_PLANTS(pca_rgb, colors = col_all[allDat_forpca$UniqueID], outline_col=col_all[allDat_forpca$UniqueID],arrow_col="black", text_col="black",choices = c(1,2))
gg_biplot_rgb
ggsave(filename = "03_EDA_and_adj/gg_biplot_rgb.png", gg_biplot_rgb
       , width=5, height=4)

pca_rgb$x %>% as.data.frame() %>%
  arrange(-PC2) %>% head()

gg_biplot_hsv <- ggbiplot_PLANTS(pca_hsv, colors = col_all[allDat_forpca$UniqueID], outline_col=col_all[allDat_forpca$UniqueID],arrow_col="black", text_col="black",choices = c(1,2)) +
  xlim(-5.5,5)
gg_biplot_hsv
ggsave(filename = "03_EDA_and_adj/gg_biplot_hsv.png", gg_biplot_hsv
       , width=5, height=4)

gg_biplot_hsl <- ggbiplot_PLANTS(pca_hsl, colors = col_all[allDat_forpca$UniqueID], outline_col=col_all[allDat_forpca$UniqueID],arrow_col="black", text_col="black",choices = c(1,2))
gg_biplot_hsl
ggsave(filename = "03_EDA_and_adj/gg_biplot_hsl.png", gg_biplot_hsl
       , width=5, height=4)

# Calculate x contrib to each axis; can we omit number of pixels?
load_coord <- t(apply(pca_rgb$rotation, 1, FUN = function(x)x*pca_rgb$sdev))^2
x_contrib <- t(apply(load_coord, 1, FUN=function(x) x/colSums(load_coord)))
# Multiply out by the variance explained by each axis
var_expl_axis <- pca_rgb$sdev^2/sum(pca_rgb$sdev^2)
total_var_expl <- apply(x_contrib, 1, function(x) x*var_expl_axis)
colSums(total_var_expl)
## 0.000087 is number pixels; or 0.0087% of variation explained by pixel count
# The first 2 PCs explain this much variation in dataset:
sum(var_expl_axis[c(1,2)])

# ###### Try without size to compare amount of variance explained by size

pca_rgb_nosize <- prcomp(allRGB[allDat_forpca$UniqueID,-1])
pca_hsv_nosize <- prcomp(allHSV[allDat_forpca$UniqueID,-1])
pca_hsl_nosize <- prcomp(allHSL[allDat_forpca$UniqueID,-1])


# PCA
# pca_all <- prcomp(alldat_mat)
gg_biplot_rgb_nosize <- ggbiplot_PLANTS(pca_rgb_nosize, colors = col_all[allDat_forpca$UniqueID], outline_col=col_all[allDat_forpca$UniqueID],arrow_col="black", text_col="black",choices = c(1,2))
gg_biplot_rgb_nosize
ggsave(filename = "03_EDA_and_adj/gg_biplot_rgb_nosize.png", gg_biplot_rgb_nosize
       , width=5, height=4)

gg_biplot_hsv_nosize <- ggbiplot_PLANTS(pca_hsv_nosize, colors = col_all[allDat_forpca$UniqueID], outline_col=col_all[allDat_forpca$UniqueID],arrow_col="black", text_col="black",choices = c(1,2))
gg_biplot_hsv_nosize
ggsave(filename = "03_EDA_and_adj/gg_biplot_hsv_nosize.png", gg_biplot_hsv_nosize
       , width=5, height=4)

gg_biplot_hsl_nosize <- ggbiplot_PLANTS(pca_hsl_nosize, colors = col_all[allDat_forpca$UniqueID], outline_col=col_all[allDat_forpca$UniqueID],arrow_col="black", text_col="black",choices = c(1,2))
gg_biplot_hsl_nosize
ggsave(filename = "03_EDA_and_adj/gg_biplot_hsl_nosize.png", gg_biplot_hsl_nosize
       , width=5, height=4)


# Note; I tried an NMDS too for rgb, but the shape looks exactly like a PCA and it's so much slowe rto run
###
# 
# 
# ### Can a separate this out into MOCK, WCS, and N2C3 treatments only?
# # unique(allDat_filt$experiment)
# toKeepIDs <- allDat_filt %>%
#   # filter(experiment == "2023-06-14_age_concentration_Pf5_PROPER") %>%
#   # select(protect_od, path_od) %>% table()
#   filter((protect=="WCS365" & path=="MOCK") | (protect == "MOCK" & path=="N2C3") | (protect=="MOCK" & path=="MOCK")) %>%
#   filter(protect_od %in% c(0,0.01), path_od %in% c(0,0.01)) %>%
#   filter(MS==0.5, MES==0.5, pH==5.8, plant=="col0", plant_age==5) %>%
#   pull(UniqueID)
# col_all_filt <- col_all[toKeepIDs]
# 
# pca_singlestrain <- prcomp(alldat_mat[toKeepIDs,])
# groups_treat <- allDat_filt %>% filter(UniqueID %in% toKeepIDs) %>% pull(StrainMix)
# # gg_biplot_singleStrains <- ggbiplot_PLANTS(pca_singlestrain, colors = col_all_filt, choices = c(1,2), groups = groups_treat, ellipse = TRUE)
# # gg_biplot_singleStrains +
# #   scale_color_manual(values=c("MOCK-MOCK"="black", "MOCK-N2C3"="darkred", "WCS365-MOCK"="blue"))
# 
# gg_pca_wcm <- ggbiplot(pca_singlestrain, choices = c(1,2), groups = groups_treat, ellipse = TRUE)+
#   scale_color_manual(values=c("MOCK-MOCK"="black", "MOCK-N2C3"="darkred", "WCS365-MOCK"="blue"))+
#   scale_x_reverse()+
#   scale_y_reverse()
# gg_pca_wcm
# ggsave("03_EDA_and_adj/gg_pca_wcm.png",gg_pca_wcm, width=5, height=4)

######## Loadings for healthy/stressed, dead/alive #########
pca_rgb_x <- data.frame(pca_rgb$x) %>%
  mutate(PC1_adj = (PC1-min(PC1))/max(PC1), PC2_adj = ((-PC2)-min(-PC2))/max(-PC2)) %>%
  mutate(combPC = PC1_adj * PC2_adj) %>%
  rownames_to_column(var="UniqueID")

pca_hsv_x_nosize <- data.frame(pca_hsv_nosize$x) %>%
  mutate(PC1_adj = (PC1-min(PC1))/max(PC1), PC2_adj = ((-PC2)-min(-PC2))/max(-PC2)) %>%
  mutate(combPC = PC1_adj * PC2_adj) %>%
  rownames_to_column(var="UniqueID")
round(var_expl_axis[c(1)],3)*100
              
gg_allplants_colsize <- pca_rgb_x %>%
  left_join(allDat_filt) %>%
ggplot() +
  geom_point(aes(x=PC1, y=-PC2, fill=UniqueID, size = all_plant_pixels), pch=21,  show.legend = FALSE)+
  scale_fill_manual(values=col_all) +
  ylab(paste0("PC2\n",round(var_expl_axis[c(2)],3)*100,"% of variation"))+
  xlab(paste0("PC1\n",round(var_expl_axis[c(1)],3)*100,"% of variation"))+
  scale_y_reverse()
gg_allplants_colsize
ggsave(filename = "03_EDA_and_adj/gg_allplants_colsize.png",gg_allplants_colsize
       ,height=4, width=5)
# Try colouring by ratios?

gg_allplants_hsv_nosize <- pca_hsv_x_nosize %>%
  left_join(allDat_filt) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=-PC2, fill=UniqueID, size = all_plant_pixels), pch=21,  show.legend = FALSE)+
  scale_fill_manual(values=col_all) +
  ylab(paste0("PC2\n",round(var_expl_axis[c(2)],3)*100,"% of variation"))+
  xlab(paste0("PC1\n",round(var_expl_axis[c(1)],3)*100,"% of variation"))+
  scale_y_reverse()
gg_allplants_hsv_nosize
ggsave(filename = "03_EDA_and_adj/gg_allplants_hsv_nosize.png",gg_allplants_hsv_nosize
       ,height=4, width=5)

gg_diff_gb <- pca_rgb_x %>%
  left_join(allDat_filt) %>%
  ggplot() +
  geom_point(aes(x=PC1_adj, y=PC2_adj, col=diff_gb_median), show.legend = TRUE)+
  scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Difference in\nGreen/Blue intensity")
gg_diff_gb
ggsave("03_EDA_and_adj/gg_diff_gb.png", gg_diff_gb, width=5, height=4)

gg_diff_gr <- pca_rgb_x %>%
  left_join(allDat_filt) %>%
  ggplot() +
  geom_point(aes(x=PC1_adj, y=PC2_adj, col=diff_gr_median), show.legend = TRUE)+
  scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Difference in\nGreen/Red intensity")
gg_diff_gr
ggsave("03_EDA_and_adj/gg_diff_gr.png", gg_diff_gr, width=5, height=4)

gg_pixels <- pca_rgb_x %>%
  left_join(allDat_filt) %>%
  ggplot() +
  geom_point(aes(x=PC1_adj, y=PC2_adj, col=log10(all_plant_pixels)), show.legend = TRUE)+
  scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Plant size (log pixels)")
gg_pixels
ggsave("03_EDA_and_adj/gg_pixels.png", gg_pixels, width=5, height=4)

gg_comb <- pca_rgb_x %>%
  left_join(allDat_filt) %>%
  mutate(comb = (diff_gb_median + diff_gr_median*2)*log10(all_plant_pixels) )%>%
  ggplot() +
  geom_point(aes(x=PC1_adj, y=PC2_adj, col=comb), show.legend = TRUE)+
  scale_y_reverse()
gg_comb
ggsave("03_EDA_and_adj/gg_comb.png", gg_comb, width=5, height=4)

###### Try hsv

gg_hsv_grad_h <- pca_hsv_x_nosize %>%
  left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID")) %>%
  ggplot() +
  geom_point(aes(x=PC1_adj, y=PC2_adj, col=hue), show.legend = TRUE)+
  scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Hue")
gg_hsv_grad_h
ggsave("03_EDA_and_adj/gg_hsv_grad_h.png", gg_hsv_grad_h, width=5, height=4)


gg_hsv_grad_s <- pca_hsv_x_nosize %>%
  left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID")) %>%
  ggplot() +
  geom_point(aes(x=PC1_adj, y=PC2_adj, col=saturation), show.legend = TRUE)+
  scale_y_reverse() + xlab("PC1") + ylab("PC2") + labs(col="Saturation")
gg_hsv_grad_s
ggsave("03_EDA_and_adj/gg_hsv_grad_s.png", gg_hsv_grad_s, width=5, height=4)


# Can this be used as a single variable?
pca_rgb_x %>%
  left_join(allDat_filt) %>%
  ggplot() +
  geom_jitter(aes(x=diff_gb_median, y=1, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
  scale_color_manual(values=col_all)

pca_rgb_x %>%
  left_join(allDat_filt) %>%
  mutate(comb=(diff_gr_median+diff_gb_median)
         ,comb2 = ((diff_gr_median*2 + diff_gb_median*1)*log10(all_plant_pixels+1)) )%>%
  ggplot() +
  geom_jitter(aes(x=comb2, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
  scale_color_manual(values=col_all)

pca_rgb_x %>%
  left_join(allDat_filt) %>%
  mutate(comb=(diff_gr_median+diff_gb_median)
         ,comb2 = ((diff_gr_median*2 + diff_gb_median*1)*(diff_gr_median/10 + log10(all_plant_pixels+1))) )%>%
  ggplot() +
  geom_jitter(aes(x=comb2, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
  scale_color_manual(values=col_all)


pca_hsv_x_nosize %>%
  left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
  left_join(allDat_filt) %>%
  mutate(comb=(1*hue)*(3*saturation))%>%
  ggplot() +
  geom_jitter(aes(x=comb, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
  scale_color_manual(values=col_all)


pca_hsv_x_nosize %>%
  left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
  left_join(allDat_filt) %>%
  mutate(comb=(hue)*(2*saturation+log10(all_plant_pixels+1)/100))%>%
  ggplot() +
  geom_jitter(aes(x=comb, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
  scale_color_manual(values=col_all)


pca_hsv_x_nosize %>%
  left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
  left_join(allDat_filt) %>%
  mutate(comb=(hue+log10(all_plant_pixels+1)/100)*(2*saturation))%>%
  ggplot() +
  geom_jitter(aes(x=comb, y=1, size=all_plant_pixels, col=UniqueID), height=0.5, width=0, show.legend = FALSE)+
  scale_color_manual(values=col_all)

#### UPDATE allDAT WITH NEW METRIC #####
allDat_filt <- allDat_filt %>%
  left_join(allHSV %>%as.data.frame() %>% rownames_to_column(var="UniqueID"))%>%
  mutate(
    # Healthiness_nosize = diff_gb_median + diff_gr_median*2
         Healthiness_rgb = (diff_gb_median + diff_gr_median*2)*log10(all_plant_pixels+1)
    , Healthiness_hsv = hue*saturation*log10(all_plant_pixels+1)*1000
    # , Healthiness_hsv = (hue)*(2*saturation+log10(all_plant_pixels+1)/100)
    # , protect_cells_log = log10(protect_cells+1)
    # , path_cells_log = log10(path_cells+1)
    # , path_cells_adj_log = log10(path_cells_adj+1)
    , ratio_protect_log = log10(ratio_protect+1)
    , ratio_path_log = log10(ratio_path + 1)
    , ratio_protpath_log = (ratio_protect_log+1)/(ratio_path_log+1)-1
    # , total_cells_log = log10(total_cells+1)
    # , total_cells_adj_log = log10(total_cells_adj + 1)
    )

#### SAVE ####
allDat_final <- allDat_filt
save(allDat_final, file = "03_EDA_and_adj/allDat_final.RData")
write.table(allDat_final, file="03_EDA_and_adj/allDat_final.txt", row.names = FALSE, quote = FALSE, sep="\t")
save(col_all, file="03_EDA_and_adj/plant_colours.RData")

#### N2C3 WCS simplified only ###

 allDat_pilot <- allDat_filt %>%   
   filter(experiment=="2023-05-16_basiceffects_n2c3")
 col_all

 gg_platemap_pilot <- allDat_pilot %>%
   mutate(row = factor(row, levels=rev(LETTERS[1:8]))) %>%
   ggplot() +
   # geom_tile(aes(x=col, y=row), fill="white") +
   geom_tile(aes(x=col, y=row, fill=protect, alpha=protect_cells_log)) +
   # geom_point(aes(x=col, y=row), col="grey") +
   geom_point(aes(x=col, y=row, col=path, alpha = path_cells_log)) +
   facet_wrap(plate~., nrow=3) +
   scale_color_manual(values=c(MOCK="grey",N2C3="darkred")) +
   scale_fill_manual(values=c(MOCK="white", WCS365= "darkgreen", CHAO="goldenrod",CH267 ="blue", PF5="purple")) +
   scale_shape_manual(values=c(19, 17, 15 )) +
   labs(col="Pathogen added", fill="Commensal added", alpha="Approx cells (log)", shape="Plant age\nat inoculation") +
   xlab("Columns in 96-well plate") + ylab("Rows in 96-well plate")
 gg_platemap_pilot
 ggsave(filename = "03_EDA_and_adj/gg_platemap_pilot.png", gg_platemap_pilot, width=8, height=8)
 
gg_pilot_n2c3wcs <- allDat_pilot %>%
  filter(protect_od %in% c(0,0.001), path_od %in% c(0,0.001)) %>%
  mutate(inoc = ifelse(StrainMix=="MOCK-MOCK", "MOCK", 
                       ifelse(StrainMix == "MOCK-N2C3", "N2C3", 
                              ifelse(StrainMix == "WCS365-MOCK", "WCS365", "N2C3 +\nWCS365")))) %>%
  mutate(inoc = factor(inoc, levels=c("MOCK","WCS365","N2C3","N2C3 +\nWCS365"))) %>%
  rowwise() %>%
  mutate(single_maxcell = round(max(c(protect_cells_log, path_cells_log)))) %>% ungroup() %>%
  ggplot(aes(x=inoc, y=Healthiness_hsv)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(fill=UniqueID, cex=(all_plant_pixels)), alpha=0.8, col="black", pch=21) +
  theme_test() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  geom_hline(aes(yintercept=400)) +
  scale_fill_manual(values=col_all) +
  # scale_color_manual(values=c("darkgrey","yellow2","gold2","goldenrod","goldenrod4","coral4"))+
  # facet_grid(.~protect, drop=TRUE, scales="free_x") +
  # scale_color_viridis_c()+
  labs(col="Inoculation dose\n(log10 cells)", size='Plant size\n(pixels)') +
  xlab("Inoculant") + ylab("Health Score")+
  guides(fill="none")+
  scale_radius(range=c(0.1,5))
gg_pilot_n2c3wcs
ggsave(filename = "03_EDA_and_adj/gg_pilot_n2c3wcs.png", gg_pilot_n2c3wcs, width=4, height=3)

# Boxplot of N2C3, WCS365, MOCK
allDat_filt %>%
  filter(StrainMix %in% c("MOCK-MOCK", "MOCK-N2C3", "WCS365-MOCK")) %>%
  filter(path_od%in% c(0.01, 0), protect_od %in% c(0.01, 0)) %>%
  ggplot() +
  geom_boxplot(aes(x=StrainMix, y=Healthiness_rgb))+
  geom_jitter(aes(x=StrainMix, y=Healthiness_rgb), height=0)+
  facet_grid(plant_age~experiment) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
allDat_filt %>%
  filter(plant_age==7) %>%
  filter(StrainMix %in% c("MOCK-MOCK", "MOCK-N2C3", "WCS365-MOCK", "WCS365-N2C3")) %>%
  # filter(path_od%in% c(0.01, 0), protect_od %in% c(0.01, 0)) %>%
  ggplot() +
  geom_boxplot(aes(x=ratio, y=Healthiness_rgb))+
  geom_jitter(aes(x=ratio, y=Healthiness_rgb), height=0)+
  facet_grid(.~StrainMix, drop=TRUE, scales="free_x") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
#
# #### Looks like first set is weird. Going to remove 2023-05-03, 2023-05-10, and 2023-05-16
# unique(allDat_filt$experiment)
# allDat_filt <- allDat_filt %>%
#   filter(!experiment %in% c("2023-05-03_age_concentration","2023-05-16_basiceffects_n2c3", "2023-05-10_age_concentration"))

#### Look at histogram of "healthiness" ###
# LIST ALL EXPERIMENTS USED
unique(allDat_filt$experiment)
toKeepExp
toKeepUsed <- c("2023-06-14_age_concentration_Pf5_PROPER","2023-06-21_ac_ccpwn", "2023_06-27_ac_ccpwn"
                , "2023-07-05_plantmutants"
                ,"2023-07-26_ac_ccpwn", "2023-08-02_ac_ccpwn"
                ,"2023-09-13_fluor","2023-09-14_fluor")

ggHist <- allDat_filt %>%
  filter(experiment %in% c(toKeepExp, toKeepUsed)) %>%
  # filter(path=="N2C3", protect=="MOCK") %>%
  ggplot() +
  geom_histogram(aes(x=Healthiness_rgb), bins=100) +
  ylab("Number of plants") +
  xlab("") +
  theme(axis.text.x = element_blank())
ggCol <- allDat_filt %>%
  filter(experiment %in% c(toKeepExp, toKeepUsed)) %>%
  ggplot() +
  geom_jitter(aes(x=Healthiness_rgb, y=1, col=UniqueID, size=all_plant_pixels), height=0.1, width=0, show.legend=FALSE)+
  scale_color_manual(values=col_all)+
  ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Health score\n((Green:Blue + 2*(Green:Red))*log10(Number of pixels))")

gg_Healthiness <- plot_grid(ggHist, ggCol, ncol=1, align='v')
gg_Healthiness
ggsave(filename = "03_EDA_and_adj/gg_Healthiness_rgb.png",
       gg_Healthiness, width=12, height=6)

ggHist <- allDat_filt %>%
  filter(experiment %in% c(toKeepExp, toKeepUsed)) %>%
  # filter(path=="N2C3", protect=="MOCK") %>%
  ggplot() +
  geom_histogram(aes(x=Healthiness_hsv), bins=100) +
  ylab("Number of plants") +
  xlab("") +
  theme(axis.text.x = element_blank())
ggCol <- allDat_filt %>%
  filter(experiment %in% c(toKeepExp, toKeepUsed)) %>%
  ggplot() +
  geom_jitter(aes(x=Healthiness_hsv, y=1, col=UniqueID, size=all_plant_pixels), height=0.1, width=0, show.legend=FALSE)+
  scale_color_manual(values=col_all)+
  ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Health score\n(Hue*Saturation*log10(Plant pixels + 1))")

gg_Healthiness <- plot_grid(ggHist, ggCol, ncol=1, align='v')
gg_Healthiness
ggsave(filename = "03_EDA_and_adj/gg_Healthiness_hsv.png",
       gg_Healthiness, width=12, height=6)


ggHist_pca <- allDat_filt %>%
  filter(experiment %in% c(toKeepExp)) %>%
  # filter(path=="N2C3", protect=="MOCK") %>%
  ggplot() +
  geom_histogram(aes(x=Healthiness_hsv), bins=100) +
  ylab("Number of plants") +
  xlab("") +
  theme(axis.text.x = element_blank())
ggCol_pca <- allDat_filt %>%
  filter(experiment %in% c(toKeepExp)) %>%
  ggplot() +
  geom_jitter(aes(x=Healthiness_hsv, y=1, col=UniqueID, size=all_plant_pixels), height=0.1, width=0, show.legend=FALSE)+
  scale_color_manual(values=col_all)+
  ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Health score\n(Hue*Saturation*log10(Plant pixels + 1))")

gg_Healthiness_pca <- plot_grid(ggHist_pca, ggCol_pca, ncol=1, align='v')
gg_Healthiness_pca
ggsave(filename = "03_EDA_and_adj/gg_Healthiness_hsv_pca.png",
       gg_Healthiness_pca, width=12, height=6)

######## prelim dta #########
# 
# allDat_filt %>%
#   filter(experiment %in% c("2022-12-06_athal_06_polymyxan2c3"
#                            , "2022-12-10_ahtal_96_polymyxan2c3"
#                            ,"2023-02-01_gradient_polymyxan2c3")) %>%
#   mutate(Alive =Healthiness_hsv>400) %>%
#   filter(ratio %in% c("0-0","1-0","0-1","1-1","5-1")) %>%
#   ggplot() +
#   geom_point(aes(x=ratio,  y=Healthiness_hsv)) +
#   facet_grid(.~ StrainMix, drop = TRUE, scales="free_x")

