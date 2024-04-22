library(dplyr)

### LOAD IN DATA ###
data <- read.csv( "dat/Plasmid-Maint-2-Nov-2023.csv")


#Create a new colunm in the data fram called %survical that is survived/dotted*100

data$survivalRate <- data$survived / data$dotted * 100

## ---- 7 day plot --- ####

day7 <- data %>% filter(Day == "7")


# Calculate mean and standard error for each strain (all experiments grouped)
day7 <- day7 %>%
  group_by(strain) %>%
  summarise(meanSurvivalRate = mean(survivalRate),
            se = sd(survivalRate)/sqrt(n())) # Assuming 'survivalRate'

day7$strain <- reorder(day7$strain, -day7$meanSurvivalRate)


# Create the plot
gg_Plasmid_Maint_combined <- ggplot(day7, aes(x = strain, y = meanSurvivalRate, fill = strain)) + 
  # Background bar for potential to reach 100%
  geom_bar(aes(y = 100), stat = "identity", fill = "lightgrey") +
  # Overlay actual survival percentage
  geom_bar(stat = "identity") +
  # Error bars
  geom_errorbar(aes(ymin = meanSurvivalRate - se, ymax = meanSurvivalRate + se), width = 0.2) +
  theme_minimal() +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
  ) +
  labs(
    title = "", 
    x = "Strain", 
    y = "% Survival on Tet Plates"
  ) +
  scale_fill_manual(values = c(
    "WCS365-Neon" = "darkseagreen4", 
    "WCS365-Crim" = "darkseagreen4", 
    "N2C3-Crim" = "darkorange", 
    "N2C3-Neon" = "darkorange",
    "N2C3" = "darkorange",
    "WCS365" = "darkseagreen4"
  )) + # Assigning colors to strains
  guides(fill = FALSE) # Remove the legend

gg_Plasmid_Maint_combined
ggsave("gg_Plasmid_Maint_combined.png", gg_Plasmid_Maint_combined, height=6, width=10)


## ---- 0 day plot --- ####

day0 <- data %>% filter(Day == "0" & !strain %in% c("N2C3", "WCS365"))#exclude n2c3 and wcs


# Calculate mean and standard error for each strain (grouped by strain and experiment)
day0 <- day0%>%
  group_by(strain) %>%
  summarise(meanSurvivalRate = mean(survivalRate),
            se = sd(survivalRate)/sqrt(n())) # Assuming 'survivalRate'

#view(day7)

# Define the order manually
desired_order <- c("WCS365-Neon", "WCS365-Crim", "N2C3-Crim", "N2C3-Neon")

# Apply the order to the 'strain' factor
day0$strain <- factor(day0$strain, levels = desired_order)


gg_Plasmid_Maint_day0 <- ggplot(day0, aes(x = strain, y = meanSurvivalRate, fill = strain)) + 
  # Background bar for potential to reach 100%
  geom_bar(aes(y = 100), stat = "identity", fill = "lightgrey", width = 0.8) +
  # Overlay actual survival percentage
  geom_bar(stat = "identity", width = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = meanSurvivalRate - se, ymax = meanSurvivalRate + se), width = .2) +
  theme_minimal() +
  theme(
    plot.background = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 19),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
  ) +
  labs(
    title = "", 
    x = "Strain", 
    y = "% Survival on Tet Plates"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("WCS365-Neon" = "darkseagreen4", 
                               "WCS365-Crim" = "darkseagreen4", 
                               "N2C3-Crim" = "darkorange", 
                               "N2C3-Neon" = "darkorange"))+ # Assigning colors to strains
  guides(fill = FALSE) # Remove the legend

gg_Plasmid_Maint_day0 
ggsave("gg_Plasmid_Maint_day0 .png", gg_Plasmid_Maint_day0 , height=6, width=10)