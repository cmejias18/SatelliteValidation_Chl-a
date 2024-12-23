#This R Code is for the analysis of the data published in "Accuracy of remote sensing of chlorophyll-a in coral reef waters relative to in situ data and assessment of eutrophication threshold concentrations in Puerto Rico". This project aims to provide a map for using in situ Chl-a data to validate shallow coral reef virtual stations around Puerto Rico for management and restoration purposes. Lines 8-45 of this code define the function that allows for the processing of extracted pixels from Sentinel 3 OLCI ocean color products (Chlorophyll-a), as described in the publication and following the process flow diagram (Figure S.1.) therein. 

# Carla L. Mejias-Rivera & Travis A. Courtney
# carla.mejias@upr.edu


####F(x): Process Extracted Pixels####
process_pixEx <- function(file_path) {
  
  # Loads necessary packages
  library(readr)
  library(dplyr)
  library(tidyr)
  
  # Defines the flags to filter
  Flags_NN <- c("WQSF_lsb_LAND", 
                "WQSF_lsb_CLOUD", 
                "WQSF_lsb_CLOUD_AMBIGUOUS",
                "WQSF_lsb_CLOUD_MARGIN",
                "WQSF_lsb_INVALID",
                "WQSF_lsb_COSMETIC", 
                "WQSF_lsb_SATURATED", 
                "WQSF_lsb_SUSPECT",
                "WQSF_lsb_HISOLZEN", 
                "WQSF_lsb_HIGHGLINT",
                "WQSF_lsb_OCNN_FAIL", 
                "WQSF_lsb_MEGLINT",
                "WQSF_lsb_COASTLINE",
                "WQSF_lsb_WHITECAPS",
                "WQSF_lsb_ADJAC")
  
  # Reads and process the Extracted Pixel file
  processed_data <- read_csv(file_path, col_types = cols(`Date` = col_date(format = "%m/%d/%Y"))) %>% 
    mutate(ALL_FLAGS = rowSums(select(.,all_of(Flags_NN)), na.rm = TRUE)) %>% 
    filter(ALL_FLAGS < 1) %>% 
    rename(Site = "Name",
           Kd490 = "KD490_M07") %>% 
    group_by(Date, Site) %>% 
    subset(select = c(Site, Latitude, Longitude, Date, CHL_NN, CHL_NN_err, Kd490, ALL_FLAGS)) %>% 
    dplyr::mutate(
      Optical_Depth = 1 / Kd490) %>%
    replace_with_na_all(condition = ~.x == -999)
  
  return(processed_data)
}


####Load Packages to be used in the Data Analysis####
library(tidyverse)
library(dplyr)
library(naniar)
library(cowplot)
library(ggplot2)
library(tigris)
library(viridis)
library(patchwork)
library(readxl)
library(ggthemes)
library(ggrepel)
library(stringr)
library(plotly)
library(ggstatsplot)
library(geosphere)
library(ggpubr)
library(nlme)
library(broom)
library(emmeans)
library(ggpmisc)
library(ggbeeswarm)
library(patchwork)
library(ggforce)
library(ggspatial)


####Import & Modify Raw Data #### 
######In Situ Results######
InSituData=read_xlsx("PRCRMP_Chla.xlsx") %>% 
  rename(Chla = "Chl_a_µg_L",
         Site = "Site_Name",
         Quarter = "Sampling_Cycle") %>% 
  filter(Chla<15)
InSituData$Date <- as.Date(InSituData$Date, format = "%m/%d/%Y")

InSituData_Virtual=read_xlsx("PRCRMP_Chla.xlsx") %>%
  rename(Chla = "Chl_a_µg_L",
         Site = "Site_Name",
         Quarter = "Sampling_Cycle") %>% 
  filter(Chla<15) %>% 
  mutate(Site=str_remove_all(Site," 5m")) %>% 
  mutate(Site=str_remove_all(Site," 05m")) %>% 
  mutate(Site=str_remove_all(Site," 10m")) %>% 
  mutate(Site=str_remove_all(Site," 15m")) %>% 
  mutate(Site=str_remove_all(Site," 20m")) %>% 
  mutate(Site=str_remove_all(Site," 30m")) 
InSituData_Virtual$Date <- as.Date(InSituData_Virtual$Date, format = "%m/%d/%Y")


######Extracted Pixels Sentinel_3 OLCI L2 using the function#####
PRCRMP_Flags <- process_pixEx("pixEx_PRCRMP_OL_2_WFR_measurements.csv")
Virtual_Flags <- process_pixEx("pixEx_Virtual_OL_2_WFR_measurements.csv")


##### Merge InSitu and Satellite Dataframes ######
PRCRMP_Merged = merge(InSituData, PRCRMP_Flags, by = c('Site', 'Date'), all.x = TRUE) %>%
  replace_with_na_all(condition = ~.x == -999)
Virtual_Merged = merge(InSituData_Virtual, Virtual_Flags, by = c('Site', 'Date'), all.x = TRUE) %>%
  replace_with_na_all(condition = ~.x == -999)


####Figure 1: Map of Puerto Rico####
options(tigris_class = "sf")
PRmap=counties(state = "PR", cb = TRUE)

site_code_map=as.data.frame(cbind(Longitude=InSituData$Longitude,Latitude=InSituData$Latitude,Site_Code=substr(InSituData$Site_Code,1,4)))%>% 
  filter() %>% 
  group_by(Site_Code) %>% 
  summarize(Longitude=mean(as.numeric(Longitude)),Latitude=mean(as.numeric(Latitude)))

PRCRMP_map=
  ggplot()+
  geom_sf(data=PRmap, aes(geometry = geometry), fill = '#A9A9A9')+
  geom_point(InSituData, mapping=aes(x= Longitude, y = Latitude), size= 3,alpha=0.9, shape = 21, fill = "coral2")+
  geom_rect(aes(xmin = -67.4898, xmax = -67.4728, ymin = 18.3915, ymax = 18.3758), colour = "black", fill="transparent", linewidth = 0.5)+ #Desecheo
  geom_rect(aes(xmin = -67.0805, xmax = -66.7200, ymin = 17.8580, ymax = 18.0200), colour = "black", fill="transparent", linewidth = 0.5)+ #SW
  geom_rect(aes(xmin = -67.2990, xmax = -67.1400, ymin = 18.3858, ymax = 18.1288), colour = "black", fill="transparent", linewidth = 0.5)+ #W
  scale_y_continuous(limits = c(17.75,18.65), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-67.55,-65.16), expand = c(0, 0))+
  xlab("")+
  ylab("")+
  theme_tufte()+
  theme(text = element_text(size=14,family="sans"),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"))+
  annotation_north_arrow(location="tr", style = north_arrow_minimal(line_width =2, line_col = "black", text_col = "black", text_size = 18))
PRCRMP_map

PRCRMP_Desecheomap=
  ggplot()+
  geom_sf(data=PRmap, aes(geometry = geometry), fill = '#A9A9A9')+
  geom_point(InSituData, mapping=aes(x= Longitude, y = Latitude), 
             size=4, alpha=0.9, shape = 21, fill = "coral2")+
  coord_sf(xlim = c(-67.4898, -67.4728), ylim = c(18.3915, 18.3758), expand = FALSE) +
  xlab("")+
  ylab("")+
  theme_tufte()+
  theme(text = element_text(size=14,family="sans"),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))
PRCRMP_Desecheomap

PRCRMP_Wmap=
  ggplot()+
  geom_sf(data=PRmap, aes(geometry = geometry), fill = '#A9A9A9')+
  geom_point(InSituData, mapping=aes(x= Longitude, y = Latitude), 
             size=4,alpha=0.9, shape = 21, fill = "coral2")+
  coord_sf(xlim = c(-67.2990, -67.1400), ylim = c(18.3858, 18.1288), expand = TRUE) +
  xlab("")+
  ylab("")+
  theme_tufte()+
  theme(text = element_text(size=14,family="sans"),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
PRCRMP_Wmap

PRCRMP_SWmap=
  ggplot()+
  geom_sf(data=PRmap, aes(geometry = geometry), fill = '#A9A9A9')+
  geom_point(InSituData, mapping=aes(x= Longitude, y = Latitude), 
             size= 4 ,alpha=0.9, shape = 21, fill = "coral2")+
  coord_sf(xlim = c(-67.0805, -66.7200), ylim = c(17.8580, 18.0200), expand = FALSE) +
  xlab("")+
  ylab("")+
  theme_tufte()+
  theme(text = element_text(size=14,family="sans"),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))
PRCRMP_SWmap

ggsave('Figure_1.jpg', PRCRMP_map, width = 10, height = 6, dpi = 300)
ggsave('Figure_1_Desecheo.jpg', PRCRMP_Desecheomap, width = 10, height = 6, dpi = 300)
ggsave('Figure_1_SW.jpg', PRCRMP_SWmap, width = 10, height = 6, dpi = 300)
ggsave('Figure_1_W.jpg', PRCRMP_Wmap, width = 10, height = 6, dpi = 300)


####Figure 2: PRCRMP/Virtual Sites vs Sat Chla + Bland Altman Data Agreement & corresponding Stats ####

##PRCRMP Sites vs Sat Chla##

# Linear regression
PRCRMP_lm <- lm(CHL_NN ~ Chla, data = PRCRMP_Merged)
summary(PRCRMP_lm)


PRCRMP_Chl <- PRCRMP_Merged %>% 
  ggplot()+
  geom_errorbar(mapping=aes(x=Chla,y=CHL_NN,ymin = CHL_NN-CHL_NN_err, ymax = CHL_NN+CHL_NN_err),linewidth=1, color="coral2")+
  geom_errorbarh(mapping=aes(xmax=Chla+Chla*0.237,xmin=Chla-Chla*0.237,y=CHL_NN), linewidth=1, color = "coral2")+
  geom_point(mapping=aes(x=Chla, y=CHL_NN), size=4, shape = 21, fill = "coral2")+
  coord_cartesian(xlim=c(0,5),ylim=c(0,5))+
  geom_abline(intercept = 0, 
              slope = 1,
              lty =1,
              col="black")+
  theme_classic()+
  labs(x = "In Situ Chl-a (µg/L)", y = "Satellite PRCRMP Chl-a (µg/L)")+
  theme(
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

ggplotly(PRCRMP_Chl)


##Virtual Sites vs Sat Chla##

# Linear regression
Virtual_lm <- lm(CHL_NN ~ Chla, data = Virtual_Merged)
summary(Virtual_lm)

conf_int <- confint(Virtual_lm)

slope <- coef(Virtual_lm)["Chla"]
slope_se <- summary(Virtual_lm)$coefficients["Chla", "Std. Error"]
slope_ci <- conf_int["Chla", ]


Virtual_Chl <- Virtual_Merged %>% 
  ggplot(aes(x = Chla, y = CHL_NN)) +
  geom_errorbar(mapping=aes(x=Chla,y=CHL_NN,ymin = CHL_NN-CHL_NN_err, ymax = CHL_NN+CHL_NN_err),linewidth=1, color = "cyan3")+
  geom_errorbarh(mapping=aes(xmax=Chla+Chla*0.237,xmin=Chla-Chla*0.237,y=CHL_NN), linewidth=1, color = "cyan3")+
  geom_point(mapping=aes(x=Chla, y=CHL_NN), size= 4, shape = 21, fill = "cyan3")+
  geom_smooth(method = "lm", se = TRUE, color = "#252525", fill = "#636363") +
  geom_abline(intercept = 0, slope = 1, lty = 1, col = "black") +
  coord_cartesian(xlim=c(0,5),ylim=c(0, 5))+
  geom_abline(intercept = 0, 
              slope = 1,
              lty =1,
              col="black")+
  annotate("text", 
           x = 4.4, y = 4.8, 
           label = paste0("Slope = ", round(slope, 3), 
                          " ± ", round(slope_se, 3), 
                          "\n95% CI: (", 
                          round(slope_ci[1], 3), ", ", 
                          round(slope_ci[2], 3), ")"), 
           size = 5, 
           hjust = 1, color = "black") +
  labs(x = "In Situ Chl-a (µg/L)", 
       y = "Satellite Virtual Chl-a (µg/L)")+
  theme_classic()+
  theme(
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

Virtual_Chl

ggplotly(Virtual_Chl)


Figure_2A<-PRCRMP_Chl + Virtual_Chl +
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(size = 18))
Figure_2A



#ANOVA + Tukey HSD

#PRCRMP 
# PRCRMP_long <- data.frame(value = c(PRCRMP_Merged$Chla, PRCRMP_Merged$CHL_NN),
#                           group = rep(c("Chla", "CHL_NN"), each = nrow(PRCRMP_Merged)))
# 
# PRCRMP_aov <- aov(value ~ group, data = PRCRMP_long)
# PRCRMP_tukey <- TukeyHSD(PRCRMP_aov)
# 
# summary(PRCRMP_aov)
# PRCRMP_tukey


#Virtual
# Virtual_long <- data.frame(value = c(Virtual_Merged$Chla, Virtual_Merged$CHL_NN),
#                            group = rep(c("Chla", "CHL_NN"), each = nrow(Virtual_Merged)))
# 
# Virtual_aov <- aov(value ~ group, data = Virtual_long)
# Virtual_tukey <- TukeyHSD(Virtual_aov)
# 
# summary(Virtual_aov)
# Virtual_tukey


#Bland-Altman Plots
#PRCRMP
PRCRMP1 <- PRCRMP_Merged$Chla
PRCRMP2 <- PRCRMP_Merged$CHL_NN

PRCRMP_valid_data <- complete.cases(PRCRMP1, PRCRMP2)
PRCRMP1 <- PRCRMP1[PRCRMP_valid_data]
PRCRMP2 <- PRCRMP2[PRCRMP_valid_data]

PRCRMP_differences <- PRCRMP1 - PRCRMP2
PRCRMP_averages <- (PRCRMP1 + PRCRMP2) / 2

# Calculate the mean difference (bias) and the limits of agreement
PRCRMP_mean_diff <- mean(PRCRMP_differences)
PRCRMP_sd_diff <- sd(PRCRMP_differences)
PRCRMP_loa_upper <- PRCRMP_mean_diff + 1.96 * PRCRMP_sd_diff
PRCRMP_loa_lower <- PRCRMP_mean_diff - 1.96 * PRCRMP_sd_diff

PRCRMP_data <- data.frame(PRCRMP_averages, PRCRMP_differences)

PRCRMP_ribbon_data <- data.frame(
  x = seq(0, 2.5, length.out = 100), 
  ymin = rep(PRCRMP_loa_lower, 100), 
  ymax = rep(PRCRMP_loa_upper, 100))

PRCRMP_BA <- ggplot() +
  geom_ribbon(data = PRCRMP_ribbon_data, 
              aes(x = x, ymin = ymin, ymax = ymax), 
              fill = "lightgray", alpha = 0.5) +
  geom_point(data = PRCRMP_data, 
             aes(x = PRCRMP_averages, y = PRCRMP_differences), 
             shape = 21, fill = "coral2", size = 4) +
  geom_hline(yintercept = 0, 
             color = "black", linetype = "solid", linewidth = 1) +
  geom_hline(yintercept = PRCRMP_mean_diff, 
             color = "black", linetype = "solid", linewidth = 2) +
  geom_hline(yintercept = PRCRMP_loa_upper, 
             color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = PRCRMP_loa_lower, 
             color = "black", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(limits = c(-5, 5)) +
  scale_x_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by=0.5), expand = c(0,0)) +
  labs(x = "Mean of In Situ and Satellite PRCRMP Chl-a (µg/L)", 
       y = "In Situ - Satellite PRCRMP Chl-a (µg/L)") +
  theme_classic() +
  theme(
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16, angle = 90, hjust = 0.5), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

PRCRMP_BA

print(PRCRMP_mean_diff)
print(PRCRMP_loa_upper)
print(PRCRMP_loa_lower)

# Create grouping variable based on ranges of PRCRMP_averages
PRCRMP_data <- PRCRMP_data %>%
  mutate(group = cut(PRCRMP_averages, 
                     breaks = seq(0, 2.5, by = 0.5), 
                     labels = paste0("Group", 1:5), 
                     include.lowest = TRUE))

# Run ANOVA
PRCRMP_aov <- aov(PRCRMP_differences ~ group, data = PRCRMP_data)

# Run Tukey post hoc test
PRCRMP_tukey <- TukeyHSD(PRCRMP_aov)

# Print results
summary(PRCRMP_aov)
PRCRMP_tukey


#Virtual Sites
Virtual1 <- Virtual_Merged$Chla
Virtual2 <- Virtual_Merged$CHL_NN

Virtual_valid_data <- complete.cases(Virtual1, Virtual2)
Virtual1 <- Virtual1[Virtual_valid_data]
Virtual2 <- Virtual2[Virtual_valid_data]

Virtual_differences <- Virtual1 - Virtual2
Virtual_averages <- (Virtual1 + Virtual2) / 2

# Calculate the mean difference (bias) and the limits of agreement
Virtual_mean_diff <- mean(Virtual_differences)
Virtual_sd_diff <- sd(Virtual_differences)
Virtual_loa_upper <- Virtual_mean_diff + 1.96 * Virtual_sd_diff
Virtual_loa_lower <- Virtual_mean_diff - 1.96 * Virtual_sd_diff

Virtual_data <- data.frame(Virtual_averages, Virtual_differences)

Virtual_ribbon_data <- data.frame(
  x = seq(0, 1.0, length.out = 100), 
  ymin = rep(Virtual_loa_lower, 100), 
  ymax = rep(Virtual_loa_upper, 100))

Virtual_BA<-ggplot() +
  geom_ribbon(data = Virtual_ribbon_data, 
              aes(x = x, ymin = ymin, ymax = ymax), 
              fill = "lightgray", alpha = 0.5) +
  geom_point(data = Virtual_data, 
             aes(x = Virtual_averages, y = Virtual_differences), 
             size=4, shape = 21, fill = "cyan3") +
  geom_hline(yintercept = 0, 
             color = "black", linetype = "solid", linewidth = 1) +
  geom_hline(yintercept = Virtual_mean_diff, 
             color = "black", linetype = "solid", size=2) +
  geom_hline(yintercept = Virtual_loa_upper, 
             color = "black", linetype = "dashed", size=1) +
  geom_hline(yintercept = Virtual_loa_lower, 
             color = "black", linetype = "dashed", size=1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-5, 5))+
  # scale_y_continuous(limits = c(-5, 5))+ 
  # scale_x_continuous(limits=c(0,1.0), breaks = seq(0,1.0, by=0.5), expand = c(0, 0))+
  labs(x = "Mean of In Situ and Satellite Virtual Chl-a (µg/L)", 
       y = "In Situ - Satellite Virtual Chl-a (µg/L)") +
  theme_classic() +
  theme(
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16, angle = 90, hjust = 0.5), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

Virtual_BA

print(Virtual_mean_diff)
print(Virtual_loa_upper)
print(Virtual_loa_lower)

# Create grouping variable based on ranges of PRCRMP_averages
Virtual_data <- Virtual_data %>%
  mutate(group = cut(Virtual_averages, 
                     breaks = seq(0, 2.5, by = 0.5), 
                     labels = paste0("Group", 1:5), 
                     include.lowest = TRUE))

# Run ANOVA
Virtual_aov <- aov(Virtual_differences ~ group, data = Virtual_data)

# Run Tukey post hoc test
Virtual_tukey <- TukeyHSD(Virtual_aov)

# Print results
summary(Virtual_aov)
Virtual_tukey


Figure_2B<-PRCRMP_BA + Virtual_BA+
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(size = 18))
Figure_2B

Fig_2 <- Figure_2A/Figure_2B+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18))
Fig_2

ggsave("Figure_2.jpg", Fig_2, width = 12, height = 12, dpi = 300)


####Figure 3. Chl-a ~ Depth Factorial Analysis ####

#Add Group Variable
InSituData$Group <- "InSitu"
PRCRMP_Merged$Group <- "PRCRMP"
Virtual_Merged$Group <- "Virtual"

Factorial_columns <- c("Depth", "CHL_NN", "Chla", "Group")

InSituData$Depth <- as.factor(InSituData$Depth)
PRCRMP_Merged$Depth <- as.factor(PRCRMP_Merged$Depth)
Virtual_Merged$Depth <- as.factor(Virtual_Merged$Depth)

for (col in Factorial_columns) {
  if (!col %in% names(InSituData)) InSituData[[col]] <- NA
  if (!col %in% names(PRCRMP_Merged)) PRCRMP_Merged[[col]] <- NA
  if (!col %in% names(Virtual_Merged)) Virtual_Merged[[col]] <- NA
}

InSituData_stats <- InSituData[Factorial_columns]
PRCRMP_Merged_stats <- PRCRMP_Merged[Factorial_columns]
Virtual_Merged_stats <- Virtual_Merged[Factorial_columns]

Factorial_combined_data <- rbind(InSituData_stats, PRCRMP_Merged_stats, Virtual_Merged_stats)

Factorial_combined_data$Group <- as.factor(Factorial_combined_data$Group)

factorial_aov <- aov(CHL_NN ~ Depth * Group, data = Factorial_combined_data)
summary(factorial_aov)

InSitu_factorial_aov <- aov(Chla ~ as.factor(Depth), data = InSituData_stats)
PRCRMP_factorial_aov <- aov(CHL_NN ~ as.factor(Depth), data = PRCRMP_Merged_stats)
Virtual_factorial_aov <- aov(CHL_NN ~ as.factor(Depth), data = Virtual_Merged_stats)

summary(InSitu_factorial_aov)
summary(PRCRMP_factorial_aov)
summary(Virtual_factorial_aov)

# Post-hoc test
TukeyHSD(InSitu_factorial_aov)
TukeyHSD(PRCRMP_factorial_aov)
TukeyHSD(Virtual_factorial_aov)


#Plot interaction

## Plot for InSituData
p1 <- ggplot(
  data= InSituData_stats, aes(x= Depth, y= Chla))+
  geom_sina(size =3, color = "black", fill= "azure4", shape = 21)+
  geom_boxplot(color = "black", fill = NA, alpha = 0.5)+
  xlab("Depth")+ 
  ylab("Chl-a (µg/L)")+
  theme_classic() + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")+
  coord_cartesian(ylim = c(0, 5))
p1

#Plot for PRCRMP_Merged
p2 <- ggplot(
  data= PRCRMP_Merged_stats, aes(x= Depth, y= CHL_NN, color=))+
  geom_sina(size =3, color = "black", fill="coral2", shape = 21)+
  geom_boxplot(color = "black", fill = NA, alpha = 0.5)+
  xlab("Depth")+ 
  ylab("Chl-a (µg/L)")+
  theme_classic() + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")+
  coord_cartesian(ylim = c(0, 5))
p2

#Plot for Virtual_Merged
p3 <- ggplot(
  data= Virtual_Merged_stats, aes(x= Depth, y= CHL_NN))+
  geom_sina(size =3, color = "black", fill="cyan3", shape=21)+
  geom_boxplot(color = "black", fill = NA, alpha = 0.5)+
  xlab("Depth")+ 
  ylab("Chl-a (µg/L)")+
  theme_classic() + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")+
  coord_cartesian(ylim = c(0, 5))
p3


Figure3 <- p1 + p2 + p3 + plot_layout(ncol = 3)+
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(size = 18))
print(Figure3)

ggsave("Figure_3.jpg", Figure3, width = 10, height = 5, dpi = 300)


####Figure 4: Monthly Analysis for years 2022-2023####
Virtual_Monthly_Flags <- process_pixEx("pixEx_Monthly_22-23_Virtual_OL_2_WFR_measurements.csv")


Virtual_Monthly_Merged<- merge(InSituData_Virtual, Virtual_Monthly_Flags, by = c('Site', 'Date'), all.x=TRUE) %>%
  replace_with_na_all(condition=~.x==-999) %>%
  filter(Quarter>1)

InSituData_Msummary=InSituData_Virtual %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), list(mean = ~mean(.x, na.rm = TRUE), 
                                           sd = ~sd(.x, na.rm = TRUE), 
                                           min =~min(.x, na.rm = TRUE), 
                                           max = ~max(.x, na.rm = TRUE))))

Virtual_Msummary=Virtual_Monthly_Flags %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), list(mean = ~mean(.x, na.rm = TRUE), 
                                           sd = ~sd(.x, na.rm = TRUE), 
                                           min = ~min(.x, na.rm = TRUE), 
                                           max = ~max(.x, na.rm = TRUE))))

Fig_Extended_InSitu <-
  ggplot(InSituData_Virtual, aes(x=Site, y=Chla)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.3, fill = "#addd8e", alpha = 0.01)+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.3, ymax = 0.5, fill = "#74c476", alpha = 0.02)+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "#74c476", alpha = 0.4)+
  geom_jitter(shape=21, fill = "#636363", size = 1) + 
  geom_boxplot(aes(group=Site), fill = "transparent", color = "black", outlier.shape = NA) +
  geom_point(data=InSituData_Msummary, aes(x=Site, y=`Chla_mean`), shape = 21, fill = "azure4", size = 3) +
  scale_y_continuous(limits = c(0,2), breaks = c(0, 0.2, 0.5, 1, 1.5, 2))+
  ylab("In Situ Chl-a (µg/L)") +
  theme_classic() +
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

Fig_Extended_InSitu

Fig_Extended_Virtual <-
  ggplot(Virtual_Monthly_Flags, aes(x=Site, y=CHL_NN)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.3, fill = "#addd8e", alpha = 0.01)+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.3, ymax = 0.5, fill = "#74c476", alpha = 0.02)+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "#74c476", alpha = 0.4)+
  geom_jitter(shape=21, fill = "cyan3", size = 1) + 
  geom_boxplot(aes(group=Site), fill = "transparent", color = "black", outlier.shape = NA) + 
  geom_point(data=Virtual_Msummary, aes(x=Site, y=`CHL_NN_mean`), shape = 21, fill = "cyan3", size = 3) + 
  scale_y_continuous(limits = c(0,2), breaks = c(0, 0.2, 0.5, 1, 1.5, 2))+
  ylab("Satellite Virtual Chl-a (µg/L)") +
  theme_classic() + 
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

Fig_Extended_Virtual

Figure_4 <- Fig_Extended_InSitu / Fig_Extended_Virtual + # use / to stack vertically
  plot_layout(ncol = 1, heights = c(1.5, 1.5)) + # heights can be adjusted if you want different relative sizes
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18))
print(Figure_4)

ggsave("Figure_4.jpg", Figure_4, width = 10, height =10, dpi = 300)


####Approximate distance between PRCRMP and Virtual stations#### 

data <- read.csv("PRCRMP_Virtual_Coords2.csv", colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
if (nrow(data)>0){
  distances <- numeric(nrow(data))
  for (i in 1:nrow(data)){
    PRCRMP_Lat <- data$PRCRMP_Lat[i]
    PRCRMP_Lon <- data$PRCRMP_Lon[i]
    Virtual_Lat <- data$Virtual_Lat[i]
    Virtual_Lon <- data$Virtual_Lon[i]
    distances[i] <- distVincentyEllipsoid(c(PRCRMP_Lat, PRCRMP_Lon), c(Virtual_Lat, Virtual_Lon))
    data$Distance<-distances
    print(data)}}

write.csv(data,"PRCRMP_vs_Virtual_Distance.csv", row.names=FALSE)
