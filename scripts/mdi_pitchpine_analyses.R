# script to analyze mdi pitch pine data

# load packages
library(tidyverse)
library(emmeans)
library(lme4)
library(car)
library(circular)
library(multcompView)
library(ggthemes)
library(agricolae)
library(patchwork)

# function to get pairwise letters from Tukey's HSD for boxplots
letters <- function(df = data, x = "Site", dfy, y){
  .checkvar <- sym(y) #convert string to variable
  .groupvar <- sym(x) #convert string to variable
  
  # find highest value of y data
  abs_max <- max(dfy, na.rm = TRUE) 
  # get the highest point for each species
  max_y <- df %>% dplyr::group_by(!! .groupvar) %>% 
    dplyr::summarise(yaxis = max(!! .checkvar, na.rm = TRUE) + 0.05 * abs_max)
    #"!!" unquotes the variable to allow for evaluation within dplyr
  # get Tukey HSD results
  hsd <- HSD.test(aov(as.formula(paste(y, x, sep = "~")), df), x, group = TRUE) 
  # add Tukey HSD results to dataframe containing graphing positions
  max_y$group <- hsd$groups$groups
  return(max_y)
}

# function to get pairwise letters from Tukey's HSD with a transformed y variable for boxplots
letters_adj <- function(df = data, x = "Site", dfy, y, adjy){
  .checkvar <- sym(y)
  .groupvar <- sym(x)
  
  abs_max <- max(dfy, na.rm = TRUE)
  max_y <- df %>% dplyr::group_by(!! .groupvar) %>% 
    dplyr::summarise(yaxis = max(!! .checkvar, na.rm = TRUE) + 0.05 * abs_max) 
  hsd <- HSD.test(aov(as.formula(paste(adjy, x, sep = "~")), df), x, group = TRUE) 
  max_y$group <- hsd$groups$groups
  return(max_y)
}

#### read in cleaned data ####

data <- read.csv('data/mdi_all_clean.csv')
data$CN_foliar <- data$C_foliar/data$N_foliar
data$CN_soil <- data$C_soil/data$N_soil

data_density <- read.csv('data/mdi_stand_density.csv')
colnames(data_density)[1] <- "Name" #rename to match data

## assign fire history status to each site
data$fire[data$Name == 'CAD'] <- 'fire' 
data$fire[data$Name == 'CADCLIFFS'] <- 'fire'
data$fire[data$Name == 'STSAUV'] <- 'no fire'
data$fire[data$Name == 'WOND'] <- 'no fire'

data_density$fire[data_density$Name == 'CAD'] <- 'fire' 
data_density$fire[data_density$Name == 'CADCLIFFS'] <- 'fire'
data_density$fire[data_density$Name == 'STSAUV'] <- 'no fire'
data_density$fire[data_density$Name == 'WOND'] <- 'no fire'

## create an elevation factor
data$elevation_fac[data$Name == 'CAD' | data$Name == 'STSAUV'] <- 'high'
data$elevation_fac[data$Name == 'CADCLIFFS' | data$Name == 'WOND'] <- 'low'
data_density$elevation_fac[data_density$Name == 'CAD' | data_density$Name == 'STSAUV'] <- 'high'
data_density$elevation_fac[data_density$Name == 'CADCLIFFS' | data_density$Name == 'WOND'] <- 'low'

## reorder levels for elevation factor from "low" to "high"
data$elevation_fac <- factor(data$elevation_fac, levels = c("low", "high"))
data_density$elevation_fac <- factor(data_density$elevation_fac, levels = c("low", "high"))

## rename site labels to match manuscript
data$Site[data$Name == "CAD"] <- "SCT"
data$Site[data$Name == "CADCLIFFS"] <- "GOR"
data$Site[data$Name == "STSAUV"] <- "STS"
data$Site[data$Name == "WOND"] <- "WON"

data_density$Site[data_density$Name == "CAD"] <- "SCT"
data_density$Site[data_density$Name == "CADCLIFFS"] <- "GOR"
data_density$Site[data_density$Name == "STSAUV"] <- "STS"
data_density$Site[data_density$Name == "WOND"] <- "WON"

## reorder levels of sites to match to increase from low to high elevation and no fire to fire
data$Site <- factor(data$Site, levels = c("WON", "GOR", "STS", "SCT"))
data_density$Site <- factor(data_density$Site, levels = c("WON", "GOR", "STS", "SCT"))

## create a generic variable set to pass to formula argument
ind_variables <- c('elevation_fac', 'fire')
dep_variables <- c("log(Elevation)", "log(height)", "log(canopy)", "log(diam)",
                  "d13C", "d15N", "C_foliar", "N_foliar", "CN_foliar", "Ca_foliar", "log(P_foliar)",
                  "log(K_foliar)", "Mg_foliar", "Al_foliar", "log(Zn_foliar)", 
                  "Ca_soil", "log(P_soil)", "K_soil", "Mg_soil", "log(Al_soil)", "log(Zn_soil)", 
                  "pH", "CEC", "C_soil", "N_soil", "log(CN_soil)", "asin(sqrt(0.01 * retention))")

#### fit models and explore results ####

## topography
### elevation
Elevation_lm <- lm(as.formula(paste("log(Elevation)",
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Elevation_lm) ~ fitted(Elevation_lm))
Anova(Elevation_lm)
cld.emmGrid(emmeans(Elevation_lm, ~elevation_fac * fire))

ggplot(data = data, aes(x = Site, y = log(Elevation))) +
  geom_boxplot()

### slope
Slope_lm <- lm(as.formula(paste("Slope",
                               paste(ind_variables, collapse = "*"),
                               sep = "~")), data = data)
#plot(resid(Slope_lm) ~ fitted(Slope_lm))
Anova(Slope_lm)
cld.emmGrid(emmeans(Slope_lm, ~elevation_fac * fire))

ggplot(data = data, aes(x = Site, y = Slope)) +
  geom_boxplot()

### aspect
#### turn aspect data into circular data that maps onto a compass
aspect_SCT <- circular(data$Aspect[data$Name == 'CAD'],
                      units = "degrees", template = "geographics")
aspect_GOR <- circular(data$Aspect[data$Name == 'CADCLIFFS'],
                            units = "degrees", template = "geographics")
aspect_STS <- circular(data$Aspect[data$Name == 'STSAUV'],
                         units = "degrees", template = "geographics")
aspect_WON <- circular(data$Aspect[data$Name == 'WOND'],
                       units = "degrees", template = "geographics")

#### Watson's Two Sample Test of Homogeneity (https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/FullLesson.html)
#### tests whether the ‘north’ and ‘south’ groups orient in different directions
watson.two.test(aspect_WON, aspect_GOR) # 0.01 < P < 0.05
watson.two.test(aspect_WON, aspect_STS) # 0.001 < P < 0.01
watson.two.test(aspect_WON, aspect_SCT) # 0.01 < P < 0.05
watson.two.test(aspect_GOR, aspect_STS) # 0.001 < P < 0.01
watson.two.test(aspect_GOR, aspect_SCT) # 0.05 < P < 0.1
watson.two.test(aspect_STS, aspect_SCT) # P < 0.001

#### circular plots for each site
jpeg(filename = "analyses/plots/plots_aspect.jpeg", width = 3000, 
     height = 3000, units = 'px')
par(mfrow = c(2, 2), cex.lab = 6, cex.main = 6, mar = c(5.5, 8.5, 5.5, 2.5))
plot_aspect_CADCLIFFS <- plot.circular(aspect_CADCLIFFS, main = 'Gorham Cliffs (a)', 
                                      ylab = "Fire", 
                                      cex = 8, col = "red", pch = 16)
plot_aspect_CAD <- plot.circular(aspect_CAD, main = 'South Cadillac (a)',
                                cex = 8, col = "red", pch = 17)
plot_aspect_WOND <- plot.circular(aspect_WOND, main = 'Wonderland (b)', 
                                 ylab = "No Fire", xlab = "Low Elevation", 
                                 cex = 8, col = "blue", pch = 16)
plot_aspect_STSAUV <- plot.circular(aspect_STSAUV, main = 'St. Sauveur (c)', 
                                   xlab = "High Elevation", 
                                   cex = 8, col = "blue", pch = 17)
dev.off()

## allometry
### height
height_lm <- lm(as.formula(paste(dep_variables[2],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(height_lm) ~ fitted(height_lm))
Anova(height_lm)
cld.emmGrid(emmeans(height_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
height_letters <- letters_adj(dfy = data$height, y = "height", adjy = "log(height)")

(plot_height <- ggplot(data = data, aes(x = Site, y = height)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = height_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) +
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Height (m)") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### canopy
canopy_lm <- lm(as.formula(paste(dep_variables[3],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
# plot(resid(canopy_lm) ~ fitted(canopy_lm))
anova(canopy_lm)
cld.emmGrid(emmeans(canopy_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
canopy_letters <- letters_adj(dfy = data$canopy, y = "canopy", adjy = "log(canopy)")

(plot_canopy <- ggplot(data = data, aes(x = Site, y = canopy)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = canopy_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Canopy Spread (m)") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### diam
diam_lm <- lm(as.formula(paste(dep_variables[4],
                              paste(ind_variables, collapse = "*"),
                              sep = "~")), data = data)
#plot(resid(diam_lm) ~ fitted(diam_lm))
anova(diam_lm)
cld.emmGrid(emmeans(diam_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
diam_letters <- letters_adj(dfy = data$diam, y = "diam", adjy = "log(diam)")

(plot_diam <- ggplot(data = data, aes(x = Site, y = diam)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = diam_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "DBH (cm)") + 
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### density
density_lm <- lm(mean_distance ~ elevation_fac * fire, data = data_density)
#plot(resid(density_lm) ~ fitted(density_lm))
anova(density_lm)
cld.emmGrid(emmeans(density_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
density_letters <- letters(df = data_density, dfy = data_density$mean_distance, y = "mean_distance")

(plot_density <- ggplot(data = data_density, aes(x = Site, y = mean_distance)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = density_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) +
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Stand Density (m)") + 
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_allometry <- plot_height + plot_canopy + plot_diam + plot_density + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## foliar organics
### C_foliar
C_foliar_lm <- lm(as.formula(paste(dep_variables[7],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(C_foliar_lm) ~ fitted(C_foliar_lm))
anova(C_foliar_lm)
cld.emmGrid(emmeans(C_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
C_foliar_letters <- letters(dfy = data$C_foliar, y = "C_foliar")

(plot_C_foliar <- ggplot(data = data, aes(x = Site, y = C_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = C_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Foliar Carbon (%)") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### N_foliar
N_foliar_lm <- lm(as.formula(paste(dep_variables[8],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(N_foliar_lm) ~ fitted(N_foliar_lm))
anova(N_foliar_lm)
cld.emmGrid(emmeans(N_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
N_foliar_letters <- letters(dfy = data$N_foliar, y = "N_foliar")

(plot_N_foliar <- ggplot(data = data, aes(x = Site, y = N_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = N_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Foliar Nitrogen (%)") +
    # theme(legend.position = "bottom", legend.title = element_text(size = 16), 
    #       legend.text = element_text(size = 12)) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### CN_foliar
CN_foliar_lm <- lm(as.formula(paste(dep_variables[9],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(CN_foliar_lm) ~ fitted(CN_foliar_lm))
anova(CN_foliar_lm)
cld.emmGrid(emmeans(CN_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
CN_foliar_letters <- letters(dfy = data$CN_foliar, y = "CN_foliar")

(plot_CN_foliar <- ggplot(data = data, aes(x = Site, y = CN_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = CN_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Foliar Carbon/Nitrogen") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_foliar_organics <- plot_C_foliar + plot_N_foliar + plot_CN_foliar + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## foliar inorganics
### Ca_foliar
Ca_foliar_lm <- lm(as.formula(paste(dep_variables[10],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Ca_foliar_lm) ~ fitted(Ca_foliar_lm))
anova(Ca_foliar_lm)
cld.emmGrid(emmeans(Ca_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Ca_foliar_letters <- letters(dfy = data$Ca_foliar, y = "Ca_foliar")

(plot_Ca_foliar <- ggplot(data = data, aes(x = Site, y = Ca_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Ca_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Foliar Calcium (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### P_foliar
P_foliar_lm <- lm(as.formula(paste(dep_variables[11],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(P_foliar_lm) ~ fitted(P_foliar_lm))
anova(P_foliar_lm)
cld.emmGrid(emmeans(P_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
P_foliar_letters <- letters_adj(dfy = data$P_foliar, y = "P_foliar", adjy = "log(P_foliar)")

(plot_P_foliar <- ggplot(data = data, aes(x = Site, y = P_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = P_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Foliar Phosphorus (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### K_foliar
K_foliar_lm <- lm(as.formula(paste(dep_variables[12],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(K_foliar_lm) ~ fitted(K_foliar_lm))
anova(K_foliar_lm)
cld.emmGrid(emmeans(K_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
K_foliar_letters <- letters_adj(dfy = data$K_foliar, y = "K_foliar", adjy = "log(K_foliar)")

(plot_K_foliar <- ggplot(data = data, aes(x = Site, y = K_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = K_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Foliar Potassium (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### Mg_foliar
Mg_foliar_lm <- lm(as.formula(paste(dep_variables[13],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Mg_foliar_lm) ~ fitted(Mg_foliar_lm))
anova(Mg_foliar_lm)
cld.emmGrid(emmeans(Mg_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Mg_foliar_letters <- letters(dfy = data$Mg_foliar, y = "Mg_foliar")

(plot_Mg_foliar <- ggplot(data = data, aes(x = Site, y = Mg_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Mg_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Foliar Magnesium (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### Al_foliar
Al_foliar_lm <- lm(as.formula(paste(dep_variables[14],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Al_foliar_lm) ~ fitted(Al_foliar_lm))
anova(Al_foliar_lm)
cld.emmGrid(emmeans(Al_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Al_foliar_letters <- letters(dfy = data$Al_foliar, y = "Al_foliar")

(plot_Al_foliar <- ggplot(data = data, aes(x = Site, y = Al_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Al_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Foliar Aluminum (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### Zn_foliar
Zn_foliar_lm <- lm(as.formula(paste(dep_variables[15],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Zn_foliar_lm) ~ fitted(Zn_foliar_lm))
anova(Zn_foliar_lm)
cld.emmGrid(emmeans(Zn_foliar_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Zn_foliar_letters <- letters_adj(dfy = data$Zn_foliar, y = "Zn_foliar", adjy = "log(Zn_foliar)")

(plot_Zn_foliar <- ggplot(data = data, aes(x = Site, y = Zn_foliar)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Zn_foliar_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Foliar Zinc (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_foliar_inorganics <- plot_Ca_foliar + plot_P_foliar + plot_K_foliar + 
    plot_Mg_foliar + plot_Al_foliar + plot_Zn_foliar +
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## foliar isotopes
### d13C
d13C_lm <- lm(as.formula(paste(dep_variables[5],
                              paste(ind_variables, collapse = "*"),
                              sep = "~")), data = data)
#plot(resid(d13C_lm) ~ fitted(d13C_lm))
anova(d13C_lm)
cld.emmGrid(emmeans(d13C_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot - not using function because of negative data
d13C_abs_max <- max(data$d13C, na.rm = TRUE)
d13C_letters <- data %>% group_by(Site) %>% 
  summarise(yaxis = max(d13C, na.rm = TRUE) - 0.02 * d13C_abs_max) # get the highest point for each species
d13C_hsd <- HSD.test(aov(d13C ~ Site, data), "Site", group = TRUE) # get Tukey HSD results
d13C_letters$group <- d13C_hsd$groups$groups

(plot_d13C <- ggplot(data = data, aes(x = Site, y = d13C)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = d13C_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression(delta^{"13"}*"C (‰)")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### d15N
d15N_lm <- lm(as.formula(paste(dep_variables[6],
                              paste(ind_variables, collapse = "*"),
                              sep = "~")), data = data)
#plot(resid(d15N_lm) ~ fitted(d15N_lm))
anova(d15N_lm)
cld.emmGrid(emmeans(d15N_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
d15N_letters <- letters(dfy = data$d15N, y = "d15N")

(plot_d15N <- ggplot(data = data, aes(x = Site, y = d15N)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = d15N_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression(delta^{"15"}*"N (‰)")) + 
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_foliar_isotopes <- plot_d13C + plot_d15N +
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil organics
### C_soil
C_soil_lm <- lm(as.formula(paste(dep_variables[24],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(C_soil_lm) ~ fitted(C_soil_lm))
anova(C_soil_lm)
cld.emmGrid(emmeans(C_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
C_soil_letters <- letters(dfy = data$C_soil, y = "C_soil")

(plot_C_soil <- ggplot(data = data, aes(x = Site, y = C_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = C_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Carbon (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### N_soil
N_soil_lm <- lm(as.formula(paste(dep_variables[25],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(N_soil_lm) ~ fitted(N_soil_lm))
anova(N_soil_lm)
cld.emmGrid(emmeans(N_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
N_soil_letters <- letters(dfy = data$N_soil, y = "N_soil")

(plot_N_soil <- ggplot(data = data, aes(x = Site, y = N_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = N_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) +
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Nitrogen (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### CN_soil
CN_soil_lm <- lm(as.formula(paste(dep_variables[26],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(CN_soil_lm) ~ fitted(CN_soil_lm))
anova(CN_soil_lm)
cld.emmGrid(emmeans(CN_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
CN_soil_letters <- letters_adj(dfy = data$CN_soil, y = "CN_soil", adjy = "log(CN_soil)")

(plot_CN_soil <- ggplot(data = data, aes(x = Site, y = CN_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = CN_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Soil Carbon/Nitrogen") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_soil_organics <- plot_C_soil + plot_N_soil + plot_CN_soil + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil inorganics
### Ca_soil
Ca_soil_lm <- lm(as.formula(paste(dep_variables[16],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Ca_soil_lm) ~ fitted(Ca_soil_lm))
anova(Ca_soil_lm)
cld.emmGrid(emmeans(Ca_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Ca_soil_letters <- letters(dfy = data$Ca_soil, y = "Ca_soil")

(plot_Ca_soil <- ggplot(data = data, aes(x = Site, y = Ca_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Ca_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Calcium (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### P_soil
P_soil_lm <- lm(as.formula(paste(dep_variables[17],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(P_soil_lm) ~ fitted(P_soil_lm))
anova(P_soil_lm)
cld.emmGrid(emmeans(P_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
P_soil_letters <- letters_adj(dfy = data$P_soil, y = "P_soil", adjy = "log(P_soil)")

(plot_P_soil <- ggplot(data = data, aes(x = Site, y = P_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = P_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Phosphorus (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### K_soil
K_soil_lm <- lm(as.formula(paste(dep_variables[18],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(K_soil_lm) ~ fitted(K_soil_lm))
anova(K_soil_lm)
cld.emmGrid(emmeans(K_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
K_soil_letters <- letters(dfy = data$K_soil, y = "K_soil")

(plot_K_soil <- ggplot(data = data, aes(x = Site, y = K_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = K_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Potassium (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### Mg_soil
Mg_soil_lm <- lm(as.formula(paste(dep_variables[19],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Mg_soil_lm) ~ fitted(Mg_soil_lm))
anova(Mg_soil_lm)
cld.emmGrid(emmeans(Mg_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Mg_soil_letters <- letters(dfy = data$Mg_soil, y = "Mg_soil")

(plot_Mg_soil <- ggplot(data = data, aes(x = Site, y = Mg_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Mg_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Magnesium (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### Al_soil
Al_soil_lm <- lm(as.formula(paste(dep_variables[20],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Al_soil_lm) ~ fitted(Al_soil_lm))
anova(Al_soil_lm)
cld.emmGrid(emmeans(Al_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Al_soil_letters <- letters_adj(dfy = data$Al_soil, y = "Al_soil", adjy = "log(Al_soil)")

(plot_Al_soil <- ggplot(data = data, aes(x = Site, y = Al_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Al_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Aluminum (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### Zn_soil
Zn_soil_lm <- lm(as.formula(paste(dep_variables[21],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Zn_soil_lm) ~ fitted(Zn_soil_lm))
anova(Zn_soil_lm)
cld.emmGrid(emmeans(Zn_soil_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
Zn_soil_letters <- letters_adj(dfy = data$Zn_soil, y = "Zn_soil", adjy = "log(Zn_soil)")

(plot_Zn_soil <- ggplot(data = data, aes(x = Site, y = Zn_soil)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = Zn_soil_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil Zinc (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_soil_inorganics <- plot_Ca_soil + plot_P_soil + plot_K_soil + 
    plot_Mg_soil + plot_Al_soil + plot_Zn_soil +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil characteristics
### pH
pH_lm <- lm(as.formula(paste(dep_variables[22],
                            paste(ind_variables, collapse = "*"),
                            sep = "~")), data = data)
#plot(resid(pH_lm) ~ fitted(pH_lm))
anova(pH_lm)
cld.emmGrid(emmeans(pH_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
pH_letters <- letters(dfy = data$pH, y = "pH")

(plot_pH <- ggplot(data = data, aes(x = Site, y = pH)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = pH_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Soil pH") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

### CEC
CEC_lm <- lm(as.formula(paste(dep_variables[23],
                             paste(ind_variables, collapse = "*"),
                             sep = "~")), data = data)
#plot(resid(CEC_lm) ~ fitted(CEC_lm))
anova(CEC_lm)
cld.emmGrid(emmeans(CEC_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
CEC_letters <- letters(dfy = data$CEC, y = "CEC")

(plot_CEC <- ggplot(data = data, aes(x = Site, y = CEC)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = CEC_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    ylab(expression("Soil CEC (cmol"[c]*" kg"^{-1}*")")) +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

## soil characteristics
### retention
retention_lm <- lm(as.formula(paste(dep_variables[27],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(retention_lm) ~ fitted(retention_lm))
anova(retention_lm)
cld.emmGrid(emmeans(retention_lm, ~elevation_fac * fire))

#### get pairwise letters for boxplot
retention_letters <- letters_adj(dfy = data$retention, y = "retention", adjy = "asin(sqrt(0.01 * retention))")

(plot_retention <- ggplot(data = data, aes(x = Site, y = retention)) +
    geom_rect(data = NULL, aes(xmin = 0, xmax = 2.5, ymin = -Inf, ymax = Inf),
              fill = "grey") +
    geom_jitter(height = 0, aes(color = fire, shape = elevation_fac), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    geom_text(data = retention_letters, aes(y = yaxis, label = group)) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Site") +
    scale_y_continuous(name = "Soil Water Retention (%)") +
    guides(color = guide_legend("Fire History")) +
    guides(shape = guide_legend("Elevation")))

(plots_soil_characteristics <- plot_retention + plot_CEC + plot_pH +
    plot_layout(guides = "collect", ncol = 3, nrow = 1) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

#### tables and posthoc ####

### topography
#### create table with mean latitude, longitude, elevation, slope, and aspect for each site
topography <- data %>% group_by(Site) %>% summarise_at(vars(latitude, longitude, Elevation, Slope, Aspect), mean, na.rm = TRUE)
write.csv(topography, "analyses/tables/topography.csv")

### allometry
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(height_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(canopy_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(diam_lm)[, c(1, 4, 5)]),
                as.matrix(anova(density_lm)[, c(1, 4, 5)])),
          'analyses/tables/allometry.csv')

### foliar organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(C_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(N_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(CN_foliar_lm)[, c(1, 4, 5)])),
          'analyses/tables/foliar_cn.csv')

### foliar inorganics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(Ca_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(P_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(K_foliar_lm)[, c(1, 4, 5)]),
                as.matrix(anova(Mg_foliar_lm)[, c(1, 4, 5)]),
                as.matrix(anova(Al_foliar_lm)[, c(1, 4, 5)]),
                as.matrix(anova(Zn_foliar_lm)[, c(1, 4, 5)])),
          'analyses/tables/foliar_inorganics.csv')

### foliar isotopes
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(d13C_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(d15N_lm)[, c(4, 5)])),
          'analyses/tables/foliar_isotope.csv')

### soil organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(C_soil_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(N_soil_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(CN_soil_lm)[, c(1, 4, 5)])),
          'analyses/tables/soil_organics.csv')

(summary(emmeans(canopy_lm, ~elevation_fac))[1,2] - summary(emmeans(canopy_lm, ~elevation_fac))[2,2])/ summary(emmeans(canopy_lm, ~elevation_fac))[2,2]

(summary(emmeans(P_foliar_lm, ~fire))[1,2] - summary(emmeans(P_foliar_lm, ~fire))[2,2])/ summary(emmeans(P_foliar_lm, ~fire))[2,2]
(summary(emmeans(K_foliar_lm, ~fire))[1,2] - summary(emmeans(K_foliar_lm, ~fire))[2,2])/ summary(emmeans(K_foliar_lm, ~fire))[2,2]
(summary(emmeans(Ca_foliar_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_foliar_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_foliar_lm, ~elevation_fac))[2,2]
(summary(emmeans(Zn_foliar_lm, ~elevation_fac))[1,2] - summary(emmeans(Zn_foliar_lm, ~elevation_fac))[2,2])/ summary(emmeans(Zn_foliar_lm, ~elevation_fac))[2,2]

(summary(emmeans(d13C_lm, ~elevation_fac))[1,2] - summary(emmeans(d13C_lm, ~elevation_fac))[2,2])/ summary(emmeans(d13C_lm, ~elevation_fac))[2,2]

(summary(emmeans(C_soil_lm, ~fire))[1,2] - summary(emmeans(C_soil_lm, ~fire))[2,2])/ summary(emmeans(C_soil_lm, ~fire))[2,2]
(summary(emmeans(C_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(C_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(C_soil_lm, ~elevation_fac))[2,2]
(summary(emmeans(CN_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(CN_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(CN_soil_lm, ~elevation_fac))[2,2]

### soil inorganics
write.csv(cbind(as.matrix(anova(Ca_soil_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(P_soil_lm)[, c(4, 5)]), 
                as.matrix(anova(K_soil_lm)[, c(4, 5)]),
                as.matrix(anova(Mg_soil_lm)[, c(4, 5)]),
                as.matrix(anova(Al_soil_lm)[, c(4, 5)]),
                as.matrix(anova(Zn_soil_lm)[, c(4, 5)])),
          'analyses/tables/soil_inorganics.csv')

(summary(emmeans(K_soil_lm, ~fire))[1,2] - summary(emmeans(K_soil_lm, ~fire))[2,2])/ summary(emmeans(K_soil_lm, ~fire))[2,2]
(summary(emmeans(Ca_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2]

## soil characteristics
write.csv(cbind(as.matrix(anova(retention_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(pH_lm)[, c(4, 5)]), 
                as.matrix(anova(CEC_lm)[, c(4, 5)])),
          'analyses/tables/soil_char.csv')

#### save graphs ####

## allometry
ggsave("analyses/plots/plots_allometry.jpeg", plot = plots_allometry,
       width = 28, height = 18, units = "cm", dpi = 600) # 4 panels

## foliar organics
ggsave("analyses/plots/plots_foliar_organics.jpeg", plot = plots_foliar_organics,
       width = 42, height = 15, units = "cm", dpi = 600) # 3 panels

## foliar inorganics
ggsave("analyses/plots/plots_foliar_inorganics.jpeg", plot = plots_foliar_inorganics,
       width = 45, height = 25, units = "cm", dpi = 600) # 6 panels

## foliar isotopes
ggsave("analyses/plots/plots_foliar_isotopes.jpeg", plot = plots_foliar_isotopes,
       width = 29, height = 12, units = "cm", dpi = 600) # 2 panels

## soil organics
ggsave("analyses/plots/plots_soil_organics.jpeg", plot = plots_soil_organics,
       width = 42, height = 15, units = "cm", dpi = 600) # 3 panels

## soil inorganics
ggsave("analyses/plots/plots_soil_inorganics.jpeg", plot = plots_soil_inorganics,
       width = 45, height = 25, units = "cm") # 6 panels

## soil characteristics
ggsave("analyses/plots/plots_soil_characteristics.jpeg", plot = plots_soil_characteristics,
       width = 42, height = 15, units = "cm", dpi = 600) # 3 panels


