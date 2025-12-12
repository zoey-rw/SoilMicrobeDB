## Gorka et al. 2023
## "Beyond PLFA: Concurrent extraction of neutral and 
##  glycolipid fatty acids provides new insights into 
##  soil microbial communities"

# This script analyses the soil fatty acid profile data.
# It creates Fig. 5 and Fig. S4

# Load packages
library(tidyverse)
library(vegan)
library(ggpubr)
#

# Prepare Data ---------------------------------------------------------------
# Import FA data
fa <- read_csv("nmolC.csv") %>% 
  dplyr::select(-info) %>% 
  pivot_longer(cols = '10Me17:0':ncol(.), 
               values_to = "nmolC", 
               names_to = "FA")

# Save FA list
#fa$FA %>% as_factor %>% levels %>% as_tibble %>% write_csv("FA_list.csv")
#--> re-import as 'phylum.csv' with FA phylum assignments

# Import phylum assignments
phylum <- read_csv("phylum.csv")

# Check if FA names correspond
levels(as.factor(phylum$FA))
levels(as.factor(fa$FA))
all(levels(as.factor(phylum$FA))==levels(as.factor(fa$FA)))
#--> look OK!

# Join FA abundance data with phylum specifications
fa <- fa %>% 
  left_join(phylum)

# Calculate mol% of total
fa$molPercent <- fa$nmolC / ave(fa$nmolC, fa$full_id, FUN = sum)
# # Check it! --> sum should be 1
aggregate(molPercent ~ full_id, data = fa, FUN = sum)

# Remove FAs <0.5 mol%
fa <- 
  fa %>% 
  filter(molPercent > 0.005) %>% 
  dplyr::select(-molPercent)

# Create soil_type column
fa <- fa %>% 
  mutate(soil_type = str_replace(soil_type, "A", "Agriculture"),
         soil_type = str_replace(soil_type, "F", "Forest"),
         soil_type = str_replace(soil_type, "G", "Grassland"))

# Resort fa_type factor levels
fa <- fa %>% 
  mutate(fa_type = fct_relevel(fa_type, "PLFA", "NLFA", "GLFA"))

# Make 'id' column
fa <- fa %>% 
  separate(full_id, 
           into = c("a", "b", "c"), 
           sep = "_", 
           remove = FALSE) %>% 
  mutate(id = paste0(b, "_", c)) %>% 
  dplyr::select(full_id, id, fa_type, soil_type, phylum, FA, nmolC)
  

# Print data
fa



#
# Boxplot + NMDS (-18:0) -- Fig. 5 --------------------------------------
# Delete FA 18:0
data_fig <- fa %>% filter(FA != "C18:0")

# Beautify FA names
data_fig <- data_fig %>%
  mutate(FA = str_replace(FA, "w", "ω"),
         FA = str_replace(FA, "C", ""))

### 1. Barchart 
fig_a <- data_fig %>% 
  group_by(full_id, fa_type, soil_type, phylum) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  mutate(phylum = fct_relevel(phylum, "General", "Actino", "Gram -", "Gram +", "Fungi", "unspecified")) %>% 
  ggplot(aes(x = phylum, y = nmolC, fill = fa_type, group = soil_type)) +
  stat_summary(fun = mean, geom = "bar", width = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0) +
  scale_x_discrete(labels = c("Gen.", "Actino", "Gram⁻\n(AMF)", "Gram⁺", "Fungi", "unsp.")) +
  facet_grid(soil_type ~ fa_type, scales = "free") +
  labs(x = "", y = "nmol C g⁻¹ dw") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 7),
        text = element_text(family = "Arial"))
fig_a




### 2. NMDS
# Pivot data
data_nmds <- data_fig %>% 
  dplyr::select(-phylum) %>% 
  pivot_wider(names_from = FA,
              values_from = nmolC,
              values_fill = list(nmolC = 0)) # Set NAs zero

# Subset data for NMDS
names(data_nmds)
env <- data_nmds[, 1:4]
species <- data_nmds[, 5:ncol(data_nmds)]

# Run model
set.seed(123)
mod <- metaMDS(comm = species)
mod

# Save coordinates
site_scores <- as.data.frame(scores(mod, "sites"))
species_scores <- as.data.frame(scores(mod)$species)
species_scores$FA <- rownames(species_scores)

# Add group vars
site_scores$soil_type <- env$soil_type
site_scores$fa_type <- env$fa_type

## Plot it!
# Calculate mean and standard error per strain
site_means <- site_scores %>% 
  group_by(soil_type, fa_type) %>% 
  summarise(mean_NMDS1 = mean(NMDS1),
            mean_NMDS2 = mean(NMDS2),
            se_NMDS1 = sd(NMDS1)/sqrt(n()),
            se_NMDS2 = sd(NMDS2)/sqrt(n()))

# Calculate upper and lower bounds of error bars
site_means$error_upper_NMDS1 <- site_means$mean_NMDS1 + site_means$se_NMDS1
site_means$error_lower_NMDS1 <- site_means$mean_NMDS1 - site_means$se_NMDS1
site_means$error_upper_NMDS2 <- site_means$mean_NMDS2 + site_means$se_NMDS2
site_means$error_lower_NMDS2 <- site_means$mean_NMDS2 - site_means$se_NMDS2

# ggplot
fig_b <-
  ggplot() +
  geom_text(data = species_scores, size = 2, alpha = 0.5,
            aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_errorbar(data = site_means, 
                aes(x = mean_NMDS1, 
                    ymin = error_lower_NMDS2, ymax = error_upper_NMDS2,
                    colour = fa_type)) +
  geom_errorbarh(data = site_means, 
                 aes(y = mean_NMDS2, 
                     xmin = error_lower_NMDS1, xmax = error_upper_NMDS1,
                     colour = fa_type)) +
  geom_point(data = site_means, size = 3,
             aes(x = mean_NMDS1, y = mean_NMDS2, 
                 colour = fa_type, shape = soil_type)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Arial"))
fig_b

ggarrange(fig_a, fig_b, nrow = 2, labels = c("a", "b"))
ggsave("Fig_5.jpg", units = "cm", width = 20, height = 22, bg = "white")


# PERMANOVA
adonis2(species ~ fa_type + soil_type, data = env)

# Print stress value
mod$stress






#
# Boxplot + NMDS (+18:0) -- Fig. S4 --------------------------------------
# Keep FA 18:0
data_fig <- fa

# Beautify FA names
data_fig <- data_fig %>%
  mutate(FA = str_replace(FA, "w", "ω"),
         FA = str_replace(FA, "C", ""))

### 1. Barchart 
fig_a <- data_fig %>% 
  group_by(full_id, fa_type, soil_type, phylum) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  mutate(phylum = fct_relevel(phylum, "General", "Actino", "Gram -", "Gram +", "Fungi", "unspecified")) %>% 
  ggplot(aes(x = phylum, y = nmolC, fill = fa_type, group = soil_type)) +
  stat_summary(fun = mean, geom = "bar", width = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0) +
  scale_x_discrete(labels = c("Gen.", "Actino", "Gram⁻\n(AMF)", "Gram⁺", "Fungi", "unsp.")) +
  facet_grid(soil_type ~ fa_type, scales = "free") +
  labs(x = "", y = "nmol C g⁻¹ dw") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 7),
        text = element_text(family = "Arial"))
fig_a




### 2. NMDS
# Pivot data
data_nmds <- data_fig %>% 
  dplyr::select(-phylum) %>% 
  pivot_wider(names_from = FA,
              values_from = nmolC,
              values_fill = list(nmolC = 0)) # Set NAs zero

# Subset data for NMDS
names(data_nmds)
env <- data_nmds[, 1:4]
species <- data_nmds[, 5:ncol(data_nmds)]

# Run model
mod <- metaMDS(comm = species)
mod

# Save coordinates
site_scores <- as.data.frame(scores(mod, "sites"))
species_scores <- as.data.frame(scores(mod)$species)
species_scores$FA <- rownames(species_scores)

# Add group vars
site_scores$soil_type <- env$soil_type
site_scores$fa_type <- env$fa_type

## Plot it!
# Calculate mean and standard error per strain
site_means <- site_scores %>% 
  group_by(soil_type, fa_type) %>% 
  summarise(mean_NMDS1 = mean(NMDS1),
            mean_NMDS2 = mean(NMDS2),
            se_NMDS1 = sd(NMDS1)/sqrt(n()),
            se_NMDS2 = sd(NMDS2)/sqrt(n()))

# Calculate upper and lower bounds of error bars
site_means$error_upper_NMDS1 <- site_means$mean_NMDS1 + site_means$se_NMDS1
site_means$error_lower_NMDS1 <- site_means$mean_NMDS1 - site_means$se_NMDS1
site_means$error_upper_NMDS2 <- site_means$mean_NMDS2 + site_means$se_NMDS2
site_means$error_lower_NMDS2 <- site_means$mean_NMDS2 - site_means$se_NMDS2

# ggplot
fig_b <-
  ggplot() +
  geom_text(data = species_scores, size = 2, alpha = 0.5,
            aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_errorbar(data = site_means, 
                aes(x = mean_NMDS1, 
                    ymin = error_lower_NMDS2, ymax = error_upper_NMDS2,
                    colour = fa_type)) +
  geom_errorbarh(data = site_means, 
                 aes(y = mean_NMDS2, 
                     xmin = error_lower_NMDS1, xmax = error_upper_NMDS1,
                     colour = fa_type)) +
  geom_point(data = site_means, size = 3,
             aes(x = mean_NMDS1, y = mean_NMDS2, 
                 colour = fa_type, shape = soil_type)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Arial"))
fig_b

ggarrange(fig_a, fig_b, nrow = 2, labels = c("a", "b"))
ggsave("Fig_S4.jpg", units = "cm", width = 20, height = 22, bg = "white")


# PERMANOVA
adonis2(species ~ fa_type + soil_type, data = env)

# Print stress value
mod$stress




#
# Summarise data in tables ------------------------------------------------------------------
# Table percent of total
fa_tot <-
  fa %>% 
  filter(FA != "C18:0") %>% 
  group_by(id, soil_type) %>% 
  summarise(nmolC_tot = sum(nmolC)) %>% 
  ungroup()

# ...by lipid fraction
fa %>% 
  filter(FA != "C18:0") %>% 
  group_by(id, fa_type, soil_type) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  left_join(fa_tot) %>% 
  mutate(perc_tot = nmolC / nmolC_tot * 100) %>% 
  group_by(fa_type, soil_type) %>% 
  summarise(perc_tot_mean = mean(perc_tot),
            perc_tot_se = sd(perc_tot) / sqrt(n()))

# ...by phylum
fa %>% 
  filter(FA != "C18:0") %>% 
  group_by(id, phylum, soil_type) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  left_join(fa_tot) %>% 
  mutate(perc_tot = nmolC / nmolC_tot * 100) %>% 
  group_by(phylum, soil_type) %>% 
  summarise(perc_tot_mean = mean(perc_tot),
            perc_tot_se = sd(perc_tot) / sqrt(n()))


# Table mean absolute values
fa %>% 
  filter(FA != "C18:0") %>% 
  group_by(id, fa_type, soil_type) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  group_by(fa_type, soil_type) %>% 
  summarise(nmolC_mean = mean(nmolC),
            nmolC_se = sd(nmolC) / sqrt(n()))

### With 18:0
fa_tot <-
  fa %>% 
  group_by(id, soil_type) %>% 
  summarise(nmolC_tot = sum(nmolC)) %>% 
  ungroup()

fa %>% 
  group_by(id, fa_type, soil_type) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  left_join(fa_tot) %>% 
  mutate(perc_tot = nmolC / nmolC_tot * 100) %>% 
  group_by(fa_type, soil_type) %>%
  summarise(mean_se(perc_tot)) %>% 
  ungroup() %>% 
  mutate(perc_tot_mean = y,
         perc_tot_se = ymax - y) %>% 
  select(-y, -ymin, -ymax)


# mol% of 18:0
fa %>% 
  left_join(fa_tot) %>% 
  mutate(perc_tot = nmolC / nmolC_tot * 100) %>% 
  filter(FA == "C18:0" | FA == "C18:2w6,9") %>% 
  group_by(FA, fa_type, soil_type) %>%
  summarise(mean_molperc = mean(perc_tot) %>% round(1))

fa %>% 
  left_join(fa_tot) %>% 
  mutate(perc_tot = nmolC / nmolC_tot * 100) %>% 
  group_by(FA, fa_type, soil_type) %>%
  summarise(mean_molperc = mean(perc_tot) %>% round(1)) %>% 
  filter(fa_type == "PLFA") %>% 
  arrange(mean_molperc %>% desc())


#
