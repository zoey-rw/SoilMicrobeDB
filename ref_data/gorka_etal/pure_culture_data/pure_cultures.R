## Gorka et al. 2023
## "Beyond PLFA: Concurrent extraction of neutral and 
##  glycolipid fatty acids provides new insights into 
##  soil microbial communities"

# This script analyses the microbial pure culture fatty acid profile data.
# It creates Fig. 3, Fig. 4, Fig. S2, and Fig. S3.

# Load packages
library(tidyverse)
library(ggrepel) # Additional geoms for ggplot2
library(viridis) # Viridis colour package
library(vegan) # Multivariate statistics
library(patchwork) # To arrange ggplots
library(rstatix) # Pipe-friendly statistics
#

# Prep data ---------------------------------------------------------------
# Import TIC area data
fa <- 
  read_csv("area.csv") %>% 
  dplyr::select(-info)

# Check data for duplicate values
fa %>% 
  group_by(full_id, FA) %>% 
  filter(n() > 1) %>% 
  nrow
# --> looks ok!

# Export list of FAs for C-atom list
#write.csv(sort(levels(as_factor(fa$FA))), file = "FA_list.csv")
# --> save as "C_atoms.csv", rename column "x" to "FA", add C.atoms

# Save internal standard extra
C19_0 <-
  fa %>% 
  filter(FA == "C19:0") %>% 
  dplyr::select(full_id, area) %>% 
  mutate(C19_0 = area / 0.639938565897674) %>% # *
  dplyr::select(-area)
# *The long value is the added amount 
#  of 19:0 in umol C (1:200, 100 ul)

# Add 19:0 column to "fa"
fa <-
  fa %>% 
  filter(FA != "C19:0") %>% # Remove FA 19:0 first
  left_join(C19_0, by = "full_id")



##### Calculations
## nmol C = area / C19_0-value * 1000
fa <-
  fa %>% 
  mutate(nmolC = area / C19_0 * 1000) %>% 
  dplyr::select(-C19_0, -area)

## Load weights data
weights <- 
  read_csv("weights.csv")

# Calc. actual and potential volumes for transfer loss correction
weights <- 
  weights %>% 
  mutate(actual_vol_g = full_vial_g - empty_vial_g, # Calc actual vol
         potential_vol_g = 5.60334210526316) %>%    # Add potential vol
  dplyr::select(-c(full_vial_g, empty_vial_g))

# Check if data corresponds correctly
all(levels(as.factor(fa$id)) == levels(as.factor(weights$id)))
all(levels(as.factor(fa$rep)) == levels(as.factor(weights$rep)))
all(levels(as.factor(paste(fa$id, fa$rep))) == levels(as.factor(paste(weights$id, weights$rep))))

# Correct for losses with transfer of lower phase after phase separation
fa <-
  fa %>% 
  left_join(weights, by = c("id", "rep", "GC_run")) %>% 
  mutate(nmolC = nmolC / actual_vol_g * potential_vol_g) %>% # CALC
  dplyr::select(-actual_vol_g, -potential_vol_g)

## Subtr Blks
blks <-
  fa %>% 
  filter(id == "Blk") %>%
  dplyr::select(fa_type, GC_run, FA, nmolC) %>% 
  group_by(fa_type, GC_run, FA) %>% 
  summarise(blk_mean = mean(nmolC)) %>% # Calc mean
  ungroup()

fa <-
  fa %>% 
  filter(id != "Blk") %>% 
  left_join(blks, by = c("fa_type", "GC_run", "FA")) %>%
  mutate(blk_mean = ifelse(is.na(blk_mean), 0, blk_mean)) %>% # Set NAs zero
  mutate(nmolC = nmolC - blk_mean) %>% # Subtract blks
  dplyr::select(-blk_mean)


## Methyl-group correction and divide by (g dw)
# Load phylum data (includes C-atoms)
#--> NOTE that the 'phylum' column refers 
#    to the phylogenetic specificity of FAs
c_atoms <- read_csv("C_atoms.csv")


# Merge dfs and Calc
all(c_atoms$FA %>% as_factor %>% levels %>% sort == fa$FA %>% as_factor %>% levels %>% sort) # Check if FA names correspond
fa <-
  fa %>% 
  left_join(c_atoms, by = "FA") %>% 
  mutate(nmolC = nmolC / (C.atoms + 1) * C.atoms ) %>% # Methyl group correction
  mutate(nmolC = nmolC / dw_mg) %>%                    # Normalise by mg dry weight
  dplyr::select(-c(C.atoms, dw_mg))

# Set negative values to zero
# --> applies to straight chain FAs which were present in blank samples
fa <-
  fa %>% 
  mutate(nmolC = ifelse(nmolC < 0, 0, nmolC)) %>% 
  droplevels()

# Reorder factor levels (PLFA > NLFA > GLFA)
fa <-
  fa %>%  
  mutate(fa_type = fct_rev(fa_type))

# Check data for duplicates
fa %>% group_by(full_id, FA) %>% filter(n() > 1) %>% nrow()
# --> looks ok!

# Add taxon column
taxon <- read_csv("species_list.csv")
all(fa$species %>% as_factor %>% levels %>% sort == taxon$species %>% as_factor %>% levels %>% sort) # Check if species names correspond
fa <- fa %>% left_join(taxon)

# Print data
fa



#
# Heatmap + Boxplot (-18:0) - Fig 3 -----------------------------------------------
# Save data for section, and delete FA 18:0
heat <- fa %>% 
  filter(FA != "C18:0")

# Calculate mol%
fa_sum <-
  heat %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  group_by(full_id) %>% 
  summarise(nmolC_sum = sum(nmolC)) %>% 
  ungroup()

heat <-
  heat %>% 
  left_join(y = fa_sum, by = "full_id") %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  mutate(molPerc = nmolC / nmolC_sum * 100) %>% 
  dplyr::select(-c(nmolC_sum))

# Check it
#--> Sum should always be 100
heat %>% 
  group_by(full_id) %>% 
  summarise(molPerc = sum(molPerc)) %>% 
  arrange(molPerc) %>% 
  print(n = 114)

# Beautify fatty acid names
heat <- 
  heat %>% 
  mutate(FA = gsub("C", "", FA),
         FA = gsub("w", "ω", FA),
         FA = gsub("_", "-", FA),
         FA = gsub("-OH", "OH", FA))

# Add genus column
heat <-
  heat %>% 
  separate(col = species,
           into = c("genus", "x"),
           sep = " ",
           remove = FALSE) %>% 
  dplyr::select(-x)

# Resort factor levels
levels(as_factor(paste(heat$taxon, heat$species))) # Print all species
heat <- heat %>% 
  mutate(species = fct_relevel(species,
                               c("Solirubrobacter soli", "Streptosporangium roseum",    # Actinobacteria
                                 "Paenibacillus alginolyticus",                         # Firmicutes
                                 "Labilithrix luteola", "Chelativorans multitrophicus", # Proteobacteria
                                 "Paraburkholderia xenovorans",
                                 "Rhizophagus irregularis",                             # Glomeromycetes
                                 "Periconia macrospinosa", "Trichoderma virens",        # Ascomycetes
                                 "Armillaria gallica", "Psilocybe cyanescens",          # Basidiomycetes
                                 "Lactarius subdulcis", "Lactarius quietus")))          # Basidiomycetes (EMF)
levels(heat$species) # Check it!

levels(as_factor(heat$taxon))
heat <- heat %>% 
  mutate(taxon = fct_relevel(taxon,
                               c("Actinobacteria", "Firmicutes", "Proteobacteria",
                                 "Glomeromycota", "Ascomycota", "Basidiomycota")))

#### Plot it!
# Heatmap
fig1 <-
heat %>% 
  dplyr::select(-nmolC) %>% 
  group_by(FA, fa_type, species) %>% 
  mutate(species = as_factor(species)) %>% 
  summarise(molPerc = mean(molPerc)) %>% 
  ungroup() %>% 
  pivot_wider(values_from = molPerc, names_from = FA) %>%
  pivot_longer(cols = c("10Me17:0":ncol(.)), values_to = "molPerc", names_to = "FA") %>% 
  mutate(molPerc = replace_na(molPerc, 0)) %>% # Replace NA values with zero
  droplevels() %>% 
  # Plot it!
  ggplot(aes(x = FA, y = species, fill = molPerc)) +
  geom_tile() +
  scale_fill_viridis(name = "Mol %") +
  labs(x = "Fatty acid", y = "Species") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~fa_type, ncol = 1, strip.position = "right") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top",
        strip.text.y = element_blank(),
        text = element_text(family = "Arial"),
        axis.text.y = element_text(face = "italic"))

# nmolC -- boxplot
fig2 <-
heat %>% 
  dplyr::select(-molPerc) %>% 
  # Calc sum of all FAs
  group_by(fa_type, species, genus, taxon, rep) %>% 
  mutate(species = as_factor(species)) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  # Plot it!
  ggplot(aes(x = nmolC, y = species, fill = taxon)) +
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  labs(x = "nmol C g⁻¹ dw", y = "Species") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~fa_type, ncol = 1, strip.position = "right") +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family = "Arial"))

# Plot and export it!
fig_JPEG <- fig1 + fig2 + plot_layout(widths = c(1.4,1))
jpeg(filename = "Fig_3.jpg", units = "cm", width = 35, height = 19, res = 300)
plot(fig_JPEG)
dev.off()


#
# Heatmap + Boxplot (+18:0) - Fig S2 -----------------------------------------------
# Save df for section
heat <- fa 

# Calculate mol%
fa_sum <-
  heat %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  group_by(full_id) %>% 
  summarise(nmolC_sum = sum(nmolC)) %>% 
  ungroup()

heat <-
  heat %>% 
  left_join(y = fa_sum, by = "full_id") %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  mutate(molPerc = nmolC / nmolC_sum * 100) %>% 
  dplyr::select(-c(nmolC_sum))

# Check it
#--> Sum should always be 100
heat %>% 
  group_by(full_id) %>% 
  summarise(molPerc = sum(molPerc)) %>% 
  arrange(molPerc)

# Beautify fatty acid :)
heat <- heat %>% 
  mutate(FA = gsub("C", "", FA),
         FA = gsub("w", "ω", FA),
         FA = gsub("_", "-", FA),
         FA = gsub("-OH", "OH", FA))

# Add genus column
heat <-
  heat %>% 
  separate(col = species,
           into = c("genus", "x"),
           sep = " ",
           remove = FALSE) %>% 
  dplyr::select(-x)

# Resort factor levels
levels(as_factor(paste(heat$taxon, heat$species)))
heat <- heat %>% 
  mutate(species = fct_relevel(species,
                               c("Solirubrobacter soli", "Streptosporangium roseum",    # Actinobacteria
                                 "Paenibacillus alginolyticus",                         # Firmicutes
                                 "Labilithrix luteola", "Chelativorans multitrophicus", # Proteobacteria
                                 "Paraburkholderia xenovorans",
                                 "Rhizophagus irregularis",                             # Glomeromycetes
                                 "Periconia macrospinosa", "Trichoderma virens",        # Ascomycetes
                                 "Armillaria gallica", "Psilocybe cyanescens",          # Basidiomycetes
                                 "Lactarius subdulcis", "Lactarius quietus")))          # Basidiomycetes (EMF)
levels(heat$species)

levels(as_factor(heat$taxon))
heat <- heat %>% 
  mutate(taxon = fct_relevel(taxon,
                                   c("Actinobacteria", "Firmicutes", "Proteobacteria",
                                     "Glomeromycota", "Ascomycota", "Basidiomycota")))

#### Plot it!
# Heatmap
fig1 <-
  heat %>% 
  dplyr::select(-nmolC) %>% 
  group_by(FA, fa_type, species) %>% 
  mutate(species = as_factor(species)) %>% 
  summarise(molPerc = mean(molPerc)) %>% 
  ungroup() %>% 
  pivot_wider(values_from = molPerc, names_from = FA) %>%
  pivot_longer(cols = c("10Me17:0":ncol(.)), values_to = "molPerc", names_to = "FA") %>% 
  mutate(molPerc = replace_na(molPerc, 0)) %>% # Replace NA values with zero
  droplevels() %>% 
  # Plot it!
  ggplot(aes(x = FA, y = species, fill = molPerc)) +
  geom_tile() +
  scale_fill_viridis(name = "Mol %") +
  labs(x = "Fatty acid", y = "Species") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~fa_type, ncol = 1, strip.position = "right") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top",
        strip.text.y = element_blank(),
        text = element_text(family = "Arial"),
        axis.text.y = element_text(face = "italic"))

# nmolC -- boxplot
fig2 <-
  heat %>% 
  dplyr::select(-molPerc) %>% 
  # Calc sum of all FAs
  group_by(fa_type, species, genus, taxon, rep) %>% 
  mutate(species = as_factor(species)) %>% 
  summarise(nmolC = sum(nmolC)) %>% 
  ungroup() %>% 
  # Plot it!
  ggplot(aes(x = nmolC, y = species, fill = taxon)) +
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  labs(x = "nmol C g⁻¹ dw", y = "Species") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~fa_type, ncol = 1, strip.position = "right") +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family = "Arial"))

# Plot it!
fig_JPEG <- fig1 + fig2 + plot_layout(widths = c(1.4,1))
jpeg(filename = "Fig_S2.jpg", units = "cm", width = 35, height = 19, res = 300)
plot(fig_JPEG)
dev.off()



#
# NMDS (-18:0) - Fig 4 ------------------------------------------------------
# Save df for section
data_nmds <- fa %>% 
  filter(FA != "C18:0") 

# Calculate mol%
fa_sum <-
  data_nmds %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  group_by(full_id) %>% 
  summarise(nmolC_sum = sum(nmolC)) %>% 
  ungroup()

data_nmds <-
  data_nmds %>% 
  left_join(y = fa_sum, by = "full_id") %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  mutate(molPerc = nmolC / nmolC_sum * 100) %>% 
  dplyr::select(-c(nmolC_sum))

# Check it
#--> Sum should always be 100
data_nmds %>% 
  group_by(full_id) %>% 
  summarise(molPerc = sum(molPerc)) %>% 
  arrange(molPerc) %>% 
  print(n = 114)

# Beautify fatty acid :)
data_nmds <- data_nmds %>% 
  mutate(FA = gsub("C", "", FA),
         FA = gsub("w", "ω", FA),
         FA = gsub("_", "-", FA),
         FA = gsub("-OH", "OH", FA))

# Check data for rowsum = 0
to_delete <- 
  data_nmds %>% 
  group_by(full_id, id, fa_type, species, FA) %>% 
  summarise(sum = sum(molPerc)) %>% 
  ungroup() %>% 
  filter(sum == 0) %>% 
  mutate(d1 = paste(full_id, FA),
         d2 = "DELETE") %>% 
  dplyr::select(d1, d2)

# Remove data with rowsum = 0
data_nmds <- 
  data_nmds %>% 
  mutate(d1 = paste(full_id, FA)) %>% 
  left_join(to_delete) %>% 
  replace_na(list(d2 = "keep")) %>% 
  filter(d2 != "DELETE") %>% 
  dplyr::select(-d1, -d2)


## Spread data
names(data_nmds)
data_nmds <- data_nmds %>% 
  dplyr::select(-nmolC, -GC_run) %>% 
  pivot_wider(
    names_from = FA,
    values_from = molPerc, 
    values_fill = list(molPerc = 0)
  )
data_nmds

# Subset data for NMDS
data_nmds %>% names
group_vars <- data_nmds[, 1:7]
data_nmds <- data_nmds[, 8:ncol(data_nmds)]

# Set seed
set.seed(123)

# Run model
mod <- metaMDS(comm = data_nmds, distance = "bray", k = 2)
mod

# Save coordinates
site_scores <- as.data.frame(scores(mod)$sites)
species_scores <- as.data.frame(scores(mod)$species)
species_scores$FA <- rownames(species_scores)

# Add group vars
site_scores$fa_type <- group_vars$fa_type
site_scores$taxon <- group_vars$taxon
site_scores$species <- group_vars$species

## Plot it!
# Calculate mean and standard error per strain
site_means <- site_scores %>% 
  group_by(taxon, species, fa_type) %>% 
  summarise(mean_NMDS1 = mean(NMDS1),
            mean_NMDS2 = mean(NMDS2),
            se_NMDS1 = sd(NMDS1)/sqrt(n()),
            se_NMDS2 = sd(NMDS2)/sqrt(n()))

# Calculate upper and lower bounds of error bars
site_means$error_upper_NMDS1 <- site_means$mean_NMDS1 + site_means$se_NMDS1
site_means$error_lower_NMDS1 <- site_means$mean_NMDS1 - site_means$se_NMDS1
site_means$error_upper_NMDS2 <- site_means$mean_NMDS2 + site_means$se_NMDS2
site_means$error_lower_NMDS2 <- site_means$mean_NMDS2 - site_means$se_NMDS2

# Resort phyla factor levels
site_means <- 
  site_means %>% 
  ungroup() %>% 
  mutate(taxon = fct_relevel(taxon,
                             c("Actinobacteria", "Firmicutes", "Proteobacteria",
                               "Glomeromycota", "Ascomycota", "Basidiomycota")))


# Easy plot without geom_repel
ggplot() +
  geom_text_repel(data = species_scores, size = 2, alpha = 0.5, 
            aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_errorbar(data = site_means, 
                aes(x = mean_NMDS1, 
                    ymin = error_lower_NMDS2, ymax = error_upper_NMDS2,
                    colour = taxon)) +
  geom_errorbarh(data = site_means, 
                 aes(y = mean_NMDS2, 
                     xmin = error_lower_NMDS1, xmax = error_upper_NMDS1,
                     colour = taxon)) +
  geom_point(data = site_means, size = 3,
             aes(x = mean_NMDS1, y = mean_NMDS2, 
                 colour = taxon, shape = fa_type)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Arial"))


#---> Many species scores (i.e., fatty acid names) are overlapping!

# Filter overlapping species scores, and group by cluster
species_scores_overlap1 <-
  species_scores %>% 
  filter(FA == "19:0cy(ω9)")

species_scores_overlap2 <-
  species_scores %>% 
  filter(FA == "10Me17:0" | FA == "10Me18:0" | FA == "17:0-Me")

species_scores_overlap3 <-
  species_scores %>% 
  filter(FA == "16:0i" | FA == "18:0i")

species_scores_overlap4 <-
  species_scores %>% 
  filter(FA == "16:2ω4,7" | FA == "16:2ω6,9")

species_scores_overlap5 <-
  species_scores %>% 
  filter(FA == "15:0" | FA == "19:1ω9")

species_scores_overlap6 <-
  species_scores %>% 
  filter(FA == "18:3ω6,9,12" | FA == "16:1ω5")

species_scores_overlap7 <-
  species_scores %>% 
  filter(FA == "19:2ω6,9")


species_scores_keep <-
species_scores %>% 
  filter(FA != "19:0cy(ω9)" & 
           FA != "10Me17:0" & FA != "10Me18:0" & FA != "17:0-Me" &
           FA != "16:0i" & FA != "18:0i" & 
           FA != "16:2ω4,7" & FA != "16:2ω6,9" & 
           FA != "15:0" & FA != "19:1ω9" & 
           FA != "18:3ω6,9,12" & FA != "16:1ω5" & 
           FA != "19:2ω6,9")

# ggplot
ggplot() +
  geom_text_repel(data = species_scores_overlap1, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.2, nudge_y = 0.1, nudge_x = -0.1,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text_repel(data = species_scores_overlap2, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.4, nudge_y = 0.1, nudge_x = 0.2,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text_repel(data = species_scores_overlap3, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.4, nudge_y = -0.1,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text_repel(data = species_scores_overlap4, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.3, nudge_x = 0.2, nudge_y = 0.05,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text_repel(data = species_scores_overlap5, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.4, nudge_x = 0.1, nudge_y = -0.1,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text_repel(data = species_scores_overlap6, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.4, nudge_y = 0.1, nudge_x = 0.1,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text_repel(data = species_scores_overlap7, size = 2, alpha = 0.5, 
                  segment.colour = "black", min.segment.length = 0, force = 0.2,
                  segment.size = 0.1,
                  box.padding = 0.5, nudge_y = -0.1,
                  aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_text(data = species_scores_keep, size = 2, alpha = 0.5, 
            aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_errorbar(data = site_means, 
                aes(x = mean_NMDS1, 
                    ymin = error_lower_NMDS2, ymax = error_upper_NMDS2,
                    colour = taxon)) +
  geom_errorbarh(data = site_means, 
                 aes(y = mean_NMDS2, 
                     xmin = error_lower_NMDS1, xmax = error_upper_NMDS1,
                     colour = taxon)) +
  geom_point(data = site_means, size = 3,
             aes(x = mean_NMDS1, y = mean_NMDS2, 
                 colour = taxon, shape = fa_type)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Arial"))


# Export plot
ggsave("Fig_4.jpg", units = "cm", width = 18, height = 12)

# PERMANOVA
adonis2(data_nmds ~ fa_type + species, data = group_vars, method = "bray")

# Print stress value
mod$stress %>% round(3)

#

# NMDS (+18:0) - Fig S3 ------------------------------------------------------
# Save df for section
data_nmds <- fa

# Calculate mol%
fa_sum <-
  data_nmds %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  group_by(full_id) %>% 
  summarise(nmolC_sum = sum(nmolC)) %>% 
  ungroup()

data_nmds <-
  data_nmds %>% 
  left_join(y = fa_sum, by = "full_id") %>% 
  mutate(nmolC = ifelse(is.na(nmolC), 0, nmolC)) %>% # Set NAs zero
  mutate(molPerc = nmolC / nmolC_sum * 100) %>% 
  dplyr::select(-c(nmolC_sum))

# Check it
#--> Sum should always be 100
data_nmds %>% 
  group_by(full_id) %>% 
  summarise(molPerc = sum(molPerc)) %>% 
  arrange(molPerc)

# Beautify fatty acid :)
data_nmds <- data_nmds %>% 
  mutate(FA = gsub("C", "", FA),
         FA = gsub("w", "ω", FA),
         FA = gsub("_", "-", FA),
         FA = gsub("-OH", "OH", FA))

# Check data for rowsum = 0
to_delete <- 
  data_nmds %>% 
  group_by(full_id, id, fa_type, species, FA) %>% 
  summarise(sum = sum(molPerc)) %>% 
  ungroup() %>% 
  filter(sum == 0) %>% 
  mutate(d1 = paste(full_id, FA),
         d2 = "DELETE") %>% 
  dplyr::select(d1, d2)

# Remove data with rowsum = 0
data_nmds <- 
  data_nmds %>% 
  mutate(d1 = paste(full_id, FA)) %>% 
  left_join(to_delete) %>% 
  replace_na(list(d2 = "keep")) %>% 
  filter(d2 != "DELETE") %>% 
  dplyr::select(-d1, -d2)


## Spread data
names(data_nmds)
data_nmds <- data_nmds %>% 
  dplyr::select(-nmolC, -GC_run) %>% 
  pivot_wider(
    names_from = FA,
    values_from = molPerc, 
    values_fill = list(molPerc = 0)
  )
data_nmds

# Subset data for NMDS
data_nmds %>% names
group_vars <- data_nmds[, 1:7]
data_nmds <- data_nmds[, 8:ncol(data_nmds)]

# Set seed
set.seed(123)

# Run model
mod <- metaMDS(comm = data_nmds, distance = "bray", k = 2)
mod

# Save coordinates
site_scores <- as.data.frame(scores(mod)$sites)
species_scores <- as.data.frame(scores(mod)$species)
species_scores$FA <- rownames(species_scores)

# Add group vars
site_scores$fa_type <- group_vars$fa_type
site_scores$taxon <- group_vars$taxon
site_scores$species <- group_vars$species

## Plot it!
# Calculate mean and standard error per strain
site_means <- site_scores %>% 
  group_by(taxon, species, fa_type) %>% 
  summarise(mean_NMDS1 = mean(NMDS1),
            mean_NMDS2 = mean(NMDS2),
            se_NMDS1 = sd(NMDS1)/sqrt(n()),
            se_NMDS2 = sd(NMDS2)/sqrt(n()))

# Calculate upper and lower bounds of error bars
site_means$error_upper_NMDS1 <- site_means$mean_NMDS1 + site_means$se_NMDS1
site_means$error_lower_NMDS1 <- site_means$mean_NMDS1 - site_means$se_NMDS1
site_means$error_upper_NMDS2 <- site_means$mean_NMDS2 + site_means$se_NMDS2
site_means$error_lower_NMDS2 <- site_means$mean_NMDS2 - site_means$se_NMDS2

# Resort phyla factor levels
site_means <- 
  site_means %>% 
  ungroup() %>% 
  mutate(taxon = fct_relevel(taxon,
                             c("Actinobacteria", "Firmicutes", "Proteobacteria",
                               "Glomeromycota", "Ascomycota", "Basidiomycota")))

# ggplot
ggplot() +
  geom_text(data = species_scores, size = 2, alpha = 0.5,
            aes(x = NMDS1, y = NMDS2, label = FA)) +
  geom_errorbar(data = site_means, 
                aes(x = mean_NMDS1, 
                    ymin = error_lower_NMDS2, ymax = error_upper_NMDS2,
                    colour = taxon)) +
  geom_errorbarh(data = site_means, 
                 aes(y = mean_NMDS2, 
                     xmin = error_lower_NMDS1, xmax = error_upper_NMDS1,
                     colour = taxon)) +
  geom_point(data = site_means, size = 3,
             aes(x = mean_NMDS1, y = mean_NMDS2, 
                 colour = taxon, shape = fa_type)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        text = element_text(family = "Arial"))

# Export plot
ggsave("Fig_S3.jpg", units = "cm", width = 18, height = 12)

# PERMANOVA
adonis2(data_nmds ~ fa_type + species, data = group_vars, method = "bray")

# Print stress value
mod$stress %>% round(3)

#

