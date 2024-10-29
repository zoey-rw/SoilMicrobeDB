## Gorka et al. 2023
## "Beyond PLFA: Concurrent extraction of neutral and 
##  glycolipid fatty acids provides new insights into 
##  soil microbial communities"

# This script analyses the pure lipid solid phase fractionation data.
# It creates Fig. 2 and Fig. S1.

# Load packages
library(tidyverse)
library(ggh4x) # For additional ggplot functions
library(rstatix) # For pipe-friendly statistic functions
#

# Prepare Data ------------------------------------------------------------
area <- read_csv("area.csv") %>% select(-info)
lipid_class <- read_csv("LipidClass.csv")

# Calculate mean of blank samples
blk <- 
  area %>% 
  filter(type == "blk") %>%
  group_by(treatment, fraction, test, lipid) %>% 
  summarise(area_blk = mean(area))

# Calculate mean of external standards
# -> external standards were not fractionated on silica SPE
std <-
  area %>% 
  filter(type == "std") %>% 
  group_by(test, lipid) %>% 
  summarise(area_std = mean(area))

# Delete blank samples and external standards from area tibble
area <-
  area %>% 
  filter(type != "blk",
         type != "std")

# Put mean values of blanks and external standards into area tibble
area <-
area %>% 
  left_join(blk) %>% 
  replace_na(list(area_blk = 0)) %>% 
  left_join(std)

# Subtract blanks, and calculate recovery
rec <-
  area %>% 
  mutate(area = area - area_blk,
         recovery = area / area_std * 100) %>% 
  select(-c(area, area_blk, area_std)) %>% 
  left_join(lipid_class)



#
# Fig. 2 -------------------------------------------------------------
rec %>% 
  group_by(treatment, fraction, test, Type, rep) %>% 
  summarise(recovery = mean(recovery)) %>% # Calculate mean values of FAs from the same lipid
  mutate(test = paste(test, treatment)) %>% # Remove first test +Eth
  filter(test != "first with_ethanol") %>% 
  mutate(Type = str_replace(Type, "GL", "Glycolipid"),
         Type = str_replace(Type, "PL", "Phospholipids"),
         Type = as_factor(Type),
         Type = fct_relevel(Type, c("TAG", "Sterols", "Glycolipid", "Betaine lipid", "Phospholipids"))) %>% 
  mutate(treatment = as_factor(treatment),
         treatment = fct_recode(treatment,
                                "Classic eluents" = "old", 
                                "Adjusted eluents" = "with_ethanol",
                                "Acetone" = "acetone"),
         treatment = fct_relevel(treatment, c("Classic eluents", "Adjusted eluents", "Acetone"))) %>% 
  mutate(fraction = as_factor(fraction),
         fraction = fct_recode(fraction,
                               "Fraction 1 (NL)" = "1",
                               "Fraction 2 (GL)" = "2",
                               "Fraction 3 (PL)" = "3")) %>% 
  ggplot(aes(x = treatment, y = recovery, fill = treatment)) +
  geom_hline(yintercept = 100, linetype = 2, linewidth = 0.1) +
  stat_summary(geom = "bar", fun = mean, width = 0.8) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  facet_grid(fraction ~ Type, scales = "free") +
  labs(x =" ", y = "Recovery (%)") +
  #scale_fill_manual(values = c("grey20", "grey50", "grey80")) +
  scale_fill_manual(values = c("#12496d", "#1f78b4", "#a6cee2")) +
  #scale_fill_manual(values = c("#7865b3", "#b37865", "#65b378")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        text = element_text(family = "Arial"))
ggsave("Fig_2.jpg", units="cm", width=18, height=14)


# Summarise data in table
rec %>% 
  group_by(treatment, fraction, test, Type, rep) %>% 
  summarise(recovery = mean(recovery)) %>% 
  group_by(treatment, fraction, test, Type) %>% 
  summarise(recovery = round(mean(recovery), 2)) %>% 
  pivot_wider(names_from = c(fraction, treatment),
              values_from = recovery)



# Fig. S1 ------------------------------------------------------------------
# --> visualises if eluent volumes change recovery
f1a <- rec %>% filter(fraction == 1, test == "first") %>% mutate(eluent_vol = "2 ml")
f1b <- rec %>% filter(fraction == 1, test == "second") %>% mutate(eluent_vol = "1.5 ml")
f2a <- rec %>% filter(fraction == 2, test == "first") %>% mutate(eluent_vol = "2 ml")
f2b <- rec %>% filter(fraction == 2, test == "second") %>% mutate(eluent_vol = "1.5 ml")
f3a <- rec %>% filter(fraction == 3, test == "first") %>% mutate(eluent_vol = "0.5 ml")
f3b <- rec %>% filter(fraction == 3, test == "second") %>% mutate(eluent_vol = "1 ml")
rec_new <- rbind(f1a, f1b, f2a, f2b, f3a, f3b)

# Replace NAs with zero + complete cases
rec_new <- rec_new %>% 
  ungroup() %>% 
  filter(treatment == "with_ethanol") %>% 
  filter(lipid != "C18:0") %>% 
  select(-Compound, -Class) %>% 
  pivot_wider(values_from = "recovery", names_from = c("Type", "lipid")) %>% 
  pivot_longer(cols = c('TAG_C14:0':'GL_C16:1'), names_to = "Type", values_to = "recovery") %>% 
  mutate(recovery = replace(recovery, is.na(recovery), 0)) %>% 
  separate(col = Type, sep = "_", into = c("Type", "lipid"))

rec_new %>% 
  mutate(Type = str_replace(Type, "GL", "Glycolipid"),
         Type = str_replace(Type, "PL", "Phospholipids"),
         Type = as_factor(Type),
         Type = fct_relevel(Type, c("TAG", "Sterols", "Glycolipid", "Betaine lipid", "Phospholipids"))) %>% 
  group_by(treatment, fraction, test, eluent_vol, Type, rep) %>% 
  summarise(recovery = mean(recovery)) %>% 
  droplevels() %>% 
  mutate(treatment = as_factor(treatment),
         treatment = fct_recode(treatment,
                                "Adjusted eluents" = "with_ethanol")) %>% 
  mutate(fraction = as_factor(fraction),
         fraction = fct_recode(fraction,
                               "Fraction 1 (NL)" = "1",
                               "Fraction 2 (GL)" = "2",
                               "Fraction 3 (PL)" = "3")) %>% 
  mutate(test = as_factor(test),
         test = fct_recode(test,
                           "Eluent volumes after Buyer & Sasser (2012)" = "first",
                           "Adjusted eluent volumes" = "second")) %>% 
  mutate(eluent_vol = as_factor(eluent_vol),
         eluent_vol = fct_relevel(eluent_vol, "2 ml", "1.5 ml", "0.5 ml", "1 ml")) %>% 
  # Plot it!
  ggplot(aes(x = eluent_vol, y = recovery, fill = test)) +
  geom_hline(yintercept = 100, linetype = 2, linewidth = 0.1) +
  stat_summary(geom = "bar", fun = mean, width = 0.8) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  facet_grid2(fraction ~ Type, scales = "free", independent = "x") +
  scale_fill_manual(values = c("#12496d", "#1f78b4")) +
  labs(x = "", y = "Recovery (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        text = element_text(family = "Arial"))
ggsave("Fig_S1.jpg", units = "cm", width = 18, height = 15)



# How much does recovery improve from 0.5 to 1 ml for PLFA?
# -> only adjusted eluents (+ 2% ethanol in chloroform)
rec %>% 
  filter(treatment == "with_ethanol",
         Type == "PL",
         fraction == 3) %>% 
  mutate(fraction = as_factor(fraction),
         fraction = fct_recode(fraction,
                               "Fraction 3" = "3")) %>% 
  mutate(test = as_factor(test),
         test = fct_recode(test, 
                           "0.5 ml" = "first",
                           "1 ml" = "second")) %>% 
  group_by(treatment, type, fraction, test, lipid, Compound, Type) %>% 
  summarise(avg_recovery = mean(recovery)) %>% 
  group_by(treatment, type, fraction, test, Type) %>% 
  summarise(avg_recovery = mean(avg_recovery))


#