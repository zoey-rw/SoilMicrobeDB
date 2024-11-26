
# Figure for manuscript showing counts for each unique genome in Soil Microbe DB

library(tidyverse)
library(data.table)

# Read in Struo2 input file
main_df <- read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/soil_genome_db_struo.tsv")



main_df %>% group_by(source,superkingdom) %>% tally() %>%
    mutate(source = recode(source, 
                           "SPIRE_MAGs" = "SPIRE catalog",
                           "SMAG" = "SMAG catalog",
                           "Mycocosm" = "JGI Mycocosm",
                           "GTDB" = "GTDB r207")) %>% 
    ungroup() %>% 
    #complete(source, superkingdom) %>% 
    ggplot(aes(x= reorder(source, -n), y = n,fill=superkingdom)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    coord_flip() +
    theme_bw(base_size = 20) +
    xlab("Genome source") + scale_y_sqrt() +
    ggtitle("SoilMicrobeDB contents") +
    ylab("Genome count") #+ facet_grid(~superkingdom)

# Values for table 1
main_df %>% group_by(source,superkingdom) %>% tally() 
