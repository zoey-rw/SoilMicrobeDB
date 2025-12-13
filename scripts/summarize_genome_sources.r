
# Figure for manuscript showing counts for each unique genome in Soil Microbe DB

library(tidyverse)
library(data.table)
if(!require(GGally, quietly = TRUE)) {
    install.packages("GGally", repos = "https://cloud.r-project.org")
    library(GGally)
} else {
    library(GGally)
}
library(ggpubr)
library(ggstatsplot)
library(patchwork)

# Read in Struo2 input file with MAG classification
main_df <- read_csv("data/genome_database/soil_microbe_db_genome_table.csv") %>% 
    mutate(source = recode(source, 
                           "SPIRE_MAGs" = "SPIRE",
                           "GEM catalog" = "GEM",
                           #"SMAG" = "SMAG",
                           "Mycocosm" = "JGI Mycocosm",
                           "GTDB" = "GTDB r207")) %>% 
    
    group_by(source,superkingdom) %>% 
    add_tally(name = "total") %>% 
    mutate(is_MAG) %>% 
    mutate(source = ordered(source, levels = 
                                rev(c("GTDB r207","SPIRE","SMAG","JGI GOLD","JGI Mycocosm","GEM","RefSoil"))))

main_df = main_df  %>% 
    mutate(overall = "Total in\nSoilMicrobeDB",
           genome_type = ifelse(is_MAG, "Metagenome-assembled\ngenomes",
                                "Isolate genomes"))

mag_count = main_df %>% ungroup %>% 
    group_by(source, superkingdom, total, is_MAG) %>% 
    tally(name = "n_MAGs") %>% filter(is_MAG) %>% 
    mutate(mag_percent = round(n_MAGs/total * 100, 1)) %>% 
    rename(n = total) %>% ungroup %>% 
    tidyr::complete(source, superkingdom, fill=list(mag_percent = NA)) %>% 
    select(source, superkingdom, n, mag_percent)


no_MAGs = main_df %>% ungroup %>% 
    group_by(source, superkingdom, total, is_MAG) %>% 
    tally(name = "n_MAGs") %>% filter(!is_MAG) %>% 
    mutate(mag_percent = round(n_MAGs/total * 100, 1)) %>% 
    rename(n = total) %>% filter(mag_percent == 100) %>% ungroup %>% 
    select(source, superkingdom, n) %>% mutate(mag_percent = 0)
mag_count = rbind.data.frame(mag_count, no_MAGs)

    
ggtable(main_df, "is_MAG", c("source"), cells = "row.prop") + ggtitle("Row proportions")


ggtable(main_df, "is_MAG", c("source"), cells = "observed")

 

main_df %>% group_by(source,superkingdom) %>% 
    tally() %>% 
    ungroup() %>% 
    #complete(source, superkingdom) %>% 
    ggplot(aes(x= reorder(source, -n), y = n,fill=superkingdom)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    coord_flip() +
    theme_bw(base_size = 20) +
    xlab("Genome source") + scale_y_sqrt() +
    ggtitle("SoilMicrobeDB contents")  +
    geom_text(data = mag_count, aes(label=mag_percent, 
                                    group=superkingdom), 
              #stat = "identity",
              position = position_dodge(width=.9),
       # width = 1),
hjust=-1) +
    ylab("Genome count") #+ facet_grid(~superkingdom)

# Values for manuscript and table 1
main_df %>% group_by(source,superkingdom) %>% tally() 

main_df %>% filter(source != "GTDB r207") %>% group_by(genome_type) %>% tally() 

main_df %>% group_by(genome_type) %>% tally() 
main_df %>% group_by(source,genome_type) %>% tally() 


# Kingdom
fig_1a = main_df %>% 
    #complete(source, superkingdom) %>% 
    ggplot(aes(x= source, fill=superkingdom)) +
    geom_bar(color = "black", position=position_fill()) +
    coord_flip() +
    theme_bw(base_size = 18) +
    ylab(NULL) + 
    xlab("Genome source") + 
    #ggtitle("SoilMicrobeDB contents")  +
    theme(legend.position = "top",legend.direction = "vertical",
          legend.key.height=unit(.7, "cm"),
          
          legend.title = element_blank(),
          
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = unit(c(0,0,0,0), "pt"),
          axis.title = element_text(face = "bold"), 
          axis.title.y.right = 
              element_text(size = 20), 
          panel.border = element_blank(), 
          strip.text = element_text(face = "bold"))  +
    scale_fill_brewer(palette="Set2") +
    scale_y_continuous(
                       breaks = c(0, 0.5, 1), 
                       labels = scales::percent(c(0, 0.5, 1))) 

# Kingdom

fig_1d = main_df %>% 
    ggplot(aes(x= overall, fill=superkingdom,  y = after_stat(prop),
               group = factor(superkingdom))) +
    geom_bar(color = "black", position=position_fill(),  stat = "prop", show.legend = F) +
    coord_flip() +
    theme_bw(base_size = 18) +
    xlab("") + 
    geom_label(aes(label = paste0(round(100 * after_stat(prop)), "%")), 
              stat = "prop", 
              position = position_fill(vjust = 0.5),
              fill="white",
              size = 5) + 
    theme(plot.margin = unit(c(0,0,0,0), "pt"),
          axis.title = element_text(face = "bold"), 
          axis.title.y.right = 
              element_text(size = 20), 
          panel.border = element_blank(), 
          strip.text = element_text(face = "bold"))  +
    scale_fill_brewer(palette="Set2") +
    scale_y_continuous(name = "% genomes",
                       breaks = c(0, 0.5, 1), 
                       labels = scales::percent(c(0, 0.5, 1))) 


# MAG
fig_1b = main_df  %>%  
    #complete(source, superkingdom) %>% 
    ggplot(aes(x= source, fill=genome_type)) +
    geom_bar(color = "black", position=position_fill()) +
    coord_flip() +
    theme_bw(base_size = 18) +
    ylab(NULL) +
    theme(legend.position = "top", legend.direction = "vertical",
          legend.key.height=unit(.9, "cm"),
          legend.title = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = unit( c(0,0,0,0), "pt"),
          axis.title = element_text(face = "bold"), 
          axis.title.y.right = 
              element_text(size = 20), 
          panel.border = element_blank(), 
          strip.text = element_text(face = "bold")) +
    xlab(NULL)  +
    scale_fill_brewer(palette="Paired") +
    scale_y_continuous(
                       breaks = c(0, 0.5, 1), 
                       labels = scales::percent(c(0, 0.5, 1))) 

fig_1e = main_df %>% 
    ggplot(aes(x= overall, fill=genome_type)) +
    geom_bar(color = "black", position=position_fill(), show.legend = F) +
    coord_flip() +
    theme_bw(base_size = 18) +
    ylab("Proportion") + 
    xlab(NULL) + 
    theme(legend.position = "top",legend.direction = "vertical",
          legend.key.height=unit(.7, "cm"),
          
          legend.title = element_blank())  +
    # geom_text(aes(label = paste0(value*100,"%")), 
    #           position = position_stack(vjust = 0.5),
    geom_label(aes(label = paste0(round(100 * after_stat(count) / sum(after_stat(count)), 1), "%")), 
              stat = "count",
              position = position_fill(vjust = 0.5),
              size = 5, fill="white") + 
    theme(  axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin = unit( c(0,0,0,0), "pt"),
            axis.title = element_text(face = "bold"), 
            axis.title.y.right = 
                element_text(size = 20), 
            panel.border = element_blank(), 
            strip.text = element_text(face = "bold"))  +
    scale_fill_brewer(palette="Paired") +
    scale_y_continuous(name = "% genomes",
                       breaks = c(0, 0.5, 1), 
                       labels = scales::percent(c(0, 0.5, 1))) 



tot <- main_df %>% 
    group_by(source) %>% 
    tally(name = "total") 

fig_1c = main_df %>% 
    ggplot(aes(x= source)) +
    geom_bar(color = "black") +
    coord_flip() +
    theme_bw(base_size = 18) +
    #scale_y_sqrt()   + 
    ylim(0,75000) + 
    #ylim(0,40000) +
    xlab(NULL) +
    ylab(NULL) +
    
    geom_text(data = tot,
              aes(label=paste0("n = ", total), y = total), 
              position = position_dodge(width=.9) , size=5,
              hjust=-.1)  + 
    theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = unit( c(0,0,0,0), "pt"),
        axis.title = element_text(face = "bold"), 
        axis.title.y.right = 
            element_text(size = 20), 
        panel.border = element_blank(), 
        strip.text = element_text(face = "bold"))  

tot_overall <- main_df %>% ungroup %>% 
    tally(name = "total") %>% mutate(overall="Total in\nSoilMicrobeDB")

fig_1f = main_df %>% 
    ggplot(aes(x= overall)) +
    geom_bar(color = "black") +
    coord_flip() +
    theme_bw(base_size = 18) +
    ylab("Genome count") + 
    #scale_y_sqrt()   +
    xlab(NULL) +
    geom_text(data = tot_overall,
              aes(label=paste0("n = ", total), y = total), 
              position = position_dodge(width=.9) , size=5,
              hjust=-.1)  + 
    ylim(0,75000) + 
    theme(  axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin = unit( c(0,0,0,0), "pt"),
            axis.title = element_text(face = "bold"), 
            axis.title.y.right = 
                element_text(size = 20), 
            panel.border = element_blank(), 
            strip.text = element_text(face = "bold")) 
fig1 = (fig_1a + fig_1b + fig_1c) / (fig_1d + fig_1e + fig_1f) +  
    plot_layout(widths=c(1,.6,1.5), heights=c(7,1)) 

fig1 


ggsave("manuscript_figures/fig1.png", fig1, width = 12, height = 8, units = "in", dpi = 300)







# Exploratory code below (commented out - not needed for figure generation)
# pluspf_inspect = read_tsv("https://genome-idx.s3.amazonaws.com/kraken/pluspf_20231009/inspect.txt", skip = 7, col_names = c(""))
# table(pluspf_inspect$X4) %>% sort
# pluspf_report = read_tsv( "https://genome-idx.s3.amazonaws.com/kraken/pluspf_20231009/library_report.tsv")
# table(pluspf_report$`#Library`) %>% sort
# pluspf_names = read_tsv("/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxonomy/names.dmp")
# length(unique(pluspf_names$`1`))
