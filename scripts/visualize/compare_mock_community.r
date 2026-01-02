library(tidyverse)
library(ggrepel)
library(ggallin)
library(ggpmisc)
library(ggsignif)
library(yardstick)
library(Metrics)
library(ggpubr)
library(ggstatsplot)
library(patchwork)

# First run scripts 01-08 within the mock_community_analysis directory
# This script plots the per-taxon error rates and per-community filter results

# Read in true vs predicted abundances per taxon (output from 07)
eval_df <- read_csv("data/mock_community/evaluation_results.csv")  %>% 
    mutate(Database = recode(db_name, 
                             "GTDB 207" = "GTDB 207"))

# Read in summary of reads classified and passing the false-positive filter (output from 08)
filter_df <- read_csv("data/mock_community/filter_results.csv") %>%
    mutate(Database = recode(db_name,
                             "soil_microbe_db" = "SoilMicrobeDB",
                             "gtdb_207_unfiltered" = "GTDB 207",
                             "pluspf" = "PlusPF"))



# Calculate error rates per taxon for each database
rmse_df <- eval_df %>% 
    # subset to taxa that are either in mock community or predicted 
    filter(abundance != 0 | predicted_abundance != 0) %>% 
    group_by(Database, fungal_proportion, taxon, readdepth,  
             fungal) %>% 
    yardstick::rmse(truth = abundance, estimate = predicted_abundance) %>% 
    rename("metric" = ".metric", "value" = ".estimate") %>% 
    mutate(metric="Root mean squared error (RMSE)")

bias_df <- eval_df %>% 
    # subset to taxa that are either in mock community or predicted 
    filter(abundance != 0 | predicted_abundance != 0) %>% 
    group_by(Database, fungal_proportion, taxon, readdepth, 
             fungal) %>% 
    summarize(value = Metrics::bias(actual = abundance, 
                                    predicted = predicted_abundance)) %>% 
    mutate(metric = "Bias")

metrics_to_plot = full_join(rmse_df, bias_df)

bias_to_plot = metrics_to_plot %>% filter(metric == "Bias" &
                                              fungal_proportion==10 & readdepth=="5e+06")

rmse_to_plot = metrics_to_plot %>% filter(metric == "Root mean squared error (RMSE)" &
                                              fungal_proportion==10 & readdepth=="5e+06")

bias_to_plot_fungi = bias_to_plot %>% filter(fungal=="Fungi")
bias_to_plot_notfungi = bias_to_plot %>% filter(fungal != "Fungi")

rmse_to_plot_fungi = rmse_to_plot %>% filter(fungal=="Fungi")
rmse_to_plot_notfungi = rmse_to_plot %>% filter(fungal != "Fungi")

rmse_df_fungi <- pairwise_comparisons(rmse_to_plot_fungi, Database, value) %>%
    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
    dplyr::arrange(group1)
rmse_df_fungi$fungal="Fungi"


bias_df_fungi <- pairwise_comparisons(bias_to_plot_fungi, Database, value) %>%
    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
    dplyr::arrange(group1)
bias_df_fungi$fungal="Fungi"


# Colorblind-friendly palette (Okabe-Ito inspired)
myColors <- c("SoilMicrobeDB" = "#E69F00",  # Orange
              "GTDB 207" = "#56B4E9",       # Sky blue
              "PlusPF" = "#009E73")         # Bluish green
colScale <- scale_colour_manual(name = "Database", values = myColors)


p1 = ggbetweenstats(
    data  = rmse_to_plot_fungi,
    x     = Database,
    y     = value,
    ylab = "RMSE (abundance in simulated community)",
    xlab = "Database",
    type="np",
    ggsignif.args    = list(textsize = 3.5, tip_length = 0.06),
    ggtheme      = theme_bw(base_size = 18),
    results.subtitle=F,
    centrality.plotting = F,
    var.equal=F,
    sample.size.label = F,
    ggplot.component = list(
        theme(plot.title=element_text(size=18),
              plot.margin = margin(t = 60, r = 5.5, b = 5.5, l = 5.5, unit = "pt"),
              axis.text.y = element_text(size = 12)),
        scale_y_sqrt(),#(limits = c(10^0,10^7),
        #               breaks = scales::trans_breaks("log10", function(x) 10^x),
        #               labels = scales::trans_format("log10", scales::math_format(10^.x))),
        annotation_logticks(sides = "l"),
        ylim(c(0,.0115)),
        coord_cartesian(clip = "off")
    )
) + ggsignif::geom_signif(step_increase = .01,
                              textsize = 3.5,
    comparisons = list(rmse_df_fungi$groups[[1]]),
    annotations = as.character(rmse_df_fungi$expression)[[1]],
    test        = NULL,
    na.rm       = TRUE,
    parse       = TRUE
)  + ggtitle("Fungi") +  #ylim(c(0,.011))  + 
    colScale +
    theme(axis.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)) +
    xlab(NULL)

p2 = ggbetweenstats(
    data  = rmse_to_plot_notfungi,
    x     = Database,
    y     = value,
    xlab = "Database",
    type="np",
    ggsignif.args    = list(textsize = 3.5, 
                            tip_length = 0.02,
                            step_increase = .04),
    ggtheme      = theme_bw(base_size = 18),
    results.subtitle=F,
    centrality.plotting = F,
    var.equal=F,
    sample.size.label = F,
    ggplot.component = list(
        theme(plot.title=element_text(size=18),
              plot.margin = margin(t = 60, r = 5.5, b = 5.5, l = 5.5, unit = "pt")),
         scale_y_sqrt(),#(limits = c(10^0,10^7),
        #               breaks = scales::trans_breaks("log10", function(x) 10^x),
        #               labels = scales::trans_format("log10", scales::math_format(10^.x))),
        annotation_logticks(sides = "l"),
        ylim(c(0,.0115)),
        coord_cartesian(clip = "off")
    )
)  + ggtitle("Prokaryote") + 
    theme(axis.title = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 12)) +
    ylab(NULL) + #ylim(c(0,.011))  + 
    colScale +
    xlab(NULL) #+ scale_y_sqrt()





p3 = ggbetweenstats(
    data  = bias_to_plot_fungi,
    x     = Database,
    y     = value,
    ylab = "Bias (simulated community)",
    xlab = "Database",
    type="np",
    ggsignif.args    = list(textsize = 5, tip_length = 0.04),
    ggtheme      = theme_bw(base_size = 18),
    results.subtitle=F,
    centrality.plotting = F,
    var.equal=F,
    sample.size.label = F,
    ggplot.component = 
        list(ggtitle("Fungal genomes"),
    geom_hline(yintercept = 0, linetype=2))
    ) + 
    ylim(c(-.005,.005))  + colScale +
    theme(axis.title = element_text(face = "bold")) +
    xlab(NULL) 

p4 = ggbetweenstats(
    data  = bias_to_plot_notfungi,
    x     = Database,
    y     = value,
    xlab = "Database",
    type="np",
    ggsignif.args = list(textsize = 5, tip_length = 0.06),
    ggtheme  = theme_bw(base_size = 18),
    results.subtitle=F,
    centrality.plotting = F,
    var.equal=F,
    sample.size.label = F,
    ggplot.component = list(ggtitle("Prokaryote genomes"),
                            geom_hline(yintercept = 0, linetype=2),
        ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis(name = NULL))
    ),
    pairwise.display = "none"
)  +# ggtitle("Prokaryote") +
    ylab(NULL) + 
    xlab(NULL) + 
    ylim(c(-.005,.005))  + colScale



# Prepare data for panel C with labels
p5_data <- filter_df %>% 
    filter(readdepth != 1e7 & metric == "percent_passing" & 
           db_name %in% c("soil_microbe_db","gtdb_207_unfiltered","pluspf")) %>%
    group_by(Database, fungal_proportion) %>%
    summarize(mean_value = mean(value), .groups = "drop")

# Get label positions (right side of plot)
label_data <- p5_data %>%
    group_by(Database) %>%
    filter(fungal_proportion == max(fungal_proportion)) %>%
    ungroup()

p5 = ggplot(p5_data, aes(y = mean_value, x = fungal_proportion, color = Database)) +
    geom_point(data = filter_df %>% 
                   filter(readdepth != 1e7 & metric == "percent_passing" & 
                          db_name %in% c("soil_microbe_db","gtdb_207_unfiltered","pluspf")),
               aes(y = value, x = fungal_proportion),
               size = 3, alpha = .3,
               position = position_jitterdodge(jitter.height = 0, jitter.width = 1, dodge.width = 1)) +
    geom_line(linewidth = 1.2) +
    geom_text(data = label_data, 
              aes(x = fungal_proportion, y = mean_value, label = Database),
              inherit.aes = FALSE, color = "black", hjust = -0.1, vjust = 0.5, size = 4, fontface = "bold",
              show.legend = FALSE) +
    theme_bw(base_size = 18) + 
    ylab("% reads confidently classified") + 
    xlab("% of fungi in simulated community") + 
    theme(axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12),
          legend.position = "none",
          plot.margin = margin(t = 5.5, r = 100, b = 5.5, l = 5.5, unit = "pt")) +
    colScale +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    coord_cartesian(clip = "off")


# Figure 3
fig3 = (p1 | p2) / p5 + 
    plot_annotation(tag_levels ="A",
                    theme = 
                        theme(plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt")))
# Save figure
ggsave("manuscript_figures/fig3.png", fig3, width = 10, height = 9, units = "in", dpi = 300)


fig_s3 = p3 + p4 +
    plot_annotation(tag_levels ="A",
                    theme = 
                        theme(plot.margin = unit(c(-.5,0,0,0), 
                                                 'lines')))
    
# Save figure
ggsave("manuscript_figures/fig_s3.png", fig_s3, width = 8, height = 4, units = "in", dpi = 300)
