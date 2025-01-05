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
eval_df <- read_csv("/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/evaluation_results.csv")  %>% 
    mutate(Database = recode(db_name, 
                             "GTDB 207" = "GTDB r207"))

# Read in summary of reads classified and passing the false-positive filter (output from 08)
filter_df <- read_csv("/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/filter_results.csv")



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


library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- c("SoilMicrobeDB", "GTDB r207", "PlusPF")
colScale <- scale_colour_manual(name = "Database",values = myColors)


p1 = ggbetweenstats(
    data  = rmse_to_plot_fungi,
    x     = Database,
    y     = value,
    ylab = "RMSE (abundance in simulated community)",
    xlab = "Database",
    type="np",
    ggsignif.args    = list(textsize = 5, tip_length = 0.06),
    ggtheme      = theme_bw(base_size = 18),
    results.subtitle=F,
    centrality.plotting = F,
    var.equal=F,
    sample.size.label = F,
    ggplot.component = list(
        theme(plot.title=element_text(size=18)),
        ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis(name = NULL)),
        scale_y_sqrt(),#(limits = c(10^0,10^7),
        #               breaks = scales::trans_breaks("log10", function(x) 10^x),
        #               labels = scales::trans_format("log10", scales::math_format(10^.x))),
        annotation_logticks(sides = "l"),
        ylim(c(0,.0115))
    )
) + ggsignif::geom_signif(step_increase = .01,
                              textsize = 5,
    comparisons = list(rmse_df_fungi$groups[[1]]),
    annotations = as.character(rmse_df_fungi$expression)[[1]],
    test        = NULL,
    na.rm       = TRUE,
    parse       = TRUE
)  + ggtitle("Fungi") +  #ylim(c(0,.011))  + 
    colScale +
    theme(axis.title = element_text(face = "bold")) +
    xlab(NULL)

p2 = ggbetweenstats(
    data  = rmse_to_plot_notfungi,
    x     = Database,
    y     = value,
    xlab = "Database",
    type="np",
    ggsignif.args    = list(textsize = 5, 
                            tip_length = 0.02,
                            step_increase = .04),
    ggtheme      = theme_bw(base_size = 18),
    results.subtitle=F,
    centrality.plotting = F,
    var.equal=F,
    sample.size.label = F,
    ggplot.component = list(
        theme(plot.title=element_text(size=18)),
        ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis(name = NULL)),
         scale_y_sqrt(),#(limits = c(10^0,10^7),
        #               breaks = scales::trans_breaks("log10", function(x) 10^x),
        #               labels = scales::trans_format("log10", scales::math_format(10^.x))),
        annotation_logticks(sides = "l"),
        ylim(c(0,.0115))
    )
)  + ggtitle("Prokaryote") + 
    theme(axis.title = element_text(face = "bold")) +
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



p5 = ggplot(filter_df %>% 
           filter(readdepth != 1e7 & metric == "percent_passing" & 
                      db_name %in% c("soil_microbe_db","gtdb_207_unfiltered","pluspf")),
       aes(y = value, x=fungal_proportion, color = Database)) +
    geom_point(#aes(color=fungal),
        size=3, alpha=.3,
        #position=position_jitter(width = .1, height=0.01)) +
        position=position_jitterdodge(jitter.height = 0, jitter.width = 1, dodge.width = 1)) +
    geom_line() +
    theme_bw(base_size = 18) + 
    #facet_grid(~readdepth, drop=T, scales="free")  +
    ylab("% reads confidently classified") + 
    xlab("% of fungi in simulated community")  + 
    theme(axis.title = element_text(face = "bold")) +
    #guides(color=guide_legend(NULL)) + 
    colScale


# Figure 3
fig3 = (p1 | p2) / p5 + 
    plot_annotation(tag_levels ="A",
                    theme = 
                        theme(plot.margin = unit(c(-.5,0,0,0), 
                                                 'lines')))
fig3 %>% 
    ggexport(height = 900, width = 800, filename = "/projectnb/frpmars/soil_microbe_db/manuscript_figures/fig3.png")


fig_s3 = p3 + p4 +
    plot_annotation(tag_levels ="A",
                    theme = 
                        theme(plot.margin = unit(c(-.5,0,0,0), 
                                                 'lines')))
    
fig_s3 %>% 
        ggexport(height = 400, width = 800, filename = "/projectnb/frpmars/soil_microbe_db/manuscript_figures/fig_s3.png")
