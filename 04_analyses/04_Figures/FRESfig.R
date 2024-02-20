library(patchwork)

FRESfig <- TPNplot + ECGplot + porewater + plot_annotation(tag_levels = "A")

FRESfig <- TPNplot + BRplot2 + porewater + plot_annotation(tag_levels = "A")

FRESfig

ggsave(plot = FRESfig, filename = "figures/FRESfig.png",
       width = 10, height = 3, units = "in")
