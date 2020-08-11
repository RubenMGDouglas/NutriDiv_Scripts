# ====================================================================================================#
# ====================================================================================================#
# Bivariate Mapping
# ====================================================================================================#

# # https://cran.r-project.org/web/packages/biscale/vignettes/biscale.html
library(biscale)
library(ggplot2)
library(cowplot)

# Once data are loaded, bivariate classes can be applied with the bi_class() function
bi.richness.data <- bi_class(na.omit(SESwrld_Moll), x = Phylogenetic_Faith, y = FRic, dim = 3)

# Once breaks are created, we can use bi_scale_fill() as part of our ggplot() call:
bi.richness.map <- ggplot() +
  geom_sf(data = bi.richness.data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  labs(
    title = "(a)  FRic & Faith Phylogenetic Diversity") +
  bi_theme() + theme(plot.title = element_text(size = 10))   

bi.richness.legend <- bi_legend(pal = "DkViolet",
                    dim = 3,
                    xlab = "Higher Faith PD",
                    ylab = "Higher FRic",
                    size = 8)

# combine map with legend
final.bi.richness.Plot <- ggdraw() +
  draw_plot(bi.richness.map, 0, 0, 1.0, 1.0) +
  draw_plot(bi.richness.legend, 0, 0.1, 0.3, 0.3)

# dispersion
bi.dispersion.data <- bi_class(na.omit(SESwrld_Moll), x = Phylogenetic_MPD, y = FDis, dim = 3)

# Once breaks are created, we can use bi_scale_fill() as part of our ggplot() call:
bi.dispersion.map <- ggplot() +
  geom_sf(data = bi.dispersion.data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  labs(
    title = "(b)  FDis & Phylogenetic Mean Pairwise Distance") +
  bi_theme() + theme(plot.title = element_text(size = 10)) 

bi.dispersion.legend <- bi_legend(pal = "DkViolet",
                                dim = 3,
                                xlab = "Higher Phylogenetic MPD",
                                ylab = "Higher FDis",
                                size = 8)

# combine map with legend
final.bi.dispersion.Plot <- ggdraw() +
  draw_plot(bi.dispersion.map, 0, 0, 1, 1) +
  draw_plot(bi.dispersion.legend, 0, 0.1, 0.3, 0.3)


# Combine Richness & Dispersion Plots with bivariate scatter plots

multiplot(final.bi.richness.Plot, final.bi.dispersion.Plot,cols=1)
multiplot(final.bi.richness.Plot,final.bi.dispersion.Plot,Bivar_FRicFaith,Bivar_FDisMPD,cols=2)

# Evenness
bi.evenness.data <- bi_class(na.omit(SESwrld_Moll), x = FEve, y = Phylogenetic_sp_evennes, dim = 3)

# Once breaks are created, we can use bi_scale_fill() as part of our ggplot() call:
bi.evenness.map <- ggplot() +
  geom_sf(data = bi.evenness.data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  labs(
    title = "FEve & Phylogenetic Species Evenness") +
  bi_theme()

bi.evenness.legend <- bi_legend(pal = "DkViolet",
                                dim = 3,
                                xlab = "Higher FEve",
                                ylab = "Higher PSE",
                                size = 8)

# combine map with legend
final.bi.evenness.Plot <- ggdraw() +
  draw_plot(bi.evenness.map, 0, 0, 1, 1) +
  draw_plot(bi.evenness.legend, 0, 0.1, 0.2, 0.2)
