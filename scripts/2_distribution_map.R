# Supplementary code to the article :
# Friesová et al. Forest herb-layer species richness in the western Balkan diversity hotspot

# Author: Klára Friesová
# Date: 2026-03-30

library(terra) # version 1.8-60
library(sf) # version 1.0-19
library(ggspatial) # version 1.1.10
library(tidyterra) # version 0.7.2
library(ggnewscale) # version 0.5.2
library(rnaturalearth) # version 1.0.1
library(rnaturalearthdata) # version 1.0.0
library(cowplot) # version 1.2.0
library(tidyverse) # version 2.0.0
# dplyr 1.1.4
# forcats 1.0.0
# ggplot2 4.0.0
# readr 2.1.5
# tibble 3.3.0
# tidyr 1.3.1


# load data ---------------------------------------------------------------

head <- read_csv('data/Montenegro_species_richness_data.csv') |> 
  mutate(veg_type = fct_relevel(veg_type, 'evergreen oak', 'mixed deciduous', 'beech', 'pine', 'fir and spruce', 'riparian')) |> 
  filter(!is.na(PH) & !is.na(twi_dem) & !is.na(hli_data) & PH != 0)


# montenegro borders
MNE <- read_sf('shps/MNE_adm0.shp') |> 
  st_transform(crs = 'epsg:3035')

# Hijmans, Robert J.. University of California, Berkeley. Museum of Vertebrate Zoology. 
# Boundary, Montenegro, 2015. [Shapefile]. University of California, Berkeley. Museum of Vertebrate Zoology. 
# Retrieved from https://maps.princeton.edu/catalog/stanford-hz193rg6075

# digital elevation model from EU-DEM
dem <- rast('data_big/DEM_Cerna_Hora/DEM_Cerna_Hora_crop.tif')

# define colors for vegetation types
veg_col <- c('beech' = '#4daf4a', 
             'fir and spruce' = '#984ea3', 
             'mixed deciduous' = '#ff7f00', 
             'evergreen oak' = '#ffff33', 
             'riparian' = '#377eb8', 
             'pine' = '#e41a1c')


# map of study sites in Montenegro
map_mn <- head |> 
  arrange(richness_herb) |> 
  ggplot()+
  geom_spatraster(data = dem, alpha = 0.8)+
  scale_fill_gradient2(low = '#4d9221', mid = '#ffffbf', high = '#8c510a', 
                       name = 'Elevation [m a.s.l.]', na.value = NA, midpoint = 1100)+
  geom_sf(data = MNE, alpha = 0)+
  new_scale_fill()+
  geom_point(aes(X, Y, size = richness_herb, fill = veg_type), pch = 21, alpha = 0.8)+
  scale_fill_manual(values = veg_col)+
  scale_size(range = c(0.6, 9), breaks = c(40, 80, 120))+
  theme_minimal()+
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        axis.title = element_blank()) +
  labs(size = 'Species richness', fill = 'Forest type')+
  annotation_scale(location = 'br', pad_y = unit(1.5, 'cm')) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

ggsave('plots/Fig1_species_richness_map.png', width = 10, height = 8)
ggsave('plots/Fig1_species_richness_map.jpg', width = 10, height = 8)


# overview map with the position of Montenegro in Europe
world_map <- ne_countries(scale = "medium", returnclass = "sf")

mn <- world_map %>% 
  filter(sovereignt == "Montenegro")

map_EU <- ggplot() + 
  geom_sf(data = world_map, fill = "#FFFFCC") +
  geom_sf(data = mn, fill = "#940312") + 
  scale_x_continuous(limits = c(-10, 33)) +
  scale_y_continuous(limits = c(37, 60))+
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", colour = "white"), axis.title = element_blank())

final_map <- ggdraw() +
  draw_plot(map_mn) + 
  draw_plot(map_EU, x = 0.522, y = 0.75, width = 0.23, height = 0.23)


final_map

# save and modify in Inkscape
ggsave('plots/Fig1_species_richness_map.png', final_map, width = 10, height = 8)
ggsave('plots/Fig1_species_richness_map.jpg', final_map,  width = 10, height = 8)
