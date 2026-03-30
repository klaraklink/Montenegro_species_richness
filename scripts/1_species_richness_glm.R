# Supplementary code to the article :
# Friesová et al. Forest herb-layer species richness in the western Balkan diversity hotspot

# Author: Klára Friesová
# Date: 2026-03-30

library(ggfortify) # version 0.4.19
library(readxl) # version 1.4.5
library(MASS) # version 7.3-65
library(ecospat) # version 4.1.2
library(emmeans) # version 1.11.2-8
library(multcomp) # version 1.4-28
library(ggeffects) # version 2.3.1
library(psych) # version 2.5.6
library(see) # version 0.12.0
library(broom) # version 1.0.10
library(corrplot) # version 0.95
library(ggnewscale) # version 0.5.2
library(patchwork) # version 1.3.2
library(tidyverse) # version 2.0.0
# dplyr 1.1.4
# forcats 1.0.0
# ggplot2 4.0.0
# purrr 1.1.0
# readr 2.1.5
# stringr 1.5.2
# tibble 3.3.0
# tidyr 1.3.1


# load data ---------------------------------------------------------------
head <- read_csv('data/Montenegro_species_richness_data.csv')

# data exploration
pdf('plots/pairs.pdf', height = 11.69, width = 16.54)
pairs.panels(head, scale = T, gap = 0)
dev.off()

# variable selection
head_eda <- head |>
  mutate(veg_type = fct_relevel(veg_type, 'evergreen oak', 'mixed deciduous', 'beech', 'pine', 'fir and spruce', 'riparian'), 
         SOIL_DEPTH = ifelse(SOIL_DEPTH > 30, 30, SOIL_DEPTH)) |> 
  filter(!is.na(PH) & !is.na(twi_dem) & !is.na(hli_data) & PH != 0) |> 
  select(richness_herb, veg_type, COV_TREES, COV_SHRUBS, cov_evergreen, SOIL_DEPTH, 
         PH, hli_data, mean_annual_temp, annual_prec, prec_warm, tri) 

# summary table for Table 1 
head_sum <- head_eda |> 
  group_by(veg_type) |> 
  summarise(n = n(), across(c(richness_herb, COV_TREES:tri), 
                            list(min = min, mean = mean, max = max), na.rm = T))

# log-transformation of skewed predictors
head2 <- head_eda |> 
  mutate(tri = log(tri))

# check correlations and distributions after transformation
pdf('plots/pairs_trans.pdf', height = 11.69, width = 16.54)
pairs.panels(head2, scale = T, gap = 0)
dev.off()

# correlation plot for Appendix S1
png('plots/corrplot.png', height = 3500, width = 3500)
cor(head2 |> select(-richness_herb, -veg_type) |> 
      rename('Cover of trees' = 'COV_TREES', 'Cover of shrubs' = 'COV_SHRUBS', 
             'Cover of evergreen' = 'cov_evergreen', 'Soil pH' = 'PH', 
             'Soil depth' = 'SOIL_DEPTH',
             'Heat load index' = 'hli_data', 
             'Mean annual temperature' = 'mean_annual_temp', 
             'Annual precipitation' = 'annual_prec', 'Precipitation warmest quarter' = 'prec_warm', 
             'Terrain ruggedness index' = 'tri'), 
    use = 'pairwise.complete.obs') |> 
  corrplot(method = 'number', type = 'upper', diag = F, tl.col = 'black', 
           number.cex = 5, tl.cex = 5, tl.srt = 45, cl.cex = 4)
dev.off()


# model of herb-layer species richness ------------------------------------

m.0 <- glm.nb(richness_herb~1, data = head2)

m.full <- glm.nb(richness_herb~(.-mean_annual_temp - veg_type - prec_warm + poly(prec_warm, 2) 
                                 + poly(mean_annual_temp, 2))^2 - poly(mean_annual_temp, 2):poly(prec_warm, 2), data = head2)


anova(m.0, m.full, test = "Chisq")
form.full <- formula(m.full)
m.step <- step(m.0, scope = form.full, direction = 'both')
summary(m.step)
anova(m.step)

# bivariate models for all predictors - linear relationships
mod_all <- head2 |> 
  pivot_longer(cols = -c(richness_herb, veg_type), names_to = 'coef') |> 
  nest(data = -coef) |> 
  mutate(mod = map(data, ~glm.nb(richness_herb~value, data = .x)),
         m_tidy = map(mod, ~tidy(.x)), 
         label = c('Cover of trees [%]', 'Cover of shrubs [%]', 'Cover of evergreen species [%]',
                   'Soil depth [cm]', 'Soil pH', 'Heat load index', 'Mean annual temperature [°C]',
                   'Annual precipitation [mm]', 'Precipitation of the warmest quarter [mm]',
                   'Terrain ruggedness index (log)'),
         m_plot = map2(mod, label, ~ .x |>
                         predict_response() |>
                         plot()+
                         labs(x = .y, y = 'Species richness') +
                         theme(title = element_blank()))
         ) |> 
  unnest(m_tidy) |> 
  filter(term != '(Intercept)' & p.value < 0.05) |> 
  arrange(desc(abs(statistic)))

# quadratic relationships
mod_poly <- head2 |> 
  pivot_longer(cols = -c(richness_herb, veg_type), names_to = 'coef') |> 
  nest(data = -coef) |> 
  filter(coef %in% c('mean_annual_temp', 'prec_warm')) |>
  mutate(mod = map(data, ~glm.nb(richness_herb~poly(value, 2), data = .x)),
         m_tidy = map(mod, ~tidy(.x)), 
         label = c('Mean annual temperature [°C]', 'Precipitation of the warmest quarter [mm]'), 
         m_plot = map2(mod, label, ~ .x |>
                         predict_response() |>
                         plot()+
                         labs(x = .y, y = 'Species richness') +
                         theme(title = element_blank()))
         ) |> 
  unnest(m_tidy) |> 
  filter(term != '(Intercept)') |> 
  arrange(desc(abs(statistic)))


# plot all significant predictors in one figure for Appendix S3
mod_poly$m_plot[[1]] + mod_all$m_plot[[6]] + mod_poly$m_plot[[2]] + 
  mod_all$m_plot[[3]] + mod_all$m_plot[[4]] + mod_all$m_plot[[2]] + 
   mod_all$m_plot[[1]] +  
   plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

ggsave('plots/appendix_sp_richness_univariate_effects.png', width = 10, height = 7)


# model with vegetation types as predictor
m_veg <- glm.nb(richness_herb ~ veg_type, data = head2)
summary(m_veg)

predict_response(m.step) |> plot()


# model diagnostics
autoplot(m.step)+
  theme_bw()

ecospat.adj.D2.glm(m.step)
AIC(m.step)


# plot results ------------------------------------------------------------

# marginal effects of predictors  
pred.list <- map(.x = c('COV_TREES', 'annual_prec', 'prec_warm', 'tri'),
                 ~predict_response(m.step, terms = c(.x), margin = 'marginalmeans') |> 
                   as_tibble() |> 
                   mutate(coef = .x)) |> 
  bind_rows() |> 
  mutate(x = ifelse(coef == 'tri', exp(x), x))


# raw data for plotting, only predictors selected in the model
head.m <- head2 |> 
  mutate(tri = exp(tri)) |> 
  pivot_longer(cols = -c(richness_herb, veg_type), names_to = 'coef') |>
  semi_join(pred.list) |> 
  rename('group' = 'veg_type')

# define colors for vegetation types
veg_col <- c('evergreen oak' = '#D9BA0D', 
             'mixed deciduous' = '#ff7f00',
             'beech' = '#4daf4a', 
             'pine' = '#e41a1c', 
             'fir and spruce' = '#984ea3', 
             'riparian' = '#377eb8') 

# pH vs temperature
head2 |> 
  arrange(richness_herb) |> 
  ggplot(aes(x = mean_annual_temp, y = PH)) +
  geom_point(aes(fill = veg_type, size = richness_herb), pch = 21, alpha = 0.7)+
  scale_fill_manual(values = veg_col)+
  scale_size_continuous(breaks = c(40, 80, 120))+
  theme_bw()+
  labs(size = 'Species richness', fill = 'Forest type', x = 'Mean annual temperature [°C]', y = 'Soil pH')

ggsave('plots/Fig2_ph_temperature.png', width = 7, height = 5.5)

# plot model results  
# significant predictors apart from the interactions
plot1 <- pred.list |> 
  filter(coef == 'COV_TREES') |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, color = 0)+
  geom_line()+
  geom_point(data = head.m |> filter(coef == 'COV_TREES'), 
             aes(x = value, y = richness_herb, color = group, shape = group), 
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  scale_x_continuous(limits = c(20, 100), expand = c(0, 0))+
  theme_bw()+
  theme(legend.position = 'bottom', legend.direction = 'horizontal')+
  labs(x = 'Cover of trees [%]', y = 'Species richness')+
  guides(color = guide_legend(title.position="top", title.hjust = 0.5, byrow = T), 
         shape = guide_legend(title.position="top", title.hjust = 0.5, byrow = T))

plot2 <- pred.list |> 
  filter(coef == 'prec_warm') |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, color = 0)+
  geom_line()+
  geom_point(data = head.m |> filter(coef == 'prec_warm'), 
             aes(x = value, y = richness_herb, color = group, shape = group), 
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  scale_x_continuous(limits = c(120, 440), expand = c(0, 0))+
  theme_bw()+
  theme(legend.position = 'bottom', legend.direction = 'horizontal')+
  labs(x = 'Precipitation of the warmest quarter [mm]', y = 'Species richness')+
  guides(color = guide_legend(title.position="top", title.hjust = 0.5, byrow = T), 
         shape = guide_legend(title.position="top", title.hjust = 0.5, byrow = T))

plot3 <- pred.list |> 
  filter(coef == 'tri') |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, color = 0)+
  geom_line()+
  geom_point(data = head.m |> filter(coef == 'tri'), 
             aes(x = value, y = richness_herb, color = group, shape = group), 
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  scale_x_continuous(limits = c(0, 92), expand = c(0, 0))+
  theme_bw()+
  theme(legend.position = 'bottom', legend.direction = 'horizontal')+
  labs(x = 'Terrain ruggedness index', y = 'Species richness')+
  guides(color = guide_legend(title.position="top", title.hjust = 0.5, byrow = T), 
         shape = guide_legend(title.position="top", title.hjust = 0.5, byrow = T))

# join together and save
plot1 / plot2 / plot3 + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')+
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")  

ggsave('plots/Fig4_sp_richness_marg_effects_all.png', width = 4.5, height = 11)


# differences between vegetation types
m_means <- emmeans(m_veg, ~ veg_type) |> 
  cld(Letters = letters) |> 
  as_tibble() |> 
  mutate(.group = str_trim(.group)) |> 
  dplyr::select(x = veg_type, sig = .group)

# violin plot for vegetation types with pointrange for model results 
predict_response(m_veg, terms = 'veg_type', margin = 'marginalmeans') |> 
  as_tibble() |> 
  left_join(m_means) |> 
  mutate(x = fct_relevel(x, 'beech', 'fir and spruce', 'mixed deciduous', 'evergreen oak', 
                         'riparian', 'pine')) |> 
  ggplot()+
  geom_violinhalf(data = head.m |> distinct(richness_herb, group, value),
                  aes(x = group, y = richness_herb, fill = group, color = group),
                  alpha = 0.4)+
  geom_pointrange(aes(x, predicted, ymin = conf.low, ymax = conf.high, fill = x), 
                  pch = 21, fatten = 5, linewidth = 0.7)+
  geom_text(aes(x, label = sig, color = x), y = 160)+
  scale_fill_manual(values = veg_col, aesthetics = c('fill', 'color'))+
  scale_y_continuous(limits = c(0, 170), expand = c(0, 0), breaks = seq(0, 160, 20))+
  scale_x_discrete(limits = rev, labels = c('riparian (9)', 'fir and spruce (33)', 'pine (16)', 'beech (66)', 
                                            'mixed deciduous (209)', 'evergreen oak (9)'))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.y = element_blank(), legend.position = 'none', 
        axis.text.y = element_text(size = 11), 
        panel.grid.major.y = element_blank())+
  labs(y = 'Species richness', fill = 'Forest type')

ggsave('plots/Fig3_sp_richness_vegtypes.png', width = 5, height = 5)

# plot interaction effects
intplot1 <- predict_response(m.step, terms = c('mean_annual_temp', 'PH[5, 7]'), margin = 'marginalmeans') |> 
  as_tibble() |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = 0)+
  geom_line(aes(color = group))+
  scale_fill_manual(values = c('#F55F73', '#A47FF5'), aesthetics = c('color', 'fill'),
                    name = 'Soil pH')+
  new_scale_color()+
  geom_point(data = head2,
             aes(x = mean_annual_temp, y = richness_herb, color = veg_type, shape = veg_type),
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  theme_bw()+
  labs(x = 'Mean annual temperature [°C]', y = 'Species richness') 

intplot2 <- predict_response(m.step, terms = c('mean_annual_temp', 'tri[2, 3.5]'), margin = 'marginalmeans') |> 
  as_tibble() |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = 0)+
  geom_line(aes(color = group))+
  scale_fill_manual(values = c('#837233', '#ffdd61'), aesthetics = c('color', 'fill'),
                    name = 'Terrain ruggedness index')+
  new_scale_color()+
  geom_point(data = head2,
             aes(x = mean_annual_temp, y = richness_herb, color = veg_type, shape = veg_type),
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  theme_bw()+
  labs(x = 'Mean annual temperature [°C]', y = 'Species richness') 

intplot3 <- predict_response(m.step, terms = c('prec_warm', 'tri[2, 3.5]'), margin = 'marginalmeans') |> 
  as_tibble() |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = 0)+
  geom_line(aes(color = group))+
  scale_fill_manual(values = c('#837233', '#ffdd61'), aesthetics = c('color', 'fill'),
                    name = 'Terrain ruggedness index')+
  new_scale_color()+
  geom_point(data = head2,
             aes(x = prec_warm, y = richness_herb, color = veg_type, shape = veg_type),
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  theme_bw()+
  labs(x = 'Precipitation of the warmest quarter [mm]', y = 'Species richness') 

intplot4 <- predict_response(m.step, terms = c('PH', 'tri[2, 3.5]'), margin = 'marginalmeans') |> 
  as_tibble() |> 
  ggplot(aes(x, predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = 0)+
  geom_line(aes(color = group))+
  scale_fill_manual(values = c('#837233', '#ffdd61'), aesthetics = c('color', 'fill'),
                    name = 'Terrain ruggedness index')+
  new_scale_color()+
  geom_point(data = head2,
             aes(x = PH, y = richness_herb, color = veg_type, shape = veg_type),
             alpha = 0.7)+
  scale_color_manual(values = veg_col, name = 'Forest type')+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 8), name = 'Forest type')+
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks = seq(0, 160, 20))+
  theme_bw()+
  labs(x = 'Soil pH', y = 'Species richness') 

# join together and save, modified afterwards in Inkscape
intplot1 + intplot2 + intplot3 + intplot4 + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom", legend.direction = 'horizontal', legend.box = 'vertical')

ggsave('plots/Fig5_sp_richness_interaction.png', width = 8, height = 8)
ggsave('plots/Fig5_sp_richness_interaction.svg', width = 8, height = 8)