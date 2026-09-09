# Presentation versions of winner_loser_boxplot.pdf and winner_loser_genus.pdf
# Same data pipeline as R/03_calculate_winners_losers.R, but colonizing_winner
# is recolored orange. Outputs go to plot/figures_presentation/.
#
# Assumes wlc, wl.sel, pools, metadata are already in the workspace
# (e.g. via load("Environment_Save_060924.RData") or by sourcing
# R/00_workflow.R and running R/03_calculate_winners_losers.R up through
# the wl.sel definition).

library(tidyverse)
library(scales)
library(geomtextpath)

cat.colors = c("#FE6100", "#06D6A0", "#EF476F", "grey80")
names(cat.colors) = c("colonizing_winner", "winner", "loser", "nc")

sp.list = read_csv("data/sp.list_TNRS.csv") %>%
  dplyr::select(Name_submitted, Accepted_name, Accepted_species) %>%
  distinct()

cat.sel =
  wl.sel %>%
  left_join(pools) %>%
  mutate(cat = ifelse(pval == "ns", "nc", wl)) %>%
  mutate(cat = ifelse(cat == "winner" & pool == "colonizing", "colonizing_winner", cat)) %>%
  mutate(cat = as.factor(cat)) %>%
  mutate(cat = factor(cat, levels = c("colonizing_winner", "winner", "loser", "nc")))

cat.sel =
  cat.sel %>%
  rename(Name_submitted = SpeciesName) %>%
  left_join(sp.list) %>%
  filter(!is.na(Accepted_name)) %>%
  rename(species = Accepted_name) %>%
  dplyr::select(-Accepted_species, -Name_submitted)

gf = c(
  "Achillea" = "forb", "Agrostis" = "grass", "Alchemilla" = "forb", "Alopecurus" = "grass",
  "Amaranthus" = "forb", "Androsace" = "forb", "Anthoxanthum" = "grass", "Anthyllis" = "forb",
  "Arrhenatherum" = "grass", "Aster" = "forb", "Astragalus" = "forb", "Avenella" = "grass",
  "Bistorta" = "forb", "Bouteloua" = "grass", "Brachypodium" = "grass", "Bromus" = "grass",
  "Calamagrostis" = "grass", "Campanula" = "forb", "Carex" = "grass", "Chaerophyllum" = "forb",
  "Crocus" = "forb", "Cyananthus" = "forb", "Dactylis" = "grass", "Deschampsia" = "grass",
  "Deyeuxia" = "grass", "Dryas" = "forb", "Drymaria" = "forb", "Elymus" = "grass",
  "Elyna" = "forb", "Equisetum" = "forb", "Erigeron" = "forb", "Eritrichium" = "forb",
  "Eurybia" = "forb", "Festuca" = "grass", "Fragaria" = "forb", "Galium" = "forb",
  "Geranium" = "forb", "Glechoma" = "forb", "Hedysarum" = "forb", "Helianthella" = "forb",
  "Helianthemum" = "forb", "Helictochloa" = "grass", "Helictotrichon" = "grass",
  "Heracleum" = "forb", "Hieracium" = "forb", "Hippocrepis" = "forb", "Homogyne" = "forb",
  "Hymenoxys" = "forb", "Hypericum" = "forb", "Knautia" = "forb", "Kobresia" = "grass",
  "Lancea" = "forb", "Laserpitium" = "forb", "Leontodon" = "forb", "Leontopodium" = "forb",
  "Ligusticum" = "forb", "Lolium" = "grass", "Lotus" = "forb", "Lupinus" = "forb",
  "Luzula" = "forb", "Microula" = "forb", "Minuartia" = "forb", "Muhlenbergia" = "grass",
  "Nardus" = "grass", "Papaver" = "forb", "Pedicularis" = "forb", "Persicaria" = "forb",
  "Phleum" = "grass", "Pimpinella" = "forb", "Plantago" = "forb", "Poa" = "grass",
  "Polygonum" = "forb", "Potentilla" = "forb", "Primula" = "forb", "Prunella" = "forb",
  "Pulsatilla" = "forb", "Pyrethrum" = "forb", "Ranunculus" = "forb", "Rosa" = "forb",
  "Salix" = "forb", "Saussurea" = "forb", "Scorzoneroides" = "forb", "Sedum" = "forb",
  "Senecio" = "forb", "Sesleria" = "grass", "Sibbaldia" = "forb", "Silene" = "forb",
  "Solidago" = "forb", "Stipa" = "grass", "Taraxacum" = "forb", "Thalictrum" = "forb",
  "Thymus" = "forb", "Trifolium" = "forb", "Trisetum" = "grass", "Trollius" = "forb",
  "Vaccinium" = "forb", "Veronica" = "forb", "Vicia" = "forb", "Viola" = "forb",
  "Youngia" = "forb", "Ligusticopsis" = "forb", "Pilosella" = "forb", "Tanacetum" = "forb",
  "Dasiphora" = "forb", "Sibbaldianthe" = "forb"
)

gf = as_tibble(gf) %>%
  rename(functional_group = value) %>%
  mutate(genus = names(gf))

per.sp =
  cat.sel %>%
  filter(!is.na(word(species, 2))) %>%
  dplyr::select(Region, originSiteID, destSiteID, species, cat) %>%
  mutate(genus = word(species, 1)) %>%
  distinct() %>%
  filter(cat != "nc") %>%
  group_by(genus) %>%
  mutate(total = n()) %>%
  group_by(genus, total, cat) %>%
  summarize(count = n()) %>%
  mutate(percent = (count / total) * 100) %>%
  left_join(gf) %>%
  mutate(gg_name = paste0(genus, " (", total, ")")) %>%
  mutate(functional_group = case_when(functional_group == "grass" ~ "Grasses",
                                      functional_group == "forb" ~ "Forbs")) %>%
  mutate(cat = factor(cat, levels = c("loser", "winner", "colonizing_winner"))) %>%
  ungroup() %>%
  arrange(functional_group, percent, cat) %>%
  mutate(genus = factor(genus, levels = unique(genus)),
         gg_name = factor(gg_name, levels = unique(gg_name)))

grayscale_color =
  per.sp %>%
  mutate(grayscale = rescale(total, to = c(0.7, 0)),
         grayscale_color = grey(grayscale)) %>%
  dplyr::select(gg_name, grayscale_color) %>%
  distinct()

label_data = data.frame(
  x = c(seq(1, 76, length.out = 100), seq(77, 100, length.out = 100)),
  y = c(rep(-10, 100), rep(-10, 100)),
  label = rep(c("Forbs", "Grasses"), each = 100)
)

p_genus =
  ggplot(per.sp, aes(x = gg_name, y = percent, fill = cat)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cat.colors) +
  coord_radial(start = 0.25 * pi, end = 1.75 * pi, inner.radius = 0.5) +
  guides(theta = guide_axis_theta(angle = 90),
         r     = guide_axis(angle = 90)) +
  geom_textpath(data = label_data, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 5, color = "black", text_smoothing = 50) +
  theme(axis.text.x = element_text(color = grayscale_color$grayscale_color),
        panel.border = element_blank(),
        legend.position = "bottom") +
  labs(x = "", y = "", fill = "")

pdf(file = "plot/figures_presentation/winner_loser_genus.pdf", width = 13, height = 10)
print(p_genus)
dev.off()

p_box =
  cat.sel %>%
  filter(!is.na(word(species, 2))) %>%
  filter(cat != "nc") %>%
  mutate(genus = word(species, 1)) %>%
  left_join(gf) %>%
  filter(!is.na(functional_group)) %>%
  mutate(functional_group = case_when(functional_group == "grass" ~ "Grasses",
                                      functional_group == "forb" ~ "Forbs")) %>%
  mutate(cat = factor(cat, levels = c("loser", "winner", "colonizing_winner"))) %>%
  ggplot(aes(cat, estimate, fill = cat, color = cat)) +
  TP_theme() +
  facet_wrap(~functional_group) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_jitter() +
  geom_text(aes(label = after_stat(sprintf("%.2f", y))),
            stat = "summary", fun = "median", vjust = -3, hjust = -0.1, size = 6) +
  scale_fill_manual(values = cat.colors) +
  scale_color_manual(values = cat.colors) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Winning/losing strength") +
  theme(axis.text.x = element_blank())

pdf(file = "plot/figures_presentation/winner_loser_boxplot.pdf", width = 8, height = 6)
print(p_box)
dev.off()

p_box_festuca =
  cat.sel %>%
  filter(!is.na(word(species, 2))) %>%
  filter(cat != "nc") %>%
  mutate(genus = word(species, 1)) %>%
  filter(genus == "Festuca") %>%
  mutate(cat = factor(cat, levels = c("loser", "winner", "colonizing_winner"))) %>%
  ggplot(aes(cat, estimate, fill = cat, color = cat)) +
  TP_theme() +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_jitter() +
  geom_text(aes(label = after_stat(sprintf("%.2f", y))),
            stat = "summary", fun = "median", vjust = -3, hjust = -0.1, size = 6) +
  scale_fill_manual(values = cat.colors) +
  scale_color_manual(values = cat.colors) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Winning/losing strength", title = "Festuca") +
  theme(axis.text.x = element_blank())

pdf(file = "plot/figures_presentation/winner_loser_boxplot_festuca.pdf", width = 5, height = 6)
print(p_box_festuca)
dev.off()

p_box_poa =
  cat.sel %>%
  filter(!is.na(word(species, 2))) %>%
  filter(cat != "nc") %>%
  mutate(genus = word(species, 1)) %>%
  filter(genus == "Poa") %>%
  mutate(cat = factor(cat, levels = c("loser", "winner", "colonizing_winner"))) %>%
  ggplot(aes(cat, estimate, fill = cat, color = cat)) +
  TP_theme() +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_jitter() +
  geom_text(aes(label = after_stat(sprintf("%.2f", y))),
            stat = "summary", fun = "median", vjust = -3, hjust = -0.1, size = 6) +
  scale_fill_manual(values = cat.colors) +
  scale_color_manual(values = cat.colors) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Winning/losing strength", title = "Poa") +
  theme(axis.text.x = element_blank())

pdf(file = "plot/figures_presentation/winner_loser_boxplot_poa.pdf", width = 5, height = 6)
print(p_box_poa)
dev.off()
