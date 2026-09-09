#Showing the method 
x = dfx %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             left_join(metadata %>% dplyr::select(Region, YearEstablished)%>%
                                         st_drop_geometry()%>%
                                         distinct())%>%
                             mutate(Year = Year-YearEstablished)))%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  mutate(comm = map(comm, ~{.} %>% 
                      filter(SpeciesName == "Silene acaulis")%>%
                      dplyr::select(ODT, destPlotID, UniqueID, Year, SpeciesName, Rel_Cover)))%>%
  unnest(comm)

pdf(file = "plot/eg.sp.1.pdf", width = 7, height = 9)
x %>%
  filter(Region == "CH_Calanda", originSiteID == "Pea", destSiteID =="Cal", ODT != "destControls")%>%
  mutate(ODT = ifelse(ODT == "warmed", "Warmed", "Origin controls"))%>%
ggplot(aes(Year, Rel_Cover, group = destPlotID, fill = ODT))+
  TP_theme()+
  facet_wrap(.~ODT, ncol = 1)+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin = 0, ymax = Rel_Cover, fill = ODT), 
               alpha = 0.1)+
  scale_fill_manual(values = c("#648FFF", "#DC267F"))+
  labs(title = "CH_Calanda - Silene acaulis", x = "Experimental years", y = "Relative abundance")+
  theme(legend.position = "none")
dev.off()

wl =
  dfx %>%
  left_join(metadata %>% 
              mutate(YearMax = YearMax-YearEstablished,
                     YearMin = YearMin - YearEstablished) %>%
              dplyr::select(Region, YearMax, YearMin) %>%
              st_drop_geometry()%>%
              distinct()) %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                            left_join(metadata %>% dplyr::select(Region, YearEstablished)%>%
                                        st_drop_geometry()%>%
                                        distinct())%>%
                             mutate(Year = Year-YearEstablished)))
wl = 
  wl %>%
  expand_grid(step = c(2, 3, 5, 7))%>%
  mutate(comm = purrr::map2(comm, step, ~{.x}%>%
                             group_by(Region, originSiteID, destSiteID, ODT, 
                                      destPlotID, SpeciesName)%>%
                             summarize(auc.trap = AUC(Year, Rel_Cover, method = "trapezoid", to = .y),
                                       auc.trap.sqrt = AUC(Year, sqrt(Rel_Cover), method = "trapezoid", to = .y))
                             ))%>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             pivot_longer(cols = auc.trap:auc.trap.sqrt, 
                                          names_to = "rel_abund", values_to = "auc")%>%
                             group_by(Region, originSiteID, destSiteID, ODT, rel_abund)%>%
                             complete(destPlotID, SpeciesName)%>%
                             #mutate(auc = ifelse(is.na(auc), 0, auc))%>%
                             ungroup()%>%
                             dplyr::select(-Region, -destSiteID, -originSiteID)))%>%
  unnest(comm)%>%
  filter(YearMax>=step)%>%
  filter(YearMin<step)%>%
  mutate(step = as.character(step))%>%
  bind_rows(wl %>%
              mutate(step = "max")%>%
              mutate(comm = purrr::map(comm, ~{.x}%>%
                                          group_by(Region, originSiteID, destSiteID, ODT, 
                                                   destPlotID, SpeciesName)%>%
                                          summarize(auc.trap = AUC(Year, Rel_Cover, method = "trapezoid"),
                                                    auc.trap.sqrt = AUC(Year, sqrt(Rel_Cover), method = "trapezoid"))))%>%
              mutate(comm = purrr::map(comm, ~{.}%>%
                                         pivot_longer(cols = auc.trap:auc.trap.sqrt, 
                                                      names_to = "rel_abund", values_to = "auc")%>%
                                         group_by(Region, originSiteID, destSiteID, ODT, rel_abund)%>%
                                         complete(destPlotID, SpeciesName)%>%
                                         #mutate(auc = ifelse(is.na(auc), 0, auc))%>%
                                         ungroup()%>%
                                         dplyr::select(-Region, -destSiteID, -originSiteID)))%>%
              unnest(comm))

pdf(file = "plot/eg.sp.2.pdf", width = 5, height = 4)
wl %>%
  filter(Region == "CH_Calanda", originSiteID == "Pea", destSiteID =="Cal", ODT != "destControls", SpeciesName == "Silene acaulis")%>%
  mutate(ODT = ifelse(ODT == "warmed", "Warmed", "Origin controls"))%>%
  filter(rel_abund == "auc.trap")%>%
  filter(step == "max")%>%
  ggplot(aes(ODT, auc, fill = ODT))+
  TP_theme()+
  geom_boxplot(alpha = 0.7)+
  geom_jitter()+
  scale_fill_manual(values = c("#648FFF", "#DC267F"))+
  labs(title = "CH_Calanda - Silene acaulis", x = "", y = "Area under relative cover")+
  theme(legend.position = "none")
dev.off()

wlx = 
  wl %>%
  nest(.by = c("Region", "originSiteID", "destSiteID", "rel_abund", "step"))%>%
  mutate(lm = purrr::map(data, ~lm(auc~ODT*SpeciesName, data = .x)))

#set contrasts
K = list(between_controls = c(1, -1, 0), to_origin = c(0, -1, 1), to_destination = c(1, 0, -1))

wlc = 
  wlx %>%
  mutate(contrast = purrr::map(lm, ~left_join(as_tibble(contrast(emmeans(.x, ~ODT | SpeciesName), K, adjust = "mvt")),
                                              as_tibble(confint(contrast(emmeans(.x, ~ODT | SpeciesName), K, adjust = "mvt")))
  )))%>%
  dplyr::select(-lm)

wlc= 
  wlc %>%
  dplyr::select(-data)%>%
  unnest(contrast)

wlc =
  wlc %>%
  mutate(pval = ifelse(p.value<0.05, "*", "ns"))%>%
  filter(contrast == "to_origin")%>%
  mutate(wl = ifelse(estimate==0, "constant", ifelse(estimate>0, "winner", "loser")))

wlc = wlc %>% filter(!is.na(p.value)) # later correct the year min for the Region NO, it is coming from that
wlc = wlc %>% mutate(step = as.factor(step),
                     rel_abund = as.factor(rel_abund))

pdf(file = "plot/compare_rel.abund_vs.sqrt.rel.abund.pdf", width = 8, height = 6)
wlc %>%
  filter(step == "max")%>%
  dplyr::select(Region, originSiteID, destSiteID, SpeciesName, 
                rel_abund, estimate)%>%
  pivot_wider(names_from = rel_abund, values_from = estimate)%>%
  ggplot(aes(auc.trap, auc.trap.sqrt))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_vline(xintercept = 0, color = "grey50")+
  geom_abline(slope = 1, intercept = 0, color = "grey50")+
  geom_point(size = 0.5)+
  geom_smooth(method = "lm", color = "red", se = FALSE)+
  stat_cor(method = "pearson", color = "red")+
  labs(x = "Winning/losing strength calculated from \nrelative abundances",
       y = "Winning/losing strength calculated from \nsquare root of relative abundances")
dev.off()

wlt=
wlc %>%
  filter(rel_abund == "auc.trap")%>%
  dplyr::select(Region, originSiteID, destSiteID, SpeciesName, 
                step, estimate)%>%
  mutate(step = paste0("step_", step))%>%
  pivot_wider(names_from = step, values_from = estimate)

p1=
wlt %>%
  dplyr::select(-step_max)%>%
  pivot_longer(cols = c("step_3", "step_5", "step_7"), names_to = "steps", values_to = "auc_value")%>%
  mutate(steps = case_when(steps == "step_3"~ "First 3 years",
                           steps == "step_5"~ "First 5 years",
                           steps == "step_7"~ "First 7 years") 
         )%>%
  ggplot(aes(step_2, auc_value, color = steps))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_vline(xintercept = 0, color = "grey50")+
  geom_abline(slope = 1, intercept = 0, color = "grey50")+
  geom_point(size = 0.5)+
  geom_smooth(method = "lm", se = FALSE)+
  
  stat_cor(method = "pearson", show.legend = FALSE)+
  scale_color_manual(values = c("grey70", "grey40", "grey5"))+
  labs(x = "Winning/losing strength calculated from first 2 years",
       y = "Winning/losing strength calculated from more than 2 years")

p2=
  wlt %>%
  dplyr::select(-step_max, - step_2)%>%
  pivot_longer(cols = c("step_5", "step_7"), names_to = "steps", values_to = "auc_value")%>%
  mutate(steps = case_when(steps == "step_3"~ "First 3 years",
                           steps == "step_5"~ "First 5 years",
                           steps == "step_7"~ "First 7 years") 
  )%>%
  ggplot(aes(step_3, auc_value, color = steps))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_vline(xintercept = 0, color = "grey50")+
  geom_abline(slope = 1, intercept = 0, color = "grey50")+
  geom_point(size = 0.5)+
  geom_smooth(method = "lm", se = FALSE)+
  
  stat_cor(method = "pearson", show.legend = FALSE)+
  scale_color_manual(values = c("grey40", "grey5"))+
  labs(x = "Winning/losing strength calculated from first 3 years",
       y = "Winning/losing strength calculated from more than 3 years")

p3 = 
  ggplot(wlt, aes(step_5, step_7))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_vline(xintercept = 0, color = "grey50")+
  geom_abline(slope = 1, intercept = 0, color = "grey50")+
  geom_point(size = 0.5)+
  geom_smooth(method = "lm", color = "black", se = FALSE)+
  stat_cor(method = "pearson", color = "black")+
  labs(x = "Winning/losing strength calculated from first 5 years",
       y = "Winning/losing strength calculated from first 7 years")

pdf(file = "plot/transient_winning_losing_strength.pdf", width =16, height = 8)
p1+p2+p3
dev.off()

wl.sel = 
  wlc %>%
  filter(rel_abund == "auc.trap")%>%
  filter(step == "max")

cat.sel =
  wl.sel %>%
  left_join(pools)%>%
  mutate(cat = ifelse(pval == "ns", "nc", wl))%>%
  mutate(cat = ifelse(cat == "winner" & pool == "colonizing", "colonizing_winner", cat))%>%
  mutate(cat = as.factor(cat))%>%
  mutate(cat = factor(cat, levels = c("colonizing_winner","winner", "loser", "nc"))) # we set nc as the last one because last level is set as reference in vgam

# Get the corrected species names
cat.sel = 
  cat.sel %>%
  rename(Name_submitted = SpeciesName)%>%
  left_join(sp.list) %>%
  filter(!is.na(Accepted_name))%>%
  rename(species = Accepted_name)%>%
  dplyr::select(-Accepted_species, -Name_submitted)

# Here getting the AUC per treatment
wl.em = 
  wlx %>%
  mutate(contrast = purrr::map(lm, ~left_join(as_tibble(emmeans(.x, ~ODT | SpeciesName)),
                                              as_tibble(confint(emmeans(.x, ~ODT | SpeciesName)))
  )))%>%
  dplyr::select(-lm)%>%
  dplyr::select(-data)%>%
  unnest(contrast)%>%
  filter(step == "max")%>%
  filter(rel_abund == "auc.trap")%>%
  dplyr::select(-SE, -lower.CL, -upper.CL, -df)%>%
  pivot_wider(names_from = "ODT", values_from = "emmean")%>% 
  left_join(wl.sel)

p1 =
wl.em %>%
  ggplot(aes(originControls, warmed, group = SpeciesName, color = pval, alpha = pval))+
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "black")+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  geom_point(size = 2, shape = 21)+
  
  scale_color_manual(values = c("black","grey50"))+
  scale_alpha_manual(values = c(1, 0.5))+
  guides(alpha = "none")+
  labs(x = "Mean area under relative cover curve (Origin controls)",
       y = "Mean area under relative cover curve (Warmed)",
       color = "Significant difference")+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )
p1


pdf(file = "plot/winners_losers.pdf", width = 9, height = 7)
p1
dev.off()

## Understanding why ----
reg = 
  cat.sel %>%
  filter(!is.na(word(species, 2)))%>% 
  left_join(metadata %>% st_drop_geometry()%>%
              dplyr::select(Region, destSiteID, YearRange, PlotSize))%>%
  left_join(clim)%>%
  left_join(ranges)%>%
  left_join(pools.traits)%>%
  left_join(traitx)%>%
  mutate(LDMC = log(LDMC*10),
         vegetative_height = log(vegetative_height*10),
         SLA = log(SLA*10),
         destP = log(destP),
         median_ap_sub = log(median_ap_sub),
         range = log(aoo_2))%>%
  mutate(dist_T = abs(median_mat_sub - destT),
         dist_P = abs(median_ap_sub - destP),
         
         dist_dest_SLA = abs(SLA - dest_SLA),
         dist_dest_LDMC = abs(LDMC - dest_LDMC),
         dist_dest_height = abs(vegetative_height - dest_vegetative_height),
         
         dist_ori_SLA = abs(SLA - ori_SLA),
         dist_ori_LDMC = abs(LDMC - ori_LDMC),
         dist_ori_height = abs(vegetative_height - ori_vegetative_height)
         )%>%
  dplyr::select(Region, originSiteID, destSiteID, estimate, SE, cat, species, YearRange, PlotSize, 
                destP, diffT, 
                SLA, LDMC, vegetative_height, 
                dist_dest_SLA, dist_dest_LDMC, dist_dest_height,
                dist_ori_SLA, dist_ori_LDMC, dist_ori_height,
                range, 
                dist_T, dist_P
                )%>% distinct()

reg = 
  reg %>%
  mutate(across(YearRange:dist_P, ~as.numeric(scale(., center = TRUE, scale = TRUE))))

reg %>%
  dplyr::select(Region, originSiteID, destSiteID, species, SLA:dist_P)%>%
  distinct()%>%
  pivot_longer(cols = SLA:dist_P)%>%
  ggplot(aes(value))+
  facet_wrap(~name, scales = "free")+
  geom_histogram()

reg %>%
  filter(cat != "nc")%>%
  pivot_longer(cols = YearRange:dist_P)%>% 
  ggplot(aes(value, estimate, color = cat))+
  geom_hline(yintercept = 0, color = "grey50")+
  TP_theme()+
  facet_wrap(~name, scales = "free", ncol = 6)+
  geom_point(size = 0.5, alpha = 0.5)+
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = cat.colors)+
  stat_cor(show.legend = FALSE)

# Multinomial model ----
reg =
  reg %>%
  mutate(cat = factor(cat, levels = c("colonizing_winner","winner", "loser", "nc")))

pvar = c("YearRange", "PlotSize", "destP", "diffT", "SLA", 
         "LDMC", "vegetative_height", "dist_dest_SLA", "dist_dest_LDMC", 
         "dist_dest_height", "dist_ori_SLA", "dist_ori_LDMC", "dist_ori_height", 
         "range", "dist_T", "dist_P")
mod.f = as.formula(paste("cat ~", paste0(pvar, collapse = "+")))

mod = vglm(mod.f, family = multinomial, data = reg)

# Summary of the model
modsum = as.data.frame(summary(mod)@coef3)%>%
  rownames_to_column(var = "terms")%>%
  as_tibble()%>%
  rename(estimate = Estimate, se = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`)%>%
  mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
  mutate(
    number = str_extract(terms, "(?<=:)[0-9]"),
    terms = str_extract(terms, "^[^:]+"),
    category = case_when(
      number == "1" ~ "colonizing_winner",
      number == "2" ~ "winner",
      number == "3" ~ "loser",
      TRUE ~ NA_character_
    )
  )%>%
  mutate(ci.low = exp(estimate - 1.96 * se),
         ci.high = exp(estimate + 1.96 * se),
         estimate = exp(estimate))%>%
  mutate(category = factor(category, levels = c("loser", "winner", "colonizing_winner")))

modsum %>%
  filter(terms != "(Intercept)")%>%
  ggplot(aes(terms, estimate, color = category, fill = pval))+
  TP_theme()+
  geom_hline(yintercept = 1, color = "grey30")+
  geom_errorbar(aes(ymax = ci.high, ymin = ci.low), position = position_dodge(width = 0.8))+
  geom_point(shape = 21, size = 2, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c("black", "grey70"))+
  labs(x= "", y = "Odds ratios", fill = "Significance", color = "")+
  scale_color_manual(values = cat.colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pvar = c("YearRange", "PlotSize", "destP", "diffT", "SLA", 
         "LDMC", "vegetative_height", "dist_dest_SLA", "dist_dest_LDMC", 
         "dist_dest_height", "dist_ori_SLA", "dist_ori_LDMC", "dist_ori_height", 
         "range", "dist_T", "dist_P")
mod.f.lmer = as.formula(paste("estimate ~", paste0(pvar, collapse = "+"), "+ (1|species)"))
mod.f = as.formula(paste("estimate ~", paste0(pvar, collapse = "+")))

regall = lm(mod.f, data = reg %>% filter(cat != "nc"))
reglmw = lm(mod.f, data = reg %>% filter(cat == "winner"))
reglml = lm(mod.f, data = reg %>% filter(cat == "loser") %>% mutate(estimate = abs(estimate)))
regcol = lm(mod.f, data = reg %>% filter(cat == "colonizing_winner"))

tidy(reglmw, conf.int = TRUE)%>%
  mutate(cat = "winner")%>%
  bind_rows(tidy(reglml, conf.int = TRUE)%>%
              mutate(cat = "loser"))%>%
  bind_rows(tidy(regcol, conf.int = TRUE)%>%
              mutate(cat = "colonizing_winner"))%>%
  bind_rows(tidy(regall, conf.int = TRUE)%>%
              mutate(cat = "all"))%>%
  mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
  filter(term != "(Intercept)")%>%
  mutate(cat = factor(cat, levels = c("all", "loser", "winner", "colonizing_winner")))%>%
  ggplot(aes(term, estimate, color = cat, shape = pval))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "black")+
  geom_point(position = position_dodge(width = 0.8), size = 3, fill = "white")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), size = 1)+
  scale_color_manual(values = c("all"="black",cat.colors))+
  scale_shape_manual(values = c(16, 21))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tidy(regall, conf.int = TRUE)%>%
  mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
  filter(term != "(Intercept)")%>%
  ggplot(aes(term, estimate, shape = pval))+
  TP_theme()+
  geom_hline(yintercept = 0, color = "black")+
  geom_point(position = position_dodge(width = 0.8), size = 3, fill = "white")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), size = 1)+
  scale_shape_manual(values = c(16, 21))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Left here!

reg %>%
  dplyr::select(estimate, SpeciesName, wl, YearRange, PlotSize, diffT, diffP)%>%
  pivot_longer(cols = YearRange:diffP)%>% 
  ggplot(aes(value, estimate, color = wl))+
  facet_grid(~name)+
  geom_point()+
  stat_smooth(method = "lm")

reglmw = lmer(estimate~diffP + diffT+ PlotSize + YearRange + (1|SpeciesName), data = reg %>% filter(wl == "winner"))
reglml = lmer(estimate~diffP + diffT+ PlotSize + YearRange + (1|SpeciesName), data = reg %>% filter(wl == "loser"))

tidy(reglmw, conf.int = TRUE, "fixed")%>%
  mutate(wl = "winner")%>%
  bind_rows(tidy(reglml, conf.int = TRUE, "fixed")%>%
              mutate(wl = "loser"))%>%
  ggplot(aes(term, estimate, color = wl))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "black")+
  geom_point(position = position_dodge(width = 0.8), size = 3)+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), size = 1)

wl.sum =
  wl.sel %>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  group_by(Region, originSiteID, destSiteID, experiment, wl, pval)%>%
  summarize(n = n())%>%
  group_by(Region, originSiteID, destSiteID, experiment)%>%
  mutate(total = sum(n))%>%
  mutate(ratio = n/total)%>%
  ungroup()

reg.sum = 
  wl.sum %>%
  filter(pval == "*")%>%
  left_join(metadata %>% st_drop_geometry()%>%
              dplyr::select(Region, destSiteID, YearRange, PlotSize) %>%
              mutate(YearRange = as.numeric(scale(YearRange, center = TRUE, scale = TRUE)))%>%
              mutate(PlotSize = as.numeric(scale(PlotSize, center = TRUE, scale = TRUE))))%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))))

reg.sum %>%
  dplyr::select(ratio, wl, YearRange, PlotSize, diffT, diffP)%>%
  pivot_longer(cols = YearRange:diffP)%>% 
  ggplot(aes(value, ratio, color = wl))+
  geom_hline(yintercept = 0, color = "grey50")+
  theme_bw()+
  facet_grid(~name)+
  geom_point()+
  stat_smooth(method = "lm")+
  labs(y = "Ratio of significant number of winners or losers to the number of species")
  
reglmw = lm(ratio~diffP + diffT+ PlotSize + YearRange , data = reg.sum %>% filter(wl == "winner"))
reglml = lm(ratio~diffP + diffT+ PlotSize + YearRange , data = reg.sum %>% filter(wl == "loser"))

tidy(reglmw, conf.int = TRUE)%>%
  mutate(wl = "winner")%>%
  bind_rows(tidy(reglml, conf.int = TRUE)%>%
              mutate(wl = "loser"))%>%
  ggplot(aes(term, estimate, color = wl))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "black")+
  geom_point(position = position_dodge(width = 0.8), size = 3)+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), size = 1)

res.pca = PCA(traitx %>% dplyr::select(C_percent:vegetative_height), graph = FALSE)
fviz_pca_var(res.pca, col.var = "black")

wl.trait = 
  wl.sel %>%
  left_join(traitx %>% 
              rename(SpeciesName = species) %>%
              mutate(across(C_percent:vegetative_height, ~as.numeric(scale(log(.), center = TRUE, scale = TRUE)))))%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))

wl.trait %>%
  pivot_longer(cols = C_percent:vegetative_height, names_to = "trait", values_to = "value")%>%
  filter(!is.na(trait))%>%
  ggplot(aes(value, estimate))+
  theme_bw()+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_smooth(aes(color = wl), method = "lm")+
  facet_grid(.~trait, scales = "free")
 
reglmw = lmerTest::lmer(estimate~C_percent + LDMC+ N_percent + seed_mass + vegetative_height + (1|SpeciesName), data = wl.trait %>% filter(wl == "winner"))
reglml = lmerTest::lmer(estimate~C_percent + LDMC+ N_percent + seed_mass + vegetative_height + (1|SpeciesName), data = wl.trait %>% filter(wl == "loser"))
reglm = lmerTest::lmer(estimate~C_percent + LDMC+ N_percent + seed_mass + vegetative_height + (1|SpeciesName), data = wl.trait)


reglmw = lmerTest::lmer(estimate~C_percent + LA+ LDMC+ N_percent + seed_mass + SLA + vegetative_height + (1|SpeciesName), data = wl.trait %>% filter(wl == "winner"))
reglml = lmerTest::lmer(estimate~C_percent + LA+ LDMC+ N_percent + seed_mass + SLA + vegetative_height + (1|SpeciesName), data = wl.trait %>% filter(wl == "loser"))
reglm = lmerTest::lmer(estimate~C_percent + LA+ LDMC+ N_percent + seed_mass + SLA + vegetative_height + (1|SpeciesName), data = wl.trait)

tidy(reglmw, conf.int = TRUE, "fixed")%>%
  mutate(wl = "winner")%>%
  bind_rows(tidy(reglml, conf.int = TRUE, "fixed")%>%
              mutate(wl = "loser"))%>%
  bind_rows(tidy(reglm, conf.int = TRUE, "fixed")%>%
              mutate(wl = "both"))%>%
  mutate(pval = ifelse(p.value<0.05, "*" , "ns"))%>%
  ggplot(aes(term, estimate, color = wl, alpha = pval))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "black")+
  geom_point(position = position_dodge(width = 0.8), size = 3)+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), size = 1)+
  scale_alpha_manual(values = c(1, 0.3))+
  scale_color_manual(values = c("#3466ffff","#f67b72ff", "#01bfc4ff"))

r.squaredGLMM(reglm)
r.squaredGLMM(reglmw)
r.squaredGLMM(reglml)

fam = read_csv("~/Desktop/TransPlant_clean_trait/data/sp.list_TNRS.csv")%>%
  dplyr::select(Accepted_name, Accepted_family)%>%
  rename(SpeciesName = Accepted_name,
         Family = Accepted_family)

wl.fam =
  wl.sel %>%
  left_join(fam)

wl.fam %>%
  filter(!is.na(Family))%>%
  ggplot(aes(Family, estimate, color = wl))+
  theme_bw()+
  geom_boxplot()+
  coord_flip()

reglmw = lmerTest::lmer(estimate~Family + (1|SpeciesName), data = wl.fam %>% filter(wl == "winner"))
reglml = lmerTest::lmer(estimate~Family + (1|SpeciesName), data = wl.fam %>% filter(wl == "loser"))
reglm = lmerTest::lmer(estimate~Family + (1|SpeciesName), data = wl.fam)

tidy(reglmw, conf.int = TRUE, "fixed")%>%
  mutate(wl = "winner")%>%
  bind_rows(tidy(reglml, conf.int = TRUE, "fixed")%>%
              mutate(wl = "loser"))%>%
  mutate(pval = ifelse(p.value<0.05, "*" , "ns"))%>%
  ggplot(aes(term, estimate, color = wl, alpha = pval))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "black")+
  geom_point(position = position_dodge(width = 0.8), size = 3)+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), size = 1)+
  scale_alpha_manual(values = c(1, 0.3))+
  scale_color_manual(values = c("#f67b72ff", "#01bfc4ff"))+
  coord_flip()


# per experiment ----
exp.tr = 
  wl.trait %>%
  filter(Region != "CN_Damxung")%>%
  filter(Region != "CN_Gongga")%>%
  nest(.by = c("Region", "originSiteID", "destSiteID", "experiment"))%>%
  mutate(lm = purrr::map(data, ~lm(estimate ~ C_percent + LDMC+ N_percent + seed_mass + vegetative_height, data  = .)))%>%
  mutate(res = purrr::map(lm, ~tidy(., confint = TRUE)))%>%
  mutate(r2 = purrr::map(lm, ~summary(.)$r.squared))

exp.tr = 
  wl.trait %>%
  filter(Region != "CN_Damxung")%>%
  filter(Region != "CN_Gongga")%>%
  filter(Region != "CN_Heibei")%>%
  nest(.by = c("Region", "originSiteID", "destSiteID", "experiment", "wl"))%>%
  mutate(lm = purrr::map(data, ~lm(estimate ~ C_percent + LDMC+ N_percent + seed_mass + vegetative_height, data  = .)))%>%
 # mutate(lm = purrr::map(data, ~rma(yi = estimate, sei = SE, mods = ~ C_percent + LDMC+ N_percent + seed_mass + vegetative_height, data  = .)))%>%
  mutate(res = purrr::map(lm, ~tidy(., confint = TRUE)))%>%
  mutate(g = purrr::map(lm, ~glance(.)))

exp.tr %>%
  dplyr::select(-lm, -data, -g)%>%
  unnest(res)%>%
  mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
  filter(term != "(Intercept)")%>%
  filter(estimate<10)%>%
  ggplot(aes(term, estimate, color = wl))+
  theme_bw()+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(aes(alpha = pval))+
  scale_alpha_manual(values = c(1, 0.2))
  
exp.tr %>%
  dplyr::select(-lm, -res, -g)%>%
  unnest(data)%>%
  pivot_longer(cols = C_percent:vegetative_height)%>%
  ggplot(aes(value, estimate))+
  theme_bw()+
  facet_grid(.~name)+
  stat_smooth(aes(group = interaction(name, experiment, wl), color = wl), 
              method = "lm", se = FALSE, fullrange = TRUE)+
  stat_smooth(color = "black", method = "lm", se = FALSE, fullrange = TRUE)+
  geom_jitter(size = 0.4, alpha = 0.1)
  
# Per category ----
cat.tr =
  wl.trait %>%
  left_join(pools)%>%
  mutate(cat = ifelse(pval == "ns", "nc", wl))%>%
  mutate(cat = ifelse(cat == "winner" & pool == "colonizing", "colonizing_winner", cat))%>%
  mutate(cat = as.factor(cat))%>%
  mutate(cat = factor(cat, levels = c("colonizing_winner","winner", "loser", "nc"))) # we set nc as the last one because last level is set as reference in vgam

levels(cat.tr$cat)

mod = vglm(cat ~ C_percent + N_percent + vegetative_height + LDMC + seed_mass, family = multinomial, data = cat.tr)
# Summary of the model
modsum = as.data.frame(summary(mod)@coef3)%>%
  rownames_to_column(var = "terms")%>%
  as_tibble()%>%
  rename(estimate = Estimate, se = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`)%>%
  mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
  mutate(
    number = str_extract(terms, "(?<=:)[0-9]"),
    terms = str_extract(terms, "^[^:]+"),
    category = case_when(
      number == "1" ~ "colonizing_winner",
      number == "2" ~ "winner",
      number == "3" ~ "loser",
      TRUE ~ NA_character_
    )
  )%>%
  mutate(ci.low = exp(estimate - 1.96 * se),
         ci.high = exp(estimate + 1.96 * se),
         estimate = exp(estimate))%>%
  mutate(category = factor(category, levels = c("loser", "winner", "colonizing_winner")))

modsum %>%
  filter(terms != "(Intercept)")%>%
ggplot(aes(terms, estimate, color = category, fill = pval))+
  theme_bw()+
  geom_hline(yintercept = 1, color = "grey30")+
  geom_errorbar(aes(ymax = ci.high, ymin = ci.low), position = position_dodge(width = 0.8))+
  geom_point(shape = 21, size = 2, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c("black", "grey50"))+
  labs(x= "", y = "Odds ratios", fill = "Significance", color = "")
  

# Per category per experiment ----
cat.tr =
  wl.trait %>%
  left_join(pools)%>%
  filter(!Region %in% c("CN_Heibei"))%>%
  filter(!experiment %in% c("CN_Damxung_HIGH_LOW", 
                            "CN_Gongga_A_M", "CN_Gongga_H_A", "CN_Gongga_H_L", "CN_Gongga_M_L",
                            "SE_Abisko_Mid_Low"))%>%
  mutate(cat = ifelse(pval == "ns", "nc", wl))%>%
  mutate(cat = ifelse(cat == "winner" & pool == "colonizing", "colonizing_winner", cat))%>%
  mutate(cat = as.factor(cat))%>%
  mutate(cat = factor(cat, levels = c("colonizing_winner","winner", "loser", "nc")))%>% # we set nc as the last one because last level is set as reference in vgam
  nest(.by = c(Region, originSiteID, destSiteID, experiment))%>%
  mutate(
    mod = purrr::map(
      data, 
      .f = function(df) {
          cat = factor(df$cat, levels = c("colonizing_winner", "winner", "loser", "nc"))
          assign("cat", cat, envir = .GlobalEnv)
        df=vglm(cat ~ C_percent + N_percent + vegetative_height + LDMC + seed_mass, 
               family = multinomial, 
               data = df)
          return(df)
      }
    )
  )%>%
  mutate(modsum = purrr::map(mod, ~as.data.frame(summary(.)@coef3)%>%
                               rownames_to_column(var = "terms")%>%
                               as_tibble()%>%
                               rename(estimate = Estimate, se = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`)%>%
                               mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
                               mutate(
                                 number = str_extract(terms, "(?<=:)[0-9]"),
                                 terms = str_extract(terms, "^[^:]+"),
                                 category = case_when(
                                   number == "1" ~ "colonizing_winner",
                                   number == "2" ~ "winner",
                                   number == "3" ~ "loser",
                                   TRUE ~ NA_character_
                                 )
                               )%>%
                               mutate(ci.low = estimate - 1.96 * se,
                                      ci.high = estimate + 1.96 * se)%>%
                               mutate(category = factor(category, levels = c("loser", "winner", "colonizing_winner")))))
  
modsum = 
  cat.tr %>%
  dplyr::select(-data, -mod)%>%
  unnest(modsum)

modsum %>%
  filter(terms != "(Intercept)")%>%
  filter(!is.na(category))%>%
  ggplot(aes(terms, estimate, color = category))+
  theme_bw()+
  geom_hline(yintercept = 1, color = "grey30")+
  geom_boxplot(outliers = FALSE)+
  geom_point(aes(fill = pval), shape = 21, size = 1, position = position_dodge(width = 0.8))+
  
  scale_fill_manual(values = c("black", "grey50"))+
  labs(x= "", y = "log-odds ratios", fill = "Significance", color = "")

regmod = 
  modsum%>%
  left_join(metadata %>% st_drop_geometry()%>%
              dplyr::select(Region, destSiteID, YearRange, PlotSize) %>%
              mutate(YearRange = as.numeric(scale(YearRange, center = TRUE, scale = TRUE)))%>%
              mutate(PlotSize = as.numeric(scale(PlotSize, center = TRUE, scale = TRUE))))%>%
  left_join(clim %>% mutate(across(destP:cumsumT, ~as.numeric(scale(., center = TRUE, scale = TRUE)))))

regmod %>%
  pivot_longer(cols = YearRange:cumsumT)%>%
  filter(terms != "(Intercept)")%>%
  filter(!is.na(category))%>%
  filter(name %in% c("diffP", "diffT", "PlotSize", "YearRange"))%>%
  ggplot(aes(value, estimate, color = category))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey50")+
  facet_grid(terms~name, scales = "free")+
  stat_smooth(method = "lm", se = FALSE)+
  labs(x = "Value of experimental co-variates (standardized)",
       y = "Log-odds ratio per experiment")

regmod =
  regmod %>%
  filter(terms != "(Intercept)")%>%
  nest(.by = c(terms, category))%>%
  mutate(data = purrr::map(data, ~tidy(lm(estimate~diffT+diffP+YearRange+PlotSize, data = .))))%>%
  unnest(data)

regmod%>%
  filter(term != "(Intercept)")%>%
  filter(!is.na(category))%>%
  mutate(pval = ifelse(p.value>0.05, "ns", "*"))%>%
  mutate(ci.low = estimate - 1.96 * std.error,
         ci.high = estimate + 1.96 * std.error)%>%
  ggplot(aes(term, estimate, color = category, fill = pval))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey30")+
  facet_grid(.~ terms)+
  geom_errorbar(aes(ymax = ci.high, ymin = ci.low), position = position_dodge(width = 0.8))+
  geom_point(position = position_dodge(width = 0.8), shape = 21, size = 2)+
  scale_fill_manual(values = c("black", "grey50"))
