#Showing the method 
x = dfx %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             left_join(metadata %>% dplyr::select(Region, YearEstablished)%>%
                                         st_drop_geometry()%>%
                                         distinct())%>%
                             mutate(Year = Year-YearEstablished)))%>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  mutate(comm = map(comm, ~{.} %>% 
                      filter(SpeciesName == "Festuca nigrescens")%>%
                      dplyr::select(ODT, destPlotID, UniqueID, Year, SpeciesName, Rel_Cover)))%>%
  unnest(comm)

x%>%
  filter(Region == "FR_Lautaret")%>%
ggplot(aes(Year, Rel_Cover, group = destPlotID, fill = ODT))+
  theme_bw()+
  facet_grid(.~ODT)+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin = 0, ymax = Rel_Cover, fill = ODT), 
              alpha = 0.1)+
  scale_fill_manual(values = odt.colors)+
  labs(title = "FR_Lautaret - Festuca nigrescens", x = "Experimental years", y = "Relative abundance")

wl =
  dfx %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                            left_join(metadata %>% dplyr::select(Region, YearEstablished)%>%
                                        st_drop_geometry()%>%
                                        distinct())%>%
                             mutate(Year = Year-YearEstablished)))%>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             group_by(Region, originSiteID, destSiteID, ODT, 
                                      destPlotID, SpeciesName)%>%
                             summarize(auc.trap = AUC(Year, Rel_Cover, method = "trapezoid"),
                                       auc.step = AUC(Year, Rel_Cover, method = "step"),
                                       auc.spline = AUC(Year, Rel_Cover, method = "spline"))
                             ))

wl %>%
  mutate(experiment = paste0(Region, "_", originSiteID, "_", destSiteID))%>%
  filter(Region == "FR_Lautaret")%>%
  mutate(comm = map(comm, ~{.} %>% 
                      filter(SpeciesName == "Festuca nigrescens")%>%
                      ungroup()%>%
                      dplyr::select(ODT, destPlotID, SpeciesName, auc.trap, auc.step, auc.spline)))%>%
  unnest(comm)%>%
  pivot_longer(cols = auc.trap:auc.spline, 
               names_to = "type_auc", values_to = "auc")%>%
  ggplot(aes(ODT, auc, color = ODT))+
  theme_bw()+
  facet_grid(.~type_auc)+
  geom_boxplot()+
  geom_jitter()+
  scale_color_manual(values = odt.colors)+
  labs(title = "FR_Lautaret - Festuca nigrescens", x = "Experimental years", y = "Relative abundance")

  

wlx = 
  wl %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                           pivot_longer(cols = auc.trap:auc.spline, 
                                        names_to = "type_auc", values_to = "auc")%>%
                           group_by(Region, originSiteID, destSiteID, ODT, type_auc)%>%
                           complete(destPlotID, SpeciesName)%>%
                           mutate(auc = ifelse(is.na(auc), 0, auc))%>%
                             ungroup()%>%
                           dplyr::select(-Region, -destSiteID, -originSiteID)))%>%
  unnest(comm)%>%
  nest(.by = c("Region", "originSiteID", "destSiteID", "type_auc"))%>%
  mutate(lm = purrr::map(data, ~lm(auc~ODT*SpeciesName, data = .x)))

#set contrasts
K = list(between_controls = c(1, -1, 0), to_origin = c(0, -1, 1), to_destination = c(1, 0, -1))

wlc = 
  wlx %>%
  mutate(contrast = purrr::map(lm, ~left_join(as_tibble(contrast(emmeans(.x, ~ODT | SpeciesName), K, adjust = "mvt")),
                                              as_tibble(confint(contrast(emmeans(.x, ~ODT | SpeciesName), K, adjust = "mvt")))
  )))%>%
  dplyr::select(-lm)

wlc = 
  wlc %>%
  dplyr::select(-data)%>%
  unnest(contrast)

wlc %>%
  filter(contrast == "to_origin")%>%
  filter(SpeciesName == "Festuca nigrescens")%>%
  filter(type_auc == "auc.trap")%>%
  mutate(pval = ifelse(p.value<0.05, "*" , "ns"))%>%
  mutate(cat = ifelse(estimate>0, "winner", "loser"))%>%
  ggplot(aes(cat, estimate, color = pval))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  geom_boxplot(position = position_dodge(width = 0.8))+
  geom_jitter(size = 3, position = position_dodge(width = 0.8))+
  scale_color_manual(values = c("black","grey50"))+
  labs(title = "Festuca nigrescens", x = "", y = "Winning/losing strength")

wlc %>%
  filter(contrast == "to_origin")%>%
  mutate(pval = ifelse(p.value<0.05, "*" , "ns"))%>%
  mutate(cat = ifelse(estimate>0, "winner", "loser"))%>%
  ggplot(aes(cat, estimate, color = pval))+
  theme_bw()+
  facet_grid(.~type_auc)+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  geom_boxplot(position = position_dodge(width = 0.8))+
  geom_jitter(size = 1, position = position_dodge(width = 0.8))+
  scale_color_manual(values = c("black","grey50"))+
  labs(title = "Festuca nigrescens", x = "", y = "Winning/losing strength")

  
# Calculate the same after 1, 3, 5 years. Is there a difference?
# Sensitivity analysis for auc calculation.
wl.trans =
  dfx %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             left_join(metadata %>% dplyr::select(Region, YearEstablished)%>%
                                         st_drop_geometry()%>%
                                         distinct())%>%
                             mutate(Year = Year-YearEstablished)))%>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             group_by(Region, originSiteID, destSiteID, ODT, 
                                      destPlotID, SpeciesName)%>%
                             summarize(auc.trap.3 = AUC(Year, Rel_Cover, method = "trapezoid", to = 3),
                                       auc.trap.5 = AUC(Year, Rel_Cover, method = "trapezoid", to = 5),
                                       auc.trap.7 = AUC(Year, Rel_Cover, method = "trapezoid", to = 7)
                                       )
  ))
wl.transx = 
  wl.trans %>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             pivot_longer(cols = auc.trap.3:auc.trap.7, 
                                          names_to = "type_auc", values_to = "auc")%>%
                             group_by(Region, originSiteID, destSiteID, ODT, type_auc)%>%
                             complete(destPlotID, SpeciesName)%>%
                             mutate(auc = ifelse(is.na(auc), 0, auc))%>%
                             ungroup()%>%
                             dplyr::select(-Region, -destSiteID, -originSiteID)))%>%
  unnest(comm)%>%
  nest(.by = c("Region", "originSiteID", "destSiteID", "type_auc"))%>%
  mutate(lm = purrr::map(data, ~lm(auc~ODT*SpeciesName, data = .x)))

wl.transc = 
  wl.transx %>%
  mutate(contrast = purrr::map(lm, ~left_join(as_tibble(contrast(emmeans(.x, ~ODT | SpeciesName), K, adjust = "mvt")),
                                              as_tibble(confint(contrast(emmeans(.x, ~ODT | SpeciesName), K, adjust = "mvt")))
  )))%>%
  dplyr::select(-lm)

wl.transc = 
  wl.transc %>%
  dplyr::select(-data)%>%
  unnest(contrast)

wl.transc %>%
  filter(contrast == "to_origin")%>%
  dplyr::select(-(SE:upper.CL))%>%
  pivot_wider(names_from = type_auc, values_from = estimate)%>%
  filter(!(auc.trap.3 == 0 & auc.trap.5 == 0))%>%
  mutate(category = ifelse(auc.trap.5 == 0, "Zero strength 5th year", "non-zero strength 5th year"))%>%
  ggplot(aes(auc.trap.3, auc.trap.5, color = category))+
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "grey80", linetype = "dashed")+
  geom_point()+
  stat_smooth(method = "lm")+
  labs(x = "Wining/losing strength until 3rd year", y = "Winning/losing strength until 5th year")+
  scale_color_manual( values = c("black", "grey50"))
  
wl.transc %>%
  filter(contrast == "to_origin")%>%
  dplyr::select(-(SE:upper.CL))%>%
  pivot_wider(names_from = type_auc, values_from = estimate)%>%
  filter(!(auc.trap.5 == 0 & auc.trap.7 == 0))%>% View
  mutate(category = ifelse(auc.trap.7 == 0, "Zero strength 7th year", "non-zero strength 7th year"))%>%
  ggplot(aes(auc.trap.5, auc.trap.7, color = category))+
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "grey80", linetype = "dashed")+
  geom_point()+
  geom_point()+
  stat_smooth(method = "lm")+
  labs(x = "Wining/losing strength until 5th year", y = "Winning/losing strength until 7th year")+
  scale_color_manual( values = c("black", "grey50"))

  
  
  
  
