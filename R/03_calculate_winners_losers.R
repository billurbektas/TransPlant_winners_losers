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
  filter(type_auc == "auc.spline")%>%
  mutate(pval = ifelse(p.value<0.05, "*" , "ns"))%>%
  mutate(cat = ifelse(estimate>0, "winner", "loser"))%>%
  ggplot(aes(cat, estimate, color = pval))+
  theme_bw()+
  geom_boxplot(position = position_dodge(width = 0.8))+
  geom_jitter(size = 0.2, position = position_dodge(width = 0.8))+
  scale_color_manual(values = c("black","grey50"))
  
# Calculate the same after 1, 3, 5 years. Is there a difference?
# Sensitivity analysis for auc calculation.

