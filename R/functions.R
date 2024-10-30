get.clim = function(sites, climdata){
  clim = 
    climdata %>%
    ungroup()%>%
    filter(var != "pet")%>%
    dplyr::select(gradient, destSiteID, year, month, var, value)%>%
    #filter(!(var == "temperature" & !month %in% c(6:8)))%>%
    #filter(!(var == "precipitation" & !month %in% c(6:8)))%>%
    pivot_wider(values_from = value, names_from = var)%>%
    group_by(gradient, destSiteID, year)%>%
    summarize(temperature = mean(temperature, na.rm = TRUE),
              precipitation = sum(precipitation, na.rm = TRUE))%>%
    rename(Year = year,
           Region = gradient,
           site = destSiteID)%>%
    ungroup()
  
  clim = 
    bind_rows(clim %>% right_join(sites %>% 
                                    ungroup() %>% 
                                    dplyr::select(Region, destSiteID, Year) %>% 
                                    rename(site = destSiteID)%>%
                                    distinct, 
                                  by = c("Region","site", "Year")),
              clim %>% right_join(sites %>% 
                                    ungroup() %>% 
                                    dplyr::select(Region, originSiteID, Year) %>% 
                                    rename(site = originSiteID)%>%
                                    distinct, 
                                  by = c("Region","site", "Year")))%>%
    distinct()%>%
    na.omit()%>% #We do not have climate data for 2020-2021
    ungroup()
  
  destclim =
    clim %>%
    rename(destSiteID = site)%>%
    right_join(sites %>% 
                 ungroup() %>% 
                 dplyr::select(Region, destSiteID, Year) %>% 
                 distinct, 
               by = c("Region","destSiteID", "Year"))%>%
    rename(dest_precipitation = precipitation,
           dest_temperature = temperature)
  
  oriclim =
    clim %>%
    rename(originSiteID = site)%>%
    right_join(sites %>% 
                 ungroup() %>% 
                 dplyr::select(Region, originSiteID, Year) %>% 
                 distinct, 
               by = c("Region","originSiteID", "Year"))%>%
    rename(ori_precipitation = precipitation,
           ori_temperature = temperature)
  
  distclim =
    sites %>% 
    ungroup() %>% 
    dplyr::select(Region, destSiteID, originSiteID, Year) %>% 
    distinct()%>%
    left_join(destclim, by = c("Region", "destSiteID", "Year"))
  
  distclim =
    distclim %>%
    left_join(oriclim %>% dplyr::select(Region, originSiteID, Year, 
                                        ori_precipitation, ori_temperature), by = c("Region", "Year","originSiteID"))
  
  distclim = 
    distclim %>%
    mutate(diff_temperature = dest_temperature-ori_temperature,
           diff_precipitation = dest_precipitation-ori_precipitation)
  
  distclim =
    distclim %>%
    group_by(Region, originSiteID, destSiteID)%>%
    mutate(cumsum_temp = cumsum(diff_temperature),
           cumsum_precipitation = cumsum(diff_precipitation))
  
  distclim =
    distclim%>%
    na.omit()%>%
    group_by(Region, originSiteID, destSiteID)%>%
    summarize(destP = mean(dest_precipitation),
              destT = mean(dest_temperature),
              oriP = mean(ori_precipitation),
              oriT = mean(ori_temperature),
              diffP = mean(diff_precipitation),
              diffT = mean(diff_temperature),
              cumsumP = max(cumsum_precipitation),
              cumsumT = max(cumsum_temp))
  return(distclim)
}