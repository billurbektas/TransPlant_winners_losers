explist = list(CH_Calanda, CH_Calanda2, CH_Lavey, #Switzerland
               CN_Damxung, CN_Gongga, CN_Heibei, #China
               DE_Grainau, DE_Susalps, DE_TransAlps, #Germany
               FR_AlpeHuez, FR_Lautaret, #France
               IN_Kashmir, #India
               IT_MatschMazia1, IT_MatschMazia2, #Italy
               NO_Gudmedalen, NO_Lavisdalen, NO_Skjellingahaugen, NO_Ulvhaugen, #Norway
               SE_Abisko, #Sweden
               US_Arizona, US_Colorado, US_Montana #United States
) 

# Get the metadata per experiment
metadata = map_dfr(explist, ~ .x$meta)
metadata = 
  metadata %>%
  mutate(YearEstablished = ifelse(Gradient == "FR_Lautaret", 2016, YearEstablished))%>% # Correction to Lautaret experiment
  mutate(Elevation = ifelse(Gradient == "CH_Calanda" & Elevation == 2800, 2600, Elevation))%>% # Correction to Calanda experiment
  rename(Region = Gradient,
         PlotSize = PlotSize_m2)%>%
  filter(!is.na(destSiteID))%>% # Correct US_Colorado metadata on cleaning scripts!
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")%>%
  mutate(Climate = extract(kgmap, .)[,2])%>%
  left_join(kglegend)

# Get the gradients to be taken out from the analysis
takeoutregion = metadata %>%
  filter(YearRange < 3) %>%
  distinct(Region) %>%
  pull(Region)

# Get the undetermined taxa to be taken out from the analysis
takeoutsp = cleanedspecies %>%
  filter(is.na(matched_name)) %>%
  distinct(original_name) %>%
  pull(original_name)%>%
  c("Equ pal/pra")

# Re-calculate the relative cover only for vegetation
df = dat %>%
  filter(!Region %in% takeoutregion)%>%
  dplyr::select(-Elevation)%>%
  group_by(Region, Year, originSiteID, originBlockID, destSiteID, destBlockID, destPlotID, Treatment, turfID, UniqueID)%>%
  mutate(Total_Cover = sum(Cover, na.rm = TRUE), Rel_Cover = Cover/Total_Cover)%>%
  ungroup()%>%
  bind_rows(dat %>% filter(Region == "US_Arizona") %>% dplyr::select(-Elevation))%>% # Only relative cover available for Arizona
  filter(!is.na(Rel_Cover)) %>% #4 NAs in DE_Susalps.
  filter(!is.na(SpeciesName))%>%
  filter(!(Region=="DE_Susalps" & Year==2016)) %>% # Not enough data
  filter(!(Region=="CH_Calanda2" & Year==2020))%>% # OriginControl not counted
  filter(!SpeciesName %in% takeoutsp) %>%
  left_join(metadata %>% dplyr::select(Region, destSiteID, YearEstablished, YearMin),
            by = c("Region","destSiteID"))%>%
  mutate(Year_0 = Year-YearEstablished)%>%
  filter(Year_0 != 0)%>% # We do not take pre-transplantation vegetation surveys
  dplyr::select(-(YearEstablished:Year_0))

# Standardize the species names 
cleanedspecies = cleanedspecies %>%
  filter(!(SpeciesName == "Ant alp" & matched_name == "Anthoxanthum odoratum subsp. nipponicum"))%>%
  filter(!(SpeciesName == "Equ pal/pra"))%>%
  filter(!(SpeciesName == "Rum ace subsp lap" & matched_name == "Rumex acetosa subsp. acetosa"))%>%
  filter(!(SpeciesName == "Sal gla" & matched_name == "Salix glauca glauca"))%>%
  filter(!(SpeciesName == "Ver alp" & matched_name == "Veronica alpina subsp. alpina"))

# Get the number of species
cleanedspecies %>%
  filter(!is.na(matched_name))%>%
  dplyr::select(matched_name)%>%
  distinct(matched_name)%>%
  filter(!grepl("_NID_", matched_name))%>%
  mutate(word_count = str_count(matched_name, "\\S+")) %>%
  group_by(word_count) %>%
  summarize(count = n())

df = df %>%
  left_join(cleanedspecies %>% 
              dplyr::select(Region, destSiteID, SpeciesName, matched_name) %>% 
              distinct(), by = c("SpeciesName", "Region", "destSiteID"))%>%
  mutate(matched_name = ifelse(is.na(matched_name), SpeciesName, matched_name))%>%
  dplyr::select(-SpeciesName)%>%
  rename(SpeciesName = matched_name)%>%
  distinct()

# Unresolved botanical issues -> we assign genus names
# CN_Damxung: Stipa capillacea, Stipa capillata -> Stipa
#             Polygonaceae -> drives the pattern too much and all of a sudden disappears to be taken out. 
# DE_Susalps (GW-FE): Festuca pratensis (in destination control plots)? Elymus repens (in origin control plots)? -> drives the pattern too much and all of a sudden disappears to be taken out.
# DE_Susalps (GW-BT): Festuca pratensis (in destination control plots)? Elymus repens (in origin control plots)? -> drives the pattern too much and all of a sudden disappears to be taken out.
df = 
  df %>%
  mutate(SpeciesName = ifelse(Region == "CN_Damxung" & SpeciesName %in% c("Stipa capillacea","Stipa capillata"), "Stipa", SpeciesName))%>%
  mutate(Rel_Cover = ifelse(Region == "CN_Damxung" & SpeciesName == "Polygonaceae", 0,Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "FE" & destSiteID == "FE" &  SpeciesName == "Festuca pratensis", 0, Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "BT" & destSiteID == "BT" & SpeciesName == "Festuca pratensis", 0, Rel_Cover))%>%
  mutate(Rel_Cover = ifelse(Region == "DE_Susalps" & originSiteID == "GW" & destSiteID == "GW" & SpeciesName == "Elymus repens", 0, Rel_Cover))%>%
  filter(Rel_Cover != 0)

# Nest the data per experiment
dfx = df %>%
  dplyr::select(Region, originSiteID, destSiteID, Treatment) %>% 
  distinct() %>% 
  filter(Treatment == "Warm") %>% 
  dplyr::select(-Treatment) %>%
  mutate(comm = pmap(.l = list(R = Region, O = originSiteID, D = destSiteID), .f = function(R, O, D){
    bind_rows(
      originControls = df %>% filter(Region == R & originSiteID == O & destSiteID == O & Treatment == "LocalControl"),
      destControls = df %>% filter(Region == R &  originSiteID == D & destSiteID == D & Treatment == "LocalControl"),
      warmed =  df %>% filter(Region == R & originSiteID == O & destSiteID == D & Treatment == "Warm"),
      .id = "ODT")
  })) 

# I take out species that occurred in communities only once over the whole experimental duration.
## This is because we cannot distinguish if there is an botanical identification error or true presence.
### However, I do not change the relative cover as the species contributed to the total cover. 
dfx = 
  dfx %>%
  mutate(comm = purrr::map(comm, ~{.} %>%
                             group_by(ODT, SpeciesName)%>%
                             mutate(n = n())%>%
                             filter(n>1)%>% 
                             ungroup()%>%
                             dplyr::select(-n)))%>%
  mutate(comm = purrr::map(comm, ~{.}%>%
                             mutate(Region = ifelse(Region == "CH_Calanda2", "CH_Calanda", Region))))%>%
  mutate(Region = case_when(Region == "CH_Calanda2"~"CH_Calanda",
                            .default = Region))

dfx = dfx %>% 
  mutate(comm = purrr::map(comm, ~{
    .x %>% 
      dplyr::select(Region, destSiteID, originSiteID, Treatment,
                    ODT, destPlotID, Year, SpeciesName, Rel_Cover)%>%
      pivot_wider(names_from = SpeciesName, values_from = Rel_Cover, values_fill = list(Rel_Cover = 0), values_fn = list(Rel_Cover = sum)) %>%
      group_by(Region, destSiteID, originSiteID, Treatment,
               ODT)%>%
      complete(Year, destPlotID)%>% # There are cases where plots are missing. Create those plots in the dataset.
      group_by(ODT, Year)%>%
      mutate(across(where(is.numeric), ~replace_na(., mean(., na.rm=TRUE))))%>% #Fill the missing plots with the mean of ODT and Year.
      ungroup()%>%
      pivot_longer(cols = c(-ODT,-Region,-Year,-originSiteID,-destSiteID,-destPlotID,-Treatment),
                   names_to = "SpeciesName",
                   values_to = "Rel_Cover")%>%
      mutate(UniqueID = paste0(Year, "_", destPlotID))
  })) 
