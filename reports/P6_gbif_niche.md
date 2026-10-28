# P6: Niche from GBIF

Species climatic optima from our own GBIF download (GBIF.org (28 October 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.7p28y9), keyed by the 903 P1 study species with their WCVP author and the reviewed ledger `prep/config/gbif_keys.csv`, cleaned with Lynn's rules plus living specimens, thinned to one value per species x CHELSA 30" cell, and summarised as T_opt = mean CHELSA 1981-2010 bio10 and P_opt = log mean bio18. Replaces Lynn's niche file (decided 28 Oct). The download record is `data/derived/prep/gbif_download.csv` and the exact request `reports/P6_gbif_download_request.json`.

## The download

| download_key | doi | total_records | size_GB | created | erase_after | sha256 |
| --- | --- | --- | --- | --- | --- | --- |
| 0011587-261023074918325 | 10.15468/dl.7p28y9 | 135,294,830 | 15.2 | 2026-10-28 | 2027-04-28 | 7e4633cd1e184f78 |

GBIF deletes downloads nobody cites after `erase_after`; the DOI goes into the methods and data availability text, so the download is kept.

## Keys

| names | class | n |
| --- | --- | --- |
| accepted name | automatic: accepted | 854 |
| accepted name | automatic: provisionally_accepted | 7 |
| accepted name | automatic: synonym | 15 |
| accepted name | ledger: GBIF files it under another species | 24 |
| accepted name | ledger: no exact match | 3 |
| other community names | ignored: no exact match | 24 |
| other community names | ledger: GBIF keeps it as another species | 21 |
| other community names | nothing to do: same GBIF species | 123 |

903 species keyed with 929 GBIF keys; 0 excluded.

Same-named usages added (7 species): GBIF holds our author's usage only provisionally or as a synonym and accepts another usage of the same name, so the records of the name are split between the two keys; both are requested:

| species | gbif_key | gbif_name | records |
| --- | --- | --- | --- |
| Anemone cathayensis | 7698129 | Anemone cathayensis | 342 |
| Artemisia tridentata | 9396703 | Artemisia tridentata | 22943 |
| Bistorta officinalis | 7925535 | Bistorta officinalis | 122527 |
| Huperzia selago | 8190643 | Huperzia selago | 44743 |
| Vaccinium membranaceum | 9060377 | Vaccinium membranaceum | 6644 |
| Valeriana edulis | 2888746 | Valeriana edulis | 3075 |
| Viscaria vulgaris | 8153621 | Viscaria vulgaris | 37677 |

### Ledger decisions (49)

| species | name | author | gbif_key | gbif_name | decision | reason |
| --- | --- | --- | --- | --- | --- | --- |
| Alopecurus gerardi | Alopecurus gerardi |  |  |  | accepted | spelling only: GBIF writes gerardii |
| Anemonoides nemorosa | Anemonoides nemorosa |  |  |  | accepted | same plant; GBIF keeps it in Anemone |
| Arctous alpina | Arctous alpina |  |  |  | accepted | same plant; GBIF keeps it in Arctostaphylos |
| Aster asteroides | Aster asteroides |  |  |  | synonym_only | GBIF has the name only as a synonym of the genus Tibetiodes: records under our name only |
| Aster bellidiastrum | Aster bellidiastrum |  |  |  | accepted | same plant; GBIF name Bellidiastrum michelii |
| Aster diplostephioides | Aster diplostephioides |  |  |  | accepted | same plant; GBIF name Tibetiodes diplostephioides |
| Berberis repens | Berberis repens |  |  |  | accepted | same plant; GBIF name Mahonia repens |
| Coluria longifolia | Coluria longifolia |  |  |  | accepted | same plant; GBIF name Geum longifolium |
| Crepis conyzifolia | Crepis conyzifolia |  |  |  | synonym_only | GBIF sinks it into Crepis blattarioides, which WCVP keeps apart: records under our name only |
| Draba muralis | Draba muralis |  |  |  | accepted | same plant; GBIF name Drabella muralis |
| Epilobium angustifolium | Epilobium angustifolium |  |  |  | accepted | same plant; GBIF name Chamaenerion angustifolium |
| Hymenidium hookeri | Hymenidium hookeri |  |  |  | accepted | same plant; GBIF name Pleurospermum hookeri |
| Hypochaeris maculata | Hypochaeris maculata |  |  |  | accepted | same plant; GBIF name Trommsdorffia maculata |
| Lithophragma glabrum | Lithophragma glabrum |  |  |  | synonym_only | GBIF: ambiguous synonym of Lithophragma tenella, which WCVP keeps apart: records under our name only |
| Microsteris gracilis | Microsteris gracilis |  |  |  | accepted | same plant; GBIF name Phlox gracilis |
| Morina chinensis | Morina chinensis |  |  |  | accepted | same plant; GBIF name Cryptothladia chinensis |
| Morina kokonorica | Morina kokonorica |  |  |  | accepted | same plant; GBIF name Cryptothladia kokonorica |
| Morina nepalensis | Morina nepalensis |  |  |  | accepted | same plant; GBIF name Acanthocalyx nepalensis |
| Mutellina purpurea | Mutellina purpurea |  |  |  | accepted | same plant (Ligusticum mutellina, recorded as such); GBIF name Mutellina adonidifolia |
| Neogaya simplex | Neogaya simplex |  |  |  | accepted | same plant (Ligusticum mutellinoides, recorded as such); GBIF name Pachypleurum mutellinoides |
| Orobanche fasciculata | Orobanche fasciculata |  |  |  | accepted | same plant; GBIF name Aphyllon fasciculatum |
| Pulsatilla alpina | Pulsatilla alpina |  | 7278071 | Pulsatilla alpina (L.) Delarbre | accepted | no match with the WCVP author (L.) Chaz.; GBIF's accepted usage (L.) Delarbre is the same plant |
| Ranunculus tuberosus | Ranunculus tuberosus |  | 7510191 | Ranunculus polyanthemos subsp. nemorosus (DC.) Schübl. & G.Martens | gbif_subtaxon | GBIF files R. tuberosus Lapeyr. under R. polyanthemos subsp. nemorosus, with R. nemorosus DC. and its other synonyms; WCVP gives all these names to R. tuberosus, and R. polyanthemos is itself a study species. The whole GBIF subspecies is taken and its records are assigned to R. tuberosus by their own name key, so they are carved out of R. polyanthemos (decided 28 Oct; see the P6 report) |
| Rhinanthus major | Rhinanthus major |  |  |  | synonym_only | GBIF: ambiguous synonym of Rhinanthus serotinus (the name has been applied to more than one plant): records under our name only |
| Rubus parviflorus | Rubus parviflorus |  |  |  | accepted | same plant (thimbleberry); GBIF name Rubus nutkanus |
| Sanguisorba minor | Sanguisorba minor |  |  |  | accepted | same plant; GBIF name Poterium sanguisorba |
| Valeriana locusta | Valeriana locusta |  |  |  | accepted | same plant; GBIF name Valerianella locusta |
| Festuca quadriflora | Festuca quadriflora | Honck. |  |  | author | the alpine grass of CH_Calanda is Honckeny's; TNRS returned Walter's American homonym |
| Pilosella piloselloides | Hieracium piloselloides |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Scorzoneroides autumnalis | Leontodon autumnalis |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Alpagrostis alpina | Agrostis alpina |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Helianthemum oelandicum | Helianthemum alpestre |  |  |  | include | WCVP includes H. alpestre in H. oelandicum; recorded under this name at the sites |
| Omalotheca sylvatica | Gnaphalium sylvaticum |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Viola biflora | Viola biflora var. rockiana |  |  |  | include | recorded at CN_Gongga; WCVP places var. rockiana in V. biflora, GBIF treats it as Viola rockiana |
| Persicaria runcinata | Polygonum runcinatum |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Pilosella hoppeana | Hieracium hoppeanum |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Noccaea fendleri | Thlaspi montanum |  |  |  | exclude | GBIF's usage is Thlaspi montanum Bory & Chaub. (Greece); the US_Montana record is N. fendleri by P1 decision |
| Omalotheca norvegica | Gnaphalium norvegicum |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Oxytropis humifusa | Oxytropis glacialis |  |  |  | include | WCVP (P1) places O. glacialis in O. humifusa, GBIF keeps O. glacialis Benth. apart; recorded under this name at CN_Damxung |
| Koeleria macrantha | Koeleria cristata |  |  |  | exclude | GBIF's Koeleria cristata Kar. & Kir. collects records of several plants labelled K. cristata (China, Russia, Ukraine, North America, France) |
| Erigeron sumatrensis | Conyza sumatrensis |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Centaurea phrygia | Centaurea pseudophrygia |  |  |  | include | WCVP includes C. pseudophrygia in C. phrygia; recorded under this name at the sites |
| Noccaea caerulescens | Noccaea brachypetala |  |  |  | include | WCVP: N. caerulescens subsp. brachypetala; recorded under this name at FR_AlpeHuez |
| Hieracium lachenalii | Hieracium vulgatum |  |  |  | include | H. vulgatum Fr. is H. lachenalii (WCVP); GBIF files the name under Sprengel's usage; recorded in the Norwegian regions |
| Empetrum nigrum | Empetrum hermaphroditum |  |  |  | include | WCVP includes E. hermaphroditum in E. nigrum; recorded under this name at the sites |
| Valeriana excelsa | Valeriana sambucifolia |  |  |  | include | WCVP: V. excelsa subsp. sambucifolia; recorded under this name at NO_Lavisdalen |
| Muhlenbergia tricholepis | Blepharoneuron tricholepis |  |  |  | include | same plant under its older genus, which GBIF keeps as a separate species |
| Dysphania incisa | Dysphania graveolens |  |  |  | include | WCVP places D. graveolens in D. incisa; recorded under this name at the sites |
| Symphyotrichum ascendens | Aster adscendens |  |  |  | include | Aster adscendens Lindl. is Symphyotrichum ascendens (P1 decision 28 Oct); GBIF keeps the name as a separate species with US records |

### Ranunculus tuberosus and Ranunculus polyanthemos (`gbif_subtaxon`, decided 28 Oct)

GBIF and WCVP cut these two plants differently. GBIF files R. tuberosus Lapeyr. as a synonym of R. polyanthemos subsp. nemorosus, under which it also files R. nemorosus DC. and the other names below; WCVP (our P1 taxonomy) gives every one of these names to R. tuberosus and keeps R. polyanthemos for the rest. Taking only the records labelled "R. tuberosus" would have given R. tuberosus a minority of its own records and left the rest in R. polyanthemos, which takes its whole GBIF species. So GBIF's subspecies nemorosus is requested, and every record filed under it, by any of its names, is assigned to R. tuberosus by the record's own name key and so carved out of R. polyanthemos. The script stops if a record under the subspecies carries a name key not listed here. Records with coordinates, no geospatial issue, 1950 onwards, before cleaning:

| species | gbif_name | gbif_key | how | records |
| --- | --- | --- | --- | --- |
| Ranunculus tuberosus | Ranunculus polyanthemos subsp. nemorosus (DC.) Schübl. & G.Martens | 7510191 | ledger: gbif_subtaxon | 59,315 |
| Ranunculus tuberosus | Ranunculus nemorosus DC. | 7993540 | ledger: gbif_subtaxon, name filed under it | 28,496 |
| Ranunculus tuberosus | Ranunculus tuberosus Lapeyr. | 3922138 | ledger: gbif_subtaxon, name filed under it | 28,335 |
| Ranunculus tuberosus | Ranunculus serpens subsp. nemorosus (DC.) G.López | 3925172 | ledger: gbif_subtaxon, name filed under it | 719 |
| Ranunculus tuberosus | Ranunculus breyninus var. nemorosus (DC.) P.Fourn. | 9436706 | ledger: gbif_subtaxon, name filed under it | 168 |
| Ranunculus tuberosus | Ranunculus nemorosus subsp. nemorosus | 7276642 | ledger: gbif_subtaxon, name filed under it | 115 |
| Ranunculus tuberosus | Ranunculus polyanthemos subsp. nemorosus (DC.) Korsh. | 6096573 | ledger: gbif_subtaxon, name filed under it | 23 |
| Ranunculus tuberosus | Ranunculastrum nemorosum (DC.) Fourr. | 3928412 | ledger: gbif_subtaxon, name filed under it | 5 |
| Ranunculus tuberosus | Ranunculus nemorosus subsp. delacourii (Labille & Gaudefroy) P.Fourn. | 9260184 | ledger: gbif_subtaxon, name filed under it | 5 |
| Ranunculus tuberosus | Ranunculus aureus Schleich. | 3925851 | ledger: gbif_subtaxon, name filed under it | 4 |
| Ranunculus tuberosus | Ranunculus mixtus Jord. | 3922012 | ledger: gbif_subtaxon, name filed under it | 4 |
| Ranunculus tuberosus | Ranunculus acris Willk. ex Freyn | 7360576 | ledger: gbif_subtaxon, name filed under it | 2 |
| Ranunculus tuberosus | Ranunculus lanuginosus Amo | 7980200 | ledger: gbif_subtaxon, name filed under it | 2 |
| Ranunculus tuberosus | Ranunculus nemorosus subsp. amansii (Jord.) Debeaux, 1898 | 12338601 | ledger: gbif_subtaxon, name filed under it | 2 |

The first row counts all names together: 59,315 records go to R. tuberosus, of which 28,335 are labelled "Ranunculus tuberosus Lapeyr." (the only ones a name-only key would have taken).

### Species whose names had more than one legitimate author in TNRS (232)

| species | author | records | in_study_countries |
| --- | --- | --- | --- |
| Achillea millefolium | L. | 1,542,438 | 596,344 |
| Agrostis capillaris | L. | 900,078 | 293,567 |
| Agrostis gigantea | Roth | 106,071 | 5,259 |
| Agrostis rupestris | All. | 24,638 | 12,369 |
| Agrostis stolonifera | L. | 850,326 | 25,632 |
| Ajuga pyramidalis | L. | 73,352 | 18,359 |
| Alchemilla monticola | Opiz | 58,300 | 10,395 |
| Alchemilla vulgaris | L. | 71,981 | 33,033 |
| Alopecurus pratensis | L. | 476,676 | 50,237 |
| Alpagrostis alpina | (Scop.) P.M.Peterson, Romasch., Soreng & Sylvester | 13,170 | 2,962 |
| Androsace chamaejasme | Wulfen | 9,549 | 4,772 |
| Androsace mariae | Kanitz | 460 | 460 |
| Androsace obtusifolia | All. | 6,268 | 2,541 |
| Anemonastrum narcissiflorum | (L.) Holub | 21,651 | 5,758 |
| Anemonastrum obtusilobum | (D.Don) Mosyakin | 131 | 88 |
| Anemonoides nemorosa | (L.) Holub | 638,349 | 35,322 |
| Antennaria microphylla | Rydb. | 2,529 | 2,216 |
| Anthoxanthum nipponicum | Honda | 39,453 | 5,719 |
| Arabidopsis thaliana | (L.) Heynh. | 224,501 | 71,376 |
| Arabis hirsuta | (L.) Scop. | 88,365 | 40,994 |
| Arctous alpina | (L.) Nied. | 32,933 | 12,915 |
| Arenaria ciliata | L. | 10,118 | 2,928 |
| Arenaria serpyllifolia | L. | 243,742 | 16,507 |
| Argentina anserina | (L.) Rydb. | 621,960 | 1,331 |
| Argentina leuconota | (D.Don) Soják | 555 | 52 |
| Argentina stenophylla | (Franch.) Soják | 246 | 245 |
| Arnica angustifolia | Vahl | 3,741 | 223 |
| Arnica montana | L. | 141,048 | 57,980 |
| Aster alpinus | L. | 41,566 | 28,644 |
| Aster asteroides | (DC.) Kuntze | 70 | 47 |
| Aster bellidiastrum | (L.) Scop. | 68,873 | 24,862 |
| Aster tataricus | L.f. | 1,849 | 575 |
| Astragalus alpinus | L. | 29,660 | 10,012 |
| Astragalus mongholicus | Bunge | 902 | 145 |
| Atriplex argentea | Nutt. | 1,059 | 1,006 |
| Avenula pubescens | (Huds.) Dumort. | 174,971 | 98,264 |
| Bartsia alpina | L. | 101,397 | 50,073 |
| Bistorta macrophylla | (D.Don) Soják | 1,285 | 1,203 |
| Bistorta officinalis | Delarbre | 197,930 | 38,642 |
| Bistorta vivipara | (L.) Delarbre | 200,941 | 143,050 |
| Bouteloua gracilis | (Kunth) Lag. ex Griffiths | 11,719 | 8,302 |
| Bromus carinatus | Hook. & Arn. | 8,969 | 3,746 |
| Bromus erectus | Huds. | 375,737 | 301,874 |
| Bromus hordeaceus | L. | 473,655 | 20,967 |
| Bromus inermis | Leyss. | 91,030 | 31,577 |
| Bromus tectorum | L. | 74,874 | 17,877 |
| Calamagrostis lapponica | (Wahlenb.) Hartm. | 13,873 | 2,751 |
| Campanula barbata | L. | 46,783 | 31,269 |
| Campanula rhomboidalis | L. | 32,925 | 24,117 |
| Campanula scheuchzeri | Vill. | 74,251 | 31,948 |
| Carduus crispus | L. | 149,198 | 7,754 |
| Carex atrata | L. | 17,578 | 6,011 |
| Carex atrofusca | Schkuhr | 8,065 | 334 |
| Carex caryophyllea | Latourr. | 160,621 | 82,341 |
| Carex curvula | All. | 22,355 | 1 |
| Carex echinata | Murray | 233,424 | 18,915 |
| Carex ferruginea | Scop. | 34,011 | 5,450 |
| Carex lehmannii | Drejer | 230 | 215 |
| Carex leporina | L. | 248,112 | 11,565 |
| Carex montana | L. | 87,304 | 27,528 |
| Carex myosuroides | Vill. | 28,045 | 20,496 |
| Carex nigra | (L.) Reichard | 498,795 | 30,114 |
| Carex norvegica | Retz. | 6,939 | 4,087 |
| Carex panicea | L. | 376,747 | 17,646 |
| Carex parviflora | Host | 11,975 | 2,604 |
| Carex pilulifera | L. | 265,162 | 13,027 |
| Carex rupestris | All. | 13,767 | 6,235 |
| Carex saxatilis | L. | 13,302 | 5,586 |
| Carex sempervirens | Vill. | 99,326 | 79,188 |
| Carex spicata | Huds. | 90,704 | 33,463 |
| Carex sylvatica | Huds. | 386,866 | 197,997 |
| Carlina acaulis | L. | 138,687 | 70,249 |
| Centaurea jacea | L. | 532,441 | 116,393 |
| Centaurea montana | L. | 77,698 | 14,175 |
| Centaurea nervosa | Willd. | 8,327 | 3,865 |
| Centaurea nigra | L. | 278,627 | 1,278 |
| Centaurea phrygia | L. | 22,109 | 55 |
| Centaurea scabiosa | L. | 270,306 | 104,894 |
| Centaurea uniflora | Turra | 10,266 | 10,038 |
| Cerastium alpinum | L. | 24,356 | 8,757 |
| Cerastium arvense | L. | 161,836 | 49,890 |
| Cerastium holosteoides | Fr. | 433,193 | 56,058 |
| Cherleria biflora | (L.) A.J.Moore & Dillenb. | 7,108 | 2,536 |
| Cherleria sedoides | L. | 18,258 | 3,495 |
| Cirsium arvense | (L.) Scop. | 1,363,883 | 48,928 |
| Cirsium palustre | (L.) Scop. | 795,451 | 16,697 |
| Clinopodium alpinum | (L.) Kuntze | 51,009 | 26,780 |
| Collomia linearis | Nutt. | 6,177 | 4,342 |
| Crepis biennis | L. | 146,477 | 26,186 |
| Crocus vernus | (L.) Hill | 51,024 | 34,967 |
| Dichodon cerastoides | (L.) Rchb. | 15,426 | 4,125 |
| Draba aizoides | L. | 24,514 | 5,513 |
| Draba muralis | L. | 18,058 | 9,683 |
| Draba nemorosa | L. | 12,514 | 773 |
| Dryas octopetala | L. | 110,258 | 30,734 |
| Drymocallis glandulosa | (Lindl.) Rydb. | 5,610 | 5,144 |
| Dysphania incisa | (Poir.) ined. | 1,719 | 918 |
| Elymus caninus | (L.) L. | 111,243 | 8,151 |
| Elymus glaucus | Buckley | 6,975 | 5,765 |
| Epilobium lactiflorum | Hausskn. | 4,389 | 656 |
| Epilobium obscurum | Schreb. | 72,016 | 2,225 |
| Eremogone kansuensis | (Maxim.) Dillenb. & Kadereit | 320 | 315 |
| Eremogone roborowskii | (Maxim.) Rabeler & W.L.Wagner | 23 | 23 |
| Erigeron annuus | (L.) Desf | 395,023 | 61,639 |
| Erigeron canadensis | L. | 463,181 | 7 |
| Erigeron divergens | Torr. & A.Gray | 5,339 | 4,797 |
| Erigeron uniflorus | L. | 19,126 | 3,586 |
| Euphorbia maculata | L. | 62,762 | 24,287 |
| Euphorbia prolifera | Buch.-Ham. ex D.Don | 80 | 60 |
| Euphorbia serpillifolia | Pers. | 4,679 | 4,002 |
| Euphrasia minima | Jacq. ex DC. | 18,746 | 14,765 |
| Euphrasia stricta | J.P.Wolff ex J.F.Lehm. | 82,990 | 4,104 |
| Eurybia conspicua | (Lindl.) G.L.Nesom | 2,916 | 400 |
| Festuca alpina | Suter | 3,017 | 1,749 |
| Festuca halleri | All. | 9,345 | 9,127 |
| Festuca laevigata | Gaudin | 29,705 | 29,507 |
| Festuca nigrescens | Lam. | 30,374 | 23,193 |
| Festuca quadriflora | Honck. | 19,362 | 3,658 |
| Galium album | Mill. | 285,016 | 32,212 |
| Galium pumilum | Murray | 77,740 | 63,593 |
| Galium uliginosum | L. | 194,907 | 6,510 |
| Gentiana crassuloides | Bureau & Franch. | 245 | 232 |
| Gentiana lutea | L. | 116,263 | 104,127 |
| Gentiana nivalis | L. | 26,687 | 12,789 |
| Gentiana prostrata | Haenke | 1,968 | 245 |
| Gentiana punctata | L. | 21,696 | 6,307 |
| Gentiana purpurea | L. | 25,648 | 14,754 |
| Gentianella amarella | (L.) Börner | 33,380 | 2,270 |
| Gentianella campestris | (L.) Börner | 74,485 | 34,774 |
| Helianthella uniflora | Torr. & A.Gray | 1,495 | 1,355 |
| Helianthus petiolaris | Nutt. | 5,070 | 0 |
| Helictochloa versicolor | (Vill.) Romero Zarco | 26,843 | 22,924 |
| Hieracium alpinum | L. | 12,861 | 1,954 |
| Hieracium triste | Willd. ex Spreng. | 1,850 | 1,304 |
| Hieracium umbellatum | L. | 191,934 | 9,223 |
| Hippocrepis comosa | L. | 204,844 | 158,838 |
| Hypericum maculatum | Crantz | 206,702 | 75,996 |
| Hypericum richeri | Vill. | 17,373 | 15,653 |
| Jacobaea alpina | (L.) Moench | 13,109 | 4,347 |
| Jacobaea incana | (L.) Veldkamp | 13,893 | 9,664 |
| Jacobaea vulgaris | Gaertn. | 748,252 | 21,403 |
| Juniperus communis | L. | 647,947 | 50,553 |
| Knautia dipsacifolia | Kreutzer | 89,235 | 24,902 |
| Knorringia sibirica | (Laxm.) Tzvelev | 713 | 260 |
| Koenigia cyanandra | (Diels) Měsíček & Soják | 177 | 177 |
| Leontodon hispidus | L. | 325,116 | 184,305 |
| Leontodon hyoseroides | Welw. ex Rchb. | 7,192 | 1,517 |
| Linaria vulgaris | Mill. | 474,103 | 105,512 |
| Linum lewisii | Pursh | 12,370 | 9,016 |
| Lolium arundinaceum | (Schreb.) Darbysh. | 343,908 | 13,990 |
| Lolium multiflorum | Lam. | 111,778 | 7,311 |
| Lolium pratense | (Huds.) Darbysh. | 272,638 | 51,313 |
| Lycoris radiata | (L'Hér.) Herb. | 7,432 | 802 |
| Medicago ruthenica | (L.) Trautv. | 1,145 | 865 |
| Melica spectabilis | Scribn. | 524 | 493 |
| Milium effusum | L. | 211,493 | 16,946 |
| Noccaea caerulescens | (J.Presl & C.Presl) F.K.Mey. | 42,492 | 14,559 |
| Omalotheca supina | (L.) DC. | 30,561 | 7,838 |
| Orchis mascula | (L.) L. | 236,843 | 104,515 |
| Ornithogalum umbellatum | L. | 123,734 | 18,944 |
| Osmorhiza berteroi | DC. | 5,779 | 4,068 |
| Oxytropis campestris | (L.) DC. | 21,764 | 4,236 |
| Oxytropis humifusa | Kar. & Kir. | 310 | 77 |
| Pedicularis comosa | L. | 5,888 | 5,071 |
| Persicaria runcinata | (Buch.-Ham. ex D.Don) H.Gross | 2,564 | 2,073 |
| Phyteuma scheuchzeri | All. | 3,352 | 9 |
| Pilosella aurantiaca | (L.) F.W.Schultz & Sch.Bip. | 152,717 | 11,598 |
| Pilosella glacialis | (Reyn. ex Lachen.) F.W.Schultz & Sch.Bip. | 4,453 | 2,622 |
| Pilosella hoppeana | (Schult.) F.W.Schultz & Sch.Bip. | 5,401 | 3,772 |
| Pilosella lactucella | (Wallr.) P.D.Sell & C.West | 70,083 | 16,582 |
| Pilosella piloselloides | (Vill.) Soják | 28,430 | 11,585 |
| Pimpinella major | (L.) Huds. | 141,997 | 40,204 |
| Pinus sylvestris | L. | 980,721 | 54,207 |
| Plantago maritima | L. | 137,793 | 8,555 |
| Plantago patagonica | Jacq. | 8,527 | 7,140 |
| Poa alpina | L. | 117,548 | 87,891 |
| Poa annua | L. | 754,782 | 15,810 |
| Poa compressa | L. | 107,355 | 4,125 |
| Poa reflexa | Vasey & Scribn. | 496 | 494 |
| Poa trivialis | L. | 879,285 | 309,716 |
| Potentilla aurea | L. | 57,647 | 45,004 |
| Potentilla clandestina | Soják | 30 | 30 |
| Potentilla erecta | (L.) Raeusch. | 866,222 | 352,163 |
| Potentilla grandiflora | L. | 21,120 | 17,201 |
| Potentilla nivea | L. | 3,834 | 110 |
| Primula serratifolia | Franch. | 59 | 59 |
| Pulsatilla alpina | (L.) Chaz. | 64,055 | 49,340 |
| Quercus robur | L. | 1,516,952 | 222,023 |
| Ranunculus tuberosus | Lapeyr. | 59,315 | 31,138 |
| Rhinanthus minor | L. | 276,296 | 41,611 |
| Rhodiola integrifolia | Raf. | 6,831 | 5,068 |
| Rubus parviflorus | Nutt. | 2,558 | 2,117 |
| Rumex arifolius | All. | 22,278 | 1,725 |
| Sabulina verna | (L.) Rchb. | 27,886 | 4,415 |
| Salix retusa | L. | 26,259 | 6,476 |
| Sanguisorba minor | Scop. | 512,768 | 51,696 |
| Saxifraga paniculata | Mill. | 65,166 | 13,472 |
| Scorzoneroides pyrenaica | (Gouan) Holub | 23,907 | 23,097 |
| Senecio doronicum | (L.) L. | 32,447 | 24,871 |
| Senecio integerrimus | Nutt. | 8,425 | 7,355 |
| Sibbaldia procumbens | L. | 34,727 | 21,067 |
| Silene latifolia | Poir. | 454,729 | 35,856 |
| Silene vulgaris | (Moench) Garcke | 440,105 | 212,338 |
| Solidago missouriensis | Nutt. | 3,944 | 2,973 |
| Solidago multiradiata | Aiton | 5,771 | 3,149 |
| Sonchus asper | (L.) Hill | 554,601 | 21,475 |
| Sonchus oleraceus | L. | 523,046 | 126,693 |
| Sorbus aucuparia | L. | 1,074,636 | 99,117 |
| Stellaria alsine | Grimm | 175,379 | 180 |
| Stellaria irrigua | Bunge | 991 | 684 |
| Symphyotrichum ascendens | (Lindl.) G.L.Nesom | 2,601 | 2,559 |
| Taraxacum officinale | F.H.Wigg. | 523,802 | 171,471 |
| Tragopogon orientalis | L. | 50,083 | 4,425 |
| Trichophorum cespitosum | (L.) Hartm. | 137,454 | 22,721 |
| Trifolium alpestre | L. | 33,520 | 10,724 |
| Trifolium alpinum | L. | 55,568 | 51,096 |
| Trifolium campestre | Schreb. | 284,781 | 6,574 |
| Vaccinium membranaceum | Douglas ex Torr. | 7,606 | 5,182 |
| Vaccinium myrtillus | L. | 885,035 | 411,846 |
| Valeriana excelsa | Poir. | 158,225 | 28,561 |
| Veronica allionii | Vill. | 12,743 | 12,564 |
| Veronica alpina | L. | 35,857 | 15,384 |
| Veronica ciliata | Fisch. | 661 | 561 |
| Veronica filiformis | Sm. | 79,400 | 6,762 |
| Veronica montana | L. | 122,755 | 7,159 |
| Veronica serpyllifolia | L. | 307,838 | 24,596 |
| Vicia americana | Muhl. ex Willd. | 19,017 | 13,062 |
| Vicia tetrasperma | (L.) Schreb. | 104,221 | 7,210 |
| Vicia villosa | Roth | 92,994 | 7,810 |
| Viola reichenbachiana | Jord. ex Boreau | 219,177 | 19,069 |
| Viscaria alpina | (L.) G.Don | 19,624 | 4,291 |
| Viscaria vulgaris | Bernh. | 83,143 | 9,729 |

Geographic check: 889 of 903 keyed species have at least 10 GBIF records in their study countries; 14 confirmed in the ledger.

## Cleaning and thinning

| step | records |
| --- | --- |
| read | 135,294,830 |
| unassigned | 0 |
| fossil or living specimen | 274,484 |
| absent | 4,846,577 |
| uncertainty > 10 km | 918,417 |
| individual count 0 | 627 |
| equal coordinates | 63 |
| zero coordinates | 756 |
| sea | 5,350,967 |
| kept | 123,902,939 |

44,661,838 species x cell values kept; 0 cells without a CHELSA value dropped.

## Niche

Study species with a niche: 901 of 903; with n_obs_sub >= 30: 880. Median n_obs_sub 7,031.

## Comparison with Lynn's niche

Lynn's names joined through the community taxonomy (774 names; 70 hand-mapped Lynn-only names left out). T_opt against Lynn's mean_sum_sub on 704 shared species: r = 0.996. Species with a niche now and none in Lynn's file: 197.

Largest differences:

| species | lynn_name | lynn_cells | our_cells | lynn_T | T_opt | difference | cause |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Helianthemum oelandicum | Helianthemum alpestre | 1251 | 4845 | 10.51 | 14.92 | 4.41 | Lynn's niche from another GBIF name |
| Oxytropis humifusa | Oxytropis glacialis | 12 | 137 | 9.67 | 6.31 | -3.36 | Lynn's niche from another GBIF name |
| Alpagrostis alpina | Agrostis alpina | 34 | 2899 | 12 | 9.63 | -2.37 | Lynn's niche from another GBIF name |
| Persicaria runcinata | Polygonum runcinatum | 231 | 595 | 21.83 | 19.75 | -2.08 | Lynn's niche from another GBIF name |
| Koeleria macrantha | Koeleria cristata | 443 | 22126 | 19.04 | 17.32 | -1.71 | Lynn's niche from another GBIF name |
| Aster diplostephioides | Aster diplostephioides | 102 | 52 | 12.67 | 11.13 | -1.55 | same name: new records, keys and cleaning |
| Anemone cathayensis | Anemone cathayensis | 44 | 56 | 21.25 | 19.76 | -1.49 | same name: new records, keys and cleaning |
| Erigeron sumatrensis | Conyza sumatrensis | 652 | 46195 | 21.17 | 19.76 | -1.4 | Lynn's niche from another GBIF name |
| Muhlenbergia tricholepis | Blepharoneuron tricholepis | 24 | 812 | 15.3 | 16.54 | 1.24 | Lynn's niche from another GBIF name |
| Empetrum nigrum | Empetrum hermaphroditum | 28472 | 120068 | 11.54 | 12.77 | 1.23 | Lynn's niche from another GBIF name |
| Pilosella hoppeana | Hieracium hoppeanum | 34 | 1176 | 11.94 | 13.06 | 1.12 | Lynn's niche from another GBIF name |
| Trollius chinensis | Trollius chinensis | 195 | 273 | 19.44 | 18.33 | -1.11 | same name: new records, keys and cleaning |
| Polypogon viridis | Polypogon viridis | 6357 | 9946 | 20 | 18.97 | -1.02 | same name: new records, keys and cleaning |
| Pedicularis rhynchodonta | Pedicularis rhynchodonta | 13 | 17 | 12.37 | 11.39 | -0.97 | same name: new records, keys and cleaning |
| Veronica szechuanica | Veronica szechuanica | 78 | 89 | 17.93 | 17.09 | -0.84 | same name: new records, keys and cleaning |

Koeleria macrantha: n_obs_sub = 22,126, T_opt = 17.32.

## Warnings captured

- name_usage() works by default with the out-of-date GBIF Backbone Taxonomy. Consider using rcol::col_usage() instead.
- Numeric taxonomic keys detected (taxonKey). These are legacy GBIF Backbone identifiers. Switching to Backbone checklistKey (d7dddbf4-2cf0-4f39-9b2a-bb099caae36c). Please consider migrating to COL XR identifiers - use gbif_to_col() to look up the COL equivalents.
- [vect] guessed crs

Expected: rgbif notes that the GBIF Backbone (chosen deliberately, as for Lynn's download) is no longer its default; "[vect] guessed crs" comes from CoordinateCleaner's zero-coordinate test, which builds its 0/0 point without a coordinate system, and terra's guess (longitude/latitude) is the right one.

## Deviations from PLAN.md

Against Lynn's script, all decided 28 Oct: keys from our species with author and ledger; records grouped by our species; bio18 for P_opt; living specimens dropped; thinning on the CHELSA grid (Lynn et al. 2021: "we created a raster grid with cell length equal to 0.0083 deg and sampled a single representative point from the occurrence data per grid cell. The grid resolution approximates the 30 s resolution of the climate data"); n_obs_sub counts cells with a climate value.
