### MACRO Morphological Data Wrangling
library(readr)
library(tidyverse)
library(ape)
library(readxl)
# Data wrangling
EG_Witt <- read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/EG_Witt.csv")

Linck_et_al_data <- read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/Linck_et-al_data.csv")

elevs <- read_csv("~/Dropbox/Research/Macro_SSDI/data/elevs.csv")

avonet <- read_csv("~/Dropbox/Research/Avian_Evo_Trans_Rev/analysis/data/Birdtree_AVONET.csv")

Quinter_Jetz_bird_elevs <- read_csv("data/Quintero_Jetz_bird_elevations.csv")

# read in internet amniote data
Amniote_2015 <- read.csv("data/Amniote_Database_Aug_2015.csv")
animal_traits <- read.csv("data/animaltraits.csv")
avian_Lislevand_2007 <- read.delim("data/avian_body_size.txt", sep="\t")

# read in Dunning excel file
# accessing all the sheets 
sheet = excel_sheets("data/Dunning.xls")

# applying sheet names to dataframe names
data_frame = lapply(setNames(sheet, sheet), 
                    function(x) read_excel("data/Dunning.xls", sheet=x))

# attaching all dataframes together
Dunning = bind_rows(data_frame, .id="Sheet")

#read in tree file
tree <- read.tree("data/birds_mcc.tre")

Dunning2 <- Dunning%>%
  dplyr::select(!Sheet)%>%
  distinct()%>%
  fill(Species)%>%
  mutate(family = str_extract(Species, "(?<=Family )\\w+"))%>%
  fill(family)%>%
  filter(!grepl('Family', Species))%>%
  filter_all(any_vars(!is.na(.)))%>%
  mutate(species = sub("^\\s*(\\S+\\s+\\S+).*", "\\1", Species))%>%
  dplyr::select(1,16:17, 3:11)%>%
  filter(Sex %in% c("M", "F"))%>%
  group_by(Species)%>%
  mutate(count = n())%>%
  pivot_wider(.,id_cols =c(Species, species), values_from = c(Mean, Min, Max, S.D.), names_from = Sex, values_fn = ~mean(.x, na.rm = TRUE))%>% #taking mean when multiple values per species
  mutate(male_mass = Mean_M, 
         female_mass = Mean_F,
         taxon = as.factor(species))%>%
  dplyr::select(taxon, male_mass, female_mass)


avonet <- avonet%>%
  mutate(taxon = factor(Species3),
         family = factor(Family3))%>%
  dplyr::select(taxon, family, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth, Tarsus.Length, Kipps.Distance, `Hand-Wing.Index`, Tail.Length, Habitat, Habitat.Density, Migration, Trophic.Level, Trophic.Niche, Min.Latitude, Max.Latitude, Centroid.Latitude, Centroid.Longitude, Range.Size)

Amniote_2015_ALL <- Amniote_2015%>%
  filter(male_body_mass_g >= -500,
         female_body_mass_g >= -500)%>%
  na.omit()%>%
  mutate(taxon = factor(paste(genus, species, sep=" ")),
         male = male_body_mass_g,
         female = female_body_mass_g,
         sex_diff_male_female = male - female,
         SDI = (male-female)/female)%>%
  dplyr::select(taxon, male, female, sex_diff_male_female, SDI)
Amniote_2015_ALL$taxon <- reorder(Amniote_2015_ALL$taxon, Amniote_2015_ALL$sex_diff_male_female)


animal_traits_all <- animal_traits%>%
  dplyr::select(genus, species, sex, body.mass)%>%
  mutate(body.mass = body.mass *1000, #change from kg to g
         taxon = factor(species),
         sex =factor(sex))%>%
  filter(!is.na(body.mass) & !is.na(sex) & !sex == "")%>%
  group_by(taxon, sex)%>%
  dplyr::summarise(body.mass = mean(body.mass, na.rm = T))%>%
  pivot_wider(.,id_cols =c(taxon), values_from = body.mass, names_from = sex)%>%
  dplyr::select(!both)%>%
  na.omit()%>%
  mutate(sex_diff_male_female = male - female,
         SDI = (male-female)/female)%>%
  dplyr::select(taxon, male, female, sex_diff_male_female, SDI)

Amniote_2015_birds <- Amniote_2015%>%
  filter(class == "Aves",
         male_body_mass_g >= -500,
         female_body_mass_g >= -500)%>%
  na.omit()%>%
  mutate(male_mass = male_body_mass_g,
         female_mass = female_body_mass_g,
         taxon = factor(paste(genus, species, sep=" ")))%>%
  dplyr::select(taxon, male_mass, female_mass)


animal_traits_birds <- animal_traits%>%
  filter(class=="Aves")%>%
  dplyr::select(genus, species, sex, body.mass)%>%
  mutate(body.mass = body.mass *1000, #change from kg to g
         taxon = factor(species))%>%
  filter(!is.na(body.mass) & !is.na(sex) & !sex == "")%>%
  group_by(taxon, sex)%>%
  dplyr::summarise(body.mass = mean(body.mass, na.rm = T))%>%
  pivot_wider(.,id_cols =c(taxon), values_from = body.mass, names_from = sex)%>%
  na.omit()%>%
  mutate(male_mass = male,
         female_mass = female)%>%
  dplyr::select(taxon, male_mass, female_mass)

Livesand_birds <- avian_Lislevand_2007%>%
  dplyr::select(Species_name, M_mass, F_mass, M_wing, F_wing, M_tarsus, F_tarsus, M_bill, F_bill)%>%
  mutate(taxon = factor(Species_name))%>%
  filter(M_mass >= -500,
         F_mass >= -500,
         M_wing >= -500,
         F_wing >= -500,
         M_tarsus >= -500,
         F_tarsus >= -500,
         M_bill >= -500,
         F_bill >= -500)%>%
  na.omit()%>%
  group_by(taxon)%>%
  dplyr::summarise(male_mass = mean(M_mass, na.rm = T),
            female_mass = mean(F_mass, na.rm = T),
            male_wing = mean(M_wing, na.rm=T),
            female_wing = mean(F_wing, na.rm=T),
            male_tarsus = mean(M_tarsus, na.rm=T),
            female_tarsus = mean(F_tarsus, na.rm=T),
            male_bill = mean(M_bill, na.rm=T),
            female_bill = mean(F_bill, na.rm=T))%>%
  dplyr::select(taxon, male_mass, female_mass, male_wing, female_wing, male_tarsus, female_tarsus, male_bill, female_bill)

df.Linck <-Linck_et_al_data%>%
  filter(!is.na(mass),
         sex %in% c("male", "female"))%>%
  mutate(elev_bin = cut(elevation, breaks=7),
         sex = factor(sex),
         taxon = factor(species)
  )%>%
  group_by(taxon) %>% 
  filter(all(levels(sex) %in% sex))%>%
  group_by(taxon, sex) %>% 
  dplyr::summarise(mean_mass = mean(mass),
  )%>%
  pivot_wider(id_cols = c(taxon), names_from = c(sex), values_from = c(mean_mass))

Quinter_Jetz_bird_elevs <-Quinter_Jetz_bird_elevs%>% #take min and max of species elev ranges hopefully this is not oversimplification
  group_by(Species)%>%
  filter(`Maximum elevation` == max(`Maximum elevation`),
         `Minimum elevation`== min(`Minimum elevation`))%>% #now delete duplicates, but first extract only rows you need
  dplyr::select(Species, `Minimum elevation`, `Maximum elevation`)%>%
  distinct()%>%
  dplyr::rename(taxon = Species,
                min_elev = `Minimum elevation`,
                max_elev = `Maximum elevation`)

#combine datasets
all_animals <- rbind(Amniote_2015_ALL, animal_traits_all)
all_birds <- Dunning2%>%
  bind_rows(.,Amniote_2015_birds, Livesand_birds, animal_traits_birds, df.Linck)%>%
  left_join(., avonet)%>%
  left_join(., Quinter_Jetz_bird_elevs)%>%
  distinct()%>%
  mutate(SSDI_mass = ((male_mass-female_mass)/female_mass),
         SSDI_bill = ((male_bill-female_bill)/female_bill),
         SSDI_wing = ((male_wing-female_wing)/female_wing),
         SSDI_tarsus = ((male_tarsus-female_tarsus)/female_tarsus))
  
  

#remove duplicates
# all_animals <- all_animals%>%
#   distinct(., taxon, .keep_all = TRUE)%>%
#   droplevels()
# all_birds <- all_birds%>%
#   distinct(., taxon, .keep_all = TRUE)%>%
#   droplevels()

#add crane no crane column to birds
# all_birds <- all_birds%>%
#   mutate(crane_y_n = factor(if_else(taxon %in% c("Grus canadensis canadensis", "Grus canadensis tabida"), 1,0)))

#Reorder by sex diff for plotting
all_animals$taxon <- reorder(all_animals$taxon, all_animals$sex_diff_male_female)
all_birds$taxon <- reorder(all_birds$taxon, all_birds$sex_diff_male_female)

all_birds <- all_birds|>
  mutate(species = factor(taxon))|>
  dplyr::select(species, male, female)



elevs <- elevs%>%
  mutate(species = str_replace(species, "_", " "),
         species = factor(species))|>
  distinct()


df <- rbind(df, all_birds)

L.elevs<- Linck_et_al_data%>%
  group_by(species)%>%
  dplyr::summarise(max_sample_elev = max(elevation, na.rm=T),
            mean_sample_elev = mean(elevation, na.rm=T),
            min_sample_elev = min(elevation, na.rm=T))

L.blood<- Linck_et_al_data%>%
  group_by(species)%>%
  dplyr::summarise(mean_hct = mean(hct, na.rm=T),
            sd_hct = sd(hct, na.rm=T),
            mean_hb = mean(hb, na.rm=T),
            sd_hb = sd(hb, na.rm=T))

L.blood.sex<- Linck_et_al_data%>%
  filter(sex %in% c("male", "female"))%>%
  group_by(species, sex)%>%
  dplyr::summarise(mean_hct = mean(hct, na.rm=T),
            sd_hct = sd(hct, na.rm=T),
            mean_hb = mean(hb, na.rm=T),
            sd_hb = sd(hb, na.rm=T))%>%
  pivot_wider(id_cols = c(species), names_from = c(sex), values_from = c(mean_hct, sd_hct, mean_hb, sd_hb))


df<-df%>%
  left_join(., L.elevs, by="species")

df<-df%>%
  left_join(., elevs, by = "species")

df<- df%>%
  left_join(., avonet, by = "species")%>%
  mutate(amplitude = max_elev - min_elev)

#generate Lovich's simple SSDI Metric

df <- df%>%
  mutate(lov.SSDI = if_else(male > female, (male/female), (female/male)),
         SSDI = (male-female)/female,
         sex_diff = male-female,
         mean_mass_mf = (male + female)/2,
         species = factor(species))

df<-df%>%
  left_join(., L.blood, by="species")

df<-df%>%
  left_join(., L.blood.sex, by="species")

# elev.d <- elevs%>%   #find duplicate names
#   filter(duplicated(species))

df <- df|>
  filter(duplicated(species) == FALSE)

#elevation range is 5000 m

df.all.birds <- all_birds%>%
  left_join(., Quinter_Jetz_bird_elevs, by = "species")%>%
  left_join(., avonet, by = "species")%>%
  mutate(amplitude = max_elev - min_elev,
         SSDI = (male-female)/female,
         mdpt.elev = (max_elev - min_elev)/2)

#get vector of species for tree THIS HAS ALL BEEN DONE ALREADY

df.all.birds2 <-df.all.birds %>% 
  group_by(family) %>% 
  filter(n() >= 10)%>%
  na.omit()

species_in_tree <- tree$tip.label

species <- df.all.birds|>
  dplyr::select(species)|>
  mutate(species = as.character(str_replace(species, " ", "_")))%>%
  filter(species %in% species_in_tree)%>%
  pull()


Livesand_birds <- Livesand_birds%>%
  mutate(species = as.character(taxon))%>%
  left_join(., Quinter_Jetz_bird_elevs)

Dunning2 <- Dunning%>%
  filter(Sex %in% c("M", "F"))%>%
  group_by(Species)%>%
  mutate(count = n())%>%
  pivot_wider(.,id_cols =c(Species, species), values_from = c(Mean, Min, Max, S.D.), names_from = Sex, values_fn = ~mean(.x, na.rm = TRUE))%>%
  left_join(., Quinter_Jetz_bird_elevs, by="species")%>%
  left_join(., avonet, by ="species")%>%
  mutate(SSDI_con = ((Mean_M -Mean_F)/Mean_F),
         amp= max_elev-min_elev,
         mdpt = (max_elev + min_elev)/2)

# 
tree <- ape::keep.tip(tree, tip = as.character(species))
# write.tree(tree, file="data/tree_4398spp.tre")
# write.tree(tree, file="data/tree_395spp.tre")
write.csv(df, "data/All_elev_macro_morpho.csv")
write.csv(Livesand_birds, "data/Livesand_sex_morph_elev.csv")
write.csv(df.all.birds, "data/All_bird_SSDI_macro_morpho.csv")
write.csv(Dunning, "data/Dunning_Quintero__avonet_mass_morph_eco_elevs.csv")
