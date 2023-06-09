### MACRO Morphological Data Wrangling
pacman::p_load(
readr,
tidyverse,
ape,
readxl
)
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
  group_by(Species, family)%>%
  mutate(count = n())%>%
  pivot_wider(.,id_cols =c(Species, species, family), values_from = c(Mean, Min, Max, S.D.), names_from = Sex, values_fn = ~mean(.x, na.rm = TRUE))%>% #taking mean when multiple values per species
  mutate(mass_male = Mean_M, 
         mass_female = Mean_F,
         taxon = as.factor(species))%>%
  dplyr::select(taxon, family, mass_male, mass_female)


avonet <- avonet%>%
  mutate(taxon = factor(Species3),
         family = Family3)%>%
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
  mutate(mass_male = male_body_mass_g,
         mass_female = female_body_mass_g,
         taxon = factor(paste(genus, species, sep=" ")))%>%
  dplyr::select(taxon, mass_male, mass_female)


animal_traits_birds <- animal_traits%>%
  filter(class=="Aves")%>%
  dplyr::select(genus, species, family, sex, body.mass)%>%
  mutate(body.mass = body.mass *1000, #change from kg to g
         taxon = factor(species))%>%
  filter(!is.na(body.mass) & !is.na(sex) & !sex == "")%>%
  group_by(taxon, family, sex)%>%
  dplyr::summarise(body.mass = mean(body.mass, na.rm = T))%>%
  pivot_wider(.,id_cols =c(taxon, family), values_from = body.mass, names_from = sex)%>%
  na.omit()%>%
  mutate(mass_male = male,
         mass_female = female)%>%
  dplyr::select(taxon, family, mass_male, mass_female)

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
  dplyr::summarise(mass_male = mean(M_mass, na.rm = T),
            mass_female = mean(F_mass, na.rm = T),
            male_wing = mean(M_wing, na.rm=T),
            female_wing = mean(F_wing, na.rm=T),
            male_tarsus = mean(M_tarsus, na.rm=T),
            female_tarsus = mean(F_tarsus, na.rm=T),
            male_bill = mean(M_bill, na.rm=T),
            female_bill = mean(F_bill, na.rm=T))%>%
  dplyr::select(taxon, mass_male, mass_female, male_wing, female_wing, male_tarsus, female_tarsus, male_bill, female_bill)


#Linck dataset to get dimoprhism in blood metrics
L.elevs<- Linck_et_al_data%>%
  group_by(species, family)%>%
  dplyr::summarise(max_sample_elev = max(elevation, na.rm=T),
                   mean_sample_elev = mean(elevation, na.rm=T),
                   min_sample_elev = min(elevation, na.rm=T))%>%
  mutate(taxon=species)
  
df.Linck <-Linck_et_al_data%>%
  filter(!is.na(mass),
         sex %in% c("male", "female"))%>%
  mutate(elev_bin = cut(elevation, breaks=7),
         sex = factor(sex),
         taxon = factor(species)
  )%>%
  group_by(taxon) %>% 
  filter(all(levels(sex) %in% sex))%>%
  group_by(taxon, sex, family) %>% 
  dplyr::summarise(mean_mass = mean(mass),
                   mean_hb = mean(hb, na.rm=T),
                   mean_hct = mean(hct, na.rm=T)
  )%>%
  pivot_wider(id_cols = c(taxon, family), names_from = c(sex), values_from = c(mean_mass, mean_hb, mean_hct))%>%
  mutate(
         SSDI_hct = ((mean_hct_male- mean_hct_female)/mean_hct_female),
         SSDI_hb = ((mean_hb_male- mean_hb_female)/mean_hb_female),
         SSDI_mass = ((mean_mass_male - mean_mass_female)/mean_mass_female))%>%
  left_join(., L.elevs)

Quinter_Jetz_bird_elevs <-Quinter_Jetz_bird_elevs%>% #take min and max of species elev ranges hopefully this is not oversimplification
  group_by(Species)%>%
  filter(`Maximum elevation` == max(`Maximum elevation`),
         `Minimum elevation`== min(`Minimum elevation`))%>% #now delete duplicates, but first extract only rows you need
  dplyr::select(Species, `Minimum elevation`, `Maximum elevation`)%>%
  distinct()%>%
  dplyr::rename(taxon = Species,
                min_elev = `Minimum elevation`,
                max_elev = `Maximum elevation`)%>%
  dplyr::select(taxon, min_elev, max_elev)



#combine datasets
all_animals <- rbind(Amniote_2015_ALL, animal_traits_all)

all_birds <- Dunning2%>%
  dplyr::select(!Species)%>%
  left_join(., Livesand_birds, by = "taxon")%>%
  mutate(taxon = as.character(taxon))%>%
  mutate_all(~ifelse(is.nan(.), NA, .))%>% #convert NAN to NA so we can ge rwo means for certain columns
  mutate(
  mean_mass_male = mean(c_across(c(mass_male.x, mass_male.y)), na.rm = TRUE),
  mean_mass_female = mean(c_across(c(mass_female.x, mass_female.y)), na.rm = TRUE))%>%
  dplyr::select(taxon, family, mean_mass_male, mean_mass_female, male_wing, female_wing, male_tarsus, female_tarsus, male_bill, female_bill)%>%
  mutate_all(~ifelse(is.nan(.), NA, .))%>% #convert NAN to NA again because better for analyses
  full_join(., df.Linck)%>%
  left_join(., avonet, by="taxon")%>%
  left_join(., Quinter_Jetz_bird_elevs)%>%
  distinct()%>%
  mutate(SSDI_mass = if_else(is.na(mean_mass_male) | is.na(mean_mass_female), NA, ((mean_mass_male-mean_mass_female)/mean_mass_female)),
         SSDI_bill = if_else(is.na(male_bill) | is.na(female_bill), NA, ((male_bill-female_bill)/female_bill)),
         SSDI_wing = if_else(is.na(male_wing) | is.na(female_wing), NA, ((male_wing-female_wing)/female_wing)),
         SSDI_tarsus = if_else(is.na(male_tarsus) | is.na(female_tarsus), NA, ((male_tarsus-female_tarsus)/female_tarsus)),
         sex_diff_male_female = if_else(is.na(mean_mass_male) | is.na(mean_mass_female), NA, (mean_mass_male-mean_mass_female)),
         family = family.x)%>%
  mutate_all(~ifelse(is.nan(.), NA, .))%>%#convert NAN to NA again because better for analyses
  filter(!is.na(mean_mass_female), !is.na(mean_mass_male))%>%
  dplyr::select(-family.x, -family.y, -Species, -species)
  
  

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

# all_birds <- all_birds|>
#   mutate(species = factor(taxon))|>
#   dplyr::select(species, male, female)



# elevs <- elevs%>%
#   mutate(species = str_replace(species, "_", " "),
#          species = factor(species))|>
#   distinct()


#df <- rbind(df, all_birds)




#generate Lovich's simple SSDI Metric


#get vector of species for tree THIS HAS ALL BEEN DONE ALREADY


species_in_tree <- tree$tip.label

species <- all_birds|>
  dplyr::select(taxon)|>
  mutate(species_tree = as.character(str_replace(taxon, " ", "_")))%>%
  filter(species_tree %in% species_in_tree)%>%
  pull()





# 
tree <- ape::keep.tip(tree, tip = as.character(species))

#Write out tree and dataframe Tree will need fruther pruning depedning on analyses
write.tree(tree, file=paste0("data/tree_", length(tree$tip.label), "spp.tre"))
write.csv(all_birds, "data/All_birds_elev_macro_morpho.csv")

