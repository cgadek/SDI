### SDI (blood) Data Wrangle###

pacman::p_load(
  tidyverse,
  reshape,
  reshape2,
  dplyr,
  tidyr,
  car,
  picante,
  rcompanion,
  GGally,
  Hmisc,
  gridExtra,
  stats,
  gplots,
  ggplot2,
  ggExtra,
  cowplot,
  colorspace,
  stats4, # Forces knitr to work when it's being wonky
  PMCMR, #Allows Kruskal-Wallis post-hocs
  effects,
  gridExtra,
  readxl
)

#Do Jessie's Hummingbird data first so we can delet from Lincks...
#actually combine Linck and Williamsons then process

avonet <- read_csv("~/Dropbox/Research/Avian_Evo_Trans_Rev/analysis/data/Birdtree_AVONET.csv")

Quinter_Jetz_bird_elevs <- read_csv("data/Quintero_Jetz_bird_elevations.csv")

avian_Lislevand_2007 <- read.delim("data/avian_body_size.txt", sep="\t")

williamson.chile<- read_csv("data/ChileBloodCatalog_NoPatagona_Wrangled_ForChauncey_2023-12-23.csv")

CRG_repro_3n <- read_csv("data/n_3_samples_for_repro_data.csv")


#Patagona data####
#pull pgc nks to filter form Linck dataset
pgc.nk <- read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%pull(nk)
pgcl.elevs <- read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  dplyr::select(species, sex, elev, hb, hct)%>% #select relevant columns
  mutate(
         family = "Trochilidae", #create family name
         species=factor(gsub("_", " ", species)))%>% #add underscore for phylo tree later
  filter(
    sex %in% c("male", "female"))%>% #keep only birds identified to sex
  dplyr::rename(elevation = elev)%>% #chnage name of elev to elevation
  dplyr::select(species,family, sex, elevation, hb, hct)%>% #select relevant columns
  bind_rows(.,
  read_csv("data/blood_data_final.csv")%>% #merge with other blood file
  dplyr::select(nk,species,family, sex, elevation, hb, hct)%>% #select relevant columns
  filter(
    !nk %in% pgc.nk, #pull only NKs not in Jessie's file
    sex %in% c("male", "female"), #only species identified to sex
    !species %in% c("Accipiter striatus", "Falco sparverius", "Passer domesticus"))%>% #get rid of a few (WHY AM I DOING THIS? because we have very few samples and I add better samples later?)
  dplyr::select(species,family, sex, elevation, hb, hct))%>% #select relevant columns
  mutate(taxon=factor(gsub(" ", "_", species)))%>% #add underscore for phylo tree later
  bind_rows(., read_csv("data/ChileBloodCatalog_NoPatagona_Wrangled_ForChauncey_2023-12-23.csv")%>% #now bind with Jessie's blood data
              dplyr::select(nk,species,family, sex, elev, hb, hct, pull_too_young)%>%#select relevant columns
              dplyr::rename(elevation = elev)%>% #change name
              filter(pull_too_young<2, #Check with Jessie here
                !nk %in% pgc.nk,
                sex %in% c("male", "female"),
                !species %in% c("Accipiter striatus", "Falco sparverius", "Passer domesticus"))%>%
              dplyr::select(species,family, sex, elevation, hb, hct))%>%
  bind_rows(., read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/blood_final_bioclim_12-9-2022.csv")%>% #bind with marshbird blood file
              dplyr::select(species, family, sex, elev, hb_final, hct_final)%>%
              dplyr::rename(elevation = elev,
                     hb = hb_final,
                     hct = hct_final)%>%
              filter(
                sex %in% c("male", "female"),
                species %in% c("Tachuris rubrigastra", "Phleocryptes melanops"))%>%
              dplyr::select(species,family, sex, elevation, hb, hct))%>%
  group_by(species, family)%>% #group by species and sex to get mean, min, and max_sample elevs
  dplyr::summarise(max_sample_elev = max(elevation, na.rm=T),
                   mean_sample_elev = mean(elevation, na.rm=T),
                   min_sample_elev = min(elevation, na.rm=T))%>%
  mutate(taxon=factor(gsub(" ", "_", species))) # a lot of work just for summary sample elevations...

pgcl <- read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
filter(
  sex %in% c("male", "female"))%>%
  mutate(
         sex = factor(sex),
         family = "Trochilidae",
         taxon = species)%>%
  group_by(taxon) %>% 
  filter(all(levels(sex) %in% sex))%>%
  dplyr::rename(elevation = elev)%>%
  dplyr::select(taxon, species,family, sex, elevation, mass, hb, hct)%>%
  bind_rows(.,
            read_csv("data/blood_data_final.csv")%>%
              dplyr::select(nk, species,family, sex, elevation, mass, hb, hct)%>%
              filter(
                !nk %in% pgc.nk,
                sex %in% c("male", "female"),
                !species %in% c("Accipiter striatus", "Falco sparverius", "Passer domesticus"))%>%
              mutate(
                     sex = factor(sex),
                     taxon = factor(gsub(" ", "_", species)))%>%
              dplyr::select(nk, taxon, species, family, sex, elevation,  mass, hb, hct))%>%
  bind_rows(., read_csv("data/ChileBloodCatalog_NoPatagona_Wrangled_ForChauncey_2023-12-23.csv")%>%
              dplyr::select(nk,species,family, sex, elev,mass, hb, hct, pull_too_young)%>%
              dplyr::rename(elevation = elev)%>%
              filter(pull_too_young==1,
                !nk %in% pgc.nk,
                sex %in% c("male", "female"))%>%
              mutate(taxon=factor(gsub(" ", "_", species)))%>%
            dplyr::select(taxon, species, family, sex, elevation,  mass, hb, hct))%>%
  bind_rows(.,read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/blood_final_bioclim_12-9-2022.csv")%>%
              filter(
                sex %in% c("male", "female"),
                species %in% c("Tachuris rubrigastra", "Phleocryptes melanops"))%>%
              mutate(elev_bin = cut(elev, breaks=7),
                     sex = factor(sex),
                     taxon = factor(gsub(" ", "_", species)),
                     hb = hb_final, 
                     hct = hct_final
              ))%>%
  ungroup()%>%
  group_by(taxon)%>%
  mutate(mean_taxon_hb = mean(hb, na.rm=T),
         mean_taxon_hct = mean(hct, na.rm=T))%>%
  group_by(taxon, sex, family, mean_taxon_hb, mean_taxon_hct) %>% 
  dplyr::summarise(mean_mass = mean(mass, na.rm=T),
                   mean_hb = mean(hb, na.rm=T),
                   mean_hct = mean(hct, na.rm=T),
                   n_hct=sum(!is.na(hct)),
                   n_hb=sum(!is.na(hb)))%>%
  pivot_wider(id_cols = c(taxon, family, mean_taxon_hb, mean_taxon_hct), names_from = c(sex), values_from = c(mean_mass, mean_hb, mean_hct, n_hct, n_hb))%>%
  mutate(
    SSDI_hct = ((mean_hct_male- mean_hct_female)/mean_hct_female),
    SSDI_hb = ((mean_hb_male- mean_hb_female)/mean_hb_female),
    SSDI_mass = ((mean_mass_male - mean_mass_female)/mean_mass_female))%>%
  left_join(., pgcl.elevs%>%dplyr::select(-c(species, family)), by="taxon")
  



# use this function found on stack exchange here:https://stackoverflow.com/questions/18509527/first-letter-to-upper-case
#to turn genera name into uppercase 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Load Stoddard et al.2017 egg data ####
stoddard <- readxl::read_excel("~/Dropbox/Data_repo/Stoddard_egg_shape/aaj1945_datas1_egg_shape_by_species_v2.xlsx")%>%
  dplyr::select(Family, Species, `AvgLength (cm)`)%>%
  dplyr::rename(family=Family,
         species = Species,
         avg.length = `AvgLength (cm)`)%>%
  mutate(taxon = gsub(" ", "_", species))%>%
  dplyr::select(family, taxon, avg.length)%>%
  dplyr::rename(egg_length = avg.length)


#Get others from Myhrvold et al. 2015
Myhrvold<- read_csv("~/Dropbox/Data_repo/Amniote_Database_Aug_2015.csv")%>%
  filter(class =="Aves")%>%
  mutate(species = paste(genus, species))%>%
  mutate(across(where(is.numeric), ~na_if(., -999)))%>%
  dplyr::select(genus, species, female_maturity_d, litter_or_clutch_size_n, litters_or_clutches_per_y, egg_mass_g, egg_length_mm, egg_width_mm)%>%
  dplyr::rename(clutch_size = litter_or_clutch_size_n,
                clutch_per_y = litters_or_clutches_per_y,
                egg_mass = egg_mass_g,
                egg_length = egg_length_mm,
                egg_width = egg_width_mm)

Myhrvold.genus.egg.traits<-Myhrvold%>%
  group_by(genus)%>%
  dplyr::summarise(gen_egg_mass = mean(egg_mass, na.rm=T),
                   gen_egg_length =mean(egg_length, na.rm=T),
                   gen_egg_width =mean(egg_length, na.rm=T),
                   gen_clutch_size = median(clutch_size, na.rm=T))


#

development<- readxl::read_excel("data/41467_2020_16257_MOESM3_ESM.xlsx", sheet=4)%>%
  dplyr::rename(taxon =  binomial)


#Minias et al. blood data####
#first elev

minias.elev <-readxl::read_excel("data/Minias_RAW_DATA.xlsx")%>%
              filter(Captivity==0,
                     Sex%in% c("M", "F"))%>%
              dplyr::select(`Scientific name`, Family, `Haemoglobin concentration (g/l)`, `Sample size`, Sex, `Elevation (m a.s.l.)`, `Clutch size`)%>%
              dplyr::rename(species = `Scientific name`,
                     hb = `Haemoglobin concentration (g/l)`,
                     sex = Sex,
                     family = Family,
                     elevation = `Elevation (m a.s.l.)`,
                     clutch_size = `Clutch size`)%>%
              dplyr::select(species,family, sex, elevation,  hb)%>%
  group_by(species, family)%>%
  dplyr::summarise(max_sample_elev = max(elevation, na.rm=T),
                   mean_sample_elev = mean(elevation, na.rm=T),
                   min_sample_elev = min(elevation, na.rm=T))%>%
  mutate(taxon=factor(gsub(" ", "_", species)))


minias.hb <- readxl::read_excel("data/Minias_RAW_DATA.xlsx")%>%
    filter(Captivity==0,
           Sex%in% c("M", "F"))%>%
  dplyr::select(`Scientific name`, `Haemoglobin concentration (g/l)`, `Sample size`, Sex, `Elevation (m a.s.l.)`, `Clutch size`)%>%
    dplyr::rename(species = `Scientific name`,
           hb = `Haemoglobin concentration (g/l)`,
           n = `Sample size`,
           elev = `Elevation (m a.s.l.)`,
           clutch_size = `Clutch size`)%>%
  ungroup()%>%
  group_by(species)%>%
  mutate(mean_taxon_hb = mean((hb/10), na.rm=T))%>%
  ungroup()%>%
    group_by(species, Sex,mean_taxon_hb)%>%
    dplyr::summarise(mean_hb = mean(hb/10, na.rm=T),
                     n = sum(n),
                     clutch_size = mean(clutch_size, na.rm=T))%>%
    pivot_wider( names_from = "Sex", values_from =c("mean_hb", "n"))%>%
    na.omit()%>%
    dplyr::rename(mean_hb_male = mean_hb_M,
                  mean_hb_female = mean_hb_F,
                  n_hb_female = n_F, 
                  n_hb_male =n_M)%>%
    mutate(SSDI_hb = (mean_hb_male-mean_hb_female)/mean_hb_female)


minias <- readxl::read_excel("data/Minias_RAW_DATA.xlsx", sheet=2)%>%
  filter(Captivity==0,
         Sex%in% c("M", "F"))%>%
  dplyr::select(`Scientific name`, `Haematocrit (%)`, `Sample size`, Sex, `Elevation (m a.s.l.)`, `Clutch size`)%>%
  dplyr::rename(species = `Scientific name`,
         hct = `Haematocrit (%)`,
         n = `Sample size`,
         elev = `Elevation (m a.s.l.)`,
         clutch_size = `Clutch size`)%>%
  ungroup()%>%
  group_by(species)%>%
  mutate(mean_taxon_hct = mean(hct, na.rm=T))%>%
  ungroup()%>%
  group_by(species, Sex, mean_taxon_hct)%>%
  dplyr::summarise(mean_hct = mean(hct/100, na.rm=T),
                   n = sum(n),
                   clutch_size = mean(clutch_size, na.rm=T))%>%
  pivot_wider( names_from = "Sex", values_from =c("mean_hct", "n"))%>%
  na.omit()%>%
  dplyr::rename(mean_hct_male = mean_hct_M,
                mean_hct_female = mean_hct_F,
                n_hct_female = n_F, 
                n_hct_male =n_M)%>%
  mutate(SSDI_hct = (mean_hct_male-mean_hct_female)/mean_hct_female)%>%
  left_join(minias.hb)%>%
  left_join(minias.elev)

#Santema et al. blood data ####
#filter out minias birds
#NOTE Masses added to Santema et al. from Lislevand et al. 2007
santema <- readxl::read_excel("data/data_comp.xlsx")%>%
  dplyr::rename(mean_hct_male= hema_m,
         mean_hct_female= hema_f,
         n_male= samplesize_m,
         n_female= samplesize_f,
         mean_sample_elev= elevation)%>%
  mutate(species=firstup(scinam),
         taxon = gsub(" ", "_", species))%>%
  filter(!species %in% c(minias$species, minias.hct$species))%>%# remove minias species because redundant
  dplyr::select(taxon, species, family, n_male, n_female, mean_sample_elev, mean_hct_female, mean_hct_male, mean_mass_male, mean_mass_female)%>%
  group_by(species, taxon, family)%>%
  dplyr::summarise(mean_hct_male = mean(mean_hct_male)/100,
                   mean_hct_female = mean(mean_hct_female)/100,
                   mean_taxon_hct = (mean_hct_female+ mean_hct_male)/2,
                   mean_mass_male = mean(mean_mass_male),
                   mean_mass_female = mean(mean_mass_female),
                   mean_sample_elev = mean(mean_sample_elev),
                   n_hct_male = sum(n_male),
                   n_hct_female = sum(n_female),
                   SSDI_hct = (mean_hct_male-mean_hct_female)/mean_hct_female,
                   SSDI_mass = (mean_mass_male-mean_mass_female)/mean_mass_female)

#combine all datastes
df.blood <- bind_rows(pgcl, santema)%>%
  left_join(., 
            avonet%>%
              mutate(species = Species3)%>%
              dplyr::select(species, `Hand-Wing.Index`, Migration, Trophic.Niche))%>%
  left_join(., stoddard)%>%
  left_join(., Myhrvold)%>%
  left_join(., minias)%>%
  left_join(., minias.hct)%>%
  left_join(development)%>%
  mutate(genus = stringr::word(species))%>%
  left_join(Myhrvold.genus.egg.traits)

# match up tree labels
tree <- read.tree("data/birds_mcc.tre")

# created crosswalk excel file (load in)
blood.crosswalk <-  read_csv("data/blood_birdtree_crosswalk.csv")
    
#load dimorphism data
dimorph<- readxl::read_excel("data/Taxon_names_for_Jessie_WithFamily_DimorphismDone.xlsx")%>%
  dplyr::rename(taxon= taxon...1)%>%
  dplyr::select(taxon, family, dimorphic_01)

df.blood%>%
  filter(!taxon %in% c("Elaenia_sp.", "Spinus_sp.", "Scytalopus_sp."))%>%
  left_join(., blood.crosswalk%>%dplyr::select(2:3), by ="taxon")%>%
  left_join(dimorph)%>%
  left_join(read_excel("data/Brown_Witt_Wright_2012.xls")%>% # add in and left join the Brown Witt, Wright et al data
              dplyr::select(4:13)%>%dplyr::rename(species=Species))%>%
  mutate(birdtree = if_else(is.na(birdtree), taxon, birdtree),
         spec_gen_egg_mass = if_else(is.na(egg_mass), gen_egg_mass, egg_mass), #these create a column with species level egg measurements if available, genus level if next best, and NA if not avaiable at all.
         spec_gen_egg_length = if_else(is.na(egg_length), gen_egg_length, egg_length),
         spec_gen_egg_width = if_else(is.na(egg_width), gen_egg_width, egg_width),
         spec_gen_clutch_size = if_else(is.na(clutch_size), gen_clutch_size, clutch_size))%>%
 write_csv(., "data/species_means_blood_final.csv")

#get mean hb value across all birds to calibrate
all.blood%>%
  ungroup()%>%
  dplyr::summarise(meanhb = mean(hb, na.rm=T),
                   meanhct = mean(hct, na.rm=T))

df.blood%>%
  ungroup()%>%
  dplyr::summarise(meanhbSDI = mean(SSDI_hb, na.rm=T),
                   meanhctSDI = mean(SSDI_hct, na.rm=T))

threshold <- 3
df.cont <- df.blood %>%
  filter(n_hb_male >= threshold,
         n_hb_female >= threshold)%>%
  left_join(dimorph)%>%
  group_by(family) %>%
  filter(!is.na(SSDI_hb)) %>%
  mutate(
    mean_family_elev = mean(mean_sample_elev, na.rm = T),
    closest_to_0_hb = min(abs(SSDI_hb), na.rm = T),
    species_0_hb = species[which.min(abs(SSDI_hb))][1],
    male_0_SDI_hb = mean_hb_male[which.min(abs(SSDI_hb))],
    female_0_SDI_hb = mean_hb_female[which.min(abs(SSDI_hb))],
    diff_from_expected_hb_male = mean_hb_male - male_0_SDI_hb,
    #Now we create contribution based on mean SSDI --assuming that low male bias is the ancestral baseline condition
    species_base_hb = species[which.min(abs(SSDI_hb - 0.0297))][1],
    male_base_SDI_hb = mean_hb_male[which.min(abs(SSDI_hb - 0.0297))],
    female_base_SDI_hb = mean_hb_female[which.min(abs(SSDI_hb - 0.0297))],
    diff_from_expected_hb_male_base = mean_hb_male - male_base_SDI_hb,
    diff_from_expected_hb_female_base = mean_hb_female - female_base_SDI_hb,
    #Now we create contribution based off global mean hb whihc is probably no good
    diff_from_expected_hb_male_global_mean = mean_hb_male - 18,
    diff_from_expected_hb_female = mean_hb_female - female_0_SDI_hb,
    diff_from_expected_hb_female_global_mean = mean_hb_female - 18,
    diff_from_expected_elev = mean_sample_elev - mean_family_elev,
    diff_from_exp_hb_elev_corr_male = diff_from_expected_hb_male + (diff_from_expected_elev /
                                                                      1000),
    diff_from_exp_hb_elev_corr_female = diff_from_expected_hb_female + (diff_from_expected_elev /
                                                                          1000),
    diff_from_exp_hb_elev_corr_male_global_mean = diff_from_expected_hb_male_global_mean  + (diff_from_expected_elev /
                                                                                               1000),
    diff_from_exp_hb_elev_corr_female_global_mean  = diff_from_expected_hb_female_global_mean  + (diff_from_expected_elev /
                                                                                                    1000),
    diff_from_exp_hb_elev_corr_male_base = diff_from_expected_hb_male_base  + (diff_from_expected_elev /
                                                                                               1000),
    diff_from_exp_hb_elev_corr_female_base  = diff_from_expected_hb_female_base  + (diff_from_expected_elev /
                                                                                                    1000),
    hb_cont = if_else(
      diff_from_exp_hb_elev_corr_female > 0 &
        abs(diff_from_exp_hb_elev_corr_female) > abs(diff_from_exp_hb_elev_corr_male),
      "female driven increase",
      if_else(
        diff_from_exp_hb_elev_corr_male > 0 &
          abs(diff_from_exp_hb_elev_corr_male) > abs(diff_from_exp_hb_elev_corr_female),
        "male driven increase",
        if_else(
          diff_from_exp_hb_elev_corr_female < 0 &
            abs(diff_from_exp_hb_elev_corr_female) > abs(diff_from_exp_hb_elev_corr_male),
          "female driven reduction",
          if_else(
            diff_from_exp_hb_elev_corr_male < 0 &
              abs(diff_from_exp_hb_elev_corr_male) > abs(diff_from_exp_hb_elev_corr_female),
            "male driven reduction",
            "at reference"
          )
        )
      )
    ),
    hb_exp_diff = diff_from_exp_hb_elev_corr_male - diff_from_exp_hb_elev_corr_female,
    hb_exp_diff_global_mean = diff_from_exp_hb_elev_corr_male_global_mean - diff_from_exp_hb_elev_corr_female_global_mean,
    hb_exp_diff_base = diff_from_exp_hb_elev_corr_male_base - diff_from_exp_hb_elev_corr_female_base
  ) %>%
  full_join(
    .,
    df.blood %>%
      group_by(family) %>%
      filter(!is.na(SSDI_hct)) %>%
      mutate(
        mean_family_elev = mean(mean_sample_elev, na.rm = T),
        closest_to_0_hct = min(abs(SSDI_hct), na.rm = T),
        species_0_hct = species[which.min(abs(SSDI_hct))][1],
        male_0_SDI_hct = mean_hct_male[which.min(abs(SSDI_hct))],
        female_0_SDI_hct = mean_hct_female[which.min(abs(SSDI_hct))],
        diff_from_expected_hct_male = mean_hct_male - male_0_SDI_hct,
        diff_from_expected_hct_female = mean_hct_female - female_0_SDI_hct,
        #Now we create contribution based on mean SSDI --assuming that low male bias is the ancestral baseline condition
        species_base_hct = species[which.min(abs(SSDI_hct - 0.028))][1],
        male_base_SDI_hct = mean_hct_male[which.min(abs(SSDI_hct - 0.028))],
        female_base_SDI_hct = mean_hct_female[which.min(abs(SSDI_hct - 0.028))],
        diff_from_expected_hct_male_base = mean_hct_male - male_base_SDI_hct,
        diff_from_expected_hct_female_base = mean_hct_female - female_base_SDI_hct,
        #Now create contribution based on global mean Hb this is probably not very useful
        diff_from_expected_hct_female_global_mean = mean_hct_female -
          0.525,
        diff_from_expected_hct_male_global_mean = mean_hct_male -
          0.525,
        diff_from_expected_elev = mean_sample_elev - mean_family_elev,
        diff_from_exp_hct_elev_corr_male = diff_from_expected_hct_male + (diff_from_expected_elev /
                                                                            10000),
        diff_from_exp_hct_elev_corr_female = diff_from_expected_hct_female + (diff_from_expected_elev /
                                                                                10000),
        diff_from_exp_hct_elev_corr_male_global_mean = diff_from_expected_hct_male_global_mean  + (diff_from_expected_elev /
                                                                                                     10000),
        diff_from_exp_hct_elev_corr_female_global_mean  = diff_from_expected_hct_female_global_mean  + (diff_from_expected_elev /
                                                                                                          10000),
        diff_from_exp_hct_elev_corr_male_base = diff_from_expected_hct_male_base  + (diff_from_expected_elev /
                                                                                                     10000),
        diff_from_exp_hct_elev_corr_female_base  = diff_from_expected_hct_female_base  + (diff_from_expected_elev /
                                                                                                          10000),
        hct_cont = if_else(
          diff_from_exp_hct_elev_corr_female > 0 &
            abs(diff_from_exp_hct_elev_corr_female) > abs(diff_from_exp_hct_elev_corr_male),
          "female driven increase",
          if_else(
            diff_from_exp_hct_elev_corr_male > 0 &
              abs(diff_from_exp_hct_elev_corr_male) > abs(diff_from_exp_hct_elev_corr_female),
            "male driven increase",
            if_else(
              diff_from_exp_hct_elev_corr_female < 0 &
                abs(diff_from_exp_hct_elev_corr_female) > abs(diff_from_exp_hct_elev_corr_male),
              "female driven reduction",
              if_else(
                diff_from_exp_hct_elev_corr_male < 0 &
                  abs(diff_from_exp_hct_elev_corr_male) > abs(diff_from_exp_hct_elev_corr_female),
                "male driven reduction",
                "at reference"
              )
            )
          )
        ),
        hct_exp_diff = diff_from_exp_hct_elev_corr_male - diff_from_exp_hct_elev_corr_female,
        hct_exp_diff_global_mean = diff_from_exp_hct_elev_corr_male_global_mean - diff_from_exp_hct_elev_corr_female_global_mean,
        hct_exp_diff_base = diff_from_exp_hct_elev_corr_male_base - diff_from_exp_hct_elev_corr_female_base
      ) %>%
      dplyr::select(
        taxon,
        closest_to_0_hct,
        species_0_hct,
        male_0_SDI_hct,
        female_0_SDI_hct,
        diff_from_expected_hct_male,
        diff_from_expected_hct_female,
        diff_from_exp_hct_elev_corr_male,
        diff_from_exp_hct_elev_corr_female,
        species_base_hct,
        hct_cont,
        hct_exp_diff,
        hct_exp_diff_global_mean,
        hct_exp_diff_base
      )
  )




#Next big task is getting some sort of species expected value for Hb and Hct (mean?) then asking if individual points adjusted for elevation fall to above or below value by sex
all.blood <- read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  filter(
    sex %in% c("male", "female"))%>%
  mutate(
    sex = factor(sex),
    family = "Trochilidae",
    taxon = species)%>%
  group_by(taxon) %>% 
  filter(all(levels(sex) %in% sex))%>%
  dplyr::rename(elevation = elev)%>%
  dplyr::select(taxon, species,family, sex, elevation, mass, hb, hct)%>%
  bind_rows(.,
            read_csv("data/blood_data_final.csv")%>%
              dplyr::select(nk, species,family, sex, elevation, mass, hb, hct)%>%
              filter(
                !nk %in% pgc.nk,
                sex %in% c("male", "female"),
                !species %in% c("Accipiter striatus", "Falco sparverius", "Passer domesticus"))%>%
              mutate(
                sex = factor(sex),
                taxon = factor(gsub(" ", "_", species)))%>%
              dplyr::select(taxon, species, family, sex, elevation,  mass, hb, hct))%>%
  bind_rows(., read_csv("data/ChileBloodCatalog_NoPatagona_Wrangled_ForChauncey_2023-12-23.csv")%>%
              dplyr::select(nk,species,family, sex, elev,mass, hb, hct, pull_too_young)%>%
              dplyr::rename(elevation = elev)%>%
              filter(pull_too_young==1,
                     !nk %in% pgc.nk,
                     sex %in% c("male", "female"))%>%
              mutate(taxon=factor(gsub(" ", "_", species)))%>%
              dplyr::select(taxon, species, family, sex, elevation,  mass, hb, hct))%>%
  bind_rows(., read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/blood_final_bioclim_12-9-2022.csv")%>%
              dplyr::select(species, family, sex, elev, mass, hb_final, hct_final)%>%
              dplyr::rename(hb = hb_final,
                     hct = hct_final,
                     elevation = elev)%>%
              filter(
                sex %in% c("male", "female"),
                species %in% c("Tachuris rubrigastra", "Phleocryptes melanops"))%>%
              mutate(
                sex = factor(sex),
                taxon = factor(gsub(" ", "_", species))
              )%>%
              group_by(taxon) %>% 
              filter(all(levels(sex) %in% sex))%>%
              dplyr::select(taxon, species, family, sex, elevation,  mass, hb, hct))%>%
  bind_rows(., readxl::read_excel("data/Minias_RAW_DATA.xlsx")%>%
              filter(Captivity==0,
                     Sex%in% c("M", "F"))%>%
              mutate(Sex= if_else(Sex =="M", "male", if_else(Sex == "F", "female", NA)))%>%
              dplyr::select(`Scientific name`, Family, `Haemoglobin concentration (g/l)`, Sex, `Elevation (m a.s.l.)`, `Clutch size`)%>%
              dplyr::rename(species = `Scientific name`,
                     hb = `Haemoglobin concentration (g/l)`,
                     elevation = `Elevation (m a.s.l.)`,
                     sex=Sex,
                     family=Family)%>%
              mutate(
                sex = factor(sex),
                taxon = factor(gsub(" ", "_", species)),
                hb= hb/10)%>%
              dplyr::select(taxon, species, family, sex, elevation, hb))%>%
                  bind_rows(., readxl::read_excel("data/Minias_RAW_DATA.xlsx", sheet=2)%>%
                              filter(Captivity==0,
                                     Sex%in% c("M", "F"))%>%
                              mutate(Sex= if_else(Sex =="M", "male", if_else(Sex == "F", "female", NA)))%>%
                              dplyr::select(`Scientific name`, Family, `Haematocrit (%)`, Sex, `Elevation (m a.s.l.)`, `Clutch size`)%>%
                              dplyr::rename(species = `Scientific name`,
                                     hct = `Haematocrit (%)`,
                                     elevation = `Elevation (m a.s.l.)`,
                                     sex=Sex,
                                     family=Family)%>%
                              mutate(
                                sex = factor(sex),
                                taxon = factor(gsub(" ", "_", species)),
                                hct = hct/100)%>%
                              dplyr::select(taxon, species, family, sex, elevation,  hct))%>%
  bind_rows(.,readxl::read_excel("data/data_comp.xlsx")%>%
              dplyr::select(scinam, family, hema_m, hema_f)%>%
              ungroup()%>%
              group_by(scinam)%>%
              mutate(SSDI_hct = (hema_m -hema_f)/hema_f)%>%
              ungroup()%>%
              pivot_longer(cols =c(3:4), names_to = "sex", values_to = "hct")%>%
              mutate(sex =  firstup(gsub("hema_", '', sex)),
                    hct = hct/100,
                    species = sub("(.)", "\\U\\1", scinam, perl=TRUE),
                    taxon = gsub(" ", "_", species),
            sex= if_else(sex =="M", "male", if_else(sex == "F", "female", NA)))%>%
              dplyr::select(taxon , species, sex, hct))%>%
  mutate(genus = str_split(species, fixed("_"))[[1]][1])%>%
  left_join(., df.blood%>%
              dplyr::select(taxon,  mean_hb_male, mean_hb_female, mean_hct_male, mean_hct_female, SSDI_hct, SSDI_hb))%>%
  left_join(., blood.crosswalk%>%dplyr::select(2:3), by ="taxon")%>%
  mutate(birdtree = if_else(is.na(birdtree), taxon, birdtree))%>%
  ungroup()%>%
  group_by(genus)%>%
  mutate(mean_genus_hb = mean(hb, na.rm=T),
         mean_genus_hct = mean(hct, na.rm =T),
         mean_genus_elev = mean(elevation, na.rm=T))%>%
  ungroup()%>%
  group_by(family)%>%
  mutate(closest_to_0_hb = min(abs(SSDI_hb), na.rm=T),
         species_0_hb = species[which.min(abs(SSDI_hb))][1],
         male_0_SDI_hb = mean(mean_hb_male[which.min(abs(SSDI_hb))], na.rm=T),
         female_0_SDI_hb = mean(mean_hb_female[which.min(abs(SSDI_hb))], na.rm=T),
         closest_to_0_hct = min(abs(SSDI_hct), na.rm=T),
         species_0_hct = species[which.min(abs(SSDI_hct))][1],
         diff_from_expected_hb_male = mean_hb_male-male_0_SDI_hb,
         diff_from_expected_hb_female = mean_hb_female-female_0_SDI_hb,
         diff_from_expected_hct = hct-closest_to_0_hct,
         diff_from_expected_elev = elevation-mean_genus_elev,
         diff_from_exp_hb_elev_corr_male = diff_from_expected_hb_male + (diff_from_expected_elev/1000),
         diff_from_exp_hb_elev_corr_female = diff_from_expected_hb_female + (diff_from_expected_elev/1000))%>%#need to correct for ~ 1 unit increase in per 1000 m increase in elevation
  left_join(dimorph)%>%
  write_csv(., "data/all_blood_final.csv")


ggplot(all.blood, aes(diff_from_expected_hb_male))+
  geom_density()+
  geom_density(data=all.blood, aes(diff_from_expected_hb_female))
  

# Gonads####

#Get gonad sizes first. This will be a pain because they are irregular character strings
gonads <- readxl::read_excel("~/Dropbox/Data_repo/blood-data-for-ethans exploration-fun.xlsx")%>%
  dplyr::select(NK,Gonads, Family, DAY, MONTH, YEAR)%>%
  separate(Gonads, c("length", "width"), "x")%>%
  filter(!is.na(DAY),
         !is.na(MONTH), 
         !is.na(YEAR))%>%
  mutate(length = strex::str_nth_number(length, n = 1),
         width = strex::str_nth_number(width, n = 1),
         log10length = log10(length),
         MONTH = match(MONTH, month.name),
         DATE = paste(DAY, MONTH, YEAR, sep="/"),
         DATE = as.Date(DATE, format = "%d/%m/%Y"),
         doy = lubridate::yday(DATE))%>%
  dplyr::rename(nk = NK)%>%
  right_join(read_csv("data/blood_data_final.csv")%>%
               dplyr::select(nk, species,family, sex, elevation, mass, hb, hct))%>%
  mutate(log10mass = log10(mass),
         log10length = if_else(is.na(log10length) | log10length=="-Inf",NA, log10length),
         log10mass = if_else(is.na(log10mass) | log10mass=="-Inf",NA, log10mass))%>%
  filter(sex %in% c("male", "female"))%>%
  group_by(sex)%>%
  mutate(
    resids.gonad.length = resid(lm(log10length ~ log10mass, dplyr::cur_data() %>% dplyr::select(log10length, log10mass), na.action = na.exclude))
  )%>%
  write_csv(., "data/gonads.csv")

# Within species
# Williamson
pgcl.elevs <- read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  dplyr::select(species, sex, elev,mass, hb, hct)%>%
  mutate(
    family = "Trochilidae",
    species=gsub("_", " ", species))%>%
  filter(
    sex %in% c("male", "female"))%>%
  rename(elevation = elev)%>%
  mutate(elev_bin = if_else(elevation <=1000, "low", if_else(elevation>=3000, "high", "mid")))%>%
  bind_rows(.,
            read_csv("data/blood_data_final.csv")%>%
              dplyr::select(species, family, sex, elevation, hb, hct)%>%
              filter(
                sex %in% c("male", "female"))%>%
              mutate(elev_bin = if_else(elevation <=1000, "low", if_else(elevation>=3000, "high", "mid"))))%>%
  group_by(species, family, elev_bin)%>%
  dplyr::summarise(max_sample_elev = max(elevation, na.rm=T),
                   mean_sample_elev = mean(elevation, na.rm=T),
                   min_sample_elev = min(elevation, na.rm=T))%>%
  mutate(taxon=gsub(" ", "_", species))%>%
  na.omit()

pgcl <- read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  filter(
    sex %in% c("male", "female"))%>%
  mutate(elev_bin = if_else(elev <=1000, "low", if_else(elev>=3000, "high", "mid")),
         sex = factor(sex),
         family = "Trochilidae",
         taxon = species)%>%
  rename(elevation = elev)%>%
  group_by(taxon) %>% 
  filter(all(levels(sex) %in% sex))%>%
  group_by(taxon, sex, family) %>% 
  dplyr::select(taxon, species,family, sex, elevation, elev_bin, mass, hb, hct)%>%
  bind_rows(.,
            read_csv("data/blood_data_final.csv")%>%
              filter(
                sex %in% c("male", "female"))%>%
              mutate(elev_bin = if_else(elevation <=1000, "low", if_else(elevation>=3000, "high", "mid")),
                     sex = factor(sex),
                     taxon = factor(gsub(" ", "_", species))
              )%>%
              group_by(taxon) %>% 
              filter(all(levels(sex) %in% sex)))%>%
  dplyr::select(taxon, species,family, sex, elevation, elev_bin, mass, hb, hct)%>%
  group_by(taxon, sex, family, elev_bin) %>% 
  dplyr::summarise(mean_mass = mean(mass, na.rm=T),
                   mean_hb = mean(hb, na.rm=T),
                   mean_hct = mean(hct, na.rm=T),
                   n=n())%>%
  pivot_wider(id_cols = c(taxon, family, elev_bin), names_from = c(sex), values_from = c(mean_mass, mean_hb, mean_hct, n))%>%
  mutate(
    SSDI_hct = ((mean_hct_male- mean_hct_female)/mean_hct_female),
    SSDI_hb = ((mean_hb_male- mean_hb_female)/mean_hb_female),
    SSDI_mass = ((mean_mass_male - mean_mass_female)/mean_mass_female))%>%
  left_join(., pgcl.elevs%>%dplyr::select(-c(species, family)), by=c("taxon",  "elev_bin"))


# #Linck
# Linck.elevs <-read_csv("data/blood_data_final.csv")%>%
#   dplyr::select(species, family, sex, elevation, hb, hct)%>%
#   filter(
#     sex %in% c("male", "female"))%>%
#   mutate(elev_bin = if_else(elevation <=1000, "low", if_else(elevation>=3000, "high", "mid")))%>%
#   group_by(species, family, elev_bin)%>%
#   dplyr::summarise(max_sample_elev = max(elevation, na.rm=T),
#                    mean_sample_elev = mean(elevation, na.rm=T),
#                    min_sample_elev = min(elevation, na.rm=T))%>%
#   mutate(taxon=gsub(" ", "_", species))%>%
#   na.omit()
# 
# df.Linck <-read_csv("data/blood_data_final.csv")%>%
#   filter(
#     sex %in% c("male", "female"))%>%
#   mutate(elev_bin = if_else(elevation <=1000, "low", if_else(elevation>=3000, "high", "mid")),
#          sex = factor(sex),
#          taxon = factor(gsub(" ", "_", species))
#   )%>%
#   group_by(taxon) %>% 
#   filter(all(levels(sex) %in% sex))%>%
#   group_by(taxon, sex, family, elev_bin) %>% 
#   dplyr::summarise(mean_mass = mean(mass, na.rm=T),
#                    mean_hb = mean(hb, na.rm=T),
#                    mean_hct = mean(hct, na.rm=T),
#                    n=n())%>%
#   pivot_wider(id_cols = c(taxon, family, elev_bin), names_from = c(sex), values_from = c(mean_mass, mean_hb, mean_hct, n))%>%
#   mutate(
#     SSDI_hct = ((mean_hct_male- mean_hct_female)/mean_hct_female),
#     SSDI_hb = ((mean_hb_male- mean_hb_female)/mean_hb_female),
#     SSDI_mass = ((mean_mass_male - mean_mass_female)/mean_mass_female))%>%
#   left_join(., Linck.elevs%>%dplyr::select(-c(species, family)), by=c("taxon", "elev_bin"))


#add marshbirds
mbb.elevs <- read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/blood_final_bioclim_12-9-2022.csv")%>%
  filter(
    sex %in% c("male", "female"),
    species %in% c("Tachuris rubrigastra", "Phleocryptes melanops"))%>%
  dplyr::select(species, family, sex, elev, hb_final, hct_final)%>%
  filter(
    sex %in% c("male", "female"))%>%
  mutate(elev_bin = if_else(elev <=1000, "low", if_else(elev>=3000, "high", "mid")))%>%
  group_by(species, family, elev_bin)%>%
  dplyr::summarise(max_sample_elev = max(elev, na.rm=T),
                   mean_sample_elev = mean(elev, na.rm=T),
                   min_sample_elev = min(elev, na.rm=T))%>%
  mutate(taxon=gsub(" ", "_", species))%>%
  na.omit()


mbb <- read_csv("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/blood_final_bioclim_12-9-2022.csv")%>%
  filter(
    sex %in% c("male", "female"),
    species %in% c("Tachuris rubrigastra", "Phleocryptes melanops"))%>%
  mutate(elev_bin = if_else(elev <=1000, "low", if_else(elev>=3000, "high", "mid")),
         sex = factor(sex),
         taxon = factor(gsub(" ", "_", species))
  )%>%
  group_by(taxon) %>% 
  filter(all(levels(sex) %in% sex))%>%
  group_by(taxon, sex, family, elev_bin) %>% 
  dplyr::summarise(mean_mass = mean(mass, na.rm=T),
                   mean_hb = mean(hb_final, na.rm=T),
                   mean_hct = mean(hct_final, na.rm=T),
                   n=n())%>%
  pivot_wider(id_cols = c(taxon, family, elev_bin), names_from = c(sex), values_from = c(mean_mass, mean_hb, mean_hct, n))%>%
  mutate(
    SSDI_hct = ((mean_hct_male- mean_hct_female)/mean_hct_female),
    SSDI_hb = ((mean_hb_male- mean_hb_female)/mean_hb_female),
    SSDI_mass = ((mean_mass_male - mean_mass_female)/mean_mass_female))%>%
  left_join(., mbb.elevs%>%dplyr::select(-c(species, family)),  by=c("taxon", "elev_bin"))

df.Linck <- rbind(pgcl, mbb)

# match up tree labels
tree <- read.tree("data/birds_mcc.tre")

# tip.labs <- as.data.frame(tree$tip.label)
# 
# not.in.tree <- df.Linck%>%dplyr::select(taxon)%>%
#   filter(!taxon %in% tree$tip.label)%>%
#   write.csv(.,"data/not_in_tree_blood.csv")

#reset_names in df to match tips in tree
# created crosswalk excel file (load in)
blood.crosswalk <- read_csv("data/blood_birdtree_crosswalk.csv")



df.Linck%>%
  filter(!taxon %in% c("Elaenia_sp.", "Spinus_sp.", "Scytalopus_sp."))%>%
  left_join(., blood.crosswalk%>%dplyr::select(2:3), by ="taxon")%>%
  mutate(birdtree = if_else(is.na(birdtree), taxon, birdtree))%>%
  write_csv(., "data/blood_final_within.csv")



# # Use `anti_join()` to determine which observations are missing
# all <- fruits %>% expand(type, size, year)
# all
# all %>% dplyr::anti_join(fruits)
# 
# # Use with `right_join()` to fill in missing rows (like `complete()`)
# fruits %>% dplyr::right_join(all)
# 
# x = 1:10
# y = rep(letters[1:5],each=2)
# z = rep(1:2,length.out =10)
# df = data.frame(x,y, z)
# df = rbind(df,c(11,"e",3))
# df$verif = paste0(df$y,df$z)
# df$x = as.numeric(df$x)
# 
# 
# 
# devtools::source_gist("https://gist.github.com/brshallo/f92a5820030e21cfed8f823a6e1d56e1")

# 
# Linck_data_within <- read_csv("data/blood_data_final.csv") %>%
#   filter(sex %in% c("male", "female")) %>%
#   mutate(elev_bin = if_else(
#     elevation < 500,
#     "0-500",
#     if_else(
#       elevation >= 500 &
#         elevation < 1000,
#       "500-1000",
#       if_else(
#         elevation >= 1000 &
#           elevation < 1500,
#         "1000-1500",
#         if_else(
#           elevation >= 1500 &
#             elevation < 2000,
#           "1500-2000",
#           if_else(
#             elevation >= 2000 &
#               elevation < 2500,
#             "2000-2500",
#             if_else(
#               elevation >= 2500 &
#                 elevation < 3000,
#               "2500-3000",
#               if_else(
#                 elevation >= 3000 &
#                   elevation < 3500,
#                 "3000-3500",
#                 if_else(
#                   elevation >= 3500 &
#                     elevation < 4000,
#                   "3500-4000",
#                   if_else(elevation >= 4000 , "4000+", paste(as.character(elevation)))
#                 )
#               )
#             )
#           )
#         )
#       )
#     )
#   ),
#   elev_bin_coarse = if_else(elevation <= 1500, "low", if_else(elevation >1500 & elevation <=2200, "mid", "high")))%>% 
#   ungroup()%>%
#   group_by(species, elev_bin_coarse)%>%
#   #mutate(mean_elev = mean(elevation, na.rm=T))%>%
#   ungroup()%>%
#   group_by(species, sex, elev_bin_coarse) %>%
#   mutate(n_sp_sex = n())%>%
#   filter(n_sp_sex >= 3) %>%
#   dplyr::select(species, family, n_sp_sex, elev_bin_coarse, mass, hb, hct)%>%
#   ungroup()%>%
#   group_by(species, family, sex, elev_bin_coarse, n_sp_sex) %>%
#   dplyr::summarize(mean_hb = mean(hb, na.rm=T),
#                    mean_hct = mean(hct, na.rm=T),
#                    mean_mass = mean(mass, na.rm=T))%>%
#   pivot_wider(
#     names_from = c("sex"),
#     values_from = c("n_sp_sex", "mean_hb", "mean_hct", "mean_mass")
#   )%>% 
#   ungroup()%>%
#   group_by(species, elev_bin_coarse)%>%
#   mutate(mean.SDI.hb = (mean_hb_male - mean_hb_female)/mean_hb_female,
#          mean.SDI.hct = (mean_hct_male - mean_hct_female)/mean_hct_female,
#          mean.SDI.mass = (mean_mass_male - mean_mass_female)/mean_mass_female,
#          elev_bin_coarse = factor(elev_bin_coarse, levels=c("low","mid", "high")))
# 

#Pairwise differences dataframe test
df<-structure(list(Levelname = c("B 1", "B 2", "B 3", 
                             "B 4", "D 1", "D 2", "D 3", "D 4"), y = c(0.679428655093332, 
                                                                       1.07554328679719, 0.883000346050764, 0.791772867506205, 0.538143790501689, 
                                                                       0.805122127560562, 0.591353204313314, 0.795225886492002), fill = c("midnightblue", 
                                                                                                                                          "dodgerblue4", "steelblue3", "lightskyblue", "midnightblue", 
                                                                                                                                          "dodgerblue4", "steelblue3", "lightskyblue"), species = c("White Grunt", 
                                                                                                                                                                                                    "White Grunt", "White Grunt", "White Grunt", "White Grunt", "White Grunt", 
                                                                                                                                                                                                    "White Grunt", "White Grunt")), row.names = c(NA, -8L), class = "data.frame")
library(tidyverse)

df %>%
  group_by(grp = str_extract(Levelname, "\\w+"))%>%
  summarise(pair = combn(Levelname, 2, str_c, collapse = " - "),
            perc_diff = combn(y, 2, function(x) 200*abs(diff(x))/sum(x)),
            .groups = 'drop')

t<-read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  mutate(elev_bin = if_else(elev <=2000, "low", if_else(elev>=2000, "high", "mid")),
    sexspp = paste0(species, sex, elev_bin))%>%
  group_by(grp = str_extract(sexspp, ""))%>%
  summarise(pair = combn(sexspp, 2, str_c, collapse = " - "),
            perc_diff = combn(hb, 2, function(x) (diff(x))/x[2]),
            .groups = 'drop')%>%
  filter(pair %in% c("Patagona_chaskimalehigh - Patagona_chaskifemalehigh", "Patagona_gigasmalehigh - Patagona_gigasfemalehigh",  "Patagona_gigasmalelow - Patagona_gigasfemalelow"),
         perc_diff <0.2 &perc_diff > -0.2)%>%
  ggplot(aes(pair, perc_diff))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_discrete(labels=c("P.chaski", "P.gigas (high)", "P.gigas (low)"))
t


p<-read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  mutate(elev_bin = if_else(elev <=2000, "low", if_else(elev>=2000, "high", "mid")))%>%
  filter(sex %in% c("male", "female"),
         species =="Patagona_gigas")%>%
  ggplot(aes(elev_bin, hb, fill=sex))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitterdodge(dodge.width=0.9))+
  scale_fill_manual(values=c("navajowhite3", "grey"))+
  coord_cartesian(ylim=c(10,22))

p1<-read_csv("data/PatagonaData_JLW_HbAndHct_ForBloodSexNElev_2023-12-17.csv")%>%
  mutate(elev_bin = if_else(elev <=2000, "low", if_else(elev>=2000, "high", "mid")))%>%
  filter(sex %in% c("male", "female"),
         species =="Patagona_chaski")%>%
  ggplot(aes(elev_bin, hb, fill=sex))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitterdodge(dodge.width=0.9))+
  scale_fill_manual(values=c("navajowhite3", "grey"))+
  coord_cartesian(ylim=c(10,22))
p1
  
  
g <- grid.arrange(p, p1, nrow=1)

  
  