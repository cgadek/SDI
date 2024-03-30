# SDI (blood) Data Wrangle####

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
              dplyr::select(species, family, sex, elevation, hb, hct))%>%
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


#Myhrvold et al. 2015####
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

#John T Rotenberry, Priya Balasubramaniam egg dataset extrapolated####
#https://academic-oup-com.libproxy.unm.edu/auk/article/137/3/ukaa019/5834541?login=true&token=#supplementary-data
egg_mass_extrap<- read_excel("data/ukaa019_suppl_supplemental_material_tables_s4-s6.xlsx")%>%
  dplyr::select(Species, egg_mass_Grams)%>%
  dplyr::rename(species = Species,
         egg_mass_extrap = egg_mass_Grams)

 
#Development dataset ####

development<- readxl::read_excel("data/41467_2020_16257_MOESM3_ESM.xlsx", sheet=4)%>%
  dplyr::rename(taxon =  binomial)


#read in generation time/demgrpahic dataset from Bird et al. 2020

demo <-readxl::read_excel("data/Bird_et_al_2020.xlsx", sheet=1)


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
  left_join(demo)%>%
  left_join(egg_mass_extrap)%>%
  mutate(genus = stringr::word(species))%>%
  left_join(Myhrvold.genus.egg.traits)


#test for normality between blood parameters and SDIs before we can use residuals
# hb.f.mod <- lm(SSDI_hb ~ mean_hb_female, data = df.blood)
# 
# shapiro.test(rstandard(hb.f.mod)) #Not normal
# 
# 
# hb.m.mod <- lm(SSDI_hb ~ mean_hb_male, data = df.blood)
# 
# shapiro.test(rstandard(hb.m.mod)) #Not normal
# 
# 
# hb.f.mass.mod <- lm(SSDI_hb ~ mean_mass_female, data = df.blood)
# 
# shapiro.test(rstandard(hb.f.mass.mod)) #Not normal
# 
# 
# hb.m.mass.mod <- lm(SSDI_hb ~ mean_mass_male, data = df.blood)
# 
# shapiro.test(rstandard(hb.m.mass.mod)) #Not normal

# 
# ggplot() +
#   geom_qq(aes(sample = rstandard(hb.f.mod))) +
#   geom_abline(color = "red") +
#   coord_fixed()
# 
# ggplot() +
#   geom_qq(aes(sample = rstandard(hb.m.mod))) +
#   geom_abline(color = "red") +
#   coord_fixed()

#both have fat tails....

#Use student t in brms
# hb.f.mod <- brm(
#   formula = SSDI_hb ~ mean_hb_female,
#   data = df.blood,
#   family = student(),
#   cores = 8,
#   chains = 4,
#   thin = 10,
#   warmup = 5000,
#   # default is iter/2; shouldn't ever be larger than iter
#   iter = 10000,
#   control = list(adapt_delta = 0.98, max_treedepth = 18),
#   sample_prior = TRUE # default priors
# )
# 
# summary(hb.f.mod)
# 
# hb.f.mod.res <- residuals(hb.f.mod)

df.blood <- df.blood %>%
  dplyr::select(species, SSDI_hb, mean_hb_female) %>%
  na.omit() %>%
  bind_cols(., hb.f.mod.res[, 1]) %>%
  dplyr::rename(female_SDI_hb_residuals = length(.)) %>%
  ungroup() %>%
  dplyr::select(species, female_SDI_hb_residuals) %>%
  right_join(df.blood)


# hb.m.mod <- brm(
#   formula = SSDI_hb ~ mean_hb_male,
#   data = df.blood,
#   family = student(),
#   cores = 8,
#   chains = 4,
#   thin = 10,
#   warmup = 5000,
#   # default is iter/2; shouldn't ever be larger than iter
#   iter = 10000,
#   control = list(adapt_delta = 0.98, max_treedepth = 18),
#   sample_prior = TRUE # default priors
# )
# 
# summary(hb.m.mod)
# hb.m.mod.res <- residuals(hb.m.mod)

df.blood <- df.blood %>%
  dplyr::select(species, SSDI_hb, mean_hb_female) %>%
  na.omit() %>%
  bind_cols(., hb.m.mod.res[, 1]) %>%
  dplyr::rename(male_SDI_hb_residuals = length(.)) %>%
  ungroup() %>%
  dplyr::select(species, male_SDI_hb_residuals) %>%
  right_join(df.blood)

#mass ~ Hb

# hb.f.mass.mod <- brm(
#   formula = mean_hb_female ~ mean_mass_female,
#   data = df.blood,
#   family = student(),
#   cores = 8,
#   chains = 4,
#   thin = 10,
#   warmup = 5000,
#   # default is iter/2; shouldn't ever be larger than iter
#   iter = 10000,
#   control = list(adapt_delta = 0.98, max_treedepth = 18),
#   sample_prior = TRUE # default priors
# )

summary(hb.f.mass.mod)

hb.f.mass.mod.res <- residuals(hb.f.mass.mod)

df.blood <- df.blood %>%
  dplyr::select(species, mean_hb_female, mean_mass_female) %>%
  na.omit() %>%
  bind_cols(., hb.f.mass.mod.res[, 1]) %>%
  dplyr::rename(female_mass_residuals = length(.)) %>%
  ungroup() %>%
  dplyr::select(species, female_mass_residuals) %>%
  right_join(df.blood)


# hb.m.mass.mod <- brm(
#   formula = mean_hb_male ~ mean_mass_male,
#   data = df.blood,
#   family = student(),
#   cores = 8,
#   chains = 4,
#   thin = 10,
#   warmup = 5000,
#   # default is iter/2; shouldn't ever be larger than iter
#   iter = 10000,
#   control = list(adapt_delta = 0.98, max_treedepth = 18),
#   sample_prior = TRUE # default priors
# )

summary(hb.m.mass.mod)
hb.m.mass.mod.res <- residuals(hb.m.mass.mod)

df.blood <- df.blood %>%
  dplyr::select(species, mean_hb_male, mean_mass_male) %>%
  na.omit() %>%
  bind_cols(., hb.m.mass.mod.res[, 1]) %>%
  dplyr::rename(male_mass_residuals = length(.)) %>%
  ungroup() %>%
  dplyr::select(species, male_mass_residuals) %>%
  right_join(df.blood)

# match up tree labels
tree <- read.tree("data/birds_mcc.tre")

# created crosswalk excel file (load in)
blood.crosswalk <-  read_csv("data/blood_birdtree_crosswalk.csv")
    
#load dimorphism data
dimorph<- readxl::read_excel("data/Taxon_names_for_Jessie_WithFamily_DimorphismDone.xlsx")%>%
  dplyr::rename(taxon= taxon...1)%>%
  dplyr::select(taxon, family, dimorphic_01)

df.blood<-df.blood%>%
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
    diff_from_exp_hb_elev_corr_female_base  = diff_from_expected_hb_female_base  + (diff_from_expected_elev /1000),
     # this hb_cont is generated from the species in each group closest to zero SDI                                                                                               1000),
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
    hb_cont_global = if_else(
      diff_from_exp_hb_elev_corr_female_global_mean > 0 &
        abs(diff_from_exp_hb_elev_corr_female_global_mean) > abs(diff_from_exp_hb_elev_corr_male_global_mean),
      "female driven increase",
      if_else(
        diff_from_exp_hb_elev_corr_male_global_mean > 0 &
          abs(diff_from_exp_hb_elev_corr_male_global_mean) > abs(diff_from_exp_hb_elev_corr_female_global_mean),
        "male driven increase",
        if_else(
          diff_from_exp_hb_elev_corr_female_global_mean < 0 &
            abs(diff_from_exp_hb_elev_corr_female_global_mean) > abs(diff_from_exp_hb_elev_corr_male_global_mean),
          "female driven reduction",
          if_else(
            diff_from_exp_hb_elev_corr_male_global_mean < 0 &
              abs(diff_from_exp_hb_elev_corr_male_global_mean) > abs(diff_from_exp_hb_elev_corr_female_global_mean),
            "male driven reduction",
            "at reference"
          )
        )
      )
    ),
    hb_cont_base = if_else(
      diff_from_exp_hb_elev_corr_female_base > 0 &
        abs(diff_from_exp_hb_elev_corr_female_base) > abs(diff_from_exp_hb_elev_corr_female_base),
      "female driven increase",
      if_else(
        diff_from_exp_hb_elev_corr_male_base > 0 &
          abs(diff_from_exp_hb_elev_corr_male_base) > abs(diff_from_exp_hb_elev_corr_female_base),
        "male driven increase",
        if_else(
          diff_from_exp_hb_elev_corr_female_base < 0 &
            abs(diff_from_exp_hb_elev_corr_female_base) > abs(diff_from_exp_hb_elev_corr_male_base),
          "female driven reduction",
          if_else(
            diff_from_exp_hb_elev_corr_male_base < 0 &
              abs(diff_from_exp_hb_elev_corr_male_base) > abs(diff_from_exp_hb_elev_corr_female_base),
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
  )%>% write_csv("data/all_blood_cont.csv")


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


  