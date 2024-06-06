#Impute Trait values using random forest####
library(TrEvol)
library(funspace)

#Read in data ####

df <- read_csv("data/species_means_blood_final.csv")

#make trimmed datatset with each variable we want to impute
df.trimmed <- df%>%
  dplyr::select(species, mean_mass_male, mean_mass_female, egg_mass, egg_length, log_generation_length, log_clutch_size)%>%
  mutate(species = gsub(" ", "_", species),
         mean_mass_male = log10(mean_mass_male),
         mean_mass_female = log10(mean_mass_female),
         egg_mass = log10(egg_mass),
         egg_length = log10(egg_length))%>%
  column_to_rownames(var="species")


#Tree wrangle ####
#read in tree
mtree <- read.nexus("data/birdtree_samples_blood.nex")[[1]]

imputed.df.list <- impute(df.trimmed, phylo = mtree)

ggplot(imputed.df.list$imputed, aes(mean_mass_female, egg_mass))+
  geom_point()







#Now prune trees to match each dataframe
for(i in 1:length(df.list)){
  
 tt <- keep.tip(mtree, df.list[[i]]%>%dplyr::pull(animal))
 
 if(length(tt$tip.label) < nrow(df.list[[i]])){
   df.list[[i]] <- df.list[[i]]%>%
     filter(species %in% tt$tip.label)
 }
 
 variance_covariance_results.list <- computeVarianceCovariancePartition(
   traits =  colnames(df.list[[i]])[-1],
   dataset = df.list[[i]],
   terminal_taxa = "animal",
   phylogeny = mtree
   
 )
 
  
  
}




variance_covariance_results.list <- computeVarianceCovariancePartition(
  traits =  c("nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "phylo_G1_trait1", "phylo_G1_trait2"),
  dataset = simulated.data$data,
  terminal_taxa = "animal",
  phylogeny = simulated.data$phylogeny
)