## There are too many data points using pairwise data and the server couldn't handle
# I'll use axes of the ordination as the variables instead of the pairwise distances

library(phyloseq)
library(tidyverse)
library(raster)
library(ncdf4)

# making a list with all the BIOM files
otu_tables <- list.files(path = "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/2ndStudy/BIOM/", pattern=".biom", full.names = T)
otu_tables

# Loop to read all the files
all_otu_tables <- list()

for (i in 1:length(otu_tables)) {
  all_otu_tables[i] <- import_biom(otu_tables[i])
}

# merging BIOM files
all_phyloseq <- merge_phyloseq(all_otu_tables[[1]], all_otu_tables[[2]], all_otu_tables[[3]], all_otu_tables[[4]],
                               all_otu_tables[[5]], all_otu_tables[[6]], all_otu_tables[[7]], all_otu_tables[[8]],
                               all_otu_tables[[9]])
all_phyloseq

# Metadata
sample.data <- read.csv("/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/2ndStudy/metadata_csv.csv", h=T)
head(sample.data)

# checking a few of the sample names
head(sort(sample.data$sample_name))
head(sort(sample_names(all_phyloseq)))

# are they factors or characters?
class(head(sample.data$sample_name))
class(head(sample_names(all_phyloseq)))

# as.character
sample.data$sample_name <- as.character(sample.data$sample_name)

# keeping only the samples with metadata information available
pruned <- prune_samples(sample.data$sample_name, all_phyloseq)
all_phyloseq
pruned

# standardizing the names

sort(unique(sample.data$host_scientific_name))

sample.data <- mutate(sample.data, std_host_name=paste0(gsub("([A-Za-z]+).*", "\\1", sample.data$host_scientific_name),"_",
                                                        gsub("^.* ([[:print:]]+)$", "\\1", sample.data$host_scientific_name)))




# read names for phylogeny
phylo_names <- read.csv(file="~/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/sponges_for_phylogeny.csv", h=T)
head(phylo_names)

to_join <- data_frame(std_host_name=phylo_names$std_host_name, keep="yes")

sample.data <- left_join(sample.data, to_join, by="std_host_name")

metadata <- sample_data(sample.data)     
head(metadata)

sample_names(metadata) <- metadata$sample_name

merged <- merge_phyloseq(pruned, metadata)

merged
sample_variables(merged)
head(sample_data(merged))


final_sponge <- prune_samples(sample_data(merged)$keep=="yes", merged)
final_sponge

sample_data(final_sponge)$std_host_name


# NOAA

coord <- data.frame(lon = get_variable(final_sponge, "longitude"), lat = get_variable(final_sponge, "latitude"))
head(coord)


paths_to_files <- list.files(path = "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/noaa", pattern=".nc", full.names = T)
paths_to_files

names <- list.files(path = "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/noaa", pattern=".nc")

rasters <- list()

for (i in 1:length(paths_to_files)) {
  name <- names[i]
  rasters[[name]] <- raster(paths_to_files[i])
}

names(rasters)


points <- SpatialPoints(coord)
par(mfrow=c(2,2))
plot(rasters[[1]], main=names(rasters[1]))
plot(points, add=T)

plot(rasters[[4]], main=names(rasters[4]))
plot(points, add=T)

plot(rasters[[7]], main=names(rasters[7]))
plot(points, add=T)

plot(rasters[[10]], main=names(rasters[10]))
plot(points, add=T)

extracts <- list()

for (i in 1:length(rasters)) {
  name <- names[i]
  extracts[[name]] <- raster::extract(rasters[[i]], coord, buffer=200000, fun=median)
}

variables <- as.data.frame(extracts)
dim(variables)
head(variables)
final_sponge

sample_data(final_sponge) <- cbind(sample_data(final_sponge), variables)
sample_variables(final_sponge)
final_sponge

final_sponge_no_na <- prune_samples(!is.na(sample_data(final_sponge)$apparent_oxygen_utilization_annual_1deg.nc), final_sponge)
sample_data(final_sponge_no_na)
final_sponge_no_na

write.csv(unique(sample_data(no_na)$std_host_name), file = "names_for_mesquite.csv")

library(ape)
library(ggtree)
library(phytools)
tree_text <- "((Plakortis_halichondrioides,Plakortis_sp.),(((Petrosia_ficiformis,(((Xestospongia_proxima,Xestospongia_sp.,Xestospongia_bocatorensis,(Xestospongia_muta,Xestospongia_testudinaria)),((Haliclona_oculata,Haliclona_mediterranea,Haliclona_mucosa,Haliclona_fascigera,Haliclona_fulva,Haliclona_maravillosa,Haliclona_vansoesti,Haliclona_sp.,Haliclona_tubifera),(Chalinula_molitba,Chalinula_sp.))),((Niphates_digitalis,Niphates_erecta),((Amphimedon_erina,Amphimedon_sp.),Amphimedon_compressa)))),((((Geodia_sp.,Geodia_barretti),((Erylus_sp.,Erylus_formosus),Ecionemia_alata)),(Cinachyrella_sp.,Cinachyrella_alloclada,Cinachyrella_levantinensis)),((Ectyoplasia_ferox,Phakellia_ventilabrum),((((Cliona_viridis,Cliona_delitrix,Cliona_celata,Spheciospongia_vagabunda,Placospongia_intermedia),((Suberites_diversicolor,Halichondria_panicea),Halichondria_melanadocia,Halichondria_magniconulosa,Halichondria_phakellioides)),((((Coelocarteria_singaporensis,(Phorbas_sp.,Phorbas_fictitius,Phorbas_tenacior)),((Tedania_ignis,Tedania_sp.,Tedania_klausi),(Lissodendoryx_colombiensis,Lissodendoryx_sp.))),((Mycale_laevis,Mycale_grandis,(Mycale_lingua,Mycale_laxissima)),(Iotrochota_sp.,Iotrochota_birotulata))),Crambe_crambe)),((Cymbastela_coralliophila,Cymbastela_sp.),((Stylissa_sp.,Stylissa_flabelliformis,Stylissa_massa),(Axinella_damicornis,Axinella_corrugata,Axinella_verrucosa,Axinella_rugosa,Axinella_sp.,Axinella_infundibuliformis,Axinella_polypoides))))))),(((Dysidea_avara,Dysidea_etheria,Dysidea_fragilis,Dysidea_sp.),((Rhopaloeides_odorabile,(Hyrtios_sp.,Hyrtios_altum,Hyrtios_violaceus)),(((((Sarcotragus_fasciculatus,Sarcotragus_sp.,Sarcotragus_spinosulus),Ircinia_variabilis),Ircinia_strobilina),(Ircinia_felix,Ircinia_oros),Ircinia_sp.),Carteriospongia_foliascens))),((Chondrilla_caribensis,Chondrilla_nucula),(Ianthella_basta,((Aplysina_fistularis,Aplysina_sp.,((Aplysina_cavernicola,Aplysina_aerophoba),(Aplysina_archeri,Aplysina_cauliformis))),(Aiolochroia_sp.,Aiolochroia_crassa)))))));"
#removed special characters
tree <- phytools::read.newick(text=tree_text)
par(mfrow=c(1,1))
plot(tree)
cbind(sort(tree$tip.label),sort(unique(sample_data(final_sponge_no_na)$std_host_name)))


####################
#removing singletons
#just selecting taxa I want to keep
final_sponge_no_na_ns <- prune_taxa(taxa_sums(final_sponge_no_na)>1, final_sponge_no_na)

#counting how many reads per sample there are
sample_sums(final_sponge_no_na_ns)

# Now it's time to remove the Archaea and the OTUs that were not assigned to anything
# A simple way to do that is to select only the Bacteria
# I'll also remove samples that had less than 10000 reads
head(tax_table(final_sponge_no_na_ns))
final_sponge_no_na_ns_10k <- subset_taxa(final_sponge_no_na_ns, Rank1=="k__Bacteria")
final_sponge_no_na_ns_10k <- prune_samples(sample_sums(final_sponge_no_na_ns_10k)>=10000, final_sponge_no_na_ns_10k)
final_sponge_no_na_ns_10k

#######RAREFIED DATA
#Below I'll run the analyses with rarefied data

final_sponge_no_na_ns_10k_raref <- rarefy_even_depth(final_sponge_no_na_ns_10k, rngseed = T)
final_sponge_no_na_ns_10k_raref
sample_sums(final_sponge_no_na_ns_10k_raref)

final_sponge_no_na_ns_10k_raref_relab <- transform_sample_counts(final_sponge_no_na_ns_10k_raref, function(x) x / sum(x) )

head(sample_data(final_sponge_no_na_ns_10k_raref_relab))
unique(sample_data(final_sponge_no_na_ns_10k_raref_relab)$sample_type)

#Network
plot_network(make_network(final_sponge_no_na_ns_10k_raref_relab, type="samples", distance="bray", max.dist = 0.8),
             final_sponge_no_na_ns_10k_raref_relab, color="std_host_name", line_weight = 0.4, label=NULL)


#Ordinations 
sponge_pcoa_relab <- ordinate(final_sponge_no_na_ns_10k_raref_relab, "MDS", "bray")
plot_ordination(final_sponge_no_na_ns_10k_raref_relab, sponge_pcoa_relab, type = "samples", color = "std_host_name") + geom_point(size = 5)

sponge_pcoa_raref <- ordinate(final_sponge_no_na_ns_10k_raref, "MDS", "bray")
plot_ordination(final_sponge_no_na_ns_10k_raref, sponge_pcoa_raref, type = "samples", color = "std_host_name") + geom_point(size = 5)

## I'll not use the relative abundance data because I can't calculate diversity measures with it

str(sponge_pcoa_raref)
dim(sponge_pcoa_raref$vectors[,1:2])
pcoa_coords_raref <- data.frame(sponge_pcoa_raref$vectors[,1:2])
pcoa_coords_raref <- mutate(pcoa_coords_raref, sample_name=rownames(pcoa_coords_raref))
head(pcoa_coords_raref)

diversity_raref <- estimate_richness(final_sponge_no_na_ns_10k_raref)
names(diversity_raref)

names(sample_data(final_sponge_no_na_ns_10k_raref))

complete_sponge_raref <- left_join(sample_data(final_sponge_no_na_ns_10k_raref), pcoa_coords_raref, by="sample_name")
complete_sponge_raref <- cbind(complete_sponge_raref, diversity_raref)
dim(complete_sponge_raref)
names(complete_sponge_raref)
class(complete_sponge_raref)
head(complete_sponge_raref)


# Now I need to include spatial variables (PCNMs) and phylogenetic variables

## Calculating node distance
library(adephylo)
node_dist <- distTips(tree)
class(node_dist)
pcoa_phylo <- ape::pcoa(node_dist)
biplot(pcoa_phylo)

pcoa_phylo_vectors <- pcoa_phylo$vectors[,1:2]
colnames(pcoa_phylo_vectors) <- c("phylo1", "phylo2")
pcoa_phylo_vectors <- mutate(as.data.frame(pcoa_phylo_vectors), std_host_name = rownames(pcoa_phylo_vectors))

complete_sponge_raref <- left_join(complete_sponge_raref, pcoa_phylo_vectors, by="std_host_name")

# Next I need the PCNMs
library(vegan)
lat_long_sponge <- dplyr::select(complete_sponge_raref, sample_name, latitude, longitude)
head(lat_long_sponge)

library(fossil)
dists <- earth.dist(data.frame(longitude=lat_long_sponge$longitude, latitute=lat_long_sponge$latitude), dist=F)
class(dists)

pcnm_sponges <- pcnm(dists)
pcnm_sponges_scores <- data.frame(scores(pcnm_sponges), sample_name=lat_long_sponge$sample_name)

ordisurf(lat_long_sponge[,2:3], scores(pcnm_sponges, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(lat_long_sponge[,2:3], scores(pcnm_sponges, choi=2), bubble = 4, main = "PCNM 2")
ordisurf(lat_long_sponge[,2:3], scores(pcnm_sponges, choi=3), bubble = 4, main = "PCNM 3")

## forward selection

# A function to make the OTU table ready to vegan
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

sponge_vegan <- veganotu(final_sponge_no_na_ns_10k_raref)

rda.1 <- rda(sponge_vegan ~ ., data=dplyr::select(pcnm_sponges_scores, -sample_name))
plot(rda.1)
anova.sponges.cca.all <- anova.cca(rda.1)
anova.sponges.cca.all

formula(rda.1)

#rda_ordistep <- rda(skin_1000_raref_vegan ~ Axis.1 + Axis.2 + Axis.3 + Axis.4 + Axis.5 + Axis.6 + Axis.7 + Axis.8
#                    + Axis.9 + Axis.10 + Axis.11 + Axis.12 + Axis.13 + Axis.14 + Axis.15 + Axis.16, data=pcoa_axes)

rda.0 <- rda(sponge_vegan ~ 1, data=dplyr::select(pcnm_sponges_scores, -sample_name))

forward_sel <- ordistep(rda.0, scope= formula(rda.1))
#taking forever to run
forward_sel

only_ordistep <- ordistep(rda_ordistep, pstep=1000)

mod0 <- rda(skin_1000_raref_vegan ~ 1, pcoa_axes)
mod1 <- rda(skin_1000_raref_vegan ~., pcoa_axes)

or2step_skin <- ordiR2step(mod0, mod1, pstep=1000)
or2step_skin$anova
rownames(or2step_skin$anova)



# calculate pairwise bray curtis distances
#vegdist(otu_table(final_sponge_no_na_ns_10k_raref), method="bray")

# A function to make the OTU table ready to vegan
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

sponge_vegan <- veganotu(final_sponge_no_na_ns_10k_raref)

sponge_vegdist <- vegdist(sponge_vegan, method="bray")

library(reshape2)
pairwise_vegdist <- melt(as.matrix(sponge_vegdist), value.name="dist.bray")
head(pairwise_vegdist)
dim(pairwise_vegdist)

# Next I'll calculate pairwise geographic distance
lat_long_sponge <- dplyr::select(complete_sponge_raref, sample_name, latitude, longitude)
head(lat_long_sponge)

library(fossil)
dists <- earth.dist(data.frame(lat_long_sponge[,3,drop=F], latitute=lat_long_sponge$latitude), dist=F)
class(dists)
pairwise_dists <- melt(dists, value.name="dist_km")
head(pairwise_dists)
dim(pairwise_dists)



## Calculating node distance
library(adephylo)
node_dist <- distTips(tree)
class(node_dist)
as.matrix(node_dist)
species_dist <- subset(melt(as.matrix(node_dist)), value!=0)
head(species_dist)

head(pairwise_vegdist)

first_join <-complete_sponge_raref %>% dplyr::select(sample_name, std_host_name) %>% mutate(Var1=sample_name, sp_var1=std_host_name)
second_join <-complete_sponge_raref %>% dplyr::select(sample_name, std_host_name) %>% mutate(Var2=sample_name, sp_var2=std_host_name)

pairwise_vegdist_species <- left_join(pairwise_vegdist, first_join, by="Var1")
head(pairwise_vegdist_species)
pairwise_vegdist_species <- left_join(pairwise_vegdist_species, second_join, by="Var2")
head(pairwise_vegdist_species)

pairwise_vegdist_species_filtered <- pairwise_vegdist_species %>% 
  mutate(sp_sp=paste0(pairwise_vegdist_species$sp_var1, pairwise_vegdist_species$sp_var2)) %>%
  dplyr::select(Var1, Var2, dist.bray, sp_sp)
head(pairwise_vegdist_species_filtered)

#now I have to merge based on the distance between species
head(species_dist)

species_dist_to_merge <- species_dist %>% mutate(sp_sp=paste0(species_dist$Var1, species_dist$Var2), node_dist=species_dist$'value') %>%
  dplyr::select(sp_sp, node_dist)
head(species_dist_to_merge)

# and merge them all

bray_and_node <- left_join(pairwise_vegdist_species_filtered, species_dist_to_merge, by="sp_sp")
head(bray_and_node)

#check if NAs are for same species distance
head(filter(bray_and_node, is.na(node_dist)==T))
bray_and_node[is.na(bray_and_node$node_dist)==T,]$node_dist <- 0

head(bray_and_node)


# Lastly, I need to calculate the distance based on all the other variables

### Now I'll select the variables I'm going to use forward and calculate pairwise distances between them
names(complete_sponge_raref)
sample_info <- c("sample_name", "std_host_name")
numeric <-  c("apparent_oxygen_utilization_annual_1deg.nc","apparent_oxygen_utilization_monthly_1deg.nc",
              "apparent_oxygen_utilization_seasonal_1deg.nc","dissolved_oxygen_annual_1deg.nc", "dissolved_oxygen_monthly_1deg.nc",
              "dissolved_oxygen_seasonal_1deg.nc","nitrate_annual_1deg.nc", "nitrate_monthly_1deg.nc","nitrate_seasonal_1deg.nc",
              "oxygen_saturation_annual_1deg.nc", "oxygen_saturation_monthly_1deg.nc","oxygen_saturation_seasonal_1deg.nc",
              "phosphate_annual_1deg.nc", "phosphate_monthly_1deg.nc", "phosphate_seasonal_1deg.nc", "salinity_annual_1deg.nc",
              "salinity_monthly_1deg.nc", "salinity_seasonal_1deg.nc", "silicate_annual_1deg.nc", "silicate_monthly_1deg.nc",
              "silicate_seasonal_1deg.nc", "temperature_annual_1deg.nc", "temperature_monthly_1deg.nc", "temperature_seasonal_1deg.nc",
              "Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

df_for_dists <- complete_sponge_raref %>% mutate_at(numeric, funs(as.character)) %>% mutate_at(numeric, funs(as.numeric))
head(df_for_dists)
rownames(df_for_dists) <- rownames(complete_sponge_raref)

#I'm also going to select sample_name
df_for_dists %>% dplyr::select(sample_name, apparent_oxygen_utilization_annual_1deg.nc,apparent_oxygen_utilization_monthly_1deg.nc,
                               apparent_oxygen_utilization_seasonal_1deg.nc,dissolved_oxygen_annual_1deg.nc, dissolved_oxygen_monthly_1deg.nc,
                               dissolved_oxygen_seasonal_1deg.nc,nitrate_annual_1deg.nc, nitrate_monthly_1deg.nc,nitrate_seasonal_1deg.nc,
                               oxygen_saturation_annual_1deg.nc, oxygen_saturation_monthly_1deg.nc,oxygen_saturation_seasonal_1deg.nc,
                               phosphate_annual_1deg.nc, phosphate_monthly_1deg.nc, phosphate_seasonal_1deg.nc, salinity_annual_1deg.nc,
                               salinity_monthly_1deg.nc, salinity_seasonal_1deg.nc, silicate_annual_1deg.nc, silicate_monthly_1deg.nc,
                               silicate_seasonal_1deg.nc, temperature_annual_1deg.nc, temperature_monthly_1deg.nc, temperature_seasonal_1deg.nc,
                               Observed, Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher) -> df_for_dists.numeric
head(df_for_dists.numeric)
rownames(df_for_dists.numeric) <- df_for_dists.numeric$sample_name
df_for_dists.numeric <- dplyr::select(df_for_dists.numeric, -sample_name)

melt(as.matrix(dist(df_for_dists.numeric[,1, drop=F])), value.name=colnames(df_for_dists.numeric)[1])

#this one works!!
dist.apply <- lapply(1:ncol(df_for_dists.numeric), function(i) melt(as.matrix(dist(df_for_dists.numeric[,i, drop=F])), value.name=colnames(df_for_dists.numeric)[i]))

data.frame.dist <- as.data.frame(dist.apply)
head(data.frame.dist)
data.frame.dist.numbers <- select_if(data.frame.dist, is.numeric)
head(data.frame.dist.numbers)
data.frame.dist.numbers <- cbind(data.frame.dist.numbers, data.frame.dist[,1:2])
head(data.frame.dist.numbers)

head(bray_and_node)

data.frame.dist.numbers <- data.frame.dist.numbers %>% mutate(unique_ident = paste0(data.frame.dist.numbers$Var1, data.frame.dist.numbers$Var2))
bray_and_node <- bray_and_node %>% mutate(unique_ident = paste0(bray_and_node$Var1, bray_and_node$Var2)) %>% dplyr::select(-Var1, -Var2)

final_sponge_dists <- left_join(data.frame.dist.numbers, bray_and_node, by="unique_ident")
head(final_sponge_dists)
head(pairwise_dists)

dim(final_sponge_dists)
dim(pairwise_dists)

final_sponge_dists <- mutate(final_sponge_dists, dist_km =pairwise_dists$dist_km)
head(final_sponge_dists)
str(final_sponge_dists)

final_sponge_dists_std <- final_sponge_dists %>% select_if(is.numeric) %>% mutate_each(funs(scale))
final_sponge_dists_std <- cbind(final_sponge_dists_std, dplyr::select(final_sponge_dists, Var1, Var2, unique_ident, sp_sp))
head(final_sponge_dists_std)
dim(final_sponge_dists_std)

boxplot(select_if(final_sponge_dists_std, is.numeric))

write.table(final_sponge_dists_std, file = "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/final_sponge_distances.txt")


library(MplusAutomation)
names(final_sponge_dists_std)
final_sponge_dists_std_toMplus <- final_sponge_dists_std[,1:36]
head(final_sponge_dists_std_toMplus)

prepareMplusData(final_sponge_dists_std_toMplus, "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/final_sponge_dists_std.dat")


plot(final_sponge_dists_std_toMplus$node_dist, final_sponge_dists_std_toMplus$dist.bray)
plot(final_sponge_dists_std_toMplus$dist_km, final_sponge_dists_std_toMplus$dist.bray)



