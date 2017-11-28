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


# Now I need to include phylogenetic variables

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


write.csv(complete_sponge_raref, file = "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/complete_sponge_raref.csv")

#Preaparing data to MPlus
library(MplusAutomation)
names(complete_sponge_raref)

names(complete_sponge_raref)



complete_sponge_raref_toMplus <- dplyr::select(complete_sponge_raref, latitude, longitude, apparent_oxygen_utilization_annual_1deg.nc,apparent_oxygen_utilization_monthly_1deg.nc,
                                               apparent_oxygen_utilization_seasonal_1deg.nc,dissolved_oxygen_annual_1deg.nc, dissolved_oxygen_monthly_1deg.nc,
                                               dissolved_oxygen_seasonal_1deg.nc,nitrate_annual_1deg.nc, nitrate_monthly_1deg.nc,nitrate_seasonal_1deg.nc,
                                               oxygen_saturation_annual_1deg.nc, oxygen_saturation_monthly_1deg.nc,oxygen_saturation_seasonal_1deg.nc,
                                               phosphate_annual_1deg.nc, phosphate_monthly_1deg.nc, phosphate_seasonal_1deg.nc, salinity_annual_1deg.nc,
                                               salinity_monthly_1deg.nc, salinity_seasonal_1deg.nc, silicate_annual_1deg.nc, silicate_monthly_1deg.nc,
                                               silicate_seasonal_1deg.nc, temperature_annual_1deg.nc, temperature_monthly_1deg.nc, temperature_seasonal_1deg.nc,
                                               Axis.1, Axis.2, Observed, Shannon, Fisher, phylo1, phylo2)
head(complete_sponge_raref_toMplus)

prepareMplusData(complete_sponge_raref_toMplus, "/Users/decio/Documents/UT-AUSTIN/13_Fall_2017/SEM course/sponges/complete_sponge_raref_toMplus.dat")


