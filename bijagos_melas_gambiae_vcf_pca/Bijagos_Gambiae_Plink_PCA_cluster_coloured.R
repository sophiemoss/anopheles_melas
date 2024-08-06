library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()
library(viridis)
library(scales)
library(stringr)

workdir <- "/mnt/storage11/sophie/bijagos_mosq_wgs/bijagos_melas_gambiae_vcf/pca_plink_files" # Working directory with plink files
prefix <- "mito_only_bijagos_gambiae_melasplusglobal" # Prefix for plink files
metadata <- "metadata_gambiae_melas.csv" # File path to metadata

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "sample"))

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA #
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling - might take a while
# saveRDS(cmd, paste0(prefix, ".dist.rds") # save to RDS format
#cmd <- readRDS(file.path(workdir, paste0(prefix, ".dist.rds"))
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Overlay region, country info
# Convert cmdscale output to a dataframe and set row names as a 'sample' column
df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
df$sample <- rownames(df)  # Set rownames as the 'sample' column

# Add metadata columns
df$country <- gsub("_", " ", desc$country)
df$island <- gsub("_", " ", desc$island)
df$species <- gsub("_", " ", desc$species)
df$cluster <- gsub("_", " ", desc$cluster)

# Correctly reorder the dataframe to make 'sample' the first column
df <- df[, c("sample", setdiff(names(df), "sample"))]

# Rename PCA columns to start with 'PC', if they start with 'V' 
# (assuming all PCA columns currently start with 'V', adjust if different)
colnames(df)[-c(1, ncol(df)-1, ncol(df))] <- gsub("^V", "PC", colnames(df)[-c(1, ncol(df)-1, ncol(df))])

color_by <- "cluster" # specify if colored by region or country

### Plot PCA old cameroon colour #7678ed
#
#my_colours <- c("An.melas.Cameroon" = "#d55e00", "An.melas.The Gambia" = "#cc79a7", 
#                "An.melas.groupA" = "#0072b2", "An.melas.groupB" = "#f0e442", 
#                "An.gambiae" = "#009e73")
##
#png("3R_only_gambiae_melas_labelled_countryspecies.png") 
#ggplot(data = df, aes(x = PC1, y = PC2, color = !!sym(color_by), shape = species)) +
#    geom_point() +
#    geom_text(aes(label = sample), vjust = 1.5, size = 2, angle = 20, check_overlap = TRUE) +
#    labs(x = paste0("PC1", " (", vars["PC1"], "%)"), y = paste0("PC2", " (", vars["PC2"], "%)"), title = "Chromosome X") +
#    scale_color_manual(values = my_colours) + 
#    scale_shape_manual(values = my_shapes) + 
#    theme_classic() +
#    theme(legend.position = "bottom", 
#          legend.box = "vertical", # Explicitly set legends to stack vertically
#          plot.title = element_text(hjust = 0.5),
#          legend.box.margin = margin(0, 0, 0, 0)) # Adjust legend box margins if needed
#dev.off()
#
# plot by cluster
color_by <- "cluster"

# Define colors
my_colours <- c("An.melas.Cameroon" = "#7678ed", "An.melas.The Gambia" = "darkorange2", 
                "An.melas.groupA" = "deeppink", "An.melas.groupB" = "#f0e442",
                "An.gambiae" = "#009e73")

# The plot command
png("mito_only_gambiae_melas_clusters.png", width = 800, height = 800) 
ggplot(data = df, aes(x = PC1, y = PC2, color = cluster)) + 
    geom_point() +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"), 
         y = paste0("PC2", " (", vars["PC2"], "%)"), 
         title = "Mitochondria") +
    scale_color_manual(values = my_colours) +  # Wrap legend text
    theme_classic() +
    theme(legend.position = "bottom", 
          legend.direction = "vertical",  # Change legend orientation
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 10, r = 80, b = 80, l = 40, unit = "pt"),  # Increase the bottom margin
          legend.text = element_text(size = 8))
dev.off()