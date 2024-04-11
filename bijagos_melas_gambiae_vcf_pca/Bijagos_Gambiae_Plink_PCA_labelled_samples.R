library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()
library(viridis)
library(scales)

workdir <- "/mnt/storage11/sophie/bijagos_mosq_wgs/bijagos_melas_gambiae_vcf/pca_plink_files" # Working directory with plink files
prefix <- "X_only_bijagos_gambiae_melasplusglobal" # Prefix for plink files
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

# Correctly reorder the dataframe to make 'sample' the first column
df <- df[, c("sample", setdiff(names(df), "sample"))]

# Rename PCA columns to start with 'PC', if they start with 'V' 
# (assuming all PCA columns currently start with 'V', adjust if different)
colnames(df)[-c(1, ncol(df)-1, ncol(df))] <- gsub("^V", "PC", colnames(df)[-c(1, ncol(df)-1, ncol(df))])

color_by <- "country" # specify if colored by region or country

# Plot PCA old cameroon colour #7678ed
#
my_colours <- c("Cameroon" = "#E74C3C", "Guinea-Bissau" = "#3d348b", "The Gambia" = "#f7b801")
my_shapes <- c("anopheles melas" = 16 , "anopheles gambiae ss" = 4)

png("X_only_gambiae_melas_labelled_countryspecies.png") 
ggplot(data = df, aes(x = PC1, y = PC2, color = !!sym(color_by), shape = species)) +
    geom_point() +
    geom_text(aes(label = sample), vjust = 1.5, size = 2, angle = 20, check_overlap = TRUE) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"), y = paste0("PC2", " (", vars["PC2"], "%)"), title = "Chromosome X") +
    scale_color_manual(values = my_colours) + 
    scale_shape_manual(values = my_shapes) + 
    theme_classic() +
    theme(legend.position = "bottom", 
          legend.box = "vertical", # Explicitly set legends to stack vertically
          plot.title = element_text(hjust = 0.5),
          legend.box.margin = margin(0, 0, 0, 0)) # Adjust legend box margins if needed
dev.off()
#
