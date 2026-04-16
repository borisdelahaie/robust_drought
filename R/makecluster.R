library(rlang)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(cluster)
library(patchwork)


### Load growth data
growth <- read.csv("data_raw/Growth_data1.csv",sep=";",h=T)

# ensure identifiers are factors
growth <- growth %>%
  mutate(
    Geno = as.factor(Geno),
    PlantID = paste0(Range, "_", Block)
  )

# select growth traits
growth_traits <- c(
  "Plant.Height",
  "Canopy.Height",
  "Canopy.Diameter",
  "Stem.Internodes",
  "Stem.Girth",
  "Number.Primaries",
  "Primary.Internodes",
  "Longest.Primary",
  "Number.Leaves",
  "Leaf.Length",
  "Leaf.Width",
  "Leaf.Area",
  "Total.Area"
)


### Genotype-level clustering
geno_growth <- growth %>%
  group_by(Geno) %>%
  summarise(across(all_of(growth_traits), mean, na.rm = TRUE))

geno_traits <- geno_growth %>%
  select(all_of(growth_traits))

scaled_geno <- scale(geno_traits)

set.seed(123)
km_geno <- kmeans(scaled_geno, centers = 2, nstart = 25)

geno_growth$Cluster_geno <- factor(km_geno$cluster)

# cluster size
table(geno_growth$Cluster_geno)



# Elbow plot
p_elbow <- fviz_nbclust(
  scaled_geno,
  kmeans,
  method = "wss",
  k.max = 10,
  nstart = 50
) +
geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
labs(
  title = "Elbow Method (WSS)",
  x = "Number of Clusters",
  y = "WSS"
) +
theme_classic(base_size = 8) +
theme(
  plot.title = element_text(size = 7),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 7)
)


### PCA visualization (genotype clusters)
pca <- prcomp(scaled_geno)

# PCA scores
scores <- as.data.frame(pca$x)
scores$Cluster <- factor(geno_growth$Cluster_geno)
scores$Cluster <- factor(
  scores$Cluster,
  levels = c(1, 2),
  labels = c("Small Biomass Cluster", "Big Biomass Cluster")
)


# PCA loadings
loadings <- as.data.frame(pca$rotation)
loadings$Trait <- rownames(loadings)

# scale arrows for plotting
arrow_scale <- 5
loadings$PC1 <- loadings$PC1 * arrow_scale
loadings$PC2 <- loadings$PC2 * arrow_scale

p_pca <- ggplot(scores, aes(PC1, PC2)) +

# points
geom_point(
  aes(color = Cluster, shape = Cluster),
  size = 1.5
) +

# loadings arrows
geom_segment(
  data = loadings,
  aes(x = 0, y = 0, xend = PC1, yend = PC2),
  arrow = arrow(length = unit(0.12, "cm")),
  color = "black"
) +

# trait labels
geom_text(
  data = loadings,
  aes(x = PC1, y = PC2, label = Trait),
  size = 2,
  hjust = 0.5
) +

# dashed axes
geom_hline(yintercept = 0, linetype = "dashed") +
geom_vline(xintercept = 0, linetype = "dashed") +

# colors similar to manuscript
scale_color_manual(
  values = c(
    "Big Biomass Cluster" = "red",
    "Small Biomass Cluster" = "blue"
  )
) +
scale_shape_manual(
  values = c(
    "Big Biomass Cluster" = 16,   # red circle
    "Small Biomass Cluster" = 17  # blue triangle
  )
) +

labs(
  title = "PCA Biplot: Biomass-Based Clustering",
  x = paste0("Dim1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
  y = paste0("Dim2 (", round(summary(pca)$importance[2,2]*100,1), "%)"),
  color = "Biomass Cluster",
  shape = "Biomass Cluster"
  ) +
  theme_classic(base_size = 7) + 
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7)
  )


final_plot <- p_elbow + p_pca +
  plot_annotation(tag_levels = "A")

tiff("outputs/S1 Figure 2.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

print(final_plot)

dev.off()




### Export cluster table for merging with physiology data
geno_growth %>%
  mutate(Cluster = geno_growth$Cluster_geno) %>%
  group_by(Cluster) %>%
  summarise(across(all_of(growth_traits), mean))

cluster_table <- geno_growth %>%
  select(Geno, Cluster_geno)

write.csv(cluster_table,
          "data_processed/Genotype_clusters.csv",
          row.names = FALSE)