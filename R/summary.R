library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(splines)
library(emmeans)
library(ggplot2)
library(here)
library(ordinal)
library(patchwork)

df <- read.csv(here("data_raw", "physio_data2.csv"), sep=";", stringsAsFactors = FALSE)

df$Day   <- as.numeric(df$Day)
df$Batch <- factor(df$Batch)
df$Geno  <- factor(df$Geno)

# Unique plant identifier
df$PlantID <- factor(paste0(df$Row, "_", df$New_Col))

# Load genotype cluster information
clusters <- read.csv(here("data_processed","Genotype_clusters.csv"))

# Check structure
str(clusters)
head(clusters)

df <- df %>%
  left_join(clusters, by = "Geno")

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP","Leafincl")
traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP","Leafincl")

# Prepare data
df <- df %>%
  mutate(
    Cluster = ifelse(Cluster_geno == 1, "Small Biomass", "Big Biomass")
  )

df_sub <- df %>% filter(Day %in% c(1,15))


# Long format
df_long <- df_sub %>%
  mutate(Leafincl = as.numeric(Leafincl)) %>%
  pivot_longer(
    cols = all_of(traits),
    names_to = "Trait",
    values_to = "value"
  ) %>%
  mutate(
    Regime = ifelse(Day == 1, "WW", "DS")
  )

df_long$Regime <- factor(df_long$Regime, levels = c("WW", "DS"))
# -----------------------------
# 3. RAW STATS
# -----------------------------
raw_stats <- df_long %>%
  group_by(Trait, Day, Cluster) %>%
  summarise(
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    CV = 100 * SD / mean(value, na.rm = TRUE),
    .groups = "drop"
  )


# Model
emm_list <- list()
pct_list <- list()

for(tr in traits){

  df_tr <- df_sub %>%
  	mutate(Leafincl = as.numeric(Leafincl))

  mod <- lmer(
    as.formula(paste(tr, "~ Batch + factor(Day)*Cluster + (1|Geno)")),
    data = df_tr
  )

  emm <- emmeans(mod, ~ Cluster | Day)
  emm_df <- as.data.frame(emm)
  emm_df$Trait <- tr

  # add Trt
  emm_df <- emm_df %>%
    mutate(Trt = ifelse(Day == 1, "WW", "DS"))

  emm_list[[tr]] <- emm_df

  # percent change
  pct <- emm_df %>%
    select(Cluster, Day, emmean) %>%
    pivot_wider(names_from = Day, values_from = emmean) %>%
    mutate(
      pct = 100 * (`15` - `1`) / abs(`1`),
      Trait = tr
    )

  pct_list[[tr]] <- pct
}

emm_all <- bind_rows(emm_list)
pct_all <- bind_rows(pct_list)

# Merge supp table

table_S1 <- emm_all %>%
  left_join(raw_stats, by = c("Trait","Day","Cluster")) %>%
  mutate(
    Mean_SE = paste0(round(emmean, 2), " ± ", round(SE, 2))
  ) %>%
  select(
    Trait, Trt, Cluster,
    Min, Max, Median,
    Mean_SE,
    SD, CV
  )

write.csv(table_S1, "outputs/summaryTable.csv", row.names = FALSE)


# Labels for plot
label_df <- pct_all %>%
  mutate(
    label = paste0(ifelse(pct > 0, "+",""), round(pct,1), "%"),
    Regime = "DS"
  )



# y position per trait
y_pos <- df_long %>%
  group_by(Trait) %>%
  summarise(y = max(value, na.rm = TRUE) * 0.85)

label_df <- label_df %>%
  left_join(y_pos, by = "Trait")

label_df$Regime <- factor(label_df$Regime, levels = c("WW", "DS"))


# Plot

plot_list <- list()

for(tr in traits){

  df_plot <- df_long %>% filter(Trait == tr)

  lab <- label_df %>% filter(Trait == tr)

  p <- ggplot(df_plot, aes(Regime, value, fill = Regime)) +

    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = NA) +

    geom_boxplot(
      width = 0.25,
      fill = "white",
      color = "black",
      linewidth = 0.3,
      outlier.shape = NA
    ) +

    facet_wrap(~Cluster, nrow = 1) +

    scale_fill_manual(values = c(WW = "#00BFC4", DS = "#F564E3")) +

    labs(title = tr, x = NULL, y = NULL) +

    theme_classic(base_size = 6) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 6),
      plot.title = element_text(size = 7, face = "bold")
    ) +

    geom_label(
      data = lab,
      aes(x = 1.5, y = y, label = label),
      fill = "white",
      color = "darkgreen",
      size = 1.5,
      inherit.aes = FALSE
    )

  plot_list[[tr]] <- p
}

final_plot <- wrap_plots(plot_list, ncol = 3)

tiff("outputs/Figure 3.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")
final_plot
dev.off()
