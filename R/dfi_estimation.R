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
library(ggpubr)
library(ggpmisc)

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
key_days <- c(1, 7, 15)


df_sub <- df %>%
  filter(Day %in% c(1, 7, 15))

mod_PI <- lmer(
  PI ~ Batch + factor(Day) +
    (0 + factor(Day) | Geno) +
    (1 | Geno:PlantID),
  data = df_sub
)

re <- ranef(mod_PI)$Geno

df_blup <- as.data.frame(re) %>%
  tibble::rownames_to_column("Geno")

fix <- fixef(mod_PI)

# extract day effects
day_effects <- fix[grep("factor\\(Day\\)", names(fix))]

intercept <- fix["(Intercept)"]

df_blup <- df_blup %>%
  mutate(
    PI_1  = intercept + `factor(Day)1`,
    PI_7  = intercept + day_effects["factor(Day)7"]  + `factor(Day)7`,
    PI_15 = intercept + day_effects["factor(Day)15"] + `factor(Day)15`
  )

df_dfi <- df_blup %>%
  mutate(
    rel_mod = PI_7 / PI_1,
    rel_sev = PI_15 / PI_1,
    DFI = log(rel_mod) + 2 * log(rel_sev)
  )


hist(df_dfi$DFI)
df_dfi[order(df_dfi$DFI, decreasing=T),]

# DFI vector
x <- df_dfi$DFI

# remove NA if needed
x <- na.omit(x)

# run kmeans (example with 3 clusters)
set.seed(123)
k <- 3
km <- kmeans(x, centers = k)

# cluster assignment
df_dfi$cluster <- NA
df_dfi$cluster[!is.na(df_dfi$DFI)] <- km$cluster

plot(x, col = km$cluster, pch = 19,
     main = "K-means clustering (1D)",
     ylab = "DFI")

wss <- sapply(1:10, function(k){
  kmeans(x, centers = k, nstart = 10)$tot.withinss
})

p1 <- wrap_elements(
  full = ~{
    par(cex.lab = 0.7, cex.axis = 0.6)
    plot(1:10, wss, type = "b", pch = 16,
               xlab = "Number of clusters",
               ylab = "Within-cluster sum of squares"
              )
              }
)

table(km$cluster)

# counts
counts <- table(km$cluster)

# create labels with n
labels <- c(
  paste0("Tolerant\n(n=", counts[2], ")"),
  paste0("Intermediate\n(n=", counts[1], ")"),
  paste0("Sensitive\n(n=", counts[3], ")")
)

# factor with correct order
df_dfi$group <- factor(km$cluster,
                       levels = c(2, 1, 3),
                       labels = labels)

# plot
p2 <- ggplot(df_dfi, aes(x = group, y = DFI, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#D55E00")) +
  labs(title = "DFI Distribution by Cluster",
       x = NULL,
       y = "DFI") +
  theme_minimal(base_size = 7) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold",size=7),
        axis.text.x = element_text(face = "bold",size=6),
        axis.title = element_text(size = 6))


tiff("outputs/DFIclustering.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

print(p1+p2)

dev.off()


df_plot <- df_dfi %>%
  left_join(df %>% select(Geno, Cluster_geno), by = "Geno") %>%
  mutate(
    group = factor(cluster,
                   levels = c(2, 1, 3),
                   labels = c("Tolerant", "Intermediate", "Sensitive")),
    
    PhenoCluster = factor(Cluster_geno,
                          levels = c(1,2),
                          labels = c("Small Biomass", "Big Biomass"))
  )

df_plot2 <- df_plot %>%
  distinct(Geno, .keep_all = TRUE)


fig4 <- ggplot(df_plot2, aes(x = PhenoCluster, y = DFI, fill = PhenoCluster)) +
  
  # smooth violin (no harsh tails)
  geom_violin(trim = TRUE, color = NA, alpha = 0.6) +
  
  # subtle boxplot
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7) +
  
  # ONE point per genotype
  geom_jitter(aes(color = group),
              width = 0.12,
              size = 1,
              alpha = 0.9) +
  
  # soft colors like your reference
  scale_fill_manual(values = c(
    "Small Biomass" = "#F4A6A6",   # soft pink
    "Big Biomass"   = "#8FD0D2"    # soft blue
  )) +
  
  scale_color_manual(values = c(
    "Sensitive" = "#E41A1C",
    "Intermediate" = "#FF9900",
    "Tolerant" = "#00C853"
  )) +
  
  labs(
    x = "Phenotypic Cluster",
    y = "Drought Factor Index (DFI)",
    color = "Ranks"
  ) +
  
  theme_classic(base_size = 8) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(face = "bold"),
    panel.border = element_blank()
  ) +
  
  # significance (clean position)
  stat_compare_means(
    comparisons = list(c("Small Biomass", "Big Biomass")),
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01
  )

tiff("outputs/Figure 4.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

fig4

dev.off()


table(df_plot2$group, df_plot2$PhenoCluster)

write.table(df_plot2,"outputs/dfi_table.csv", row.names=F)


df_sub2 <- df %>%
  filter(Day %in% c(1, 7, 15))


df_sub2 <- df %>%
  filter(Day %in% c(1, 7, 15)) %>%
  left_join(
    df_plot2 %>% select(Geno, cluster),
    by = "Geno"
  ) %>% mutate(
    group = factor(cluster,
                   levels = c(2, 1, 3),
                   labels = c("Tolerant", "Intermediate", "Sensitive")),
    
    PhenoCluster = factor(Cluster_geno,
                          levels = c(1,2),
                          labels = c("Small Biomass", "Big Biomass"))
  )

traits <- c("YII","ETR","gs","Fv.Fm","PI","PDLWP")



results <- lapply(traits, function(tr){

  formula <- as.formula(
    paste(tr, "~ Batch + group * factor(Day) + (1 | Geno)")
  )

  mod <- lmer(formula, data = df_sub2)

  emm <- emmeans(mod, ~ group | Day, at = list(Day = c(1,15)))
  contr <- summary(contrast(emm, method = "pairwise"))

  list(
    model = mod,
    emm = as.data.frame(emm),
    contrast = as.data.frame(contr)
  )
})

names(results) <- traits

emm_df <- bind_rows(lapply(names(results), function(tr){
  df <- results[[tr]]$emm
  df$Trait <- tr
  df
}))

emm_df2 <- emm_df %>%
  mutate(group = factor(group,
                        levels = c("Sensitive", "Tolerant")))

emm_df2 <- emm_df2 %>%
  filter(group %in% c("Sensitive", "Tolerant")) %>%
  mutate(
    group = factor(group, levels = c("Sensitive", "Tolerant")),
    Condition = ifelse(Day == 1, "WW", "DS"),
    Condition = factor(Condition, levels = c("WW", "DS"))
  )

contr_df <- bind_rows(lapply(names(results), function(tr){
  df <- results[[tr]]$contrast
  df$Trait <- tr
  df
}))

trait_labels <- c(
  PI = "Photosynthetic Performance Index (PI)",
  YII = "Effective Quantum Yield of PSII (YII)",
  Fv.Fm = "Maximum Quantum Efficiency (Fv/Fm)",
  ETR = "Electron Transport Rate (ETR)",
  gs = "Stomatal Conductance (gs)",
  PDLWP = "Predawn Leaf Water Potential (PDLWP)"
)

theme_fig6 <- theme_classic() +
  theme(
    axis.line = element_line(size = 0.6, color = "black"),
    axis.ticks = element_line(size = 0.4),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 4, face = "bold"),
    strip.background = element_blank(),
    legend.position = "right",
    panel.spacing = unit(1.4, "lines")
  )

y_labels <- c(
  PI = "PI",
  YII = "YII",
  Fv.Fm = "Fv/Fm",
  ETR = expression(ETR~(mu*mol~m^{-2}~s^{-1})),
  gs = expression(gs~(mol~m^{-2}~s^{-1})),
  PDLWP = "PDLWP (MPa)"
)

p <- ggplot(emm_df2, aes(x = Condition, y = emmean, fill = group)) +
  
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    color = "black",
    linewidth = 0.25
  ) +
  
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL), 
    position = position_dodge(width = 0.7),
    width = 0.2,
    size = 0.7,
    linewidth = 0.3
  ) +
  
  facet_wrap(~ Trait, scales = "free_y", ncol = 3,
           labeller = labeller(Trait = trait_labels),
           strip.position = "top") +
  
  scale_fill_manual(values = c(
    "Sensitive" = "#d73027",
    "Tolerant" = "#2e7d32"
  )) +
  
  labs(
    x = "Treatment",
    y = "Trait value",
    fill = NULL
  ) +
  
  theme_fig6




tiff("outputs/Figure 6.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

p

dev.off()

## Correlation - Regression

# Day 1
df_DS <- df %>%
  filter(Day == 1)

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP","Leafincl")

emm_list <- lapply(traits, function(tr){

  formula <- as.formula(paste(tr, "~ Batch + Geno"))

  mod <- lm(formula, data = df_DS)

  emm <- emmeans(mod, ~ Geno)

  df_out <- as.data.frame(emm)
  df_out$Trait <- tr

  return(df_out)
})

emm_all <- bind_rows(emm_list)

traits_wide <- emm_all %>%
  select(Geno, Trait, emmean) %>%
  pivot_wider(names_from = Trait, values_from = emmean)

df_final <- traits_wide %>%
  left_join(df_dfi, by = "Geno")

cor_mat <- cor(df_final[,c("DFI",traits)], use = "complete")

tiff("outputs/correlation.Day1.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

corrplot(
  cor_mat,
  method = "circle",        # circles
  type = "lower",           # lower triangle
  col = colorRampPalette(c("red", "white", "blue"))(200),
  tl.col = "red",           # variable names in red
  tl.srt = 90,
  tl.cex = 0.8,              # vertical labels
  addCoef.col = "black",    # numbers
  number.cex = 0.4,
  diag = FALSE
)

dev.off()

#Day 7
df_DS <- df %>%
  filter(Day == 7)

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP","Leafincl")

emm_list <- lapply(traits, function(tr){

  formula <- as.formula(paste(tr, "~ Batch + Geno"))

  mod <- lm(formula, data = df_DS)

  emm <- emmeans(mod, ~ Geno)

  df_out <- as.data.frame(emm)
  df_out$Trait <- tr

  return(df_out)
})

emm_all <- bind_rows(emm_list)

traits_wide <- emm_all %>%
  select(Geno, Trait, emmean) %>%
  pivot_wider(names_from = Trait, values_from = emmean)

df_final <- traits_wide %>%
  left_join(df_dfi, by = "Geno")

cor_mat <- cor(df_final[,c("DFI",traits)], use = "complete")

tiff("outputs/correlation.Day7.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

corrplot(
  cor_mat,
  method = "circle",        # circles
  type = "lower",           # lower triangle
  col = colorRampPalette(c("red", "white", "blue"))(200),
  tl.col = "red",           # variable names in red
  tl.srt = 90,
  tl.cex = 0.8,              # vertical labels
  addCoef.col = "black",    # numbers
  number.cex = 0.4,
  diag = FALSE
)

dev.off()





#Day 15
df_DS <- df %>%
  filter(Day == 15)

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP","Leafincl")

emm_list <- lapply(traits, function(tr){

  formula <- as.formula(paste(tr, "~ Batch + Geno"))

  mod <- lm(formula, data = df_DS)

  emm <- emmeans(mod, ~ Geno)

  df_out <- as.data.frame(emm)
  df_out$Trait <- tr

  return(df_out)
})

emm_all <- bind_rows(emm_list)

traits_wide <- emm_all %>%
  select(Geno, Trait, emmean) %>%
  pivot_wider(names_from = Trait, values_from = emmean)

df_final <- traits_wide %>%
  left_join(df_dfi, by = "Geno")

cor_mat <- cor(df_final[,c("DFI",traits)], use = "complete")

tiff("outputs/correlation.Day15.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

corrplot(
  cor_mat,
  method = "circle",        # circles
  type = "lower",           # lower triangle
  col = colorRampPalette(c("red", "white", "blue"))(200),
  tl.col = "red",           # variable names in red
  tl.srt = 90,
  tl.cex = 0.8,              # vertical labels
  addCoef.col = "black",    # numbers
  number.cex = 0.4,
  diag = FALSE
)

dev.off()




# Regression
#Day 1
df_DS <- df %>%
  filter(Day == 1)

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP")

emm_list <- lapply(traits, function(tr){

  formula <- as.formula(paste(tr, "~ Batch + Geno"))

  mod <- lm(formula, data = df_DS)

  emm <- emmeans(mod, ~ Geno)

  df_out <- as.data.frame(emm)
  df_out$Trait <- tr

  return(df_out)
})

emm_all <- bind_rows(emm_list)

traits_wide <- emm_all %>%
  select(Geno, Trait, emmean) %>%
  pivot_wider(names_from = Trait, values_from = emmean)

df_final <- traits_wide %>%
  left_join(df_dfi, by = "Geno")

df_final <- df_final %>%
  left_join(
    df_plot2 %>%
      select(Geno, PhenoCluster),
    by = "Geno"
  ) %>%
  rename(Cluster = PhenoCluster)

plot_list <- list()

for (tr in traits) {

  p <- ggplot(df_final, aes_string(x = tr, y = "DFI")) +
    
    geom_point(aes(color = Cluster), size = 2, alpha = 0.8) +
    
    geom_smooth(method = "lm", se = TRUE,
                color = "black", fill = "grey70", alpha = 0.3, size = 1) +
    
    stat_poly_eq(
      aes_string(
        x = tr, y = "DFI",
        label = "after_stat(paste(..eq.label.., ..adj.rr.label.., sep = '~~~'))"
      ),
      formula = y ~ x,
      parse = TRUE,
      size = 2,
      geom = "text_npc",
      npcx = "left",
      npcy = 0.95
    ) +
    
    scale_color_manual(values = c(
      "Small Biomass" = "#2c5aa0",
      "Big Biomass" = "#c43c35"
    )) +
    
    labs(x = tr, y = "DFI") +
    
    theme_classic() +
    theme(
      legend.position = "none",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      axis.title = element_text(face = "bold", size = 6),
      axis.text = element_text(color = "black", size = 6),
    )
    

  plot_list[[tr]] <- p
}

plot_list

combined_plot <- wrap_plots(plot_list, ncol = 3)

tiff("outputs/regression.Day1.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

print(combined_plot)

dev.off()

#Day 7
df_DS <- df %>%
  filter(Day == 7)

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP")

emm_list <- lapply(traits, function(tr){

  formula <- as.formula(paste(tr, "~ Batch + Geno"))

  mod <- lm(formula, data = df_DS)

  emm <- emmeans(mod, ~ Geno)

  df_out <- as.data.frame(emm)
  df_out$Trait <- tr

  return(df_out)
})

emm_all <- bind_rows(emm_list)

traits_wide <- emm_all %>%
  select(Geno, Trait, emmean) %>%
  pivot_wider(names_from = Trait, values_from = emmean)

df_final <- traits_wide %>%
  left_join(df_dfi, by = "Geno")

df_final <- df_final %>%
  left_join(
    df_plot2 %>%
      select(Geno, PhenoCluster),
    by = "Geno"
  ) %>%
  rename(Cluster = PhenoCluster)

plot_list <- list()

for (tr in traits) {

  p <- ggplot(df_final, aes_string(x = tr, y = "DFI")) +
    
    geom_point(aes(color = Cluster), size = 2, alpha = 0.8) +
    
    geom_smooth(method = "lm", se = TRUE,
                color = "black", fill = "grey70", alpha = 0.3, size = 1) +
    
    stat_poly_eq(
      aes_string(
        x = tr, y = "DFI",
        label = "after_stat(paste(..eq.label.., ..adj.rr.label.., sep = '~~~'))"
      ),
      formula = y ~ x,
      parse = TRUE,
      size = 2,
      geom = "text_npc",
      npcx = "left",
      npcy = 0.95
    ) +
    
    scale_color_manual(values = c(
      "Small Biomass" = "#2c5aa0",
      "Big Biomass" = "#c43c35"
    )) +
    
    labs(x = tr, y = "DFI") +
    
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black", size = 6)
    )

  plot_list[[tr]] <- p
}

plot_list

combined_plot <- wrap_plots(plot_list, ncol = 3)

tiff("outputs/regression.Day7.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

print(combined_plot)

dev.off()



#Day 15
df_DS <- df %>%
  filter(Day == 15)

traits <- c("YII","ETR","gs","E","Fo.Fm","Fv.Fm","PI","Fv.Fo","Chlo","PDLWP")

emm_list <- lapply(traits, function(tr){

  formula <- as.formula(paste(tr, "~ Batch + Geno"))

  mod <- lm(formula, data = df_DS)

  emm <- emmeans(mod, ~ Geno)

  df_out <- as.data.frame(emm)
  df_out$Trait <- tr

  return(df_out)
})

emm_all <- bind_rows(emm_list)

traits_wide <- emm_all %>%
  select(Geno, Trait, emmean) %>%
  pivot_wider(names_from = Trait, values_from = emmean)

df_final <- traits_wide %>%
  left_join(df_dfi, by = "Geno")

df_final <- df_final %>%
  left_join(
    df_plot2 %>%
      select(Geno, PhenoCluster),
    by = "Geno"
  ) %>%
  rename(Cluster = PhenoCluster)

plot_list <- list()

for (tr in traits) {

  p <- ggplot(df_final, aes_string(x = tr, y = "DFI")) +
    
    geom_point(aes(color = Cluster), size = 2, alpha = 0.8) +
    
    geom_smooth(method = "lm", se = TRUE,
                color = "black", fill = "grey70", alpha = 0.3, size = 1) +
    
    stat_poly_eq(
      aes_string(
        x = tr, y = "DFI",
        label = "after_stat(paste(..eq.label.., ..adj.rr.label.., sep = '~~~'))"
      ),
      formula = y ~ x,
      parse = TRUE,
      size = 2,
      geom = "text_npc",
      npcx = "left",
      npcy = 0.95
    ) +
    
    scale_color_manual(values = c(
      "Small Biomass" = "#2c5aa0",
      "Big Biomass" = "#c43c35"
    )) +
    
    labs(x = tr, y = "DFI") +
    
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black", size = 6)
    )

  plot_list[[tr]] <- p
}

plot_list

combined_plot <- wrap_plots(plot_list, ncol = 3)

tiff("outputs/regression.Day15.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

print(combined_plot)

dev.off()
