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
key_days <- c(1, 7, 15)

fit_spline_model <- function(trait, df) {
  formula <- as.formula(
    paste0(trait, " ~ ns(Day, df = 2)*Cluster_geno + Batch + ",
           "(1 | Geno) + (1|Geno:PlantID)")
  )
  lmer(formula, data = df)
}

extract_emmeans <- function(model, trait) {
  emm <- emmeans(model, ~ Cluster_geno | Day, at = list(Day = key_days))
  as.data.frame(emm) %>% mutate(Trait = trait)
}

compute_percent_change <- function(df_trait) {
  
  df_trait %>%
    select(Trait, Cluster_geno, Day, emmean) %>%   # keep only needed columns
    pivot_wider(names_from = Day, values_from = emmean) %>%
    mutate(
      pct_mid = (`7` - `1`) / abs(`1`) * 100,
      pct_end = (`15` - `1`) / abs(`1`) * 100
    )
}


all_emm <- list()
all_pct <- list()
models <- list()

for(tr in traits){

  cat("Running:", tr, "\n")

  mod <- fit_spline_model(tr, df)
  models[[tr]] <- mod

  # model checks


  cat("  Checking model...\n")

  # Residual summary
  res <- residuals(mod)
  cat("   Residual SD:", round(sd(res, na.rm = TRUE), 3), "\n")

  # Singular fit check
  if(lme4::isSingular(mod)){
    cat("   ⚠ Singular fit detected\n")
  }


  # visual checks


try({

    set.seed(1)

    geno_sample <- df %>%
      distinct(Geno, Cluster_geno) %>%
      group_split(Cluster_geno) %>%
      lapply(function(x) x[sample(nrow(x), min(3, nrow(x))), ]) %>%
      bind_rows() %>%
      pull(Geno)

    df_sub <- df %>% filter(Geno %in% geno_sample)

    df_sub$pred <- predict(mod, newdata = df_sub, re.form = NULL)

    df_plot <- df_sub %>%
      group_by(Geno, Cluster_geno, Day) %>%
      summarise(
        pred = mean(pred, na.rm = TRUE),
        obs  = mean(.data[[tr]], na.rm = TRUE),
        .groups = "drop"
      )

    p_traj <- ggplot(df_plot, aes(Day, color = factor(Geno))) +
      geom_point(aes(y = obs), alpha = 0.6) +
      geom_line(aes(y = pred), linewidth = 1.2) +
      facet_wrap(~Cluster_geno) +
      theme_bw() +
      labs(title = paste("trajectories -", tr), y = tr)

    ggsave(
      filename = paste0("outputs/", tr, "_trajectories.png"),
      plot = p_traj,
      width = 8,
      height = 5,
      dpi = 300
    )

    df_res <- data.frame(
      fitted = fitted(mod),
      resid = residuals(mod)
    )

    p_res <- ggplot(df_res, aes(fitted, resid)) +
      geom_point(alpha = 0.3) +
      geom_smooth(se = FALSE, color = "red") +
      theme_bw() +
      labs(title = paste("residuals -", tr))

    ggsave(
      filename = paste0("outputs/", tr, "_residuals.png"),
      plot = p_res,
      width = 6,
      height = 4,
      dpi = 300
    )

  }, silent = TRUE)

  # main analysis


  emm <- extract_emmeans(mod, tr)
  all_emm[[tr]] <- emm

  pct <- compute_percent_change(emm)
  pct$Trait <- tr
  all_pct[[tr]] <- pct
}

df$Cluster_geno <- factor(df$Cluster_geno)


plot_emm_points <- function(model, df, trait){

  df_points <- as.data.frame(
    emmeans(
      model,
      ~ Cluster_geno | Day,
      at = list(Day = sort(unique(df$Day)))
    )
  )
  df_points$Cluster_geno <- factor(df_points$Cluster_geno)
  
  ggplot(df_points, aes(Day, emmean, color = Cluster_geno)) +
    geom_point(size = 1) +
    geom_line(linewidth = 0.4) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.25) +
    theme_classic(base_size=6) + 
    theme(
      legend.position = "none",
      plot.title = element_text(size = 7), 
      axis.title = element_text(size = 6), 
      axis.text = element_text(size = 5),     
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(
      title = trait,
      x = "Day",
      y = NULL,
      color = "Biomass cluster"
    ) +
    scale_color_manual(
      values = c("steelblue", "firebrick"),
    labels = c("Small", "Big")
)
}

# --- create all plots ---
plots <- lapply(names(models), function(tr){
  plot_emm_points(models[[tr]], df, tr)
})

# --- combine into multi-panel figure ---
final_plot <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )

tiff("outputs/S1 Fig3.tiff",
     width = 1800, height = 1200, units = "px",
     res = 300, compression = "lzw")

print(final_plot)

dev.off()

# Results table
for (tr in traits[traits!="Leafincl"]) print(summary(models[[tr]]))

anova_list <- lapply(models[traits!="Leafincl"], anova)

anova_table <- bind_rows(
  lapply(names(anova_list), function(tr){
    df <- as.data.frame(anova_list[[tr]])
    df$Effect <- rownames(df)
    df$Trait <- tr
    df
  })
)

anova_clean <- anova_table %>%
  filter(Effect %in% c("ns(Day, df = 2)", 
                      "Cluster_geno", 
                      "ns(Day, df = 2):Cluster_geno"))





write.csv(anova_table, "outputs/anova_results.csv", row.names = FALSE)
write.csv(anova_clean, "outputs/anova_main_effects.csv", row.names = FALSE)



# Leaf inclination special case
library(ordinal)

df$Leafincl <- factor(df$Leafincl, ordered = TRUE)
df$Geno <- factor(df$Geno)

# full model
mod_full <- clmm(
  Leafincl ~ ns(Day, df = 2)*Cluster_geno + Batch +
    (1 | Geno) + (1 | Geno:PlantID),
  data = df
)

# no interaction
mod_no_inter <- clmm(
  Leafincl ~ ns(Day, df = 2) + Cluster_geno + Batch +
    (1 | Geno) + (1 | Geno:PlantID),
  data = df
)

# no cluster
mod_no_cluster <- clmm(
  Leafincl ~ ns(Day, df = 2) + Batch +
    (1 | Geno) + (1 | Geno:PlantID),
  data = df
)

# no day
mod_no_day <- clmm(
  Leafincl ~ Cluster_geno + Batch +
    (1 | Geno) + (1 | Geno:PlantID),
  data = df
)

mod_no_batch <- clmm(
  Leafincl ~ ns(Day, df = 2)*Cluster_geno +
    (1 | Geno) + (1 | Geno:PlantID),
  data = df
)

anova(mod_full, mod_no_inter)   # interaction
anova(mod_no_inter, mod_no_cluster)  # cluster
anova(mod_no_inter, mod_no_day)      # day
anova(mod_full, mod_no_batch)      # batch
