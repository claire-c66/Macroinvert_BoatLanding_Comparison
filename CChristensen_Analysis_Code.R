# =========================================================
# Natural Boat Landing Impacts on Benthic Invertebrates
# Author: Claire Christensen
# Output: Statistical analysis comparing boat landing 
# vs undisturbed shoreline sites and figure generation
# =========================================================

# ---------------------------
# Libraries
# ---------------------------
library(readxl)
library(tidyverse)
library(vegan)
library(lme4)
library(grid)
library(gridExtra)

# ---------------------------
# Data Import
# ---------------------------
setwd("#insert your working directory#")

invertebrate_counts <- read_excel("invertebrate_counts.xlsx") %>%
  column_to_rownames("site_ID")

site_info <- read_excel("site_info.xlsx") %>%
  select(site_ID, zone_type)

# Ensure consistent factor order
site_info$zone_type <- factor(site_info$zone_type,
                              levels = unique(site_info$zone_type))

# ---------------------------
# Species Richness
# ---------------------------
sppr <- specnumber(invertebrate_counts)

sppr_df <- enframe(sppr, name = "site_ID", value = "richness") %>%
  left_join(site_info, by = "site_ID")

shapiro.test(sppr_df$richness)
bartlett.test(richness ~ zone_type, data = sppr_df)

model_sppr <- lm(richness ~ zone_type, data = sppr_df)
plot(model_sppr, which = 1)
plot(model_sppr, which = 2)

wil_sppr <- wilcox.test(richness ~ zone_type,
                        data = sppr_df,
                        paired = TRUE,
                        alternative = "less")
wil_sppr

sppr_plot <- ggplot(sppr_df, aes(x = zone_type, y = richness)) +
  geom_violin(fill = "grey85", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x = NULL, y = "Species Richness")

# ---------------------------
# Shannon Diversity
# ---------------------------
shannon <- diversity(invertebrate_counts, index = "shannon")

shannon_df <- site_info %>%
  mutate(shannon = shannon)

shapiro.test(shannon_df$shannon)
bartlett.test(shannon ~ zone_type, data = shannon_df)

t_shannon <- t.test(shannon ~ zone_type,
                    paired = TRUE,
                    data = shannon_df)
t_shannon

sppdiv_plot <- ggplot(shannon_df, aes(x = zone_type, y = shannon)) +
  geom_violin(fill = "grey85", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x = NULL, y = "Shannon Diversity Index")

# ---------------------------
# Simpson Diversity
# ---------------------------
simpson <- diversity(invertebrate_counts, index = "simpson")

simpson_df <- site_info %>%
  mutate(simpson = simpson)

shapiro.test(simpson_df$simpson)
bartlett.test(simpson ~ zone_type, data = simpson_df)

model_simpson <- lm(simpson ~ zone_type, data = simpson_df)
plot(model_simpson, which = 1)
plot(model_simpson, which = 2)

wil_simpson <- wilcox.test(simpson ~ zone_type,
                           data = simpson_df,
                           paired = TRUE,
                           alternative = "less")
wil_simpson

simp_plot <- ggplot(simpson_df, aes(x = zone_type, y = simpson)) +
  geom_violin(fill = "grey85", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x = NULL, y = "Simpson Diversity Index")

# ---------------------------
# HBI
# ---------------------------
HBI_data <- read_excel("HBI.xlsx") %>%
  mutate(HBI = as.numeric(HBI))

shapiro.test(HBI_data$HBI)
bartlett.test(HBI ~ zone_type, data = HBI_data)

wil_HBI <- wilcox.test(HBI ~ zone_type,
                       data = HBI_data,
                       paired = TRUE,
                       alternative = "less")
wil_HBI

HBI_plot <- ggplot(HBI_data, aes(x = zone_type, y = HBI)) +
  geom_violin(fill = "grey85", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x = "Zone Type", y = "HBI Value")

# ---------------------------
# EPT Index
# ---------------------------
EPT_data <- read_excel("EPT.xlsx")

shapiro.test(EPT_data$EPT_index_value)
bartlett.test(EPT_index_value ~ zone_type, data = EPT_data)

wil_EPT <- wilcox.test(EPT_index_value ~ zone_type,
                       data = EPT_data,
                       paired = TRUE,
                       alternative = "less")
wil_EPT

EPT_plot <- ggplot(EPT_data, aes(x = zone_type, y = EPT_index_value)) +
  geom_violin(fill = "grey85", color = "black") +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x = "Zone Type", y = "EPT Index Value")

# ---------------------------
# Community Composition (PERMANOVA)
# ---------------------------
community_data <- read_excel("PERMANOVA.xlsx")

community_matrix <- community_data %>%
  select(-zone_type, -total) %>%
  column_to_rownames("site_ID")

log_comm <- log(community_matrix + 1)

diss_matrix <- vegdist(log_comm, method = "bray")

env_data <- community_data %>%
  select(site_ID, zone_type) %>%
  mutate(zone_type = factor(zone_type)) %>%
  filter(site_ID %in% rownames(log_comm))

permanova <- adonis2(diss_matrix ~ zone_type, data = env_data)

print(permanova)

# ---------------------------
# Multi-panel Figure
# ---------------------------

final_plot <- grid.arrange(
  arrangeGrob(sppr_plot, top = textGrob("a)", x = 0,
                                        gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(sppdiv_plot, top = textGrob("b)", x = 0,
                                          gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(simp_plot, top = textGrob("c)", x = 0,
                                        gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(HBI_plot, top = textGrob("d)", x = 0,
                                       gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(EPT_plot, top = textGrob("e)", x = 0,
                                       gp = gpar(fontsize = 16, fontface = "bold"))),
  ncol = 2
)

final_plot