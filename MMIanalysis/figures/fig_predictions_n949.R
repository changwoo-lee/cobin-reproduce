rm(list = ls())
library(ggplot2)
library(maps)
library(usmap)
library(sf)
library(dplyr)
library(spNNGP)
# check version
packageVersion("spNNGP") #### spNNGP version is important: this demo is based on 1.0.1


# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))





df = read.csv("../mmi_lakecat.csv")
df = df[df$MMI_BENT!=0,]
dim(df) # nrow = 949

drawdf4_beta = readRDS("../n949/res/predict_beta.rds")

p11_beta = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4_beta, aes(x = x, y = y, color = munew_hat), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "pred_MMI (beta)")
p11_beta

# ggsave(p11_beta, width = 7, height = 4, filename = "fig_beta_muhat.pdf")
# ggsave(p11_beta, width = 7, height = 4, filename = "fig_beta_muhat.png")
# ggsave(p11_beta, width = 7, height = 4, filename = "fig_beta_muhat.jpg")



drawdf4_cobin = readRDS("../n949/res/predict_cobin.rds")

p11_cobin = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4_cobin, aes(x = x, y = y, color = munew_hat), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "pred_MMI (cobin)")
p11_cobin

# ggsave(p11_cobin, width = 7, height = 4, filename = "fig_cobin_muhat.pdf")
# ggsave(p11_cobin, width = 7, height = 4, filename = "fig_cobin_muhat.png")
# ggsave(p11_cobin, width = 7, height = 4, filename = "fig_cobin_muhat.jpg")

drawdf4_betarec = readRDS("../n949/res/predict_betarec.rds")
p11_betarec = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4_betarec, aes(x = x, y = y, color = munew_hat), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "pred_MMI (betarec)")
p11_betarec

# ggsave(p11_betarec, width = 7, height = 4, filename = "fig_betarec_muhat.pdf")
# ggsave(p11_betarec, width = 7, height = 4, filename = "fig_betarec_muhat.png")
# ggsave(p11_betarec, width = 7, height = 4, filename = "fig_betarec_muhat.jpg")

drawdf4_micobin = readRDS("../n949/res/predict_micobin.rds")
p11_micobin = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4_micobin, aes(x = x, y = y, color = munew_hat), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "predicted MMI")
p11_micobin

p22_micobin = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4_micobin, aes(x = x, y = y, color = munew_sd), size = 0.5, alpha = 0.5) +
  scale_colour_distiller(palette = "Spectral")+
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  ) + labs(color = "p.p.d. stdev")
p22_micobin

library(ggpubr)
p11p22 = ggarrange(p11_micobin, p22_micobin)
p11p22
# ggsave(p11p22, width = 12, height = 4, filename = "fig_micobin_muhatsd.pdf")
# ggsave(p11p22, width = 12, height = 4, filename = "fig_micobin_muhatsd.png", bg = "white")
# ggsave(p11p22, width = 12, height = 4, filename = "fig_micobin_muhatsd.jpg", bg = "white")




####################
# contrast plots
df_contrast = data.frame(
  x = drawdf4_cobin$x,
  y = drawdf4_cobin$y,
  wnew_hat_cobin_micobin = drawdf4_cobin$wnew_hat - drawdf4_micobin$wnew_hat,
  munew_hat_cobin_micobin = drawdf4_cobin$munew_hat - drawdf4_micobin$munew_hat,
  munew_sd_beta_micobin = drawdf4_beta$munew_sd - drawdf4_micobin$munew_sd,
  wnew_hat_beta_micobin = drawdf4_beta$wnew_hat - drawdf4_micobin$wnew_hat,
  munew_hat_beta_micobin = drawdf4_beta$munew_hat - drawdf4_micobin$munew_hat,
  wnew_hat_betarec_micobin = drawdf4_betarec$wnew_hat - drawdf4_micobin$wnew_hat,
  munew_hat_betarec_micobin = drawdf4_betarec$munew_hat - drawdf4_micobin$munew_hat
  #  wnew_hat_betarec_micobin = drawdf4_t$wnew_hat - drawdf4_micobin$wnew_hat,
  #  munew_hat_betarec_micobin = drawdf4_t$munew_hat - drawdf4_micobin$munew_hat
)

idx1 = which(df$comid==9201925)
idx2 = which(df$comid==22845861)

transformed_coords_outlier <- sf::st_coordinates(usmap::usmap_transform(df[c(idx1, idx2),c("lon","lat")]))


df_outlier <- data.frame(
  x = transformed_coords_outlier[,1],
  y = transformed_coords_outlier[,2],
  shape_label = c("Square", "Triangle")  # for legend
)
shape_values <- c("Square" = 15, "Triangle" = 17)

p111 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = df_contrast, aes(x = x, y = y, color = munew_hat_beta_micobin), size = 0.5, alpha = 0.5) +
  geom_point(data = df_outlier, aes(x = x, y = y, shape = shape_label),
             color = "black", size = 2,
             show.legend = FALSE) +  # <- new layer
  scale_shape_manual(values = shape_values) +  # <- link shape names to shape codes
  scale_color_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.15,0.04)) + 
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "beta-micobin")
p111
p222 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = df_contrast, aes(x = x, y = y, color = munew_hat_cobin_micobin), size = 0.5, alpha = 0.5) +
  geom_point(data = df_outlier, aes(x = x, y = y, shape = shape_label),
             color = "black", size = 2,
             show.legend = FALSE) +  # <- new layer
  scale_shape_manual(values = shape_values) +  # <- link shape names to shape codes
  scale_color_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.15,0.04)) + 
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+  labs(color = "cobin-micobin")
p222
p333 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = df_contrast, aes(x = x, y = y, color = munew_hat_betarec_micobin), size = 0.5, alpha = 0.5) +
  geom_point(data = df_outlier, aes(x = x, y = y, shape = shape_label),
             color = "black", size = 2,
             show.legend = FALSE) +  # <- new layer
  scale_shape_manual(values = shape_values) +  # <- link shape names to shape codes
  scale_color_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.15,0.04)) + 
  theme(
    legend.position = c(0.9, 0.05),
    legend.justification = c("right", "center"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+  labs(color = "betarec-micobin")
p333

# ggsave(p111, width = 7, height = 4, filename = "fig_contrast_beta_micobin.pdf")
# ggsave(p222, width = 7, height = 4, filename = "fig_contrast_cobin_micobin.pdf")
# ggsave(p333, width = 7, height = 4, filename = "fig_contrast_betarec_micobin.pdf")

library(patchwork)
gtogether <- (p11_beta  | p111)  / 
  (p11_cobin | p222)  / 
  (p11_betarec | p333)
gtogether
# ggsave(gtogether, width = 13, height = 12, filename = "fig_all.pdf")
# ggsave(gtogether, width = 13, height = 12, filename = "fig_all.png")
# ggsave(gtogether, width = 13, height = 12, filename = "fig_all.jpg")



