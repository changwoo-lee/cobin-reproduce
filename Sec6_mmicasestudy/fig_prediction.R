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


df = read.csv("mmi_lakecat.csv")
df = df[df$MMI_BENT!=0,]
dim(df) # nrow = 949

# step 2: transform covariates
df = df %>% mutate(agkffact_log = log2(1+agkffact),
                   bfi_log = log2(1+bfi),
                   cbnf_log = log2(1+cbnf),
                   conif_log = log2(1+conif),
                   crophay_log = log2(1+crophay),
                   fert_log = log2(1+fert),
                   manure_log = log2(1+manure),
                   pestic1997_log = log2(1+pestic1997),
                   urbmdhi_log = log2(1+urbmdhi)
)

# step 3: standardize covariates
df = df %>% mutate(agkffact_logstd = scale(agkffact_log),
                   bfi_logstd = scale(bfi_log),
                   cbnf_logstd = scale(cbnf_log),
                   conif_logstd = scale(conif_log),
                   crophay_logstd = scale(crophay_log),
                   fert_logstd = scale(fert_log),
                   manure_logstd = scale(manure_log),
                   pestic1997_logstd = scale(pestic1997_log),
                   urbmdhi_logstd = scale(urbmdhi_log)
)

# find mean and sd of each covariate
x_center = df %>% select(agkffact_log,
                         bfi_log,
                         cbnf_log,
                         conif_log,
                         crophay_log,
                         fert_log,
                         manure_log,
                         pestic1997_log,
                         urbmdhi_log) %>% summarise_all(mean)
x_sd = df %>% select(agkffact_log,
                     bfi_log,
                     cbnf_log,
                     conif_log,
                     crophay_log,
                     fert_log,
                     manure_log,
                     pestic1997_log,
                     urbmdhi_log) %>% summarise_all(sd)
# load previously fitted ata

fit_micobin1 = readRDS("results_main_n949/fit_micobin1_n949.rds")

# load prediction data

df4 = read.csv("lakecat-over40000m2.csv")

################################

scaleback <- function(x, centers, scales){
  as.numeric(c(x[1] - sum((x[-1] * centers) / scales),
               x[-1] / scales))
}

x_center = as.numeric(x_center)
x_sd = as.numeric(x_sd)


# plug-in posterior samples
fit_micobin1$spNNGPfit$p.beta.samples = fit_micobin1$post_save[,1:ncol(fit_micobin1$X)]
# order important
fit_micobin1$spNNGPfit$p.theta.samples = fit_micobin1$post_save[,c("sigma.sq", "phi")]
fit_micobin1$spNNGPfit$p.w.samples = t(as.matrix(fit_micobin1$post_u_save))


df4 = df4 %>% mutate(agkffact_logstd = scale(log2(1+agkffact), center = x_center[1], scale = x_sd[1]),
                     bfi_logstd = scale(log2(1+bfi), center = x_center[2], scale = x_sd[2]),
                     cbnf_logstd = scale(log2(1+cbnf), center = x_center[3], scale = x_sd[3]),
                     conif_logstd = scale(log2(1+conif), center = x_center[4], scale = x_sd[4]),
                     crophay_logstd = scale(log2(1+crophay), center = x_center[5], scale = x_sd[5]),
                     fert_logstd = scale(log2(1+fert), center = x_center[6], scale = x_sd[6]),
                     manure_logstd = scale(log2(1+manure), center = x_center[7], scale = x_sd[7]),
                     pestic1997_logstd = scale(log2(1+pestic1997), center = x_center[8], scale = x_sd[8]),
                     urbmdhi_logstd = scale(log2(1+urbmdhi), center = x_center[9], scale = x_sd[9])
)
# select
Xnew = data.frame(
  df4$agkffact_logstd,
  df4$bfi_logstd,
  df4$cbnf_logstd,
  df4$conif_logstd,
  df4$crophay_logstd,
  df4$fert_logstd,
  df4$manure_logstd,
  df4$pestic1997_logstd,
  df4$urbmdhi_logstd
)
Xnew = as.matrix(cbind(1, Xnew))
str(Xnew)

# predict using spNNGP; takes time
pred_micobin1 <- predict(fit_micobin1$spNNGPfit, 
                         X.0 = Xnew, 
                         coords.0 = cbind(df4$easting,df4$northing),
                         sub.sample = list(start = 1, end = 5000, thin = 1))


###################################################

pred_micobin1 = readRDS("../real2/pred_micobin1.rds")

transformed_coords4 <- sf::st_coordinates(usmap::usmap_transform(df4[,c("lon","lat")]))

# replace NaN to NA: likely due to numerical issue with too close coordinates
pred_micobin1$p.w.0[is.nan(pred_micobin1$p.w.0)] = NA
pred_micobin1$p.y.0[is.nan(pred_micobin1$p.y.0)] = NA

wnew_hat = rowMeans(pred_micobin1$p.w.0, na.rm = TRUE)
wnew_sd = apply(pred_micobin1$p.w.0, 1, sd, na.rm = TRUE)

munew_hat = rowMeans(cobin::bftprime(pred_micobin1$p.y.0), na.rm = TRUE)
munew_sd = apply(cobin::bftprime(pred_micobin1$p.y.0), 1, sd, na.rm = TRUE)



drawdf4 = data.frame(
  x = transformed_coords4[,1],
  y = transformed_coords4[,2],
  urbmdhi = log2(1+df4$urbmdhi),
  wnew_hat = wnew_hat,
  wnew_sd = wnew_sd,
  munew_hat = munew_hat,
  munew_sd = munew_sd
)


p11 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4, aes(x = x, y = y, color = munew_hat), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
  #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
  # guides(
  #   color = guide_colorbar(order = 1),
  #   shape = guide_legend(order = 2)
  # )+
  theme(
    # Position the legend using relative coordinates (x, y)
    legend.position = c(0.9, 0.05),
    # Optionally adjust the justification so that the legend is anchored correctly
    legend.justification = c("right", "center"),
    # Reduce the left margin of the legend box to nudge it further left
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "predicted MMI")
p11


p22 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4, aes(x = x, y = y, color = munew_sd), size = 0.5, alpha = 0.5) +
  #scale_color_viridis_c(option = "plasma", direction = +1, limit = c(0,0.1)) +
  scale_colour_distiller(palette = "Spectral")+
  #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
  # guides(
  #   color = guide_colorbar(order = 1),
  #   shape = guide_legend(order = 2)
  # )+
  theme(
    # Position the legend using relative coordinates (x, y)
    legend.position = c(0.9, 0.05),
    # Optionally adjust the justification so that the legend is anchored correctly
    legend.justification = c("right", "center"),
    # Reduce the left margin of the legend box to nudge it further left
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  ) + labs(color = "stdev")

p22

p11



library(ggpubr)
p11p22 = ggarrange(p11, p22)
ggsave(p11p22, width = 12, height = 4, filename = "fig_micobin_muhat_01.pdf")
ggsave(p11p22, width = 12, height = 4, filename = "fig_micobin_muhat_01.png", bg = "white")
ggsave(p11p22, width = 12, height = 4, filename = "fig_micobin_muhat_01.jpg", bg = "white")



p33 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4, aes(x = x, y = y, color = wnew_hat), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "mako") +
  #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
  # guides(
  #   color = guide_colorbar(order = 1),
  #   shape = guide_legend(order = 2)
  # )+
  theme(
    # Position the legend using relative coordinates (x, y)
    legend.position = c(0.9, 0.05),
    # Optionally adjust the justification so that the legend is anchored correctly
    legend.justification = c("right", "center"),
    # Reduce the left margin of the legend box to nudge it further left
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )+ labs(color = "latent NNGP")
p33


p44 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4, aes(x = x, y = y, color = wnew_sd), size = 0.5, alpha = 0.5) +
  #scale_color_viridis_c(option = "plasma", direction = +1, limit = c(0,0.1)) +
  scale_colour_distiller(palette = "Spectral")+
  #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
  # guides(
  #   color = guide_colorbar(order = 1),
  #   shape = guide_legend(order = 2)
  # )+
  theme(
    # Position the legend using relative coordinates (x, y)
    legend.position = c(0.9, 0.05),
    # Optionally adjust the justification so that the legend is anchored correctly
    legend.justification = c("right", "center"),
    # Reduce the left margin of the legend box to nudge it further left
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  ) + labs(color = "p.p.d. stdev")

p44



library(ggpubr)
p33p44 = ggarrange(p33, p44)
ggsave(p33p44, width = 12, height = 4, filename = "fig_micobin_uhat.pdf")
ggsave(p33p44, width = 12, height = 4, filename = "fig_micobin_uhat.png", bg = "white")
ggsave(p33p44, width = 12, height = 4, filename = "fig_micobin_uhat.jpg", bg = "white")

# 
# ##################################################################################
# 
# 
# pred_beta1 = readRDS("../real2/pred_beta1.rds")
# 
# transformed_coords4 <- sf::st_coordinates(usmap::usmap_transform(df4[,c("lon","lat")]))
# 
# drawdf4 = data.frame(
#   x = transformed_coords4[,1],
#   y = transformed_coords4[,2],
#   urbmdhi = df4$logurbmdhi,
#   logconif = df4$logconif,
#   logconif01 = factor(df4$conif==0),
#   wnew_hat = pred_beta1$wnew_hat,
#   wnew_sd = pred_beta1$wnew_sd,
#   munew_hat = pred_beta1$munew_hat,
#   munew_sd = pred_beta1$munew_sd
# )
# 
# 
# p11 = plot_usmap("states", exclude = c("HI","AK")) +
#   geom_point(data = drawdf4, aes(x = x, y = y, color = munew_hat), size = 0.5, alpha = 0.5) +
#   scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
#   #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
#   # guides(
#   #   color = guide_colorbar(order = 1),
#   #   shape = guide_legend(order = 2)
#   # )+
#   theme(
#     # Position the legend using relative coordinates (x, y)
#     legend.position = c(0.9, 0.05),
#     # Optionally adjust the justification so that the legend is anchored correctly
#     legend.justification = c("right", "center"),
#     # Reduce the left margin of the legend box to nudge it further left
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
#   )+ labs(color = "predicted MMI")
# p11
# 
# 
# p22 = plot_usmap("states", exclude = c("HI","AK")) +
#   geom_point(data = drawdf4, aes(x = x, y = y, color = munew_sd), size = 0.5, alpha = 0.5) +
#   #scale_color_viridis_c(option = "plasma", direction = +1, limit = c(0,0.1)) +
#   scale_colour_distiller(palette = "Spectral")+
#   #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
#   # guides(
#   #   color = guide_colorbar(order = 1),
#   #   shape = guide_legend(order = 2)
#   # )+
#   theme(
#     # Position the legend using relative coordinates (x, y)
#     legend.position = c(0.9, 0.05),
#     # Optionally adjust the justification so that the legend is anchored correctly
#     legend.justification = c("right", "center"),
#     # Reduce the left margin of the legend box to nudge it further left
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
#   ) + labs(color = "p.p.d. stdev")
# 
# p22
# 
# p11
# 
# 
# 
# library(ggpubr)
# p11p22 = ggarrange(p11, p22)
# p11p22
# ggsave(p11p22, width = 12, height = 4, filename = "fig_beta_muhat_01.pdf")
# ggsave(p11p22, width = 12, height = 4, filename = "fig_beta_muhat_01.png", bg = "white")
# ggsave(p11p22, width = 12, height = 4, filename = "fig_beta_muhat_01.jpg", bg = "white")
# 
# 
# 
# p33 = plot_usmap("states", exclude = c("HI","AK")) +
#   geom_point(data = drawdf4, aes(x = x, y = y, color = wnew_hat), size = 0.5, alpha = 0.5) +
#   scale_color_viridis_c(option = "mako") +
#   #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
#   # guides(
#   #   color = guide_colorbar(order = 1),
#   #   shape = guide_legend(order = 2)
#   # )+
#   theme(
#     # Position the legend using relative coordinates (x, y)
#     legend.position = c(0.9, 0.05),
#     # Optionally adjust the justification so that the legend is anchored correctly
#     legend.justification = c("right", "center"),
#     # Reduce the left margin of the legend box to nudge it further left
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
#   )+ labs(color = "latent NNGP")
# p33
# 
# 
# p44 = plot_usmap("states", exclude = c("HI","AK")) +
#   geom_point(data = drawdf4, aes(x = x, y = y, color = wnew_sd), size = 0.5, alpha = 0.5) +
#   #scale_color_viridis_c(option = "plasma", direction = +1, limit = c(0,0.1)) +
#   scale_colour_distiller(palette = "Spectral")+
#   #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
#   # guides(
#   #   color = guide_colorbar(order = 1),
#   #   shape = guide_legend(order = 2)
#   # )+
#   theme(
#     # Position the legend using relative coordinates (x, y)
#     legend.position = c(0.9, 0.05),
#     # Optionally adjust the justification so that the legend is anchored correctly
#     legend.justification = c("right", "center"),
#     # Reduce the left margin of the legend box to nudge it further left
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
#   ) + labs(color = "p.p.d. stdev")
# 
# p44
# 
# 
# 
# library(ggpubr)
# p33p44 = ggarrange(p33, p44)
# ggsave(p33p44, width = 12, height = 4, filename = "fig_beta_uhat.pdf")
# ggsave(p33p44, width = 12, height = 4, filename = "fig_beta_uhat.png", bg = "white")
# ggsave(p33p44, width = 12, height = 4, filename = "fig_beta_uhat.jpg", bg = "white")
# 
# 
# 
# 
# pcompare = plot_usmap("states", exclude = c("HI","AK")) +
#   geom_point(data = drawdf4, aes(x = x, y = y, color = (pred_micobin1$wnew_hat - pred_beta1$wnew_hat) ), size = 0.5, alpha = 0.5) +
#   scale_colour_distiller(palette = "Spectral")+
#   #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
#   # guides(
#   #   color = guide_colorbar(order = 1),
#   #   shape = guide_legend(order = 2)
#   # )+
#   theme(
#     # Position the legend using relative coordinates (x, y)
#     legend.position = c(0.9, 0.05),
#     # Optionally adjust the justification so that the legend is anchored correctly
#     legend.justification = c("right", "center"),
#     # Reduce the left margin of the legend box to nudge it further left
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
#   )+ labs(color = "micobin - beta (latent NNGP)")
# pcompare
# ggsave(pcompare, width = 6, height = 4, filename = "fig_beta_micobin_uhat_compare.pdf")
# ggsave(pcompare, width = 6, height = 4, filename = "fig_beta_micobin_uhat_compare.png")
#
# df$state[which(df$conif==0)]
# summary(df$conif)
#
# df$conif_zeroone = df$conif<5
# which(table(df$conif_zeroone, df$state)[1,]==0)
#
#
#
# fields::quilt.plot(cbind(df$lon, df$lat), df$conif==0)
#
#
# length(unique(df$state))
# # Load necessary package
# library(ggplot2)
#
# # Get map data for the 48 contiguous US states
# states_map <- map_data("state")
#
# # Define states to be highlighted (in lowercase)
# highlighted_states <- c("connecticut", "iowa", "kansas", "kentucky",
#                         "nebraska", "north dakota", "ohio",
#                         "pennsylvania", "west virginia")
#
# # Create a new column indicating whether the state should be highlighted
# states_map$highlight <- ifelse(states_map$region %in% highlighted_states,
#                                "highlight", "other")
#
# # Plot the map
# ggplot() +
#   geom_polygon(data = states_map,
#                aes(x = long, y = lat, group = group, fill = highlight),
#                color = "black") +
#   scale_fill_manual(values = c("highlight" = "red", "other" = "grey90")) +
#   coord_fixed(1.3) +
#   theme_void() +
#   ggtitle("US 48 States with Highlighted States")
