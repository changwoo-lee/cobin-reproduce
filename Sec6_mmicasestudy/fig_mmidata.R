rm(list = ls())

library(ggplot2)
library(maps)
library(usmap)
##library(ggrepel) # if need to repel labels


# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(usmap)
library(sf)
library(ggplot2)

df = read.csv("mmi_lakecat.csv")
df = df[df$MMI_BENT!=0,]
df4 = read.csv("lakecat-over40000m2.csv")
df4$logurbmdhi = log2(1+df4$urbmdhi)
df4$logconif = log2(1+df4$conif)
df4$logconif01 = factor(df4$conif==0)

# Lat/Lon of Sioux Falls, SD
test_data <- data.frame(lon = -96.70, lat = 43.55)
df$MMI_BENT
transformed_coords <- sf::st_coordinates(usmap::usmap_transform(df[,c("lon","lat")]))

drawdf = data.frame(
  x = transformed_coords[,1],
  y = transformed_coords[,2],
  MMI = df$MMI_BENT,
  conif_is_0 = factor(df$conif==0)
)

p1 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf, aes(x = x, y = y, color = MMI), size = 1) +
  scale_color_viridis_c(option = "turbo", direction = -1, limits = c(0,1)) +
  #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
  guides(
    color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)
  )+
  theme(
    # Position the legend using relative coordinates (x, y)
    legend.position = c(0.9, 0.05),
    # Optionally adjust the justification so that the legend is anchored correctly
    legend.justification = c("right", "center"),
    # Reduce the left margin of the legend box to nudge it further left
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )

p1

df$state[which(df$conif==0)]
summary(df$conif)

df$conif_zeroone = df$conif<5
which(table(df$conif_zeroone, df$state)[1,]==0)







transformed_coords4 <- sf::st_coordinates(usmap::usmap_transform(df4[,c("lon","lat")]))

drawdf4 = data.frame(
  x = transformed_coords4[,1],
  y = transformed_coords4[,2],
  urbmdhi = df4$logurbmdhi,
  logconif = df4$logconif,
  logconif01 = factor(df4$conif==0)
)

p2 = plot_usmap("states", exclude = c("HI","AK")) +
  geom_point(data = drawdf4, aes(x = x, y = y, color = urbmdhi), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(option = "A", direction = -1) +
  #scale_shape_manual(name = NULL, values = c("circle","triangle"), label = c("conif>0","conif=0"))+
  guides(
    color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)
  )+
  theme(
    # Position the legend using relative coordinates (x, y)
    legend.position = c(0.9, 0.05),
    # Optionally adjust the justification so that the legend is anchored correctly
    legend.justification = c("right", "center"),
    # Reduce the left margin of the legend box to nudge it further left
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10)
  )
p2

library(ggpubr)
p1p2 = ggarrange(p1, p2)
p1p2
#ggsave(p1p2, width = 12, height = 4, filename = "fig_mmidata.pdf")
#ggsave(p1p2, width = 12, height = 4, filename = "fig_mmidata.png", bg = "white")
