# NHMS projects 09589.
# Analysis of 2019-2022 samples

setwd("N:\\PAPER NHMS Jenny 2023\\")
######################################################################
# Fig 2 -  plots of sample locations

# spatial
library(raster)
library(ggmap)
library(tmap)# https://r-spatialecology.github.io/ialena-2020/#9
library(exactextractr) #  extract summed raster values from an overlay polygon
library(mapview)  #  simple visualisatin of sample location

#-------------------------------------------------------------------

# You will need to request the UKCEH landcover map from https://www.ceh.ac.uk/data/ukceh-land-cover-maps
#in this case We are using the 25m raster version for 2021 (LCM2021_25m.tif).

# Note if you do not do this and  do not have the UKCEH Land Cover Map you can
# plot the location of the samples using the mapview package.  However, this will only
# be the spatial locations and not the underlying land cover classifications

LCM_2021 <- raster(paste0("LCM2021_25m.tif"))
crs(LCM_2021) #  coordinate reference system
#  +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs 
# for each dataset there is a Great Britain version in the

#------------------------------------------------------------------------------------------------------
# 21 LCM classes -  see help pdf file lcm2017-2019product_documentation_v1_5_1.pdf provided as meta data with UKCEH map 
# requests

# reclassify the values into 3 groups relating to urban/suburban, arable, and other
# all values > 0 and <= 1 become 1, etc.
m <- c(
  0, 2, 0,
  2, 3, 1,
  3, 19, 0,
  19,  21, 2)
rclmat <- matrix(m, ncol=3, byrow=TRUE) #  a 3 column matrix shoing c1=start value, c2=end value, c3 =  new descriptor)
LCM_reclasify21 <- reclassify(LCM_2021, rclmat)

plot(LCM_reclasify21)




#-----------------------------------------------------------------------------
# Locations of all sites where eDNA extractions of hive samples were undertkaen.
# For General Data Protection Regulations  in the UK we have resolved these to the nearest 1km
# Original maps presented in the paper were produced to a higher resolution


# 2018 sample year
location_2018<-read.csv("2018sites_1km.csv")
sp_2018<- SpatialPoints(location_2018[, c("X", "Y")])
crs(sp_2018) <-"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs "
mapview(sp_2018) #  simple visualization of locations


# 2019 sample year
location_2019<-read.csv("2019sites_1km.csv")
sp_2019<- SpatialPoints(location_2019[, c("X", "Y")])
crs(sp_2019) <-"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs "
mapview(sp_2019) #  simple visualization of locations

#2020 sample year
location_2020<-read.csv("2020sites_1km.csv")
sp_2020<- SpatialPoints(location_2020[, c("X", "Y")])
crs(sp_2020) <-"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs "
mapview(sp_2020) #  simple visualization of locations


#2021 sample year
location_2021<-read.csv("2021sites_1km.csv")
sp_2021<- SpatialPoints(location_2021[, c("X", "Y")])
crs(sp_2021) <-"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs "
mapview(sp_2021) #  simple visualization of locations

# 2022 sample year
location_2022<-read.csv("2022sites.csv")
sp_2022<- SpatialPoints(location_2022[, c("X", "Y")])
crs(sp_2022) <-"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs "
mapview(sp_2022) #  simple visualization of locations


read.csv("C:\\Users\\bawood\\Downloads\\pd10000_top40plant_maysept.csv")

#-------------------------------------------
# blank map plot showing arable and  urban land use

# create color object with nice new colors!
chm_colors <- c( "gray80", "chartreuse4",  "black")

# plot reclassified data
plot(LCM_reclasify21  ,
     legend = FALSE,
     col = chm_colors,
     axes = FALSE,
     # remove the box around the plot
     box = FALSE,
     main = "")


legend("topleft",
       legend = c("Other",	"Arable",
                  "Urban"),
       fill = chm_colors,
       border = FALSE,
       bty = "n")   #  https://www.datasciencemadesimple.com/add-legend-plot-legend-function-r/

#-------------------------------------------
# plot 2018



# plot reclassified data
plot(LCM_reclasify21  ,
     legend = FALSE,
     col = chm_colors,
     axes = FALSE,
     # remove the box around the plot
     box = FALSE,
     main = "")



points(sp_2018,   pch = 19, col = "blue", cex = .75)  


legend("topleft",
       legend = c("Other",	"Arable",
                  "Urban"),
       fill = chm_colors,
       border = FALSE,
       bty = "n")   #  https://www.datasciencemadesimple.com/add-legend-plot-legend-function-r/


#-------------------------------------------------
#2019

# plot reclassified data
plot(LCM_reclasify21  ,
     legend = FALSE,
     col = chm_colors,
     axes = FALSE,
     # remove the box around the plot
     box = FALSE,
     main = "")


points(sp_2019,   pch = 19, col = "blue", cex = .75)  



#-------------------------------------------------
#2020

# plot reclassified data
plot(LCM_reclasify21 ,
     legend = FALSE,
     col = chm_colors,
     axes = FALSE,
     # remove the box around the plot
     box = FALSE,
     main = "")


points(sp_2020,   pch = 19, col = "blue", cex = .75)  



#-------------------------------------------------
#2021

# plot reclassified data
plot(LCM_reclasify21  ,
     legend = FALSE,
     col = chm_colors,
     axes = FALSE,
     # remove the box around the plot
     box = FALSE,
     main = "")


points(sp_2021,   pch = 19, col = "blue", cex = .75)  

#-------------------------------------------------
#2022

# plot reclassified data
plot(LCM_reclasify21 ,
     legend = FALSE,
     col = chm_colors,
     axes = FALSE,
     # remove the box around the plot
     box = FALSE,
     main = "")


points(sp_2022,   pch = 19, col = "blue", cex = .75)  

#========================================================================================
# Fig. 3.  Regional proportions of OTU reads (Abundance) from the  different genera.  Note the final figure was modified in photoshop to add 
# Pictures of plants and the map -  but not to change the graphs.


library(ggplot2)

summary_genus_reads<-read.csv("Pivoted_summed_OTUreads_by_genus.csv")

# Colour palette
color_values <- c(
  "Brassica" = "#F7DC6F",   # yellow
  "Trifolium" = "#8E44AD",  # purple
  "Impatiens" = "#FF69B4",  # pink
  "Rubus" = "#87CEFA",      # light blue
  "Myosotis" = "#4682B4",   # steel blue
  "Other species" = "grey80"
)


# reordering for grpahs
summary_genus_reads$month <- factor(summary_genus_reads$month,
                                           levels = c("May", "June", "July", "August", "September"))

summary_genus_reads$genus_grouped  <- factor(summary_genus_reads$genus_grouped ,
                                                    levels = c("Other species", "Myosotis", "Rubus","Impatiens",  "Trifolium", "Brassica"))


month_ukband <- ggplot(summary_genus_reads, aes(x = month, y = mean_reads, fill = genus_grouped)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~region, scales = "free_x") +
  scale_fill_manual(values = color_values, name = "Plant genus") +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  labs(x = "Month sampled", y = "Proportion of DNA sequences") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold")
  ) 
print(month_ukband)

#======================================================================================
# Fig. 4 -  ladn use surrounding hives pie chart

landuse<-read.csv("LandUse.csv")


# Load required packages
# Load required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)


# --- Step 1: Define colour palette ---
cats <- unique(landuse$LandUse)
cols <- setNames(
  rep(brewer.pal(min(8, length(cats)), "Dark2"), length.out = length(cats)),
  cats
)

# --- Step 2: Create pie chart ---
ggplot(landuse, aes(x = "", y = Average, fill = LandUse)) +
  geom_bar(width = 1, stat = "identity", color = "black", linewidth = 0.4) +  # black borders
  coord_polar(theta = "y") +  # convert bar to pie
  scale_fill_manual(values = cols) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # wrap legend into 2 rows

#-------------------------------------------
#  Crop cover

landuse_crop<-read.csv("CropCover.csv")



# Define colour palette ---
cats <- unique(landuse_crop$LandUse)
cols <- setNames(
  rep(brewer.pal(min(8, length(cats)), "Dark2"), length.out = length(cats)),
  cats
)

# Create pie chart ---
ggplot(landuse_crop, aes(x = "", y = Average, fill = LandUse)) +
  geom_bar(width = 1, stat = "identity", color = "black", linewidth = 0.4) +  # black borders
  coord_polar(theta = "y") +  # convert bar to pie
  scale_fill_manual(values = cols) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # wrap legend into 2 rows