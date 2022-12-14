##########################################   Start Functions    ###############################################
# This function takes input presence data points (data.frame with lon,lat) and creates
# a buffer (of specified distance), then rasterizes it for a given mask. Finally, it
# saves the buffered raster in a specified file. It also returns the raster.
# @param (data.frame) xy: presence locations
# @param (raster) mask: mask that gives the native area (for geographic extent and resolution)
# @param (numeric) buff_dist: distance of buffer in degree or km (as needed)
# @param (string) format: output file format
# @param (string) filename: output file name
# @return (raster): rasterized buffer

create_buffers <- function(xy, msk, buff_dist = 0.5, format = "GTiff", filename) {
  
  
  msk <- raster(msk)
  
  if(!file.exists(paste0(filename, "/buffer_radius_to_omit.tif"))){
    
    # Ensure msk has a CRS assigned to
    proj4string(msk) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
    
    # Making points spatial object with coordinates
    xy_coords <- as.data.frame(cbind(xy$lon, xy$lat))
    xy_coords <- xy_coords[complete.cases(xy_coords),]
    xy_coords <- unique(xy_coords)
    colnames(xy_coords) <- c("x","y")
    coordinates(xy_coords) <- ~x+y
    proj4string(xy_coords) <- CRS("+proj=longlat +datum=WGS84")
    
    # Buffering
    buffer <- gBuffer(xy_coords, width = buff_dist)
    
    # Rasterizing and making it into a mask
    buffer_rs <- rasterize(buffer, msk)
    buffer_rs[which(!is.na(buffer_rs[]))] <- 1
    buffer_rs[which(is.na(buffer_rs[]) & msk[] == 1)] <- 0
    buffer_rs[which(is.na(msk[]))] <- NA
    
    # Writing raster
    writeRaster(buffer_rs, paste0(filename, "/buffer_radius_to_omit.tif"), format = format, overwrite = T)
    writeOGR(obj = as(buffer, "SpatialPolygonsDataFrame"), dsn = filename, layer = "buffer_radius_to_omit_shp", driver = "ESRI Shapefile", overwrite_layer = T)
    
  } else {
    # Load raster in case it exists
    buffer_rs <- raster(paste0(filename, "/buffer_radius_to_omit.tif"))
  }
  
  return(buffer_rs)
  
}

# testing the function
# base_dir <- "~/nfs"
# source("~/Repositories/aichi13/src/config.R")
# config(dirs=T, )
# xy <- read.csv(paste(gap_dir,"/2686262/",run_version,"/occurrences/2686262.csv",sep=""))[,c("lon","lat")]
# load(paste(gap_dir,"/2686262/",run_version,"/bioclim/crop_narea.RDS",sep=""))
# msk <- biolayers_cropc[[1]]; rm(biolayers_cropc)
# msk[which(!is.na(msk[]))] <- 1
# x <- create_buffers(xy, msk, buff_dist=0.5, format="GTiff", filename=paste(gap_dir,"/2686262/",run_version,"/modeling/alternatives/ca50_total_narea.tif",sep=""))
