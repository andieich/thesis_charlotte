#27.05.2024

library(httr)
library(here)
library(tidyverse)
library(pbapply) #with progress bar, or use "normal" lapply
library(raster) #DON'T USE terra PACKAGE!
library(lubridate)


make_url <- function(lat1,
                     lat2,
                     lon1,
                     lon2,
                     day1,
                     day2,
                     parameter){

  if (parameter == "CRW_SST"){
    
    base_url <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.nc?CRW_SST"
    
  } else if (parameter == "CRW_DHW"){
    
    base_url <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.nc?CRW_DHW"
    
  } 

  url <- paste0(base_url,
                "[(",
                day1, "):1:(",
                day2, ")][(",

                lat1, "):1:(",
                lat2, ")][(",

                lon1, "):1:(",
                lon2, ")]"
  )

  return(url)

}






download_dhw <- function(lat1,
                         lat2,
                         lon1,
                         lon2,
                         day1,
                         day2,
                         url,
                         folder,
                         parameter){
  
  timeout = 600


  filename <- paste0(paste(lat1, lat2, day1, day2,parameter, sep = "_"),
                     ".nc")

  path <- paste(folder, filename, sep = "/")

  cat("Will download between",
      as.character(day1),
      "and",
      as.character(day2),
      "from",
      url,
      "to",
      path)


  httr::GET(url,
            httr::write_disk(path,
                             overwrite = T),
            httr::add_headers('Accept-Encoding' = 'gzip, deflate'),
            httr::progress("down"),
            httr::timeout(timeout))
}



download_all <- function(lat1,
                         lat2,
                         lon1,
                         lon2,
                         day1,
                         day2,
                         tresh = 365*3,
                         folder,
                         parameter = "CRW_SST"){
  
  if (!parameter%in% c("CRW_SST", "CRW_DHW")){
    
    stop(paste("Parameter has to be CRW_SST or CRW_DHW but is", parameter))
    
  }
  day1 <- as.Date(day1)
  day2 <- as.Date(day2)

  day1_new <- day1

  diff <-  as.numeric(day2 - day1)
  
  

  while (diff > tresh) {

    print("New while: day1_new")
    print(day1_new)

    day2_new <- day1_new + tresh

    print("day2_new")
    print(day2_new)


    day1_new <- as.character(day1_new)
    day2_new <- as.character(day2_new)

    url <- make_url(lat1 = lat1,
                    lat2 = lat2,
                    lon1 = lon1,
                    lon2 = lon2,
                    day1 = day1_new,
                    day2 = day2_new,
                    parameter = parameter)

    download_dhw(lat1 = lat1,
                 lat2 = lat2,
                 lon1 = lon1,
                 lon2 = lon2,
                 day1 = day1_new,
                 day2 = day2_new,
                 url = url,
                 folder = folder,
                 parameter = parameter)

    day1_new <- as.Date(day1_new)
    day2_new <- as.Date(day2_new)

    day1_new <- day2_new + 1

    print("end while, new d1_new")
    print(day1_new)

    diff <-  as.numeric(day2 - day1_new)

  }

  url <- make_url(lat1 = lat1,
                  lat2 = lat2,
                  lon1 = lon1,
                  lon2 = lon2,
                  day1 = day1_new,
                  day2 = day2,
                  parameter = parameter)

  download_dhw(lat1 = lat1,
               lat2 = lat2,
               lon1 = lon1,
               lon2 = lon2,
               day1 = day1_new,
               day2 = day2,
               url = url,
               folder = folder,
               parameter = parameter)

}





#Read .nc files


# Reads separate files, if lat/lon borders are given, filters the data
# The raster package is used to read .nc files and transforms them into a data.frame
# Inputs
# path: path to .nc file
# lat_min: lower latitude border. If set, will be included
# lat_max: upper latitude border. If set, will be included
# lon_min: lower longitude border. If set, will be included
# lon_max: upper longitude border. If set, will be included
# varname used in .nc file.
# values_to:  column name for varname

read_and_format <- function(path,
                            lat_min = NA,
                            lat_max = NA,
                            lon_min = NA,
                            lon_max = NA,
                            varname = varname,
                            values_to = values_to){


  dat <- raster::brick(path,
                       varname = varname) %>%
    raster::as.data.frame(xy = T) %>%
    pivot_longer(3:ncol(.), names_to = "date", values_to = values_to) %>%
    mutate(date = substr(date,2,11)) %>%
    mutate(date = as.Date(lubridate::fast_strptime(date, "%Y.%m.%d"))) %>%
    mutate(lat = round(y,3),
           lon = round(x, 3)) %>%
    dplyr::select(-x, -y)


  if(!is.na(lat_min) & !is.na(lat_max)){

    if(lat_min > lat_max){
      lat_tmp <- lat_min
      lat_min <- lat_max
      lat_max <- lat_tmp}

    dat <- dat %>%
      filter(lat >= lat_min, lat <= lat_max)
  }

  if(!is.na(lon_min) & !is.na(lon_max)){

    if(lon_min > lon_max){
      lon_tmp <- lon_min
      lon_min <- lon_max
      lon_max <- lon_tmp}

    dat <- dat %>%
      filter(lon >= lon_min, lon <= lon_max)
  }

  return(dat)
}


# Helper function to get file extension. Used to check if input file is really .nc

getExtension <- function(path){
  ext <- strsplit(path, "\\.")[[1]]
  ext <-tail(ext, 1)
  return(ext)
}


# Reads .nc files
# Can handle folders or separate files
# if lat/lon _min/_max is set, data will be filtered (including borders)
# if path is folder, all .nc files will be combined

read_nc <- function(path,
                    lat_min = NA,
                    lat_max = NA,
                    lon_min = NA,
                    lon_max = NA,
                    varname = "CRW_SST",
                    values_to = "sst"){



  if(dir.exists(path)){

    files <- list.files(path = path,
                        pattern = "\\.nc$",
                        full.names = T)

    if(length(files) == 0){
      stop("Folder does not contain .nc files")
    } else {
      files <- list.files(path = path,
                          pattern = "\\.nc$",
                          full.names = T)

      print(paste("Will read", length(files), "files from folder"))

      dat <- do.call(rbind,
                     pbapply::pblapply(files,
                                       read_and_format,
                                       lat_min,
                                       lat_max,
                                       lon_min,
                                       lon_max,
                                       varname,
                                       values_to))

    }

  } else if(file.exists(path)){

    if(getExtension(path) != "nc"){
      stop("File is no .nc file")
    } else {

      print("Will read one file")
      dat <-  read_and_format(path = path,
                              lat_min = lat_min,
                              lat_max = lat_max,
                              lon_min = lon_min,
                              lon_max = lon_max,
                              varname = varname,
                              values_to = values_to)
    }

  } else {
    stop("File or folder not found")
  }

  return(dat)

}


bounding_box <- function(center_lat, center_lon, width_km, length_km) {
  # Earth's radius in km
  earth_radius <- 6371.0
  
  # Convert width and length to radians
  width_rad <- width_km / earth_radius
  length_rad <- length_km / earth_radius
  
  # Convert center coordinates to radians
  center_lat_rad <- center_lat * (pi / 180)
  center_lon_rad <- center_lon * (pi / 180)
  
  # Calculate bounding box
  lat_min <- center_lat_rad - width_rad
  lat_max <- center_lat_rad + width_rad
  lon_min <- center_lon_rad - length_rad
  lon_max <- center_lon_rad + length_rad
  
  # Convert back to degrees
  lat_min <- lat_min * (180 / pi)
  lat_max <- lat_max * (180 / pi)
  lon_min <- lon_min * (180 / pi)
  lon_max <- lon_max * (180 / pi)
  
  return(list(lat_min = lat_min, lat_max = lat_max, lon_min = lon_min, lon_max = lon_max))
}


subset_ncdf <- function(ncdf_file, lat_min, lat_max, lon_min, lon_max) {
  # Create an extent object with the provided coordinates
  extent_obj <- extent(c(lon_min, lon_max, lat_min, lat_max))
  
  # Crop the NetCDF4 file with the extent
  cropped_ncdf <- crop(ncdf_file, extent_obj)
  
  #to df
  cropped_ncdf_df <- raster::as.data.frame(cropped_ncdf, xy = T) %>% 
    mutate(lat = round(y,3),
           lon = round(x, 3)) %>%
    dplyr::select(-x, -y) %>% 
    rename_with(.cols = 1, ~"mmm")
  
  
  
  return(cropped_ncdf_df)
}




find_closest_point <- function(data, lat, lon) {
  data <- as.data.frame(data)
  # Calculate the distance between each point in the data and the given point
  distances <- distm(data[, c("lon", "lat")], c(lon, lat), fun = distVincentySphere)
  
  # Find the index of the minimum distance
  closest_point_index <- which.min(distances)

  
  # Get the closest lat and lon values
  closest_lat <- data[closest_point_index, "lat"]
  closest_lon <- data[closest_point_index, "lon"]

  

  
  # Return all rows that match the closest lat and lon
  return(data %>% 
           filter(lat == closest_lat, lon == closest_lon))

}
