read_and_clean <- function(file){
  dat <- read_csv(file, show_col_types = FALSE)[,2:7]#get rid of type and blank corrected values (calculate bc below)
  dat <- dat %>% 
    clean_names() %>% 
    mutate(position = as.character(position)) %>% 
    mutate_at(vars(2:6),
              list(as.numeric)) %>% 
    filter(!if_all(2:6, ~ is.na(.))) %>% # filter if complete column is NA
    mutate(filename = basename(file) %>% 
             substr(1, nchar(.)-4))
  
  return(dat)
}


read_all_files <- function(folder){
  # list all csv files
  csv_files <- list.files(folder, pattern = '.csv')
  
  csv_paths <- paste(folder , csv_files, sep = '/')
  
  dat <- csv_paths %>% 
    map_dfr(read_and_clean)
  
  return(dat)
  
}


get_area <- function(data){
  data %>% 
    mutate(delta_w = w2_scratched - w1_scratched) %>% 
    mutate(area = delta_w * 34.32) %>% 
    select(sample_id, bag_number, area)
}



get_counts <- function(data_counts, 
                       data_overview, 
                       data_area, 
                       v_ml_count = 0.1,
                       v_ml_sw_added = 3,
                       v_ml_sw_pipetted = 1){
  
  # all sample_ids from data_counts in data_overview and data_area?
  missing_in_ov <- data_counts$sample_id[which(!data_counts$sample_id %in% data_overview$sample_id)]
  
  if (length(missing_in_ov) > 0){
    warning(paste("Sample ID", missing_in_ov, "is in data_counts but not in data_overview \n"))
  }
  
  missing_in_area <- data_counts$sample_id[which(!data_counts$sample_id %in% data_area$sample_id)]
  
  if (length(missing_in_area) > 0){
    warning(paste("Sample ID", missing_in_area, "is in data_counts but not in data_area \n"))
  }
  

  
  data_counts %>% 
    left_join(data_overview %>% 
                select(sample_id, w1, w2), by = "sample_id") %>% #add overview data (for Volume)
    left_join(data_area %>% 
                select(sample_id, area), by = "sample_id") %>%   # add area data
    select(-mean_count_per_sample) %>% #get rid of mean
    pivot_longer(cols = c1:c6, names_to = "count_replicate", values_to = "count_subsample") %>% #make long
    mutate(dillution_factor = v_zoox / (v_zoox + v_sw)) %>%
    mutate(count_per_mL = 1000 * count_subsample / (dillution_factor * v_ml_count)) %>%
    mutate(count_per_sample = count_per_mL * (v_ml_sw_added / v_ml_sw_pipetted) * (w2 / w1)) %>%
    mutate(count_per_cm2 = count_per_sample / area) %>%
    select(sample_id, bag_number, count_per_cm2)
}

# why not multiplied with total volume in Laets sheet?



get_chlorophyll <- function(data_chlorophyll, 
                            data_chlorophyll_overview, 
                            data_area, 
                            data_overview,
                            pl = 1,
                            v_ml_sw_added = 3,
                            v_ml_sw_pipetted = 1){

  
  # remove overview for 24 h data (will be added later again, to make 24 h and 48 h data more similar)
  
  if("date_extraction" %in% names(data_chlorophyll_overview) & 
     "w1" %in% names(data_chlorophyll_overview) & 
     "w2" %in% names(data_chlorophyll_overview))
    {
    data_chlorophyll_overview <- data_chlorophyll_overview %>% 
      select(-date_extraction, -w1, -w2)
  }
  
  # add absorption data to chl overview data
    # check if file names are missing in chl overview (which filename from the combined chlorophyll data does not exist in chl overview data?)
  
  missing_fn_in_ov <- data_chlorophyll$filename[which(!data_chlorophyll$filename %in% data_chlorophyll_overview$filename)]
  
  # if a filename is missing, print a warning
  
  if (length(missing_fn_in_ov) > 0){
    warning(paste("Filename", missing_fn_in_ov, "is in data_chlorophyll but not in data_chlorophyll_overview \n"))
  }
  
  
  # make long, i.e. position of replicate measurement in 2 rows
  data_chlorophyll_overview <- data_chlorophyll_overview %>% 
    pivot_longer(cols = c('m1', 'm2'), 
                 names_to = 'measurement', 
                 values_to = 'position')
  

  # test if unique positions
  dup_pos_in_chl_ov <- data_chlorophyll_overview %>% 
    group_by(filename) %>% 
    filter(duplicated(position))
  
  if (nrow(dup_pos_in_chl_ov) > 0){
    warning("Some positions in data_chlorophyll_overview are duplicated: \n", print_and_capture(dup_pos_in_chl_ov))
  }
  
  dup_pos_in_chl_dat <- data_chlorophyll %>% 
    group_by(filename) %>% 
    filter(duplicated(position))
  
  if (nrow(dup_pos_in_chl_dat) > 0){
    warning("Some positions in data_chlorophyll_overview are duplicated: \n", print_and_capture(dup_pos_in_chl_dat))
  }
  

  # test if any position does not exist in overview data
  
  test_pos <- left_join(data_chlorophyll %>% 
                          select(filename, position), 
                        data_chlorophyll_overview,
                        by = c("filename", "position"))
  
  
  pos_not_exist<- test_pos %>% 
    filter(is.na(type)) 
  
  if (nrow(pos_not_exist) > 0){
    warning(paste("Some filnames or positions in data_chlorophyll do not exist in data_chlorophyll_overview  \n"))
    print(pos_not_exist)
  }

  # combine chl ov and chl data
  data_chlorophyll <- data_chlorophyll_overview %>% 
    left_join(data_chlorophyll, 
              by = c("filename", "position"))
  #test if any file name or position from chl overview does not exist (NA after left join)
  
  fn_not_exist<- data_chlorophyll %>% 
    filter(is.na(abs_663))

  if (nrow(fn_not_exist) > 0){
    warning("Some filnames or positions in data_chlorophyll_overview do not exist in data_chlorophyll: \n", print_and_capture(fn_not_exist))
  }
  
  
  

  # add area
  # check: All sample IDs from chl data in area data?
  missing_ID_in_area <- data_chlorophyll$sample_id[which(!data_chlorophyll$type == "blank" & !data_chlorophyll$sample_id %in% data_area$sample_id)]
  
  if (length(missing_ID_in_area) > 0){
    warning(paste("Sample ID", missing_ID_in_area, "is in data_chlorophyll but not in data_area \n"))
  }
  
  # check: All sample IDs from area data in chl data?
  missing_area_ID_in_chl <- data_area$sample_id[which(!data_area$sample_id %in% data_chlorophyll$sample_id)]
  
  if (length(missing_area_ID_in_chl) > 0){
    warning(paste("Sample ID", missing_area_ID_in_chl, "is in data_area but not in data_chlorophyll \n"))
  }

  # OK to combine
  
  data_chlorophyll <- data_chlorophyll %>% 
    left_join(data_area %>% 
                select(sample_id, area), 
              by = "sample_id")
  

 

  # add dat ov
  
  # same tests for dat overview:
  # check: All sample IDs from chl data in ov data?
  missing_ID_in_dat_ov <- data_chlorophyll$sample_id[which(!data_chlorophyll$type == "blank" & !data_chlorophyll$sample_id %in% data_overview$sample_id)]
  
  if (length(missing_ID_in_dat_ov) > 0){
    warning(paste("Sample ID", missing_ID_in_dat_ov, "is in data_chlorophyll but not in data_overview \n"))
  }
  
  # check: All sample IDs from chl ov data in chl data?
  missing_dat_ov_ID_in_chl <- data_overview$sample_id[which(!data_overview$type == "blank" & !data_overview$sample_id %in% data_chlorophyll$sample_id)]
  
  if (length(missing_dat_ov_ID_in_chl) > 0){
    warning(paste("Sample ID", missing_dat_ov_ID_in_chl, "is in data_overview but not in data_chlorophyll \n"))
  }

  
  #ok to combine
  data_chlorophyll <- data_chlorophyll %>% 
    left_join(data_overview %>% 
                filter(!is.na(sample_id)) %>% #get rid of blanks
                select(sample_id, w1, w2), by = "sample_id")
  

  
  # separate blanks and calculate mean per batch
  blanks <- data_chlorophyll %>% 
    filter(type == "blank") %>% 
    group_by(batch_number) %>% 
    summarise(mean_abs_663 = round(mean(abs_663),3),
              mean_abs_630 = round(mean(abs_630),3),
              mean_abs_750 = round(mean(abs_750),3),
              mean_abs_470 = round(mean(abs_470),3)) %>% 
    ungroup()
  
  
  # add blank means to original data
  
  data_chlorophyll <- data_chlorophyll %>% 
    left_join(blanks, by = "batch_number") %>% 
    mutate(abs_663_c = round(abs_663 - mean_abs_663,3),
           abs_630_c = round(abs_630 - mean_abs_630,3),
           abs_750_c = round(abs_750 - mean_abs_750,3),
           abs_470_c = round(abs_470 - mean_abs_470,3)) 
  


  
  
  data_chlorophyll %>% 
    group_by(type) %>% 
    summarise_at(c("abs_663_c", "abs_630_c", "abs_750_c", "abs_470_c"), 
                 list(min = min, max = max, mean = mean), 
                 na.rm = TRUE) %>% 
    pivot_longer(abs_663_c_min:abs_470_c_mean, 
                 names_to = "parameter", 
                 values_to = "value") %>% 
    mutate(wavelength = substr(parameter, 5,7),
           summary = substr(parameter, 11,14),
           value = round(value,3)) %>% 
    select(wavelength, summary, value, type) %>% 
    arrange(wavelength, summary, value) %>% 
    pivot_wider(names_from = type, values_from = value) %>% 
    print()
  
  # filter only samples
  
  data_chlorophyll <- data_chlorophyll %>% 
    filter(type == "sample")
  
  # calculate chlorophyll with formula, correct for volume
  data_chlorophyll <- data_chlorophyll %>% 
    mutate(chl_a_subsample = 11.43 * (abs_663_c - abs_750_c) / pl - 0.64 * (abs_630_c - abs_750_c)/pl) %>% #in µg per mL
    mutate(chl_a_per_sample = chl_a_subsample * (v_ml_sw_added / v_ml_sw_pipetted) * (w2 / w1)) %>% # µg in whole sample. of V_slurry, 40 mL were centrifuged
    #and filled up with 3 mL SSW, 1 mL of SSW was centrifuged and 1 mL acetone added, 0.5 mL acetone measured.
    # c in 0.5 mL is the absolute amount in the 1 mL, 3 is multiplied because after centrifucation of the 40 mL, only 1 mL is taken for chl. 
    # This is the absolute amount in the 40 mL that were centrifuged
    # (w2 / w1) is multiplied to consider that V_slurry was more than 40 mL
    mutate(chl_a_per_cm2 = chl_a_per_sample/area)
  
  # clean
  data_chlorophyll <- data_chlorophyll %>%
    select(sample_id,
           bag_number,
           chl_a_per_cm2)
  
  return(data_chlorophyll)
  
}


get_info <- function(data){
  data %>% 
    mutate(
      t = substr(sample_id, 1,1),
      site = substr(sample_id, 2,2),
      depth = substr(sample_id, 3,3),
      spec = substr(sample_id, 4,4),
      size = substr(sample_id, 5,5),
      replicate = as.numeric(substr(sample_id, 6,6))) %>% 
    mutate(site = recode(site,
                         "T" = "Temae",
                         "E" = "E2B"),
           depth = recode(depth,
                          "D" = "20 m",
                          "S" = "5 m"),
           spec = recode(spec,
                         "A" = "A. hya",
                         "P" = "P. ver",
                         "M" = "P. mea"),
           size = recode(size,
                         "L" = "Adult",
                         "S" = "Juvenile"))
}

print_and_capture <- function(x)
{
  paste(capture.output(print(as.data.frame(x))), collapse = "\n")
}


