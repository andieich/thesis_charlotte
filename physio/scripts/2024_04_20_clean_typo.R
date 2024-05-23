# I made a type in the data export setup in the photometer setup
# the second 'abs_663' columns should be 'abs_630'
# This script corrects this
# I also corrected the data export in the photometer software

library(here)
files <- list.files(here('chlorophyll/data/raw_data'),pattern = '.csv',recursive = T)

files <- paste(here('chlorophyll/data/raw_data'), files, sep = '/')

for (file_i in files) {
  dat_temp <-  read.csv(file_i)
  names(dat_temp)[4] <- "abs_630"
  write.csv(dat_temp,file = file_i, row.names = F)
  
}






