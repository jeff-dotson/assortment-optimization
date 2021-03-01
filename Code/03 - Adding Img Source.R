# Creating new table

library(tidyverse)
library(googledrive)
library(furrr)

data <- read_csv('Data/peanutbutter.csv')

# Looking up image pointers using function in file 2

data <- add_image_src(data, upc)

# We drop rows with missing values in the image db
clean_data <- drop_na(data, img_src)

# Then we resave the tbl

write_csv(clean_data, 'Data/peanutbutter.csv')

# Re-upload to google drive

