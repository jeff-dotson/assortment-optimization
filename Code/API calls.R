# This is the file used to get image URL's associated with the barcodes

library(httr)
library(jsonlite)
library(tidyverse)

key <- '3zpakfh2ispl447k2y69jsjtev5qpf'
test_key <- 'm6pf1gbkogbeys5ryld8udlywj4b14'
barcode_lookup <- "https://api.barcodelookup.com/v2/products?"

# Function used to grab the image URL's using barcodelookup.com
# Make sure to provide the barcode
barcode_image_lookup <- function(barcode) {
  img_src <- GET(url = paste(barcode_lookup, "barcode=", barcode, "&key=", key, sep = "")) %>% 
    content("parsed") %>% 
    flatten()
  unlist(img_src$products$images)
}


# Optional Google Custom Search API 
# a work In progress . . . 
g_key <- "AIzaSyA3J9hf9CZ-mOsQT0IzkAMwHpR9eXKAkfY"
search_engine <- "4aa35251bd443419c"
google_lookup <- "https://www.googleapis.com/customsearch/v1"
barcode_image_google <- function(barcode) {
  
}

# Function used to add additional columns to a lookup table

# Your Data here:

tbl <- #Obj name

tbl <- tbl %>% 
  mutate(across())
