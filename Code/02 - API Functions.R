# This is the file used to get image URL's associated with the barcodes
# NOTE: you will need an API key.

library(httr)
library(jsonlite)
library(tidyverse)

barcode_lookup <- "https://api.barcodelookup.com/v2/products?"

# Function used to grab the image URL's using barcodelookup.com
# Make sure to provide the barcode
barcode_image_lookup <- function(barcode) {
  img_src <-
    GET(url = paste(barcode_lookup, "barcode=", barcode, "&key=", key, sep = ""))
  if (img_src$status_code == 200) { 
## Only compute if the product is found in barcodelookup database.
    img <- content(img_src, "parsed") %>%
      flatten()
    unlist(img$products$images)
  }
  else {
    img <- NA
    img
  }
}

# Function used to add additional columns to a look up table

add_image_src <- function(tbl, barcode) {
  tbl %>%
    group_by({{ barcode }}) %>% 
    count({{ barcode }}) %>% 
    arrange(desc(n)) %>% 
    mutate(img_src = NA) %>% 
    rowwise() %>%
    mutate(img_src = barcode_image_lookup({{ barcode }}))
}
