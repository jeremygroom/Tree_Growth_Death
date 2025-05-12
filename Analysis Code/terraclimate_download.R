library(rvest)

links <- read_html("http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/data/catalog.html")

# Where to save the output
out_dir <- "E:/terraclimate/all/" 

# Get all links
all_links <- links |>
  html_elements("a") |>
  html_attr("href") |>
  as.character()

# Subset to only links ending with .nc
nc_links <- all_links[endsWith(all_links, ".nc")]

# Get all file names, which are used to construct download links later
nc_names <- basename(nc_links)

# Download links
download_links <- paste0("http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/", nc_names)

# Create destination filenames
out_names <- file.path(out_dir, nc_names)

for (i in seq_along(out_names)) {
  print(paste("Downloading file:", nc_names[[i]]))
  download.file(
    url = download_links[[1]],
    destfile = out_names[[1]],
    quiet = TRUE
  )
}

