### failed script to apply .points transoforms to dems

input_raster <- "data_raw/records/CP.tif"
points_file <- "data_raw/2022Nov_RR_alignmentpoints/202211_cookspath_aligned.tif.points"
temp_output <- "temp_with_gcps.tif"
final_output <- "data_raw/records/CPt.tif"

lines <- readLines(points_file)

gcps <- lapply(lines, function(line) {
  # Trim and skip empty lines
  line <- trimws(line)
  if (nchar(line) == 0) return(NULL)
  
  # Try to parse numeric values
  vals <- suppressWarnings(as.numeric(strsplit(line, ",")[[1]]))
  
  # Check that there are at least 5 values and the 5th is 1 (enabled GCP)
  if (length(vals) >= 5 && !any(is.na(vals)) && vals[5] == 1) {
    return(c("-gcp", vals[1], vals[2], vals[3], vals[4]))
  } else {
    return(NULL)
  }
})

# Flatten into a single character vector
gcp_args <- unlist(gcps)

# Check output
print(gcp_args)

cmd_translate <- c(
  "gdal_translate",
  "-of", "GTiff",
  gcp_args,
  shQuote(input_raster),
  shQuote(temp_output)
)

system(paste(cmd_translate, collapse = " "))

# Warp to georeference
cmd_warp <- paste(
  "gdalwarp",
  "-tps",                         # Thin Plate Spline
  "-r bilinear",                  # Resampling method
  "-overwrite",                   # Overwrite output if exists
  shQuote(temp_output),
  shQuote(final_output)
)

system(cmd_warp)
