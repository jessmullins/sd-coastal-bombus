# Run once to install everything sd_bombus.R needs.
pkgs <- c("readr", "dplyr", "stringr", "tidyr", "sf", "terra", "tigris",
          "ggplot2", "leaflet", "htmlwidgets", "openxlsx", "here")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)
message("All packages installed.")
