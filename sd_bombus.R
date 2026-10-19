# Coastal San Diego Bombus diet study: data prep
# 1) flag observations missing "Interaction->Visited flower of"
# 2) define "coastal" (ecoregion primary, distance-to-coast as sensitivity check)
# 3) remove urban/garden records using a county-wide urbanness mask built once from
#    NLCD impervious surface (mean % impervious in a circular window), plus a
#    protected-area check
#
# =============================================================================
# DATA DOWNLOADS: do this once before running
# =============================================================================
# Put everything in the data/ folder of this project. Unzip each download into
# its own subfolder; the script finds files by name, so folder names don't matter.
#
# 1. iNaturalist observations  ->  observations-XXXXXX.csv
#    https://www.inaturalist.org/observations/export
#    Log in as a manager of the San Diego County bees project first, so the
#    private_latitude / private_longitude columns get filled for observers who
#    trust the project.
#    - Query: taxon_id=52775 (Bombus), place = San Diego County
#    - Columns: id, url, observed_on, latitude, longitude, private_latitude,
#      private_longitude, positional_accuracy, coordinates_obscured, geoprivacy,
#      taxon_geoprivacy, captive_cultivated, quality_grade, scientific_name,
#      user_login
#    - Observation fields: tick "Interaction->Visited flower of" (required), plus
#      "nectar / pollen delivering plant", "nectar plant", "Associated species
#      with names lookup", "name of associated plant", "Plant association",
#      "feeding on"
#    Never commit this file: it contains coordinates shared in trust.
#
# 2. EPA Level IV ecoregions  ->  ca_eco_l4.shp (+ .dbf, .shx, .prj)
#    https://www.epa.gov/eco-research/ecoregion-download-files-state-region-9
#    Under California, download the Level IV shapefile (ca_eco_l4), not Level III.
#
# 3. NLCD fractional impervious surface  ->  a .tif with "FctImp" in the name
#    https://www.mrlc.gov/data
#    Choose Annual NLCD, Fractional Impervious Surface, the most recent year,
#    CONUS. The national file is several GB; the script crops it to San Diego
#    County on the first run and caches the crop in outputs/.
#    If you clip it yourself instead (e.g. with the MRLC Viewer), make sure
#    "FctImp" is still in the file name.
#
# 4. California Protected Areas Database  ->  a .shp with "Holdings" in the name
#    https://www.calands.org/cpad/
#    Download "California Protected Areas Database 2026a Release" (the shapefile
#    zip, not the ArcGIS REST services). It contains Holdings, Units and
#    SuperUnits; the script uses Holdings.
#
# 5. Vegetation, western San Diego County  ->  a .shp or .gpkg with "ds2964" or
#    "veg" in the name
#    https://catalog.data.gov/dataset/vegetation-western-san-diego-county-ds2964
#    (CDFW BIOS ds2964). The class column is CalVegType. If it downloads as a
#    geodatabase (.gdb folder), set veg_file by hand in the PATHS block below.
#
# Everything else (county boundary, coastline) downloads automatically.
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(terra)
library(tigris)   # county boundary + coastline
library(ggplot2)  # static map
library(leaflet)  # interactive map
library(htmlwidgets)
options(tigris_use_cache = TRUE)

crs_ca <- 3310    # California Albers, meters

# ---- PATHS ------------------------------------------------------------------------
# The project folder is found automatically when you open sd-coastal-bombus.Rproj
# (via the here package), so the script runs on anyone's computer without edits.
# Downloaded inputs go in data/; results are written to outputs/.
base_dir <- here::here()
data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "outputs")
if (!dir.exists(data_dir))
  stop("No data/ folder in ", base_dir, ". Open sd-coastal-bombus.Rproj first, ",
       "and see README.md for the files to download.")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Find exactly one file inside data/ (searches its subfolders)
find_file <- function(pattern, what) {
  hits <- list.files(data_dir, pattern = pattern, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0)
    stop("Missing input: ", what, "\nNo file matching '", pattern, "' in ", data_dir,
         "\nSee DATA DOWNLOADS at the top of this script for where to get it.",
         call. = FALSE)
  if (length(hits) > 1)
    stop("Several files match '", pattern, "' (", what, "):\n",
         paste(hits, collapse = "\n"), "\nKeep just one.", call. = FALSE)
  hits
}

inat_file <- find_file("^observations-.*\\.csv$",           "iNaturalist export (download 1)")
eco_file  <- find_file("^ca_eco_l4\\.shp$",                 "EPA Level IV ecoregions (download 2)")
nlcd_file <- find_file("FctImp.*\\.tif$",                   "NLCD impervious surface (download 3)")
cpad_file <- find_file("Holdings.*\\.shp$",                 "CPAD Holdings (download 4)")
veg_file  <- find_file("(ds2964|veg).*\\.(shp|gpkg)$",      "vegetation map (download 5)")
# If the vegetation download is a geodatabase (.gdb folder), set the path directly:
# veg_file <- file.path(data_dir, "FOLDER", "ds2964.gdb")

# Cached intermediate files (built once, reused on every later run)
sd_imp_file <- file.path(out_dir, "nlcd_imperv_sandiego.tif")
imp_win_file <- function(r) file.path(out_dir, paste0("imperv_mean_", r, "m_sandiego.tif"))

# ---- Settings --------------------------------------------------------------------
buffer_radii  <- c(100, 250, 500)   # m; windows for mean % impervious
imp_radius    <- 100                # m; which window the main urban rule uses
imp_thresh    <- 25                 # % impervious within imp_radius counted as "natural"
cpad_imp_cap  <- 50                 # looser cap for points inside CPAD protected areas
# Starting values only: replace them with the best rule from the calibration (4g)

# Vegetation map settings (fill in after the first run prints the column names
# and the list of classes; see section 4f)
veg_class_field     <- "CalVegType"    # column holding the vegetation class
non_natural_classes <- c(
  # clearly developed or planted
  "Urban",
  "Urban-Related Bare Soil",
  "Developed Water Feature",
  "Non-native/Ornamental Conifer",
  "Non-native/Ornamental Hardwood",
  "Non-native Ornamental Grasses",
  "Pastures and Crop Agriculture",
  # not foraging habitat on its own (edge rule still catches adjacent vegetation)
  "Water"
  # Judgment calls currently counted as NATURAL; move up here to exclude:
  # "Eucalyptus", "Giant Reed/Pampas Grass", "Tamarisk",
  # "Annual Grasses and Forbs", "Barren", "Beach Sand"
)
# Polygons with no class (NA) are always treated as not natural
veg_buffer_m        <- 30              # m; how close to natural vegetation still counts
max_pos_acc   <- 500                # m; drop points with worse positional accuracy
coast_dists   <- c(10, 15, 20)      # km; distance-to-coast sensitivity definitions
coast_pad_m   <- 1000               # m; extend coastal ecoregions seaward to catch shoreline records

# ---- County boundary (used throughout) -------------------------------------------
# cb = FALSE: full-resolution boundary (the generalized cb = TRUE version cuts off
# beaches, bluffs and points along the shoreline)
sd_county <- counties("CA", cb = FALSE) |>
  filter(NAME == "San Diego") |>
  st_transform(crs_ca)

# ---- 0. Read and clean iNat data ------------------------------------------------
obs <- read_csv(inat_file, show_col_types = FALSE)

# ---- Flower fields ---------------------------------------------------------------
# People log flower visits in several observation fields. They are listed here in
# priority order: if a record has values in more than one, the first one wins.
# Delete any you decide not to trust (e.g. "plant association" can mean resting on
# a plant, not visiting its flowers).
flower_fields <- c(
  "interaction->visited flower of",
  "nectar / pollen delivering plant",
  "nectar plant",
  "associated species with names lookup",
  "name of associated plant",
  "plant association",
  "feeding on"
)

# iNat exports these as "field:<name>"; match ignoring case and extra spaces
norm_name  <- function(x) str_squish(tolower(str_remove(x, "^field:")))
col_lookup <- setNames(names(obs), norm_name(names(obs)))

found_fields <- intersect(flower_fields, names(col_lookup))   # keeps priority order
missing_fields <- setdiff(flower_fields, found_fields)
if (!"interaction->visited flower of" %in% found_fields)
  stop("No 'Interaction->Visited flower of' column found; re-export with that field ticked")
if (length(missing_fields) > 0)
  message("Flower fields not in this export (fine if nobody used them): ",
          paste(missing_fields, collapse = "; "))

# One column per field, blanks turned into NA
flower_vals <- obs |>
  select(all_of(unname(col_lookup[found_fields]))) |>
  mutate(across(everything(), ~ na_if(str_trim(as.character(.x)), "")))
names(flower_vals) <- found_fields

first_filled <- apply(!is.na(flower_vals), 1, function(r) which(r)[1])

obs <- obs |>
  mutate(
    visited_flower  = do.call(coalesce, unname(as.list(flower_vals))),
    flower_source   = found_fields[first_filled],      # which field the plant came from
    n_flower_fields = rowSums(!is.na(flower_vals)),
    # TRUE when two fields name different plants; worth a look by hand
    flower_conflict = apply(flower_vals, 1, function(r) length(unique(tolower(na.omit(r)))) > 1),
    # Taxon-type fields can export as numeric taxon IDs instead of names
    flower_is_taxon_id = str_detect(coalesce(visited_flower, ""), "^\\d+$")
  )

print(obs |> count(flower_source, sort = TRUE))
cat(sum(obs$flower_conflict), "records name different plants in different fields\n")
cat(sum(obs$flower_is_taxon_id), "plant values are numeric taxon IDs, not names\n")

for (col in c("private_latitude", "private_longitude")) {
  if (!col %in% names(obs)) obs[[col]] <- NA_real_
}

is_obscured <- function(x) !is.na(x) & toupper(as.character(x)) == "TRUE"

obs <- obs |>
  mutate(
    has_flower     = !is.na(visited_flower),
    has_private    = !is.na(private_latitude) & !is.na(private_longitude),
    obscured       = is_obscured(coordinates_obscured),
    lat_use = if_else(has_private, private_latitude,  latitude),
    lon_use = if_else(has_private, private_longitude, longitude),
    loc_source = case_when(
      has_private & obscured ~ "recovered_private",
      has_private            ~ "private_unobscured",
      obscured               ~ "obscured_dropped",
      TRUE                   ~ "public"
    )
  ) |>
  filter(!is.na(lat_use), !is.na(lon_use))

# ---- 1. Which observations still need the plant field? -------------------------
needs_plant <- obs |>
  filter(!has_flower) |>
  select(id, url, observed_on, scientific_name, user_login, quality_grade)

write_csv(needs_plant, file.path(out_dir, "needs_visited_flower_field.csv"))
cat(nrow(needs_plant), "of", nrow(obs), "observations have no plant in any flower field\n")

# Plant is recorded, but only in a non-standard field. Optional to-do list if you
# want to copy these into "Interaction->Visited flower of" on iNat.
plant_other_field <- obs |>
  filter(has_flower, flower_source != "interaction->visited flower of") |>
  select(id, url, scientific_name, visited_flower, flower_source, user_login)
write_csv(plant_other_field, file.path(out_dir, "plant_only_in_other_field.csv"))

# Records whose fields disagree, for checking by hand
write_csv(obs |> filter(flower_conflict) |>
            select(id, url, scientific_name, user_login, flower_source, visited_flower),
          file.path(out_dir, "flower_field_conflicts.csv"))

# ---- 2. Location quality --------------------------------------------------------
loc_summary <- obs |>
  mutate(taxon_obscured = !is.na(taxon_geoprivacy) & taxon_geoprivacy != "open") |>
  count(loc_source, taxon_obscured)
print(loc_summary)

# Are dropped records biased toward particular species or observers?
print(obs |> count(scientific_name, loc_source) |>
        pivot_wider(names_from = loc_source, values_from = n, values_fill = 0))
print(obs |> filter(loc_source == "obscured_dropped") |>
        count(user_login, sort = TRUE))   # who to ask for trust

obs_good <- obs |>
  filter(loc_source != "obscured_dropped") |>
  filter(is.na(positional_accuracy) | positional_accuracy <= max_pos_acc)

obs_sf <- st_as_sf(obs_good, coords = c("lon_use", "lat_use"), crs = 4326, remove = FALSE) |>
  st_transform(crs_ca)

# ---- 3. Define "coastal" --------------------------------------------------------
# 3a. Ecoregion-based (primary)
eco <- st_read(eco_file, quiet = TRUE) |> st_transform(crs_ca) |> st_make_valid()
eco_sd <- st_intersection(eco, sd_county) |>
  st_collection_extract("POLYGON")   # drop stray line/point slivers from the clip

print(eco_sd |> st_drop_geometry() |> distinct(US_L4CODE, US_L4NAME) |> arrange(US_L4CODE))
plot(eco_sd["US_L4NAME"], key.pos = 4, key.width = lcm(6),
     main = "Level IV ecoregions, SD County")

# Coastal units, chosen a priori by Level IV code
coastal_codes <- c("85e", "85f")
stopifnot(all(coastal_codes %in% eco_sd$US_L4CODE))   # catches typos
print(eco_sd |> st_drop_geometry() |> filter(US_L4CODE %in% coastal_codes) |>
        distinct(US_L4CODE, US_L4NAME))               # confirm the names
coastal_poly <- eco_sd |> filter(US_L4CODE %in% coastal_codes) |> st_union()

# Ecoregion shorelines are generalized, so beach, bluff and headland records
# (La Jolla, Point Loma, etc.) can fall just outside the polygon. Extend the coastal
# polygon by coast_pad_m, but ONLY into areas no ecoregion covers (the shoreline gap
# and the ocean). The original coastal polygon is always kept whole.
coastal_core <- coastal_poly
any_eco      <- eco |>
  st_filter(st_buffer(coastal_core, 2 * coast_pad_m)) |>
  st_union()
shore_gap    <- st_difference(st_buffer(coastal_core, coast_pad_m), any_eco)
coastal_poly <- st_union(coastal_core, shore_gap)

n_core <- sum(lengths(st_intersects(obs_sf, coastal_core)) > 0)
n_ext  <- sum(lengths(st_intersects(obs_sf, coastal_poly)) > 0)
cat(n_core, "observations inside the coastal ecoregions;",
    n_ext - n_core, "more recovered along the shoreline\n")
if (n_ext < n_core) stop("Coastal extension lost observations; check coastal_poly")

obs_sf$in_coastal_eco <- lengths(st_intersects(obs_sf, coastal_poly)) > 0

# 3b. Distance to coast (sensitivity check / continuous covariate)
coast <- coastline(year = 2023) |>
  st_transform(crs_ca) |>
  st_crop(st_buffer(sd_county, 20000))
obs_sf$dist_coast_km <- as.numeric(st_distance(obs_sf, st_union(coast))) / 1000

# ---- 4. Urbanness mask (built once, then instant lookups) -----------------------
# 4a. Crop the national NLCD raster to the county (first run only)
if (!file.exists(sd_imp_file)) {
  message("Cropping NLCD to San Diego County (first run only)...")
  imp_full <- rast(nlcd_file)
  sd_box   <- vect(st_transform(st_buffer(sd_county, 2000), crs(imp_full)))
  writeRaster(crop(imp_full, sd_box), sd_imp_file, overwrite = TRUE)
}
imp <- rast(sd_imp_file)

# 4b. County-wide mean % impervious within a circular window, one raster per radius
#     (first run only; NLCD class breaks for reference: developed open space <20%,
#      low 20-49%, medium 50-79%, high >=80%)
for (r in buffer_radii) {
  if (!file.exists(imp_win_file(r))) {
    message("Building ", r, " m impervious window raster (first run only)...")
    w <- focalMat(imp, r, type = "circle")   # weights sum to 1
    writeRaster(focal(imp, w, fun = "sum", na.rm = TRUE), imp_win_file(r), overwrite = TRUE)
  }
}

# 4c. Point lookups
pts <- vect(st_transform(obs_sf, crs(imp)))
for (r in buffer_radii) {
  obs_sf[[paste0("imp_", r, "m")]] <- terra::extract(rast(imp_win_file(r)), pts)[, 2]
}
obs_sf$imp_main <- obs_sf[[paste0("imp_", imp_radius, "m")]]   # the window the rule uses

# 4d. Save a binary wildlands mask for mapping / checking in QGIS
wild_mask <- rast(imp_win_file(imp_radius)) < imp_thresh
writeRaster(wild_mask, file.path(out_dir, "wildlands_mask_sandiego.tif"), overwrite = TRUE)

# 4e. Inside protected open space? (positive "natural habitat" check)
cpad <- st_read(cpad_file, quiet = TRUE) |>
  st_transform(crs_ca) |>
  st_make_valid() |>
  st_filter(sd_county)
obs_sf$in_cpad <- lengths(st_intersects(obs_sf, cpad)) > 0

# ---- 4f. Natural vegetation map (positive habitat layer) -------------------------
# Impervious windows misjudge narrow canyons because the circle catches the houses
# on the rims. A mapped vegetation polygon follows the canyon itself.
veg <- st_read(veg_file, quiet = TRUE) |> st_transform(crs_ca) |> st_make_valid()

if (!veg_class_field %in% names(veg)) {
  print(names(veg))
  stop("Set veg_class_field (Settings) to the column above that holds the vegetation class")
}
print(veg |> st_drop_geometry() |> as_tibble() |>
        count(.data[[veg_class_field]], sort = TRUE), n = Inf)
if (any(non_natural_classes == "FILL IN"))
  stop("Use the class list above to fill in non_natural_classes (Settings)")
stopifnot(all(non_natural_classes %in% veg[[veg_class_field]]))   # catches typos

veg_natural <- veg |>
  filter(!is.na(.data[[veg_class_field]]),
         !.data[[veg_class_field]] %in% non_natural_classes) |>
  st_filter(st_buffer(coastal_poly, 5000))     # only the coastal area is needed

obs_sf$in_veg_poly <- lengths(st_intersects(obs_sf, veg_natural)) > 0
obs_sf$near_veg    <- lengths(st_is_within_distance(obs_sf, veg_natural,
                                                    dist = veg_buffer_m)) > 0
obs_sf$veg_class   <- veg[[veg_class_field]][st_nearest_feature(obs_sf, veg)]

# ---- 4g. Calibrate the urban rule against hand-labeled points --------------------
# Step 1 (first run): writes a labeling sheet of coastal points sampled across the
#   whole range of urbanness. Impervious values are left OFF the sheet so your
#   labels stay blind to them (label from the satellite links, not the R map).
# Step 2: open each `satellite` link and fill in `label` with one of:
#     garden  = yard, lawn, street planting, landscaped campus or business park
#     canyon  = natural vegetation in an urban canyon or small fragment
#     reserve = natural vegetation in a large open space
#     unsure  = can't tell (ignored)
#   Save it as a CSV with the same name. The script never overwrites this file.
# Step 3: re-run. Once at least 20 points are labeled, every radius x threshold
#   combination is scored and the results are written to calibration_scores.csv.
calib_file <- file.path(out_dir, "PRIVATE_calibration_sheet.csv")

if (!file.exists(calib_file)) {
  set.seed(42)
  calib_sheet <- obs_sf |>
    st_drop_geometry() |>
    filter(in_coastal_eco, !is.na(imp_100m)) |>
    mutate(imp_bin = cut(imp_100m, c(-Inf, 5, 10, 20, 35, 50, Inf))) |>
    group_by(imp_bin) |>
    slice_sample(n = 15) |>          # up to 15 points per urbanness bin (~90 total)
    ungroup() |>
    slice_sample(prop = 1) |>        # shuffle so the bins aren't in order
    transmute(id, url, scientific_name, observed_on,
              satellite = paste0("https://www.google.com/maps?q=", lat_use, ",",
                                 lon_use, "&t=k&z=19"),
              label = "", notes = "")
  write_csv(calib_sheet, calib_file)
  message("Calibration sheet written: ", calib_file,
          "\nFill in the `label` column, save, and re-run.")
}

calib_labels <- read_csv(calib_file, col_types = cols(.default = col_character())) |>
  mutate(id = as.numeric(id), label = tolower(str_trim(label))) |>
  filter(label %in% c("garden", "canyon", "reserve"))

if (nrow(calib_labels) >= 20) {
  print(count(calib_labels, label))
  
  calib <- calib_labels |>
    select(id, label) |>
    inner_join(st_drop_geometry(obs_sf) |>
                 select(id, in_cpad, in_veg_poly, near_veg, starts_with("imp_")),
               by = "id") |>
    mutate(truth_natural = label != "garden")
  
  # How well does a TRUE/FALSE "keep" vector match your labels?
  score_pass <- function(pass) tibble(
    gardens_removed   = mean(!pass[calib$label == "garden"]),
    canyons_kept      = mean(pass[calib$label == "canyon"]),
    reserves_kept     = mean(pass[calib$label == "reserve"]),
    balanced_accuracy = (mean(!pass[!calib$truth_natural]) +
                           mean(pass[calib$truth_natural])) / 2
  )
  
  window_pass <- function(radius, threshold, cpad_rule) {
    imp <- calib[[paste0("imp_", radius, "m")]]
    !is.na(imp) & (imp < threshold | (cpad_rule & calib$in_cpad & imp < cpad_imp_cap))
  }
  
  score_rule <- function(radius, threshold, cpad_rule)
    bind_cols(tibble(radius, threshold, cpad_rule),
              score_pass(window_pass(radius, threshold, cpad_rule)))
  
  # Vegetation-based rules, scored the same way
  current_window <- window_pass(imp_radius, imp_thresh, TRUE)
  veg_scores <- bind_rows(
    bind_cols(rule = "inside natural vegetation",            score_pass(calib$in_veg_poly)),
    bind_cols(rule = paste0("within ", veg_buffer_m, " m of natural vegetation"),
              score_pass(calib$near_veg)),
    bind_cols(rule = "near natural vegetation OR current window rule",
              score_pass(calib$near_veg | current_window))
  )
  cat("\nVegetation-based rules:\n")
  print(veg_scores)
  write_csv(veg_scores, file.path(out_dir, "calibration_scores_vegetation.csv"))
  
  rule_grid <- expand_grid(radius = buffer_radii,
                           threshold = seq(5, 50, by = 5),
                           cpad_rule = c(FALSE, TRUE))
  calib_scores <- bind_rows(Map(score_rule, rule_grid$radius,
                                rule_grid$threshold, rule_grid$cpad_rule)) |>
    arrange(desc(balanced_accuracy))
  
  cat("\nBest rules (balanced accuracy = average of gardens removed and natural kept):\n")
  print(head(calib_scores, 10))
  cat("\nCurrent rule:\n")
  print(filter(calib_scores, radius == imp_radius, threshold == imp_thresh, cpad_rule))
  
  write_csv(calib_scores, file.path(out_dir, "calibration_scores.csv"))
  
  p_calib <- ggplot(calib_scores,
                    aes(threshold, balanced_accuracy, colour = factor(radius),
                        linetype = cpad_rule)) +
    geom_line(linewidth = 0.8) +
    labs(x = "Impervious threshold (%)", y = "Balanced accuracy",
         colour = "Window radius (m)", linetype = "CPAD allowance") +
    theme_minimal(base_size = 11)
  print(p_calib)
  ggsave(file.path(out_dir, "calibration_curve.png"), p_calib, width = 7, height = 4.5, dpi = 300)
} else {
  message(nrow(calib_labels), " points labeled so far; calibration runs at 20+.")
}

# ---- 5. Final sets ---------------------------------------------------------------
obs_sf <- obs_sf |>
  mutate(
    window_ok     = !is.na(imp_main) &
      (imp_main < imp_thresh | (in_cpad & imp_main < cpad_imp_cap)),
    # Why each point counts as natural (or not); first match wins
    habitat_basis = case_when(
      in_veg_poly ~ "inside natural vegetation",
      near_veg    ~ "edge of natural vegetation",
      window_ok   ~ "low impervious surface",
      TRUE        ~ "urban"
    ),
    natural      = habitat_basis != "urban",
    keep_primary = in_coastal_eco & natural & has_flower
  )

print(table(coastal = obs_sf$in_coastal_eco, natural = obs_sf$natural))
print(obs_sf |> st_drop_geometry() |> filter(in_coastal_eco) |> count(habitat_basis))

# Sensitivity: drop the "edge of natural vegetation" points (possible rim gardens)
sens_edge <- obs_sf |>
  st_drop_geometry() |>
  filter(in_coastal_eco, natural, has_flower,
         habitat_basis != "edge of natural vegetation") |>
  count(scientific_name, visited_flower, name = "n")

# Sensitivity: does the diet change with a distance-based "coastal" definition?
sens <- lapply(coast_dists, function(d) {
  obs_sf |>
    st_drop_geometry() |>
    filter(dist_coast_km <= d, natural, has_flower) |>
    count(scientific_name, visited_flower, name = "n") |>
    mutate(def = paste0(d, "km"))
}) |> bind_rows()

# Sensitivity: does the diet change with the impervious threshold?
sens_imp <- lapply(c(10, 15, 20, 25, 30, 40), function(t) {
  obs_sf |>
    st_drop_geometry() |>
    filter(in_coastal_eco, has_flower, !is.na(imp_main),
           imp_main < t | (in_cpad & imp_main < cpad_imp_cap)) |>
    count(scientific_name, visited_flower, name = "n") |>
    mutate(imp_threshold = t)
}) |> bind_rows()

final <- obs_sf |> filter(keep_primary)
cat(nrow(final), "observations in the primary coastal + natural + plant-recorded set\n")

# ---- 6. Write outputs ------------------------------------------------------------
# Internal files: contain true coordinates. Keep these private.
write_csv(st_drop_geometry(final),
          file.path(out_dir, "PRIVATE_coastal_natural_bombus_with_plants.csv"))
st_write(obs_sf, file.path(out_dir, "PRIVATE_bombus_sd_flagged.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)   # inspect in QGIS
write_csv(sens,     file.path(out_dir, "sensitivity_distance_definitions.csv"))
write_csv(sens_imp, file.path(out_dir, "sensitivity_impervious_thresholds.csv"))
write_csv(sens_edge, file.path(out_dir, "sensitivity_without_veg_edge_points.csv"))

# Shareable version: no true coordinates for publicly obscured observations
final_public <- final |>
  st_drop_geometry() |>
  select(-private_latitude, -private_longitude) |>
  mutate(
    lat_use = if_else(loc_source == "recovered_private", NA_real_, lat_use),
    lon_use = if_else(loc_source == "recovered_private", NA_real_, lon_use)
  )
write_csv(final_public, file.path(out_dir, "coastal_natural_bombus_with_plants_SHAREABLE.csv"))

# ---- 7. Maps of what survived the filters ---------------------------------------
status_levels <- c("Kept", "Kept habitat, needs plant ID", "Dropped: urban", "Dropped: not coastal")
status_cols   <- c("#d95f02", "#7570b3", "grey25", "grey70")
names(status_cols) <- status_levels

obs_sf <- obs_sf |>
  mutate(status = factor(case_when(
    keep_primary                            ~ "Kept",
    in_coastal_eco & natural & !has_flower  ~ "Kept habitat, needs plant ID",
    in_coastal_eco & !natural               ~ "Dropped: urban",
    TRUE                                    ~ "Dropped: not coastal"
  ), levels = status_levels))

print(obs_sf |> st_drop_geometry() |> count(status))

# Wildland mask at 150 m for display only (30 m is too many pixels to draw)
wild_disp <- aggregate(wild_mask * 1, fact = 5, fun = "mean", na.rm = TRUE) >= 0.5
wild_disp <- ifel(wild_disp, 1, NA)

# 7a. Static map (PNG for slides / supplement)
zoom_bb  <- st_bbox(st_buffer(coastal_poly, 3000))
veg_zoom <- st_crop(veg_natural, zoom_bb)

p <- ggplot() +
  geom_sf(data = sd_county, fill = "grey97", colour = "grey60") +
  geom_sf(data = veg_zoom, fill = "#cfe3c3", colour = NA) +
  geom_sf(data = coastal_poly, fill = NA, colour = "#2b6c8f", linewidth = 0.6) +
  geom_sf(data = arrange(obs_sf, desc(status)),          # kept points drawn on top
          aes(colour = status, size = status), alpha = 0.85) +
  scale_colour_manual(values = status_cols, name = NULL) +
  scale_size_manual(values = c(1.8, 1.4, 0.8, 0.8), name = NULL) +
  coord_sf(xlim = zoom_bb[c("xmin", "xmax")], ylim = zoom_bb[c("ymin", "ymax")]) +
  labs(title = "Bombus observations after habitat and impervious-surface filters",
       subtitle = "Green = mapped natural vegetation; blue outline = coastal ecoregions") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", axis.title = element_blank())

print(p)
ggsave(file.path(out_dir, "PRIVATE_map_filtered_points.png"), p, width = 8, height = 9, dpi = 300)

# 7b. Interactive map over satellite imagery (click a point for details + iNat link)
#     Use this to spot-check: are "Kept" points really in canyons, not backyards?
obs_ll <- st_transform(obs_sf, 4326)
pal <- colorFactor(status_cols, levels = status_levels)

build_map <- function(pts, veg_tolerance = 5) {
  m <- leaflet() |>
    addProviderTiles("Esri.WorldImagery", group = "Satellite") |>
    addProviderTiles("CartoDB.Positron",  group = "Light") |>
    addPolygons(data = st_transform(st_simplify(veg_natural, dTolerance = veg_tolerance), 4326),
                fillColor = "#7fbf7b", fillOpacity = 0.35, stroke = FALSE,
                group = "Natural vegetation") |>
    addRasterImage(wild_disp, colors = "#fdb863", opacity = 0.35,
                   project = TRUE, method = "ngb", group = "Low-impervious mask") |>
    addPolylines(data = st_transform(st_cast(coastal_poly, "MULTILINESTRING"), 4326),
                 color = "#2b6c8f", weight = 2, group = "Coastal ecoregions")
  
  for (s in status_levels) {
    d <- filter(pts, status == s)
    if (nrow(d) == 0) next
    m <- m |> addCircleMarkers(
      data = d, radius = 5, stroke = TRUE, weight = 1, color = "white",
      fillColor = ~pal(status), fillOpacity = 0.9, group = s,
      popup = ~paste0("<b>", scientific_name, "</b><br>", observed_on,
                      "<br>Plant: ", coalesce(visited_flower, "none recorded"),
                      ifelse(is.na(flower_source), "", paste0(" (", flower_source, ")")),
                      "<br>Impervious: ", round(imp_100m, 1), "% (100 m), ",
                      round(imp_250m, 1), "% (250 m)",
                      "<br>In CPAD: ", in_cpad,
                      "<br>Habitat: ", habitat_basis, " (", veg_class, ")",
                      "<br><a href='", url, "' target='_blank'>Open on iNaturalist</a>"))
  }
  
  m |>
    addLayersControl(baseGroups = c("Satellite", "Light"),
                     overlayGroups = c(status_levels, "Natural vegetation",
                                       "Low-impervious mask", "Coastal ecoregions"),
                     options = layersControlOptions(collapsed = FALSE)) |>
    hideGroup(c("Dropped: not coastal", "Low-impervious mask")) |>
    addLegend("bottomright", pal = pal, values = status_levels, title = "Filter result")
}

m <- build_map(obs_ll)
m   # shows in the RStudio Viewer
saveWidget(m, file.path(out_dir, "PRIVATE_map_filtered_points.html"), selfcontained = FALSE)

# ---- 8. Shareable outputs --------------------------------------------------------
# 8a. A map you can send to someone: ONE html file (opens in any browser), with the
#     observations whose true location came from trust (recovered_private) left off.
pandoc_ok <- nzchar(Sys.getenv("RSTUDIO_PANDOC")) || nzchar(Sys.which("pandoc"))
obs_ll_share <- filter(obs_ll, loc_source != "recovered_private")
cat(nrow(obs_ll) - nrow(obs_ll_share),
    "trusted-location observations left off the shareable map\n")

m_share <- build_map(obs_ll_share, veg_tolerance = 15)   # simpler polygons = smaller file
share_html <- file.path(out_dir, "SHAREABLE_map_filtered_points.html")
saveWidget(m_share, share_html, selfcontained = pandoc_ok,
           title = "Coastal San Diego Bombus: filter results")
if (!pandoc_ok) message("Pandoc not found, so the map was saved with a companion folder; ",
                        "run this from RStudio to get a single file.")
cat("Shareable map:", share_html, "(",
    round(file.size(share_html) / 1e6, 1), "MB )\n")

# Static version of the same map
p_share <- p
p_share$layers[[length(p_share$layers)]]$data <-
  arrange(filter(obs_sf, loc_source != "recovered_private"), desc(status))
ggsave(file.path(out_dir, "SHAREABLE_map_filtered_points.png"), p_share,
       width = 8, height = 9, dpi = 300)

# 8b. Spreadsheet of observations that still need a plant ID, 2019 onward.
#     No coordinates, so it is safe to share with volunteers. Two tabs:
#     priority = coastal + natural habitat (these go straight into the study),
#     all      = every San Diego County record still missing a plant.
library(openxlsx)
start_date <- as.Date("2019-01-01")

plant_todo_cols <- function(df) df |>
  mutate(observed_on = as.Date(observed_on)) |>
  filter(!has_flower, observed_on >= start_date) |>
  arrange(desc(observed_on)) |>
  transmute(observed_on, scientific_name, observer = user_login, quality_grade,
            url, plant_id = "", notes = "")

todo_priority <- obs_sf |> st_drop_geometry() |>
  filter(status == "Kept habitat, needs plant ID") |>
  plant_todo_cols()
todo_all <- obs |> plant_todo_cols()

wb <- createWorkbook()
for (tab in list(list("Coastal natural (priority)", todo_priority),
                 list("All San Diego County", todo_all))) {
  d <- tab[[2]]
  class(d$url) <- "hyperlink"           # clickable links in Excel
  addWorksheet(wb, tab[[1]])
  writeData(wb, tab[[1]], d, headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb, tab[[1]], firstRow = TRUE)
  addFilter(wb, tab[[1]], rows = 1, cols = seq_along(d))
  setColWidths(wb, tab[[1]], cols = seq_along(d), widths = c(12, 24, 18, 12, 45, 28, 30))
}
todo_file <- file.path(out_dir, "needs_plant_ID_2019_onward.xlsx")
saveWorkbook(wb, todo_file, overwrite = TRUE)
cat(nrow(todo_priority), "priority and", nrow(todo_all),
    "total observations need a plant ID (2019 onward):", todo_file, "\n")