# Coastal San Diego bumble bee diet

Data preparation for a study of the diet of bumble bees (*Bombus*) in coastal San
Diego County, using iNaturalist observations with flower-visit records.

`sd_bombus.R` takes a raw iNaturalist export and:

1. combines the plant recorded in any of several flower observation fields, and
   lists the observations that still need a plant ID;
2. recovers true coordinates for obscured observations where observers have
   trusted the project, and drops records whose location is still obscured;
3. keeps observations in the coastal Level IV ecoregions (85e and 85f), extended a
   little seaward to catch beach and bluff records;
4. removes urban and garden records, keeping natural habitat, including narrow
   urban canyons;
5. writes the filtered dataset, sensitivity tables, maps, and a plant-ID to-do
   spreadsheet.

## Quick start

1. Clone this repo and open `sd-coastal-bombus.Rproj` in RStudio. Paths are set
   relative to the project, so nothing in the script needs editing.
2. Run `install_packages.R` once.
3. Download the inputs into `data/` (see [data/README.md](data/README.md)). They
   aren't in the repo, because they are large and the iNaturalist export can hold
   private coordinates.
4. Run `sd_bombus.R`. The first run takes a few minutes, while it crops NLCD to the
   county and builds the impervious-surface rasters. Those are cached in `outputs/`,
   so later runs take seconds.

## iNaturalist export settings

Export from the [iNaturalist export page](https://www.inaturalist.org/observations/export)
**while logged in as a manager of the San Diego County bees project**. That way, the
`private_latitude` / `private_longitude` columns are filled in for observers who
trust the project.

- **Query:** `taxon_id=52775` (*Bombus*), place = San Diego County
- **Columns:** `id`, `url`, `observed_on`, `latitude`, `longitude`,
  `private_latitude`, `private_longitude`, `positional_accuracy`,
  `coordinates_obscured`, `geoprivacy`, `taxon_geoprivacy`, `captive_cultivated`,
  `quality_grade`, `scientific_name`, `user_login`
- **Observation fields** (in the order the script prefers them):
  - Interaction->Visited flower of (**required**)
  - nectar / pollen delivering plant
  - nectar plant
  - Associated species with names lookup
  - name of associated plant
  - Plant association
  - feeding on

## How the filters work

| Step | Rule | Setting(s) in `sd_bombus.R` |
|---|---|---|
| Location | Drop records still obscured; drop positional accuracy worse than 500 m | `max_pos_acc` |
| Coastal | Inside Level IV ecoregions 85e or 85f, extended 1 km seaward into areas no ecoregion covers | `coastal_codes`, `coast_pad_m` |
| Natural habitat | Kept if any of these hold (first match is recorded in `habitat_basis`): inside mapped natural vegetation; within 30 m of it; or mean impervious surface within 100 m under 25% (under 50% inside a CPAD protected area). Otherwise dropped as urban. | `veg_class_field`, `non_natural_classes`, `veg_buffer_m`, `imp_radius`, `imp_thresh`, `cpad_imp_cap` |
| Plant | At least one flower field filled in | `flower_fields` |

Why the vegetation map: circular impervious-surface windows catch the houses on
the rims of narrow urban canyons, so canyon records were being dropped even
though the bees were in coastal sage scrub. Mapped vegetation polygons follow the
canyon itself. The impervious rule is kept as a fallback for gaps in the
vegetation map.

Vegetation classes counted as not natural: Urban, Urban-Related Bare Soil,
Developed Water Feature, the three Non-native/Ornamental classes, Pastures and
Crop Agriculture, Water, and unclassified polygons. Judgment calls currently
counted as natural: Eucalyptus, Giant Reed/Pampas Grass, Tamarisk, Annual Grasses
and Forbs, Barren, Beach Sand.

### Calibrating the urban rule

Section 4g checks the rules against hand-labeled points:

1. The first run writes `outputs/PRIVATE_calibration_sheet.csv`, with about 90
   coastal points sampled across the urbanness range. Each has a satellite link.
2. Open each link and fill in `label` as `garden`, `canyon`, `reserve`, or
   `unsure`. Label from the links rather than the R map, since the map popups show
   impervious values. The script never overwrites this file.
3. Re-run. With 20 or more labels, every radius and threshold combination, plus
   the vegetation-based rules, is scored for gardens removed, canyons kept, and
   reserves kept (`calibration_scores*.csv`, `calibration_curve.png`).

## Outputs (`outputs/`, not committed)

| File | What it is | Shareable? |
|---|---|---|
| `PRIVATE_coastal_natural_bombus_with_plants.csv` | Final dataset, true coordinates | No |
| `coastal_natural_bombus_with_plants_SHAREABLE.csv` | Final dataset, trusted coordinates blanked | Yes |
| `PRIVATE_bombus_sd_flagged.gpkg` | Every record with all flags, for QGIS | No |
| `PRIVATE_map_filtered_points.html` / `.png` | Maps of what passed each filter | No |
| `SHAREABLE_map_filtered_points.html` / `.png` | Same maps without trusted-location records | Yes |
| `needs_plant_ID_2019_onward.xlsx` | Plant-ID to-do list with iNat links, no coordinates | Yes |
| `needs_visited_flower_field.csv` | Every record with no plant in any field | Yes |
| `plant_only_in_other_field.csv` | Plant logged only in a non-standard field | Yes |
| `flower_field_conflicts.csv` | Records where the flower fields disagree | Yes |
| `sensitivity_*.csv` | Diet tallies under alternative coastal, impervious, and edge definitions | Yes |

## Privacy

Observers who trust the project share their true coordinates with project managers
only. Never commit or publicly share anything prefixed `PRIVATE_`, or the raw
iNaturalist export. `.gitignore` excludes `data/` and `outputs/` for this reason.
Before every commit, check `git status` to make sure no data files are staged.
