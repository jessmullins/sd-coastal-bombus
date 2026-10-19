# data/

Put the downloaded input files here. Subfolders are fine (for example, unzip each
download into its own folder). Nothing in this folder except this README is
committed to GitHub.

The script finds each file by name, so it needs exactly one match for each:

| Input | File the script looks for | Where to get it |
|---|---|---|
| iNaturalist export | `observations-*.csv` | iNaturalist export page (settings in the main README) |
| EPA Level IV ecoregions | `ca_eco_l4.shp` (+ .dbf, .shx, .prj) | [EPA Region 9 ecoregion downloads](https://www.epa.gov/eco-research/ecoregion-download-files-state-region-9) |
| NLCD fractional impervious surface | a `.tif` with `FctImp` in the name | [MRLC Annual NLCD](https://www.mrlc.gov/data) |
| California Protected Areas Database | a `.shp` with `Holdings` in the name | [CPAD at calands.org](https://www.calands.org/cpad/) (the "2026a Release" download) |
| Vegetation, western San Diego County | a `.shp` or `.gpkg` with `ds2964` or `veg` in the name | [CDFW ds2964 on data.gov](https://catalog.data.gov/dataset/vegetation-western-san-diego-county-ds2964) |

If the vegetation map downloads as a geodatabase (`.gdb` folder), set `veg_file`
by hand in the PATHS block of `sd_bombus.R`.
