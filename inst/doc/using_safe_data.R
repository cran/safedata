## ----knithook, echo=FALSE-----------------------------------------------------
library(knitr)
# Knit hook to truncate output.
hook_output <- knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
lines <- options$output.lines
if (is.null(lines)) {
    return(hook_output(x, options))  # pass to default hook
}
x <- unlist(strsplit(x, "\n"))
more <- "..."
if (length(lines)==1) {        # first n lines
	if (length(x) > lines) {
		# truncate the output, but add ....
		x <- c(head(x, lines), more)
	}
	} else {
		x <- c(more, x[lines], more)
	}
	# paste these lines together
	x <- paste(c(x, ""), collapse = "\n")
	hook_output(x, options)
})
# Comments as #> not ##
knitr::opts_chunk$set(
  comment = '#>'
)

## ----load_remote_data, echo=FALSE---------------------------------------------
# All the remote data used in building this vignette is downloaded
# and cached in R/sysdata.rda, so that vignette code checking and 
# building does not rely on remote resources. The vignette_objects.R
# script _does_ need to be run before fresh releases to update the file.

use_remote <- as.logical(Sys.getenv('BUILD_VIGNETTE_USING_REMOTE', FALSE))
if (! use_remote){
  soil_datasets <- safedata:::vignette_soil_search
  ants <- safedata:::vignette_ants_search
  all_taxa <- safedata:::vignette_taxon_coverage
}

## ----install_cran, eval=FALSE-------------------------------------------------
#  install.packages("safedata")

## ----install_git, eval=FALSE--------------------------------------------------
#  devtools::install_github("ImperialCollegeLondon/safedata")

## ----load_lib, collapse=TRUE--------------------------------------------------
library(safedata)

## ----create_new_dir, eval=FALSE-----------------------------------------------
#  set_safe_dir('my_safe_directory', create=TRUE)
#  ## Safe data directory created

## ----show_safe_dir, eval=FALSE------------------------------------------------
#  getOption('safedata.dir')
#  ## [1] "my_safe_directory"

## ----set_safe_dir, eval=FALSE-------------------------------------------------
#  set_safe_dir('my_safe_directory')
#  ## Checking for updates
#  ##  - Index up to date
#  ##  - Gazetteer up to date
#  ##  - Location aliases up to date
#  ## Validating directory

## ----load_example_dir---------------------------------------------------------
set_example_safe_dir()

## ----search_text, eval=use_remote---------------------------------------------
#  soil_datasets <- search_text('soil')

## ----search_text_out, collapse=TRUE-------------------------------------------
print(soil_datasets)

## ----search_taxa, eval=use_remote---------------------------------------------
#  ants <- search_taxa('Formicidae')

## ----search_taxa_out, collapse=TRUE-------------------------------------------
print(ants)

## ----search_taxa_id, eval=use_remote------------------------------------------
#  ants <- search_taxa(gbif_id=4342)

## ----search_spatial, eval=use_remote------------------------------------------
#  # Datasets that include sampling within experimental block A
#  within_a <- search_spatial(location='BL_A')
#  # Datasets that sampled within 2 km of the Maliau Basin Field Study Centre
#  near_maliau <- search_spatial(wkt='POINT(116.97394 4.73481)', distance=2000)

## ----combining_searches, collapse=TRUE, eval=use_remote-----------------------
#  # Three searches
#  fish <- search_taxa("Actinopterygii")
#  odonates <- search_taxa("Odonata")
#  ewers <- search_authors("Ewers")
#  # Logical combinations
#  aquatic <- fish | odonates
#  aquatic_ewers <- aquatic & ewers
#  all_in_one <- (fish | odonates) & ewers

## ----restricting_searches, collapse=TRUE, eval=use_remote---------------------
#  fish <- search_taxa("Actinopterygii")
#  ewers <- search_authors("Ewers", ids=fish)

## ----validate_recs, echo=TRUE-------------------------------------------------
recs <- validate_record_ids(c('https://doi.org/10.5281/zenodo.3247631',
                              '10.5281/zenodo.3266827',
                              'https://zenodo.org/record/3266821'))
print(recs)

## ----show_concepts------------------------------------------------------------
show_concepts(recs)

## ----show_record--------------------------------------------------------------
show_record(1400562)

## ----show_worksheet-----------------------------------------------------------
show_worksheet(1400562, 'EnvironVariables')

## ----show_worksheet_long, output.lines=15-------------------------------------
show_worksheet(1400562, 'EnvironVariables', extended_fields=TRUE)

## ----download_safe, eval=FALSE------------------------------------------------
#  download_safe_files(within_a)
#  ## 26 files requested from 26 records
#  ##  - 0 local (0 bytes)
#  ##  - 4 embargoed or restricted (2.2 Mb)
#  ##  - 22 to download (43.6 Mb)
#  ##
#  ## 1: Yes
#  ## 2: No
#  ##
#  ## Selection:

## ----download_safe_external, eval=FALSE---------------------------------------
#  download_safe_files(3697804, xlsx_only=FALSE)
#  # 2 files requested from 1 records
#  #  - 1 local (11.3 Kb)
#  #  - 0 embargoed or restricted (0 bytes)
#  #  - 1 to download (74.1 Kb)
#  #
#  # 1: Yes
#  # 2: No
#  #
#  # Selection: 1
#  # 2 files for record 3697804: 1 to download
#  #  - Downloaded: Sampling_area_borders.xlsx,
#  # - Downloaded: Sampling_area_borders_UTM50N_WGS84.zip

## ----load_safe_data-----------------------------------------------------------
beetle_abund <- load_safe_data(1400562, 'Ant-Psel')
str(beetle_abund)
print(beetle_abund)

## ----explore_loaded-----------------------------------------------------------
show_concepts(beetle_abund)
show_record(beetle_abund)
show_worksheet(beetle_abund)

## ----get_taxa-----------------------------------------------------------------
beetle_taxa <- get_taxa(beetle_abund)
str(beetle_taxa)

## ----add_taxa-----------------------------------------------------------------
beetle_morph <- load_safe_data(1400562, 'MorphFunctTraits')
beetle_morph <- add_taxa(beetle_morph)
str(beetle_morph)

## ----get_taxon_coverage, eval=use_remote--------------------------------------
#  all_taxa <- get_taxon_coverage()

## ----get_taxon_coverage_out---------------------------------------------------
str(all_taxa)

## ----get_phylo----------------------------------------------------------------
library(ape)
beetle_phylo <- get_phylogeny(1400562)
plot(beetle_phylo, show.node.label=TRUE, font=1, no.margin=TRUE)

## ----load_gazetteer-----------------------------------------------------------
gazetteer <- load_gazetteer()
print(gazetteer)

## ----get_locations------------------------------------------------------------
library(sf)
beetle_locs <- get_locations(1400562)
print(beetle_locs)
fragments <- subset(gazetteer, type=='SAFE forest fragment')
par(mar=c(3,3,1,1))
plot(st_geometry(fragments), col='khaki', graticule=TRUE)
plot(st_geometry(beetle_locs), add=TRUE, col='red', pch=4)

## ----add_locations------------------------------------------------------------
beetle_env <- load_safe_data(1400562, 'EnvironVariables')
beetle_env <- add_locations(beetle_env)
print(beetle_env)
plot(beetle_env['Cover'], key.pos=4, breaks=seq(0,100, by=5))

## ----insert_dataset-----------------------------------------------------------
files <- system.file('safedata_example_dir', 'template_ClareWfunctiondata.xlsx', 
                     package='safedata')
insert_dataset(1237719, files)
dat <- load_safe_data(1237719, 'Data')
str(dat)

## ----unset, echo=FALSE--------------------------------------------------------
unset_example_safe_dir()

