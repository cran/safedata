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

## ----search_text, echo=TRUE---------------------------------------------------
soil_datasets <- search_text('soil')
print(soil_datasets)

## ----search_taxa, collapse=TRUE-----------------------------------------------
print(ants <- search_taxa('Formicidae'))

## ----search_taxa_id, collapse=TRUE--------------------------------------------
ants <- search_taxa(gbif_id=4342)

## ----search_spatial, collapse=TRUE--------------------------------------------
# Datasets that include sampling within experimental block A
within_a <- search_spatial(location='BL_A')
# Datasets that sampled within 2 km of the Maliau Basin Field Study Centre
near_maliau <- search_spatial(wkt='POINT(116.97394 4.73481)', distance=2000)

## ----validate_recs, echo=TRUE-------------------------------------------------
recs <- validate_record_ids(c('https://doi.org/10.5281/zenodo.3247631',
                              '10.5281/zenodo.3266827',
                              'https://zenodo.org/record/3266821'))
print(recs)

## ----show_concepts------------------------------------------------------------
show_concepts(recs)

## ----show_record--------------------------------------------------------------
show_record(recs[3,])

## ----show_worksheet-----------------------------------------------------------
show_worksheet(recs[3,], 'Data')

## ----show_worksheet_long, output.lines=15-------------------------------------
show_worksheet(recs[3,], 'Data', extended_fields=TRUE)

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

## ----get_taxon_coverage-------------------------------------------------------
all_taxa <- get_taxon_coverage()
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

