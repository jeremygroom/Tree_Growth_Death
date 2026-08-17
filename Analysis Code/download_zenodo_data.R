## ---------------------------------------------------------------------------
## download_zenodo_data.R
##
## Purpose:
##   Ensure the large, zipped processed datasets required by the
##   Tree_Growth_Death analysis are present in the local `Data/` folder.
##   Any listed file that is missing is downloaded from the project's Zenodo
##   archive (identified by DOI) using the zen4R package. Files already on
##   disk are left untouched, so this script is safe to run repeatedly.
##
## Usage:
##   (a) Standalone, once after cloning the repository (run from the project
##       root, which is the default working directory in the RStudio project):
##         source("Analysis Code/download_zenodo_data.R")
##
##   (b) Automatically, from Data_Prep_Analysis.R. Immediately after
##         source("Global.R")          # defines DATA.LOC / CODE.LOC
##       add this line (before the data-loading block):
##         source(paste0(CODE.LOC, "download_zenodo_data.R"))
##
## ---------------------------------------------------------------------------


download_zenodo_data <- function(doi   = ZENODO_DOI,
                                 files = ZENODO_FILES,
                                 dir   = DATA.LOC) {

  ## Which requested files are not yet on disk?
  missing <- files[!file.exists(file.path(dir, files))]

  if (length(missing) == 0L) {
    message("All Zenodo data files already present in '", dir,
            "/'. Nothing to download.")
    return(invisible(character(0)))
  }

  
  message("Missing ", length(missing), " file(s): ",
          paste(missing, collapse = ", "),
          "\nDownloading from Zenodo DOI ", doi, " -> '", dir, " ...")

  ## Download ONLY the missing files. zen4R checks each file's MD5 checksum
  ## after download and reports a mismatch.
  ok <- tryCatch({
    zen4R::download_zenodo(doi   = doi,
                           path  = dir,
                           files = as.list(missing),
                           quiet = FALSE)
    TRUE
  }, error = function(e) {
    message("Zenodo download failed: ", conditionMessage(e))
    FALSE
  })

  ## Confirm every requested file now exists (catches filename mismatches).
  still_missing <- missing[!file.exists(file.path(dir, missing))]
  if (!ok || length(still_missing) > 0L) {
    stop("These file(s) were not retrieved: ",
         paste(if (length(still_missing)) still_missing else missing,
               collapse = ", "),
         ".\nCheck that (1) the DOI is correct, (2) the filenames above match ",
         "the record exactly, and (3) you have internet access. For a restricted ",
         "or embargoed record you may need a Zenodo token (see ?zen4R::ZenodoManager).",
         call. = FALSE)
  }

  message("Done. Retrieved: ", paste(missing, collapse = ", "))
  invisible(missing)
}


##### ---- Run on source() ---- #####
download_zenodo_data()
