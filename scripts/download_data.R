##############################################
# Fetch the AK data.
#
# The friendship networks and covariates are not distributed with this
# repository -- they belong to the study they come from, and are hosted
# separately.  See data/README.md for the source and the terms.
#
# Everything that needs them calls ensure_ak_data() first, which downloads
# whatever is missing into data/ and then gets out of the way.  Files already
# present are left alone, so this is safe to call repeatedly and costs nothing
# once the data are there.
#
# Run from this bundle's root (the folder holding index.qmd).
##############################################

## ---------------------------------------------------------------------------
## PLACEHOLDER ADDRESS -- replace before circulating.
##
## Expected to be a directory containing the files named in AK_FILES below, so
## that <AK_DATA_URL><filename> resolves.  Keep the trailing slash.
## ---------------------------------------------------------------------------
AK_DATA_URL <- "https://example.org/workshops/AK/"

## The three files the code actually reads.  The pupil questionnaire is at the
## same address but nothing here needs it, so it is not fetched by default.
AK_FILES <- c("AK_friendship.RData", "AK_gender.RData", "AK_primary.RData")


ensure_ak_data <- function(files = AK_FILES,
                           dir   = "data",
                           url   = AK_DATA_URL,
                           quiet = FALSE) {

    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    dest <- file.path(dir, files)
    need <- !file.exists(dest)

    if (!any(need)) return(invisible(dest))

    if (grepl("example\\.org", url, fixed = FALSE))
        warning("scripts/download_data.R still has the placeholder address. ",
                "Set AK_DATA_URL, or put the files in ", dir, "/ by hand.",
                call. = FALSE)

    for (i in which(need)) {
        if (!quiet) message("downloading ", files[i], " ...")
        ok <- tryCatch({
            utils::download.file(paste0(url, files[i]), dest[i],
                                 mode = "wb", quiet = TRUE)
            TRUE
        }, error = function(e) FALSE, warning = function(w) FALSE)

        ## A failed download can still leave a stub behind; do not let a later
        ## call mistake it for the real file.
        if (!ok || !file.exists(dest[i]) || file.size(dest[i]) == 0) {
            unlink(dest[i])
            stop("could not fetch ", files[i], " from ", url, "\n",
                 "  Download the AK files by hand and put them in ", dir, "/ -- ",
                 "see data/README.md for where to get them.", call. = FALSE)
        }
    }

    invisible(dest)
}
