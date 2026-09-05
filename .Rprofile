## Development convenience only.
##
## This bundle is meant to be self-contained: with the packages listed in
## scripts/README.md installed, nothing here is needed. If a directory named
## renv/ happens to sit one level up (e.g. this repo checked out as a sibling
## of a renv-managed RSiena checkout), R will use its library instead of the
## system one. Harmless no-op otherwise, which is the common case.
local({
  activate <- file.path("..", "renv", "activate.R")
  if (file.exists(activate)) {
    Sys.setenv(RENV_PROJECT = normalizePath(".."))
    source(activate)
  }
})
