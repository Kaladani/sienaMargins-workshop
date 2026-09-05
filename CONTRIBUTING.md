# Editing this material

The site rebuilds itself. Push to `main` and a GitHub Action re-renders the book
and publishes it, usually within a minute or two. There is no need to run Quarto
locally, or to run `quarto publish`.

## Changing words

Edit the `.qmd` and push — including straight from GitHub's web editor. That is
all. This covers the reading chapters (`01`, `02`, `99`, `index`) completely, and
covers the prose in the executable chapters (`03`, `04`, `05`) too: the results in
those chapters come from the cache in `_freeze/`, which is reused rather than
recomputed, so text can be edited without anything being refitted.

## Changing code

If you change a **code chunk** in `03`, `04` or `05`, the cache no longer matches
the code and has to be rebuilt — and that cannot happen on the build machine,
which has neither RSiena nor the data. So:

1. Render locally: `quarto render`. This needs the `sienaMargins` branch of
   RSiena and the data in `data/` (see `data/README.md`).
2. Commit the updated `_freeze/` along with your change.

If you forget, the build fails rather than quietly publishing stale results — you
will see a red mark on the commit, and the site keeps the previous version until
it is fixed.

## Why it is arranged this way

`_freeze/` holds the output of every chunk, so rendering the book does not require
re-estimating anything. That is what lets the site rebuild without RSiena, and
without the data, which are not distributed here.

Not everything in the repository is rebuilt by the site. `results/` and `output/`
come from the scripts in `scripts/`, which are run by hand when the analysis
changes; `scripts/README.md` has the run order and what each stage costs.
