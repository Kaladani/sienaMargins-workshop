##############################################
# R script analyzing AK friendship data
# for waves W -> X in all available classes;
# Script 2: Postestimation -- ALL models (0-5), all three estimands
#
# Uses the configuration-object marginalEffects() interface.  Reads the single
# file written by the estimation step, results/estimation_results.rds.
#
# Run from this bundle's root (the folder holding index.qmd):
#   Rscript scripts/Postestimation.R
#
# What gets run is set in the config block below, and can be overridden from
# the shell -- see there.  Static is seconds per model; dynamic and accumulated
# are minutes.
##############################################

## ---------------------------------------------------------------------------
## Getting the package.
##   ## install.packages("remotes")
##   ## remotes::install_github("stocnet/rsiena@sienaMargins")
##   ## renv users:
##   ## renv::install("stocnet/rsiena@sienaMargins")
library(RSiena)


## ---------------------------------------------------------------------------
## WHAT TO RUN.  Edit here, or override any of it from the shell:
##
##   MODELS=0,2 RUN_DYNAMIC=FALSE Rscript scripts/Postestimation.R
##
## The names that are read, and nothing else -- an unrecognised one is ignored
## without warning, so check the "models / estimands / conditioning" line the
## script prints on startup against what you asked for:
##
##   MODELS  ESTIMANDS  CONDITIONINGS  DESC_ESTIMANDS   comma-separated lists
##   RUN_STATIC  RUN_DYNAMIC  RUN_ACCUMULATED           drop one estimand
##   RUN_BYDENSITY  RUN_DESCRIPTIVE                     TRUE/FALSE
##
## Every combination of an ESTIMAND and a CONDITIONING is run for every model
## in MODELS, so the cost is roughly (models x estimands x conditionings).
## Static is seconds; dynamic and accumulated are minutes per model.
## ---------------------------------------------------------------------------

## Which models.  "0".."5", or a subset.
MODELS <- c("0", "1", "2", "3", "4", "5")

## Which estimands.  static = observed waves; dynamic = simulated chains;
## accumulated = summed over the ministeps of a chain (per period, not per
## opportunity).
ESTIMANDS <- c("static", "dynamic", "accumulated")

## Which conditionings.  Each entry is a name and the change statistic(s) to
## stratify by; NULL means no stratification.  "bydensity" splits tie CREATION
## from tie MAINTENANCE, which are different decisions.
CONDITIONINGS <- list(
    pooled    = NULL,
    bydensity = "density"
)

## Model-implied descriptives (the last section): the same comparisons as the
## marginal effects, but WITHOUT holding the rest of the decision fixed.  Set
## RUN_DESCRIPTIVE=FALSE to skip.
RUN_DESCRIPTIVE <- TRUE

## Which estimands to compute them for.  Static is one predict() per class per
## model and takes seconds; dynamic re-simulates the chains and takes minutes,
## which is the only reason it is not on by default.  It is the more complete
## quantity -- see the note at the head of that section.
DESC_ESTIMANDS <- c("static")

set.seed(42)

## --- shell overrides ---------------------------------------------------------
.env_chr <- function(v, default) {
    e <- Sys.getenv(v, NA_character_)
    if (is.na(e) || !nzchar(e)) default else trimws(strsplit(e, ",")[[1L]])
}
.env_lgl <- function(v, default) {
    e <- Sys.getenv(v, NA_character_)
    if (is.na(e)) default else isTRUE(as.logical(e))
}

MODELS    <- .env_chr("MODELS", MODELS)
ESTIMANDS <- .env_chr("ESTIMANDS", ESTIMANDS)
for (nm in c("static", "dynamic", "accumulated"))
    if (!.env_lgl(paste0("RUN_", toupper(nm)), TRUE))
        ESTIMANDS <- setdiff(ESTIMANDS, nm)
CONDITIONINGS <- CONDITIONINGS[.env_chr("CONDITIONINGS", names(CONDITIONINGS))]
if (!.env_lgl("RUN_BYDENSITY", TRUE)) CONDITIONINGS$bydensity <- NULL
RUN_DESCRIPTIVE <- .env_lgl("RUN_DESCRIPTIVE", RUN_DESCRIPTIVE)
DESC_ESTIMANDS  <- .env_chr("DESC_ESTIMANDS", DESC_ESTIMANDS)

stopifnot(length(MODELS) > 0L, length(ESTIMANDS) > 0L,
          length(CONDITIONINGS) > 0L)
cat("models     : ", paste(MODELS, collapse = ", "), "\n",
    "estimands  : ", paste(ESTIMANDS, collapse = ", "), "\n",
    "conditioning: ", paste(names(CONDITIONINGS), collapse = ", "), "\n",
    sep = "")

n3pointest <- 500

# DATA ------------------------------------------------------------------------

## Everything the estimation step produced, in one file.
estimation   <- readRDS("results/estimation_results.rds")
thedata      <- estimation$thedata
alg_controls <- estimation$alg_controls
selection    <- estimation$selection
## Downstream keys these by bare model number; the file names them model_N.
models <- setNames(estimation$results[paste0("model_", 0:5)], as.character(0:5))
effs   <- setNames(estimation$models [paste0("model_", 0:5)], as.character(0:5))

# TARGETS ---------------------------------------------------------------------
#
# The old scripts listed targets as `effectList` entries, each repeating the
# interaction triple (interaction1 / intEffectNames1 / modEffectNames1) with the
# moderator differing per entry.  Here the relationship between effects is
# declared ONCE per targets object and the engine works out which targets it
# reaches, on both sides of a second difference.
#
#   ENDO  unspInt  is the user-defined interaction of recip and transTrip1
#   EXO   egoXaltX is the interaction of egoX and altX
#
# `name =` reproduces the old effectList entry names exactly, so the outputs
# line up one-to-one with the stored old results.

DEP_ENDO <- unspInt  ~ recip:transTrip1
DEP_EXO  <- egoXaltX ~ egoX:altX

mk_base <- function(fit, effects, dynamic, accumulated, deps, condition = NULL,
                    level = "period")
    make_marginal_targets(fit, effects = effects, depvar = NULL,
                         type = "tieProb", level = level,
                         dynamic = dynamic, accumulated = accumulated,
                         condition = condition,
                         dependencies = deps,
                         includeDefaults = FALSE)

add_X <- function(tg)
    set_target(tg, X, contrast = c(0, 1), name = "X_primary", verbose = FALSE)

add_gender <- function(tg, second = TRUE) {
    tg <- set_target(tg, egoX, contrast = c(0, 1),
                     name = "egoX_gender_m_fd", verbose = FALSE)
    tg <- set_target(tg, altX, contrast = c(0, 1),
                     name = "altX_gender_m_fd", verbose = FALSE)
    if (second)
        tg <- set_second_diff(tg, list(egoX = list(contrast = c(0, 1)),
                                       altX = list(contrast = c(0, 1))),
                              name = "egoXaltX_gender_m_sd", verbose = FALSE)
    tg
}

add_recip      <- function(tg)
    set_target(tg, recip, contrast = c(0, 1), name = "recip_fd", verbose = FALSE)
add_transTrip1 <- function(tg)
    set_target(tg, transTrip1, diff = 1, name = "transTrip1_fd", verbose = FALSE)
add_endo_sd    <- function(tg)
    set_second_diff(tg, list(transTrip1 = list(diff = 1),
                             recip      = list(contrast = c(0, 1))),
                    name = "interaction_sd", verbose = FALSE)

## ---------------------------------------------------------------------------
## Per-model specification: which dependencies hold, which targets to compute
## ---------------------------------------------------------------------------
MODEL_SPEC <- list(
    "0" = list(deps    = list(DEP_ENDO),
               build   = function(tg) add_endo_sd(add_transTrip1(add_recip(tg))),
               order   = c("recip_fd", "transTrip1_fd", "interaction_sd")),
    "1" = list(deps    = list(DEP_EXO),
               build   = function(tg) add_gender(add_X(tg)),
               order   = c("X_primary", "egoX_gender_m_fd", "altX_gender_m_fd",
                           "egoXaltX_gender_m_sd")),
    "2" = list(deps    = list(DEP_EXO, DEP_ENDO),
               build   = function(tg)
                   add_endo_sd(add_transTrip1(add_recip(add_gender(add_X(tg))))),
               order   = c("X_primary", "egoX_gender_m_fd", "altX_gender_m_fd",
                           "egoXaltX_gender_m_sd",
                           "recip_fd", "transTrip1_fd", "interaction_sd")),
    "3" = list(deps    = list(DEP_EXO),
               build   = function(tg)
                   add_transTrip1(add_recip(add_gender(add_X(tg)))),
               order   = c("X_primary", "egoX_gender_m_fd", "altX_gender_m_fd",
                           "egoXaltX_gender_m_sd", "recip_fd",
                           "transTrip1_fd")),
    "4" = list(deps    = list(DEP_EXO),
               build   = function(tg) add_recip(add_gender(add_X(tg))),
               order   = c("X_primary", "egoX_gender_m_fd", "altX_gender_m_fd",
                           "egoXaltX_gender_m_sd", "recip_fd")),
    ## Model 5 has no egoXaltX effect, so there is no relation to declare: the
    ## second difference of egoX and altX is then a plain cross-partial with
    ## nothing moving alongside.
    "5" = list(deps    = NULL,
               build   = function(tg) add_gender(add_X(tg)),
               order   = c("X_primary", "egoX_gender_m_fd", "altX_gender_m_fd",
                           "egoXaltX_gender_m_sd"))
)

# RUNNER ----------------------------------------------------------------------

run_block <- function(model, estimand, condition = NULL) {
    spec        <- MODEL_SPEC[[model]]
    fits        <- FITS[[model]]
    effs        <- EFFS[[model]]
    dynamic     <- estimand %in% c("dynamic", "accumulated")
    accumulated <- estimand == "accumulated"
    unc <- set_postest_uncertainty_saom(mode = if (dynamic) "deltaFull" else "delta")
    out <- set_postest_output_saom(combineSameLevel = FALSE)

    res <- vector("list", length(selection))
    for (i in seq_along(selection)) {
        cat(sprintf("  %-11s %-9s model %s  class_%s\n", estimand,
                    if (is.null(condition)) "pooled" else paste(condition, collapse="+"),
                    model, selection[i]))
        utils::flush.console()
        algo <- if (dynamic)
            set_postest_algo_saom(algorithm = alg_controls[[i]],
                                  n3PointEst = n3pointest, verbose = FALSE)
        else
            set_postest_algo_saom(verbose = FALSE)

        tg <- spec$build(mk_base(fits[[i]], effs[[i]], dynamic, accumulated,
                                 spec$deps, condition))
        r  <- marginalEffects(fits[[i]], data = thedata[[i]], targets = tg,
                              control_uncertainty = unc, control_algo = algo,
                              control_out = out)
        if (is.data.frame(r)) r <- setNames(list(r), spec$order[1L])
        stopifnot(setequal(names(r), spec$order))
        res[[i]] <- r[spec$order]
    }
    setNames(res, selection)
}

# COMPILING DATA FRAMES -------------------------------------------------------
#
# Long format: ONE row per (model, class, target, period).
#
# One row per class per target per period, which is what the meta-analysis
# needs: rma() counts rows as studies, so a duplicated row is a fabricated
# class.  The check at the foot of this script enforces it.
#
# `se` is the standard error for whichever uncertainty mode was asked for --
# under deltaFull that means the one INCLUDING the path-distribution channel.
# The frozen-chain SE is still on the per-model results (`SE_conditional`) if
# it is wanted as a diagnostic; it is not what a meta-analysis should pool,
# because it treats one realisation of the simulated chains as if it were the
# only one.

## A conditioned run adds one column per conditioning statistic, named
## <depvar>_<statistic>.  Pick it out by suffix rather than by position, and
## carry it as `stratum` so every run has the same columns whether it was
## conditioned or not.
get_any_effect <- function(target, effect_name, class_id, condition) {
    est_col <- if (grepl("_sd$", effect_name)) "secondDiff" else "firstDiff"
    stratum <- if (is.null(condition)) NA_real_ else {
        col <- grep(paste0("_", condition, "$"), names(target), value = TRUE)[1L]
        if (is.na(col))
            stop("conditioned on '", condition, "' but no matching column in ",
                 "the result for '", effect_name, "': ",
                 paste(names(target), collapse = ", "), call. = FALSE)
        target[[col]]
    }
    data.frame(class_id = class_id,
               effect   = effect_name,
               period   = if (!is.null(target[["period"]])) target[["period"]]
                          else NA_integer_,
               stratum  = stratum,
               est      = target[[est_col]],
               se       = target[["SE"]],
               stringsAsFactors = FALSE)
}

process_model_list <- function(model_list, condition) {
    do.call(rbind, lapply(names(model_list), function(cl) {
        m <- model_list[[cl]]
        do.call(rbind, lapply(names(m),
                              function(nm) get_any_effect(m[[nm]], nm, cl,
                                                          condition)))
    }))
}

compile_block <- function(lst, estimand, cond_label, condition) {
    parts <- lapply(names(lst), function(k) {
        d <- process_model_list(lst[[k]], condition)
        d$model     <- paste0("model_", k)
        d$estimand  <- estimand
        d$condition <- cond_label
        d
    })
    do.call(rbind, parts)
}

# EXECUTE ---------------------------------------------------------------------

## Every estimand is computed twice: pooled over all ministeps, and split by
## density -- which is the split between deciding to CREATE a tie and deciding
## to MAINTAIN one.  Those are different questions and, in this data, often
## have different answers; a pooled average over both is an average over two
## populations rather than an estimate for one.
PARTS <- "results/postest_parts"
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

for (estimand in ESTIMANDS) {
    for (cond_label in names(CONDITIONINGS)) {
        condition <- CONDITIONINGS[[cond_label]]
        cat(sprintf("\n=== %s / %s ===\n", toupper(estimand), cond_label))
        blocks <- setNames(lapply(MODELS, run_block, estimand = estimand,
                                  condition = condition),
                           MODELS)

        for (m in MODELS)
            saveRDS(blocks[[m]],
                    sprintf("results/res_margins_%s_%s_%s.rds",
                            estimand, cond_label, m))

        ## One part file per (estimand, conditioning, MODEL).  Keying on the
        ## model rather than on the set of models is what makes a partial
        ## re-run safe: MODELS=2 overwrites exactly model 2's parts instead of
        ## depositing a second file that the assembly would rbind alongside
        ## the first, duplicating every row for that model.
        df <- compile_block(blocks, estimand, cond_label, condition)
        for (m in MODELS) {
            dm <- df[df$model == paste0("model_", m), , drop = FALSE]
            saveRDS(dm, file.path(PARTS, sprintf("%s_%s_model%s.rds",
                                                 estimand, cond_label, m)))
        }
        cat(sprintf("%s / %s: %d rows\n", estimand, cond_label, nrow(df)))
    }
}

## The combined frame across every estimand -- models 0-5 in one place with
## `estimand` as a column.
##
## Assembled from what is ON DISK rather than from this invocation, because the
## blocks are meant to be runnable a few at a time (see the config block).
## Building it from the current run alone would leave a "combined" frame
## holding only the estimand that ran last.
{
    files <- list.files(PARTS, pattern = "\\.rds$", full.names = TRUE)
    stopifnot(length(files) > 0L)
    parts <- setNames(lapply(files, readRDS), files)

    ## Reading from disk means a file could have been written by an older
    ## version of this script.  Say so, rather than failing inside rbind() with
    ## "number of columns of arguments do not match".
    cols <- lapply(parts, names)
    ok   <- vapply(cols, identical, logical(1L), cols[[1L]])
    if (!all(ok))
        stop("These saved frames have a different set of columns to ",
             names(parts)[1L], ", so they were written by an earlier run:\n  ",
             paste(names(parts)[!ok], collapse = "\n  "),
             "\nDelete them and re-run the corresponding block.", call. = FALSE)
    postest_all_df <- do.call(rbind, parts)
    rownames(postest_all_df) <- NULL
    saveRDS(postest_all_df, "results/postest_all_df.rds")
    cat(sprintf("\ncombined: %d rows -> results/postest_all_df.rds\n",
                nrow(postest_all_df)))

    ## One row per (model, class, target, period).  Cheap, and it fails loudly
    ## rather than silently inflating a downstream meta-analysis.
    key <- postest_all_df[, c("estimand", "condition", "stratum", "model",
                              "class_id", "effect", "period")]
    stopifnot(!anyDuplicated(key))
    cat("no duplicated (estimand, condition, stratum, model, class, target, ",
        "period) rows\n", sep = "")
}


# MODEL-IMPLIED DESCRIPTIVES ---------------------------------------------------
#
# Everything above asks what ONE mechanism contributes, holding the rest of the
# decision fixed.  This section asks the confounded question instead, and asks
# it of the same model at the same ministep: split the predicted tie
# probabilities by a change statistic and take the difference between strata.
#
#     descriptive   P(tie | recip = 1) - P(tie | recip = 0)
#     adjusted      the recip_fd first difference
#
# Both are probabilities attached to one ministep of one model, so they are on
# the same scale and differ in exactly one respect: whether the rest of the
# decision is held fixed.  That is what makes their difference readable as the
# confounding.
#
# Static first.  Static evaluates in the observed contexts, which keeps it close
# to the data and cheap enough to explore with.  The dynamic version follows the
# model's own simulated sequences from one wave to the next and is the fuller
# answer; set DESC_ESTIMANDS to run it.

## For each marginal target: which change statistic reproduces it, whether the
## contrast is a plain difference or the 2x2 that mirrors a second difference,
## which levels to use, and whose choice set to compare inside.
##
## `levels`  the levels the marginal target itself contrasts.  Pinning matters:
##           RSiena mean-imputes a MISSING covariate, which adds a third level
##           between 0 and 1 in classes where some pupils' gender is unrecorded.
##
## `scope`   an actor chooses among ITS OWN alternatives, so where a statistic
##           varies across those alternatives the comparison belongs inside one
##           choice set.  An ego covariate is the exception: it is the same for
##           all of an actor's alternatives -- an actor cannot choose to be a
##           different gender -- so that comparison is necessarily across actors.

DESC_SPEC <- list(
    recip_fd             = list(stats = "recip",      kind = "fd",
                                levels = list(c(0, 1)), scope = "withinEgo"),
    transTrip1_fd        = list(stats = "transTrip1", kind = "fd",
                                levels = list(NULL),    scope = "withinEgo"),
    X_primary            = list(stats = "X",          kind = "fd",
                                levels = list(c(0, 1)), scope = "withinEgo"),
    altX_gender_m_fd     = list(stats = "altX",       kind = "fd",
                                levels = list(c(0, 1)), scope = "withinEgo"),
    interaction_sd       = list(stats = c("transTrip1", "recip"), kind = "sd",
                                levels = list(NULL, c(0, 1)),
                                scope = "withinEgo"),
    egoX_gender_m_fd     = list(stats = "egoX",       kind = "fd",
                                levels = list(c(0, 1)), scope = "betweenEgo"),
    egoXaltX_gender_m_sd = list(stats = c("egoX", "altX"), kind = "sd",
                                levels = list(c(0, 1), c(0, 1)),
                                scope = "betweenEgo")
)

## Conditioning columns are named <depvar>_<statistic>, and for a covariate
## effect the covariate name is appended (friendship_egoX_sex_m), so the
## statistic is matched as a whole segment rather than as a suffix.
desc_col <- function(d, stat) {
    hit <- grep(paste0("_", stat, "($|_)"), names(d), value = TRUE)
    if (!length(hit))
        stop("no column for conditioning statistic '", stat, "' in: ",
             paste(names(d), collapse = ", "), call. = FALSE)
    hit[1L]
}

## How many decision alternatives sit at each level of each change statistic.
## These are the WEIGHTS: they come from the same change statistics predict()
## aggregates, so the descriptive is averaged over the same decisions the
## marginal effect is.
desc_weights <- function(fit, data, effects, depvar) {
    cs <- RSiena:::getStaticChangeContributions(ans = fit, data = data, effects = effects,
                                       depvar = depvar, returnWide = TRUE)
    W        <- as.data.frame(cs$changeStats$csMat)
    names(W) <- cs$changeStats$csNames
    W$period <- cs$period
    W$ego    <- cs$ego
    W
}

## The descriptive counterpart of a `diff = 1` target.  That target moves each
## decision one unit up from wherever it sits, then averages; so this compares
## strata one unit apart and averages over the same decisions:
##
##     sum_l  w_l * [ P(l_next) - P(l) ] / (l_next - l)   /   sum_l w_l
##
## Two levels one unit apart -- every binary statistic here -- reduce this to
## the plain difference P(hi) - P(lo), so one rule covers both kinds.
##
## It stays CONFOUNDED: P(l_next) and P(l) average over different sets of dyads,
## which differ in reciprocation, gender and primary school as well as in the
## statistic.  Nothing is held fixed within a step.
##
## Returned as a coefficient vector rather than a number, so the estimate stays
## a linear combination of the strata -- which is what makes the standard error
## exact.
desc_coefs <- function(lev, n) {
    o   <- order(lev)
    cf  <- numeric(length(lev))
    tot <- 0
    for (k in seq_len(length(o) - 1L)) {
        a <- o[k]; b <- o[k + 1L]
        step <- lev[b] - lev[a]
        wk   <- n[a]
        if (!is.finite(step) || step <= 0) next
        if (!is.finite(wk)   || wk   <= 0) next
        cf[b] <- cf[b] + wk / step
        cf[a] <- cf[a] - wk / step
        tot   <- tot + wk
    }
    if (tot <= 0) return(NULL)
    cf / tot
}

## One prediction frame -> one row per (period, conditioning stratum).
desc_contrast <- function(d, W, V, effect_name, dspec, condition, class_id) {
    J <- attr(d, "delta_jacobians")$J_full
    if (is.null(J))
        stop("no delta Jacobian on the prediction for '", effect_name,
             "' -- uncertainty must be enabled with mode = 'delta'.",
             call. = FALSE)
    V <- V[colnames(J), colnames(J), drop = FALSE]
    d <- as.data.frame(d)

    stat_cols <- vapply(dspec$stats, desc_col, character(1L), d = d)
    if (!all(stat_cols %in% names(W)))
        stop("no change statistic column for '", effect_name, "' in the ",
             "weights: ", paste(setdiff(stat_cols, names(W)), collapse = ", "),
             call. = FALSE)

    ## Drop levels this target does not contrast (the mean-imputed
    ## "covariate missing" level above all) before grouping or counting.
    allow <- rep(TRUE, nrow(d))
    for (k in seq_along(stat_cols)) {
        lk <- dspec$levels[[k]]
        if (!is.null(lk)) allow <- allow & d[[stat_cols[k]]] %in% lk
    }

    ## Group by period and by the conditioning stratum, so the comparison is
    ## taken WITHIN a density level rather than across it.
    within    <- identical(dspec$scope, "withinEgo")
    strat_col <- if (is.null(condition)) NULL
                 else vapply(condition, desc_col, character(1L), d = d)
    grp_cols  <- c(if (!is.null(d[["period"]])) "period", strat_col)
    gid <- if (length(grp_cols))
        do.call(paste, c(d[, grp_cols, drop = FALSE], sep = "\r"))
    else rep("1", nrow(d))

    ## Decision counts at the levels of `sc`, among the alternatives matching
    ## this group's period and density stratum (and, for a second difference,
    ## the moderator level the inner first difference is taken at).
    counts_for <- function(sc, levels_wanted, keep) {
        tb <- table(W[[sc]][keep])
        n  <- as.numeric(tb[as.character(levels_wanted)])
        n[is.na(n)] <- 0
        n
    }

    ## The coefficient vector for one block of rows sharing a choice set: within
    ## one ego for a withinEgo target, or the whole stratum for a betweenEgo one.
    block_coefs <- function(ix, keep) {
        if (identical(dspec$kind, "fd")) {
            lev <- d[[stat_cols[1L]]][ix]
            cf  <- desc_coefs(lev, counts_for(stat_cols[1L], lev, keep))
            if (is.null(cf)) return(NULL)
            out <- numeric(nrow(d)); out[ix] <- cf
            return(out)
        }
        ## First difference in the first statistic WITHIN each level of the
        ## moderator, then differenced -- the mirror of a second difference.
        ## For two binary statistics this is the 2x2 interaction contrast.
        mod  <- d[[stat_cols[2L]]][ix]
        mlev <- sort(unique(mod[!is.na(mod)]))
        if (length(mlev) != 2L) return(NULL)
        out <- numeric(nrow(d))
        for (s in seq_along(mlev)) {
            jx  <- ix[mod == mlev[s]]
            lev <- d[[stat_cols[1L]]][jx]
            cf  <- desc_coefs(lev,
                counts_for(stat_cols[1L], lev,
                           keep & W[[stat_cols[2L]]] == mlev[s]))
            if (is.null(cf)) return(NULL)
            out[jx] <- if (s == 2L) cf else -cf
        }
        out
    }

    do.call(rbind, lapply(unique(gid), function(g) {
        ix <- which(gid == g & allow)
        if (!length(ix)) return(NULL)
        ## The alternatives this group's strata were averaged over.
        keep0 <- rep(TRUE, nrow(W))
        if ("period" %in% grp_cols) keep0 <- keep0 & W$period == d$period[ix[1L]]
        for (sc in strat_col) keep0 <- keep0 & W[[sc]] == d[[sc]][ix[1L]]

        if (within) {
            ## One contrast per ego, then the egos averaged equally -- the same
            ## way the marginal effect averages them.  An ego whose choice set
            ## holds only one level of the statistic drops out; `n_egos` records
            ## how many remained.
            egos <- unique(d$ego[ix])
            cfs  <- lapply(egos, function(e)
                block_coefs(ix[d$ego[ix] == e], keep0 & W$ego == e))
            good <- !vapply(cfs, is.null, logical(1L))
            kept <- egos[good]
            if (!length(kept)) return(NULL)
            cvec <- Reduce(`+`, cfs[good]) / length(kept)
        } else {
            cvec <- block_coefs(ix, keep0)
            kept <- NULL
            if (is.null(cvec)) return(NULL)
        }

        ## Only the rows the contrast actually touches.
        nz  <- which(cvec != 0)
        est <- sum(cvec[nz] * d$tieProb[nz])
        gv  <- colSums(cvec[nz] * J[nz, , drop = FALSE])
        se  <- sqrt(max(as.numeric(t(gv) %*% V %*% gv), 0))

        data.frame(class_id = class_id,
                   effect   = effect_name,
                   period   = if ("period" %in% grp_cols) d$period[ix[1L]]
                              else NA_integer_,
                   stratum  = if (is.null(strat_col)) NA_real_
                              else d[[strat_col[1L]]][ix[1L]],
                   est      = est,
                   se       = se,
                   scope    = dspec$scope,
                   n_egos   = if (within) length(kept) else NA_integer_,
                   egos     = I(list(kept)),
                   stringsAsFactors = FALSE)
    }))
}

## The ADJUSTED effect on exactly the egos the descriptive could use.  A
## within-ego comparison needs an ego whose choice set holds both levels, so
## some egos drop out; measuring the adjusted effect on a different set of egos
## would put a second difference between the two numbers.
matched_adjusted <- function(m, V, rows) {
    J <- attr(m, "delta_jacobians")$J_full
    if (is.null(J))
        stop("no delta Jacobian on the marginal effect -- uncertainty must be ",
             "enabled with mode = 'delta'.", call. = FALSE)
    V <- V[colnames(J), colnames(J), drop = FALSE]
    m <- as.data.frame(m)
    est_col <- if ("secondDiff" %in% names(m)) "secondDiff" else "firstDiff"

    ## Egos with no usable decision are dropped; `n_adj` records how many
    ## carried the average.
    usable <- is.finite(m[[est_col]]) & is.finite(rowSums(J))

    do.call(rbind, lapply(seq_len(nrow(rows)), function(r) {
        keep <- usable
        if (!is.na(rows$period[r]) && !is.null(m[["period"]]))
            keep <- keep & m$period == rows$period[r]
        if (!is.na(rows$stratum[r])) {
            sc <- grep("_density$", names(m), value = TRUE)[1L]
            if (!is.na(sc)) keep <- keep & m[[sc]] == rows$stratum[r]
        }
        eg <- rows$egos[[r]]
        if (!is.null(eg)) keep <- keep & m$ego %in% eg
        ix <- which(keep)
        if (!length(ix))
            return(data.frame(est_adj = NA_real_, se_adj = NA_real_,
                              n_adj = 0L))
        gv <- colMeans(J[ix, , drop = FALSE])
        data.frame(est_adj = mean(m[[est_col]][ix]),
                   se_adj  = sqrt(max(as.numeric(t(gv) %*% V %*% gv), 0)),
                   n_adj   = length(ix))
    }))
}

run_descriptive_block <- function(model, condition, estimand = "static") {
    spec    <- MODEL_SPEC[[model]]
    fits    <- FITS[[model]]
    effs    <- EFFS[[model]]
    wanted  <- intersect(spec$order, names(DESC_SPEC))
    if (!length(wanted)) return(NULL)
    dynamic <- estimand == "dynamic"

    unc <- set_postest_uncertainty_saom(mode = if (dynamic) "deltaFull" else "delta")
    out <- set_postest_output_saom(combineSameLevel = FALSE)

    res <- vector("list", length(selection))
    for (i in seq_along(selection)) {
        cat(sprintf("  descriptive %-9s %-7s model %s  class_%s\n",
                    if (is.null(condition)) "pooled" else paste(condition, collapse = "+"),
                    estimand, model, selection[i]))
        utils::flush.console()

        algo <- if (dynamic)
            set_postest_algo_saom(algorithm = alg_controls[[i]],
                                  n3PointEst = n3pointest, verbose = FALSE)
        else set_postest_algo_saom(verbose = FALSE)

        tg <- make_predict_targets(fits[[i]], effects = effs[[i]], depvar = NULL,
                                   type = "tieProb", level = "period",
                                   dynamic = dynamic,
                                   includeDefaults = FALSE)
        ## The conditioning statistic goes on TOP of the density split, so the
        ## contrast is taken within creation and within maintenance.  `level`
        ## carries the scope: an ego-level prediction keeps one row per ego, so
        ## the contrast can be taken inside a single choice set.
        for (nm in wanted)
            tg <- set_condition(tg, c(condition, DESC_SPEC[[nm]]$stats),
                                level = if (identical(DESC_SPEC[[nm]]$scope,
                                                      "withinEgo")) "ego"
                                        else "period",
                                name = nm, verbose = FALSE)

        pr <- predict(fits[[i]], data = thedata[[i]], targets = tg,
                      control_uncertainty = unc, control_algo = algo,
                      control_out = out)
        if (!is.list(pr) || inherits(pr, "data.frame"))
            pr <- setNames(list(pr), wanted[1L])

        ## The adjusted effects at EGO level, so they can be averaged over
        ## whichever egos the descriptive kept.
        mtg <- spec$build(mk_base(fits[[i]], effs[[i]], dynamic, FALSE,
                                  spec$deps, condition, level = "ego"))
        mr <- marginalEffects(fits[[i]], data = thedata[[i]], targets = mtg,
                              control_uncertainty = unc, control_algo = algo,
                              control_out = out)
        if (is.data.frame(mr)) mr <- setNames(list(mr), spec$order[1L])

        ## One call per class, shared by every target: the weights are a
        ## property of the class's decision alternatives, not of the target.
        W <- desc_weights(fits[[i]], thedata[[i]], effs[[i]],
                          attr(tg, "depvar"))
        V <- vcov(fits[[i]])
        res[[i]] <- do.call(rbind, lapply(wanted, function(nm) {
            rows <- desc_contrast(pr[[nm]], W, V, nm, DESC_SPEC[[nm]], condition,
                                  selection[i])
            if (is.null(rows)) return(NULL)
            cbind(rows[, setdiff(names(rows), "egos")],
                  matched_adjusted(mr[[nm]], V, rows))
        }))
    }
    do.call(rbind, res)
}

if (RUN_DESCRIPTIVE) {
    DPARTS <- "results/postest_descr_parts"
    dir.create(DPARTS, showWarnings = FALSE, recursive = TRUE)

    for (de in DESC_ESTIMANDS) {
    for (cond_label in names(CONDITIONINGS)) {
        condition <- CONDITIONINGS[[cond_label]]
        cat(sprintf("\n=== DESCRIPTIVE / %s / %s ===\n", de, cond_label))
        for (m in MODELS) {
            dm <- run_descriptive_block(m, condition, estimand = de)
            if (is.null(dm)) next
            dm$model     <- paste0("model_", m)
            dm$estimand  <- de
            dm$condition <- cond_label
            saveRDS(dm, file.path(DPARTS, sprintf("descr_%s_%s_model%s.rds",
                                                  de, cond_label, m)))
            cat(sprintf("descriptive / %s / %s / model_%s: %d rows\n",
                        de, cond_label, m, nrow(dm)))
        }
    }
    }

    ## Assembled from disk for the same reason the marginal frame is: the blocks
    ## are meant to be runnable a few models at a time.
    files <- list.files(DPARTS, pattern = "\\.rds$", full.names = TRUE)
    if (length(files)) {
        parts <- setNames(lapply(files, readRDS), files)
        cols  <- lapply(parts, names)
        ok    <- vapply(cols, identical, logical(1L), cols[[1L]])
        if (!all(ok))
            stop("These saved frames have a different set of columns to ",
                 names(parts)[1L], ", so they were written by an earlier run:\n  ",
                 paste(names(parts)[!ok], collapse = "\n  "),
                 "\nDelete them and re-run the corresponding block.", call. = FALSE)
        descr_df <- do.call(rbind, parts)
        rownames(descr_df) <- NULL
        saveRDS(descr_df, "results/postest_descr_df.rds")
        cat(sprintf("\ndescriptive combined: %d rows -> ", nrow(descr_df)),
            "results/postest_descr_df.rds\n", sep = "")

        key <- descr_df[, c("estimand", "condition", "stratum", "model",
                            "class_id", "effect", "period")]
        stopifnot(!anyDuplicated(key))
        nmiss <- sum(is.na(descr_df$est))
        if (nmiss)
            cat(sprintf("note: %d contrast(s) undefined -- the class has no ",
                        nmiss), "ministeps at one of the two levels\n", sep = "")
    }
}
