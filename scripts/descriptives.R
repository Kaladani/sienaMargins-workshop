##############################################
# R script analyzing AK friendship data
# for waves W -> X in all available classes;
# Script 1b: Descriptives across school classes
#
# Two things, both computed from the data alone with no model involved:
#
#   (1) ordinary network descriptives per class and wave -- size, density,
#       degree, reciprocity, transitivity, and the Jaccard index between the
#       two waves.
#
#   (2) the DESCRIPTIVE DIFFERENCE for each condition that a later marginal
#       effect targets: the difference in tie proportions at t2 between dyads
#       where a t1 condition holds and dyads where it does not,
#
#           P(x_ij(t2) = 1 | c_ij = 1)  -  P(x_ij(t2) = 1 | c_ij = 0)
#
#       This needs no model and is confounded by construction -- it reflects
#       every mechanism at once.
#
# Run from this bundle's root (the folder holding index.qmd):
#   Rscript scripts/descriptives.R
##############################################

WORK <- "results"   # .rds inputs and outputs
OUT  <- "output"    # tables and plots
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

estimation <- readRDS(file.path(WORK, "estimation_results.rds"))
thedata    <- estimation$thedata
selection  <- estimation$selection
names(thedata) <- selection

# HELPERS ---------------------------------------------------------------------

## Structural zeros/ones are coded 10/11 in RSiena; the diagonal is not a dyad.
as_tie <- function(m) { m[m >= 10] <- NA; diag(m) <- NA; m }

## Two-path counts, treating missing as absent.
two_paths <- function(x) { a <- x; a[is.na(a)] <- 0; tp <- a %*% a; diag(tp) <- 0; tp }

net_descriptives <- function(net, class_id) {
    do.call(rbind, lapply(seq_len(dim(net)[3L]), function(w) {
        x  <- as_tie(net[, , w])
        a  <- x; a[is.na(a)] <- 0
        tp <- two_paths(x)
        ties <- sum(x == 1, na.rm = TRUE)
        data.frame(class_id       = class_id,
                   wave           = w,
                   actors         = nrow(x),
                   ties           = ties,
                   density        = mean(x == 1, na.rm = TRUE),
                   mean_degree    = ties / nrow(x),
                   recip_edgewise = sum(x == 1 & t(x) == 1, na.rm = TRUE) / ties,
                   transitivity   = sum(tp * a) / sum(tp),
                   missing_frac   = mean(is.na(x)) -  1 / nrow(x),
                   stringsAsFactors = FALSE)
    }))
}

jaccard <- function(a, b) {
    a <- as_tie(a); b <- as_tie(b)
    k <- !is.na(a) & !is.na(b)
    n11 <- sum(a[k] == 1 & b[k] == 1)
    n10 <- sum(a[k] == 1 & b[k] == 0)
    n01 <- sum(a[k] == 0 & b[k] == 1)
    n11 / (n11 + n10 + n01)
}

## Each descriptive has to be the observable counterpart of a PARTICULAR
## marginal-effect contrast, which means matching how that contrast is defined.
## Two kinds appear here.
##
##   "level"  the condition is binary and the marginal effect contrasts 0 with 1
##            (reciprocated or not, same primary school or not, male or not).
##            The counterpart is the plain difference in tie proportions.
##
##   "unit"   the condition is a count and the marginal effect is a +1 shift
##            from whatever the dyad actually has -- not a contrast between
##            "none" and "some".  Comparing dyads with any two-paths against
##            dyads with none would answer a different and larger question,
##            because dyads with many two-paths sit in the first group.
##
## Read the spread across classes as the measure of how firm either is.

## "level": difference in tie proportions at t2 between dyads meeting a binary
## condition at t1 and dyads not meeting it.
descriptive_difference <- function(x1, x2, cond, label, class_id) {
    x2   <- as_tie(x2)
    keep <- !is.na(x2) & !is.na(cond)
    n1   <- sum(keep &  cond); n0 <- sum(keep & !cond)
    if (n1 == 0L || n0 == 0L)
        return(data.frame(class_id = class_id, condition = label, kind = "level",
                          n_cond = n1, n_not = n0,
                          p_t2_cond = NA_real_, p_t2_not = NA_real_,
                          n_steps = NA_integer_, est = NA_real_,
                          stringsAsFactors = FALSE))
    p1 <- mean(x2[keep &  cond] == 1)
    p0 <- mean(x2[keep & !cond] == 1)
    data.frame(class_id = class_id, condition = label, kind = "level",
               n_cond = n1, n_not = n0,
               p_t2_cond = p1, p_t2_not = p0,
               n_steps = NA_integer_, est = p1 - p0,
               stringsAsFactors = FALSE)
}

## "unit": the average effect of one more, over the levels that actually occur.
## For each level k, the change in tie proportion at t2 between dyads with k+1
## and dyads with k at t1; averaged over k weighted by how many dyads sit at k,
## since k is where a +1 shift starts.  Levels with fewer than `min_n` dyads on
## either side are skipped, as are gaps in the observed counts.
descriptive_unit_difference <- function(x1, x2, cnt, label, class_id,
                                       min_n = 5L) {
    x2   <- as_tie(x2)
    keep <- !is.na(x2) & !is.na(cnt)
    y    <- x2[keep] == 1
    k    <- cnt[keep]

    empty <- data.frame(class_id = class_id, condition = label, kind = "unit",
                        n_cond = 0L, n_not = 0L,
                        p_t2_cond = NA_real_, p_t2_not = NA_real_,
                        n_steps = 0L, est = NA_real_, stringsAsFactors = FALSE)
    if (!length(k)) return(empty)

    lev <- sort(unique(k))
    f   <- factor(k, levels = lev)
    p   <- tapply(y, f, mean)
    n   <- tapply(y, f, length)

    steps <- do.call(rbind, lapply(seq_len(length(lev) - 1L), function(i) {
        if (lev[i + 1L] != lev[i] + 1L) return(NULL)      # gap in the counts
        if (n[i] < min_n || n[i + 1L] < min_n) return(NULL)
        data.frame(d = p[i + 1L] - p[i], w = n[i],
                   p_hi = p[i + 1L], p_lo = p[i], n_hi = n[i + 1L])
    }))
    if (is.null(steps)) return(empty)

    data.frame(class_id = class_id, condition = label, kind = "unit",
               n_cond = sum(steps$n_hi), n_not = sum(steps$w),
               p_t2_cond = weighted.mean(steps$p_hi, steps$w),
               p_t2_not  = weighted.mean(steps$p_lo, steps$w),
               n_steps   = nrow(steps),
               est       = weighted.mean(steps$d, steps$w),
               stringsAsFactors = FALSE)
}

# CONDITIONS ------------------------------------------------------------------
#
# Each condition is the observable, uncontrolled counterpart of a target in
# Postestimation.R.  The names on the left are the marginal-effect target
# names, so the two tables join on them.

conditions_for <- function(dat) {
    x1  <- as_tie(dat$depvars$friendship[, , 1L])
    sex <- as.vector(dat$cCovars$sex_m)
    prm <- dat$dycCovars$primary
    male <- sex == 1

    ## Third element is the kind: "level" for a binary 0-vs-1 contrast, "unit"
    ## for one-more on a count.  transTrip1 is the only count here.
    list(
        recip_fd             = list(t(x1) == 1,
                                    "j nominates i at t1", "level"),
        transTrip1_fd        = list(two_paths(x1),
                                    "one more two-path at t1", "unit"),
        X_primary            = list(prm == 1,
                                    "primary-school friend", "level"),
        egoX_gender_m_fd     = list(outer(male, rep(TRUE, length(male))),
                                    "ego is male", "level"),
        altX_gender_m_fd     = list(outer(rep(TRUE, length(male)), male),
                                    "alter is male", "level"),
        egoXaltX_gender_m_sd = list(outer(male, male),
                                    "both male", "level")
    )
}

# COMPUTE ---------------------------------------------------------------------

desc_rows <- list(); diff_rows <- list(); jac_rows <- list()

for (cl in selection) {
    dat <- thedata[[cl]]
    net <- dat$depvars$friendship

    desc_rows[[cl]] <- net_descriptives(net, cl)
    jac_rows[[cl]]  <- data.frame(class_id = cl,
                                  jaccard = jaccard(net[, , 1L], net[, , 2L]),
                                  stringsAsFactors = FALSE)

    conds <- conditions_for(dat)
    diff_rows[[cl]] <- do.call(rbind, lapply(names(conds), function(tgt) {
        cnd <- conds[[tgt]]
        fn  <- if (cnd[[3L]] == "unit") descriptive_unit_difference
               else                     descriptive_difference
        cbind(target = tgt,
              fn(net[, , 1L], net[, , 2L], cnd[[1L]], cnd[[2L]], cl))
    }))
}

descriptives  <- do.call(rbind, desc_rows)
jaccards      <- do.call(rbind, jac_rows)
descriptive_diffs <- do.call(rbind, diff_rows)
rownames(descriptives) <- rownames(descriptive_diffs) <- NULL

saveRDS(descriptives,      file.path(WORK, "descriptives_by_class.rds"))
saveRDS(jaccards,          file.path(WORK, "jaccard_by_class.rds"))
saveRDS(descriptive_diffs, file.path(WORK, "descriptive_diffs_by_class.rds"))
write.csv(descriptives,      file.path(OUT, "table_descriptives_by_class.csv"),
          row.names = FALSE)
write.csv(descriptive_diffs, file.path(OUT, "table_descriptive_diffs_by_class.csv"),
          row.names = FALSE)

# REPORT ----------------------------------------------------------------------

rng <- function(v) sprintf("%.3f [%.3f, %.3f]", mean(v, na.rm = TRUE),
                           min(v, na.rm = TRUE), max(v, na.rm = TRUE))

cat("\n--- network descriptives across", length(selection), "classes ---\n")
for (w in sort(unique(descriptives$wave))) {
    d <- descriptives[descriptives$wave == w, ]
    cat(sprintf("wave %d: actors %s | density %s | degree %s | recip %s | trans %s\n",
                w, rng(d$actors), rng(d$density), rng(d$mean_degree),
                rng(d$recip_edgewise), rng(d$transitivity)))
}
cat("jaccard w1->w2: ", rng(jaccards$jaccard), "\n", sep = "")

cat("\n--- descriptive differences (mean over classes) ---\n")
cat("  level = binary 0 vs 1; unit = one more, averaged over levels\n")
agg <- aggregate(cbind(est, p_t2_cond, p_t2_not) ~ target + condition + kind,
                 data = descriptive_diffs, FUN = mean, na.rm = TRUE)
agg <- agg[order(-abs(agg$est)), ]
for (i in seq_len(nrow(agg))) with(agg[i, ],
    cat(sprintf("  %-22s %-26s %-5s  P(hi)=%.3f  P(lo)=%.3f  diff=%+.3f\n",
                target, condition, kind, p_t2_cond, p_t2_not, est)))

cat("\nSaved: descriptives_by_class.rds, descriptive_diffs_by_class.rds\n")
