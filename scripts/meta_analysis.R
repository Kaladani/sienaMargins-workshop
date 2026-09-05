##############################################
# R script analyzing AK friendship data
# for waves W -> X in all available classes;
# Script 3: Meta-Analysis & ggplotting -- ALL models (0-5), all estimands
#
# Consumes the output of Postestimation.R:
#   results/postest_all_df.rds
#
# Run from this bundle's root (the folder holding index.qmd):
#   Rscript scripts/meta_analysis.R
#
# One row per (model, class, target, period) goes in, so `k` is the number of
# classes contributing to each pooled estimate -- 17 here.  It is printed, and
# checked at the foot of the script, because rma() counts rows as studies: a
# duplicated row is a fabricated class, and it moves more than the standard
# error.  Copies carry no new between-class information, so they bias tau^2
# downward, and since the random-effects weights are 1/(se^2 + tau^2) the
# pooled ESTIMATE shifts too.  Q, and hence the heterogeneity test, scales with
# k directly.
##############################################

library(metafor)
library(ggplot2)

WORK <- "results"   # .rds inputs and outputs
OUT  <- "output"    # tables and plots
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# LOAD DATA -------------------------------------------------------------------

## Raw siena07 coefficients: one table for all models, built by the
## estimation step from the same list it estimated.
raw_df <- readRDS(file.path(WORK, "estimation_results.rds"))$coefs

## Marginal effects: one frame, `estimand` as a column.
postest_all <- readRDS(file.path(WORK, "postest_all_df.rds"))

MODEL_LEVELS <- paste0("model_", 0:5)
MODEL_COLOURS <- c(model_0 = "#4477AA", model_1 = "#EE6677", model_2 = "#AA3377",
                   model_3 = "#228833", model_4 = "#CCBB44", model_5 = "#66CCEE")

# FUNCTIONS -------------------------------------------------------------------

calculate_meta <- function(df_long) {
  effs <- unique(as.character(df_long$effect))
  mods <- intersect(MODEL_LEVELS, unique(as.character(df_long$model)))

  do.call(rbind, lapply(effs, function(eff) {
    do.call(rbind, lapply(mods, function(mod) {
      sub_df <- subset(df_long, effect == eff & model == mod)
      if (nrow(sub_df) == 0L) return(NULL)  # effect not in this model's spec
      ## REML does not converge on a few of these pools; DerSimonian-Laird is
      ## closed-form and cannot fail, so it is the fallback rather than losing
      ## the row.  `method` records which one produced it.
      m <- tryCatch(
          rma(yi = sub_df$est, sei = sub_df$se, control = list(maxiter = 10000)),
          error = function(e)
              rma(yi = sub_df$est, sei = sub_df$se, method = "DL"))
      data.frame(est = m$b[1, 1], se = m$se, ci_low = m$ci.lb, ci_high = m$ci.ub,
                 effect = eff, model = mod, class_id = "Summary",
                 pval = m$pval, tau2 = m$tau2, QEp = m$QEp, k = m$k,
                 method = m$method,
                 stringsAsFactors = FALSE)
    }))
  }))
}

plot_meta <- function(df_long, meta_results, plot_title, file_suffix) {
  if (is.null(meta_results) || nrow(meta_results) == 0L) {
    message("plot_meta: nothing to plot for '", file_suffix, "'")
    return(invisible(NULL))
  }
  if (!"ci_low" %in% names(df_long)) {
    df_long$ci_low  <- df_long$est - 1.96 * df_long$se
    df_long$ci_high <- df_long$est + 1.96 * df_long$se
  }

  original_classes <- sort(unique(as.character(df_long$class_id)))
  visual_order     <- c("Summary", "", rev(original_classes))
  df_long$class_id       <- factor(df_long$class_id, levels = visual_order)
  meta_results$class_id  <- factor("Summary", levels = visual_order)

  p <- ggplot() +
    geom_point(data = df_long,
               aes(x = est, y = class_id, color = model),
               position = position_dodge(width = 0.5), size = 1.5) +
    geom_errorbarh(data = df_long,
                   aes(xmin = ci_low, xmax = ci_high, y = class_id, color = model),
                   height = 0, position = position_dodge(width = 0.5)) +
    geom_point(data = meta_results,
               aes(x = est, y = class_id, color = model),
               shape = 18, size = 5, position = position_dodge(width = 0.5)) +
    geom_errorbarh(data = meta_results,
                   aes(xmin = ci_low, xmax = ci_high, y = class_id, color = model),
                   height = 0.3, position = position_dodge(width = 0.5), lwd = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 2, linetype = "solid", color = "black", linewidth = 0.8) +
    facet_wrap(~effect, ncol = 4, scales = "free_x", drop = TRUE) +
    scale_y_discrete(limits = visual_order) +
    scale_color_manual(values = MODEL_COLOURS) +
    labs(x = "Estimate", y = "School class", title = plot_title) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "right")

  ggsave(file.path(OUT, paste0("plot_", file_suffix, ".png")), p,
         width = 16, height = 10, bg = "white")
  invisible(p)
}

# RUN -------------------------------------------------------------------------

ESTIMANDS <- intersect(c("static", "dynamic", "accumulated"),
                       unique(postest_all$estimand))

## `postest_all` holds each estimand twice -- pooled over all ministeps and
## split by density.  The tables and plots below use the POOLED rows; the
## density split gets its own section at the foot of the script.
postest_pooled <- subset(postest_all, condition == "pooled")

metas <- list(raw = calculate_meta(raw_df))
for (e in ESTIMANDS)
  metas[[e]] <- calculate_meta(subset(postest_pooled, estimand == e))

for (nm in names(metas)) {
  write.csv(metas[[nm]], file.path(OUT, sprintf("table_%s_meta.csv", nm)),
            row.names = FALSE)
  cat(sprintf("%-12s %3d meta-analyses | k: %s\n", nm, nrow(metas[[nm]]),
              paste(sort(unique(metas[[nm]]$k)), collapse = ", ")))
}

plot_meta(raw_df, metas$raw, "Raw Coefficients", "raw_all")
for (e in ESTIMANDS)
  plot_meta(subset(postest_pooled, estimand == e), metas[[e]],
            paste0(tools::toTitleCase(e), " Effects"), paste0(e, "_all"))

# FILTERED SUBSETS ------------------------------------------------------------
#
# The endogenous interaction is called `unspInt` among the raw coefficients and
# `interaction` among the marginal effects; the two greps below reflect that,
# not a difference in the quantity.

sub_plot <- function(df, meta, pattern, negate, title, suffix) {
  keep_d <- grepl(pattern, df$effect)
  keep_m <- grepl(pattern, meta$effect)
  if (negate) { keep_d <- !keep_d; keep_m <- !keep_m }
  plot_meta(df[keep_d, , drop = FALSE], meta[keep_m, , drop = FALSE],
            title, suffix)
}

sub_plot(raw_df, metas$raw, "unspInt", FALSE,
         "TransTrip1 x Rec Only", "raw_transTrip1Rec_only")
for (e in ESTIMANDS)
  sub_plot(subset(postest_pooled, estimand == e), metas[[e]], "interaction",
           FALSE, "TransTrip1 x Rec Only", paste0(e, "_transTrip1Rec_only"))

sub_plot(raw_df, metas$raw, "rate|rec|trans|density|unspInt", TRUE,
         "Covariate Effects Only", "raw_cov_only")
for (e in ESTIMANDS)
  sub_plot(subset(postest_pooled, estimand == e), metas[[e]],
           "rate|rec|trans|density|interaction", TRUE,
           "Covariate Effects Only", paste0(e, "_cov_only"))

sub_plot(raw_df, metas$raw, "rec|trans|density|unspInt", FALSE,
         "Network Effects Only", "raw_net_only")
sub_plot(raw_df, metas$raw, "rate", FALSE, "Rate Only", "raw_rate_only")


# CREATION VERSUS MAINTENANCE --------------------------------------------------
#
# Density is +1 when an actor is deciding whether to CREATE a tie and -1 when
# deciding whether to MAINTAIN one.  Those are different decisions and often
# have different answers, so a pooled average over both is an average over two
# populations.  Everything above is repeated here split that way.

by_density <- subset(postest_all, condition == "bydensity" & stratum != 0)
by_density$context <- factor(ifelse(by_density$stratum > 0,
                                    "creation", "maintenance"),
                             levels = c("creation", "maintenance"))

meta_context <- do.call(rbind, lapply(ESTIMANDS, function(e) {
  d <- subset(by_density, estimand == e)
  do.call(rbind, lapply(levels(d$context), function(ctx) {
    m <- calculate_meta(subset(d, context == ctx))
    if (is.null(m)) return(NULL)
    m$estimand <- e; m$context <- ctx
    m
  }))
}))

if (!is.null(meta_context) && nrow(meta_context)) {
  write.csv(meta_context, file.path(OUT, "table_context_meta.csv"),
            row.names = FALSE)
  cat(sprintf("\nby density  %3d meta-analyses | k: %s\n", nrow(meta_context),
              paste(sort(unique(meta_context$k)), collapse = ", ")))

  ## MODEL gets the colour, so the models sit next to each other within an
  ## effect: whether creation and maintenance differ, and whether that survives
  ## controlling for more, are both read ACROSS models.
  ##
  ## Estimand goes in the columns so each one gets its own estimate scale.
  ## Accumulated sums over the ministeps of a period and is an order of
  ## magnitude larger than the per-opportunity estimands; on a shared scale it
  ## would flatten static and dynamic onto the zero line.
  meta_context$model <- factor(meta_context$model, levels = MODEL_LEVELS)

  p_ctx <- ggplot(meta_context,
                  aes(x = est, xmin = ci_low, xmax = ci_high, y = effect,
                      colour = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.3,
                    orientation = "y") +
    facet_grid(context ~ estimand, scales = "free_x") +
    scale_colour_manual(values = MODEL_COLOURS, breaks = MODEL_LEVELS,
                        drop = FALSE) +
    labs(x = "pooled marginal effect (95% CI)", y = NULL, colour = NULL,
         title = "Tie creation versus tie maintenance, by model") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  ggsave(file.path(OUT, "plot_context_forest.png"), p_ctx,
         width = 16, height = 9, bg = "white")

  ## The same numbers as a DIFFERENCE -- the subtraction the plot above asks the
  ## eye to do.  Zero means the mechanism does the same work in both decisions.
  wide_ctx <- merge(
      subset(meta_context, context == "creation",
             select = c("effect", "model", "estimand", "est", "se")),
      subset(meta_context, context == "maintenance",
             select = c("effect", "model", "estimand", "est", "se")),
      by = c("effect", "model", "estimand"), suffixes = c("_cr", "_mt"))
  ## Treated as independent, which they are not -- both rest on the same fitted
  ## thetas.  The interval is for reading, not for testing.
  wide_ctx$diff    <- wide_ctx$est_cr - wide_ctx$est_mt
  wide_ctx$diff_se <- sqrt(wide_ctx$se_cr^2 + wide_ctx$se_mt^2)
  wide_ctx$model   <- factor(wide_ctx$model, levels = MODEL_LEVELS)
  write.csv(wide_ctx, file.path(OUT, "table_context_difference.csv"),
            row.names = FALSE)

  p_ctxd <- ggplot(wide_ctx,
                   aes(x = diff, xmin = diff - 1.96 * diff_se,
                       xmax = diff + 1.96 * diff_se, y = effect,
                       colour = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.3,
                    orientation = "y") +
    facet_wrap(~ estimand, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = MODEL_COLOURS, breaks = MODEL_LEVELS,
                        drop = FALSE) +
    labs(x = "creation - maintenance (95% CI)", y = NULL, colour = NULL,
         title = "How differently does each mechanism act on the two decisions?") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  ggsave(file.path(OUT, "plot_context_difference.png"), p_ctxd,
         width = 16, height = 7, bg = "white")
}


# THE MODEL-IMPLIED DESCRIPTIVE COUNTERPART -----------------------------------
#
# Postestimation.R computes, per class, the SAME comparisons as the
# marginal effects but without holding the rest of the decision fixed.  Pooled
# here with the same estimator, they are the confounded version of the adjusted
# numbers above.
#
# This is the comparison to use for "how much of the raw association is the
# mechanism": same model, same ministep, same probability scale, so only the
# adjustment differs.  The wave-to-wave descriptive at the foot of this script
# is a different thing -- it also sums over a whole period, so a gap between it
# and a marginal effect cannot be pinned on either difference.

descr_file <- file.path(WORK, "postest_descr_df.rds")
if (file.exists(descr_file)) {
  descr <- readRDS(descr_file)

  ## A contrast is undefined where a class has no decisions at one of the two
  ## levels; those rows are dropped, so `k` records how many classes
  ## contributed.
  descr <- descr[!is.na(descr$est) & !is.na(descr$se) & descr$se > 0, ]

  meta_descr <- calculate_meta(subset(descr, condition == "pooled"))
  write.csv(meta_descr, file.path(OUT, "table_descriptive_model_meta.csv"),
            row.names = FALSE)
  cat(sprintf("\nmodel-implied descriptive  %3d meta-analyses | k: %s\n",
              nrow(meta_descr), paste(sort(unique(meta_descr$k)), collapse = ", ")))

  ## Side by side with the MATCHED adjusted estimate: the marginal effect
  ## averaged over exactly the egos the descriptive could use.  It differs
  ## slightly from the headline marginal effect in the tables above, on purpose.
  ##
  ## No ratio column.  The two are on the same scale, so their DIFFERENCE is the
  ## quantity that means something; a ratio would read like a "share explained",
  ## which it is not -- either number can cross zero.
  descr_adj <- subset(descr, condition == "pooled" &
                             !is.na(est_adj) & !is.na(se_adj) & se_adj > 0)
  descr_adj$est <- descr_adj$est_adj
  descr_adj$se  <- descr_adj$se_adj
  meta_adj <- calculate_meta(descr_adj)
  if (!is.null(meta_adj)) {
    conf <- merge(
        data.frame(effect = meta_descr$effect, model = meta_descr$model,
                   descriptive = meta_descr$est, descriptive_se = meta_descr$se,
                   k_descriptive = meta_descr$k, stringsAsFactors = FALSE),
        data.frame(effect = meta_adj$effect, model = meta_adj$model,
                   adjusted = meta_adj$est, adjusted_se = meta_adj$se,
                   k_adjusted = meta_adj$k, stringsAsFactors = FALSE),
        by = c("effect", "model"))
    ## Scope carried through so a reader can see which rows are within-ego.
    sc <- unique(descr[, c("effect", "scope")])
    conf <- merge(conf, sc, by = "effect", all.x = TRUE)
    conf$confounding <- conf$descriptive - conf$adjusted
    conf <- conf[order(conf$effect, conf$model), ]
    write.csv(conf, file.path(OUT, "table_confounding.csv"),
              row.names = FALSE)

    cat("\n--- confounded vs adjusted, same model and ministep (static) ---\n")
    cat("    NOTE: most of this gap is the tie's own current state, not covariate\n",
        "   confounding.  See the density split below before quoting it.\n", sep = "")
    for (i in seq_len(nrow(conf))) with(conf[i, ],
        cat(sprintf("  %-22s %-8s %-10s descriptive %+.4f   adjusted %+.4f   difference %+.4f\n",
                    effect, model, scope, descriptive, adjusted, confounding)))

    ## Paired, so the shrinkage under adjustment is one segment per effect and
    ## model rather than two numbers to hunt for in separate panels.
    long <- rbind(
        data.frame(effect = conf$effect, model = conf$model,
                   quantity = "descriptive (confounded)", est = conf$descriptive,
                   se = conf$descriptive_se, stringsAsFactors = FALSE),
        data.frame(effect = conf$effect, model = conf$model,
                   quantity = "adjusted (marginal effect)", est = conf$adjusted,
                   se = conf$adjusted_se, stringsAsFactors = FALSE))
    long$quantity <- factor(long$quantity,
                            levels = c("descriptive (confounded)",
                                       "adjusted (marginal effect)"))
    long$model <- factor(long$model, levels = MODEL_LEVELS)

    p_conf <- ggplot(long, aes(x = est, y = model, colour = quantity)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_linerange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                     position = position_dodge(width = 0.5),
                     orientation = "y") +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      facet_wrap(~ effect, scales = "free_x", ncol = 4) +
      scale_colour_manual(values = c("descriptive (confounded)"   = "#d73027",
                                     "adjusted (marginal effect)" = "#2166ac")) +
      labs(x = "pooled estimate (95% CI)", y = NULL, colour = NULL,
           title = "What survives adjustment",
           subtitle = paste("Same model, same ministep, same probability scale;",
                            "the two differ only in whether the rest of the",
                            "decision is held fixed")) +
      theme_minimal() +
      theme(panel.grid.major.y = element_blank(), legend.position = "top")
    ggsave(file.path(OUT, "plot_confounding.png"), p_conf,
           width = 16, height = 9, bg = "white")
  }

  ## Split by creation and maintenance, so the question becomes whether the
  ## confounding itself differs between the two decisions.
  descr_ctx <- subset(descr, condition == "bydensity" & stratum != 0)
  if (nrow(descr_ctx)) {
    descr_ctx$context <- ifelse(descr_ctx$stratum > 0, "creation", "maintenance")
    meta_descr_ctx <- do.call(rbind, lapply(unique(descr_ctx$context), function(ctx) {
      m <- calculate_meta(subset(descr_ctx, context == ctx))
      if (is.null(m)) return(NULL)
      m$context <- ctx
      m
    }))
    write.csv(meta_descr_ctx,
              file.path(OUT, "table_descriptive_model_context_meta.csv"),
              row.names = FALSE)
    cat(sprintf("model-implied descriptive by density  %3d meta-analyses | k: %s\n",
                nrow(meta_descr_ctx),
                paste(sort(unique(meta_descr_ctx$k)), collapse = ", ")))

    ## READ THIS BEFORE THE POOLED TABLE ABOVE.
    ##
    ## The pooled gap is real, but it is almost entirely ONE confounder: whether
    ## the tie already exists.  Reciprocation and two-paths both go together
    ## with being tied already, and under `tieProb` an existing tie is
    ## near-certain to survive a ministep while an absent one is near-certain
    ## not to appear.  So the pooled descriptive largely measures "already
    ## tied", not the mechanism.  The marginal effect has no such problem: it
    ## moves one statistic on one dyad and leaves density where it was.
    ##
    ## Split by density, the two nearly coincide -- within creation the gap for
    ## transTrip1 is ~0.0004, against 0.196 pooled.  Where a gap DOES survive
    ## the split, that is the interesting part: recip in maintenance runs
    ## descriptive below adjusted.
    ## Matched adjusted again, per context.
    da <- subset(descr_ctx, !is.na(est_adj) & !is.na(se_adj) & se_adj > 0)
    da$est <- da$est_adj; da$se <- da$se_adj
    adj_ctx <- do.call(rbind, lapply(unique(da$context), function(ctx) {
      m <- calculate_meta(subset(da, context == ctx))
      if (is.null(m)) return(NULL)
      m$context <- ctx
      m
    }))
    if (!is.null(adj_ctx) && nrow(adj_ctx)) {
      conf_ctx <- merge(
          data.frame(effect = meta_descr_ctx$effect, model = meta_descr_ctx$model,
                     context = meta_descr_ctx$context,
                     descriptive = meta_descr_ctx$est,
                     descriptive_se = meta_descr_ctx$se,
                     k_descriptive = meta_descr_ctx$k, stringsAsFactors = FALSE),
          data.frame(effect = adj_ctx$effect, model = adj_ctx$model,
                     context = adj_ctx$context, adjusted = adj_ctx$est,
                     adjusted_se = adj_ctx$se, k_adjusted = adj_ctx$k,
                     stringsAsFactors = FALSE),
          by = c("effect", "model", "context"))
      conf_ctx$confounding <- conf_ctx$descriptive - conf_ctx$adjusted
      conf_ctx <- conf_ctx[order(conf_ctx$effect, conf_ctx$model,
                                 conf_ctx$context), ]
      write.csv(conf_ctx, file.path(OUT, "table_confounding_context.csv"),
                row.names = FALSE)

      ## How much of the pooled gap survives the density split.  Near zero
      ## means the pooled gap was the density confound and nothing else.
      if (exists("conf")) {
        shrink <- merge(conf[, c("effect", "model", "confounding")],
                        aggregate(list(confounding_split = abs(conf_ctx$confounding)),
                                  by = list(effect = conf_ctx$effect,
                                            model = conf_ctx$model), FUN = max),
                        by = c("effect", "model"))
        cat("\n--- how much of the pooled gap survives the density split ---\n")
        for (i in seq_len(nrow(shrink))) with(shrink[i, ],
            cat(sprintf("  %-22s %-8s pooled gap %+.4f   largest within-stratum gap %+.4f\n",
                        effect, model, confounding, confounding_split)))
      }

      long_ctx <- rbind(
          data.frame(effect = conf_ctx$effect, model = conf_ctx$model,
                     context = conf_ctx$context,
                     quantity = "descriptive (confounded)",
                     est = conf_ctx$descriptive, se = conf_ctx$descriptive_se,
                     stringsAsFactors = FALSE),
          data.frame(effect = conf_ctx$effect, model = conf_ctx$model,
                     context = conf_ctx$context,
                     quantity = "adjusted (marginal effect)",
                     est = conf_ctx$adjusted, se = conf_ctx$adjusted_se,
                     stringsAsFactors = FALSE))
      long_ctx$quantity <- factor(long_ctx$quantity,
                                  levels = c("descriptive (confounded)",
                                             "adjusted (marginal effect)"))
      long_ctx$model <- factor(long_ctx$model, levels = MODEL_LEVELS)

      p_conf_ctx <- ggplot(long_ctx, aes(x = est, y = model, colour = quantity)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
        geom_linerange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                       position = position_dodge(width = 0.5),
                       orientation = "y") +
        geom_point(position = position_dodge(width = 0.5), size = 2) +
        facet_grid(context ~ effect, scales = "free_x") +
        scale_colour_manual(values = c("descriptive (confounded)"   = "#d73027",
                                       "adjusted (marginal effect)" = "#2166ac")) +
        labs(x = "pooled estimate (95% CI)", y = NULL, colour = NULL,
             title = "What survives adjustment, within tie creation and tie maintenance",
             subtitle = paste("The pooled version of this comparison is dominated by",
                              "whether the tie already exists; splitting by density",
                              "removes that confounder")) +
        theme_minimal() +
        theme(panel.grid.major.y = element_blank(), legend.position = "top")
      ggsave(file.path(OUT, "plot_confounding_context.png"), p_conf_ctx,
             width = 18, height = 9, bg = "white")
    }
  }
} else {
  message("postest_descr_df.rds not found -- run Postestimation.R ",
          "with RUN_DESCRIPTIVE=TRUE for the confounded/adjusted comparison.")
}


# THE UNCONTROLLED COUNTERPART -------------------------------------------------
#
# scripts/descriptives.R computes, per class, the difference in tie proportions at
# t2 between dyads where a t1 condition holds and dyads where it does not.  It
# needs no model and confounds every mechanism at once.  Pooled here with the
# same estimator, it is the uncontrolled version of the same comparison, on the
# same probability scale as the accumulated marginal effect.
#
# The descriptive carries no standard error -- dyads within a class are not
# independent -- so it is summarised as the mean across classes with its range.
# The comparison below is for reading, not for testing, and there is
# deliberately no ratio column; see the note at the merge.

desc_file <- file.path(WORK, "descriptive_diffs_by_class.rds")
if (file.exists(desc_file)) {
  desc <- readRDS(desc_file)

  ## No pooling model for the descriptive: dyads within a class are not
  ## independent, so there is no honest per-class standard error to weight by.
  ## The simple mean across classes, with the range, is what it supports.
  meta_desc <- do.call(rbind, lapply(split(desc, desc$target), function(d) {
    d <- d[!is.na(d$est), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    data.frame(effect = d$target[1L], condition = d$condition[1L],
               descriptive = mean(d$est),
               descriptive_min = min(d$est), descriptive_max = max(d$est),
               n_classes = nrow(d), stringsAsFactors = FALSE)
  }))
  write.csv(meta_desc, file.path(OUT, "table_descriptive_meta.csv"),
            row.names = FALSE)

  adjusted <- do.call(rbind, lapply(ESTIMANDS, function(e) {
    m <- metas[[e]]
    if (is.null(m)) return(NULL)
    data.frame(effect = m$effect, model = m$model, estimand = e,
               adjusted = m$est, adjusted_se = m$se, stringsAsFactors = FALSE)
  }))

  ## NO RATIO COLUMN, deliberately.  The accumulated effect is a SUM over the
  ## ministeps of a period and is very nearly lambda times the per-opportunity
  ## effect (verified at ~1% in this model), so accumulated / descriptive
  ## mostly measures the rate parameter.  The descriptive is a difference of
  ## two proportions, bounded in [-1, 1]; the accumulated is a cumulative
  ## score with no such bound.  Dividing one by the other does not yield a
  ## "share of the raw association explained", and reporting it as though it
  ## did invites exactly that reading.
  comparison <- merge(adjusted, meta_desc, by = "effect")
  comparison <- comparison[order(comparison$effect, comparison$model,
                                 comparison$estimand), ]
  write.csv(comparison, file.path(OUT, "table_descriptive_vs_adjusted.csv"),
            row.names = FALSE)
  cat("\n--- uncontrolled vs adjusted (model_2) ---\n")
  show <- subset(comparison, estimand == "accumulated" & model == "model_2")
  if (!nrow(show)) {
    show <- subset(comparison, model == "model_2")
    cat("  NOTE: accumulated results not available; showing ",
        paste(unique(show$estimand), collapse = "/"),
        ", which is per decision opportunity rather than per period.\n", sep = "")
  }
  if (nrow(show)) for (i in seq_len(nrow(show))) with(show[i, ],
      cat(sprintf("  %-22s %-28s descriptive %+.3f [%+.3f,%+.3f]   %-11s %+.3f\n",
                  effect, condition, descriptive, descriptive_min,
                  descriptive_max, estimand, adjusted)))
} else {
  message("descriptive_diffs_by_class.rds not found -- ",
          "run scripts/descriptives.R first for the uncontrolled comparison.")
}


# SANITY ----------------------------------------------------------------------
#
# Every pooled estimate must rest on one observation per class.  A k above
# the class count means the input frame carries duplicated rows, which would
# shrink every standard error and bias tau^2.

n_classes <- length(unique(postest_all$class_id))
all_metas <- do.call(rbind, lapply(metas[ESTIMANDS], function(m) m[, "k", drop = FALSE]))
## The descriptive pools can legitimately sit BELOW the class count -- a class
## with no decisions at one of the two levels drops out -- but never above it.
if (exists("meta_descr") && !is.null(meta_descr))
  all_metas <- rbind(all_metas, meta_descr[, "k", drop = FALSE])
if (!is.null(meta_context))
  all_metas <- rbind(all_metas, meta_context[, "k", drop = FALSE])
if (any(all_metas$k > n_classes))
  stop("k exceeds the number of classes (", n_classes,
       ") -- the input frame still contains duplicated rows.", call. = FALSE)
cat(sprintf("\nOK: %d classes; no meta-analysis pools more than %d studies.\n",
            n_classes, n_classes))
