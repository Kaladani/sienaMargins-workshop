# Figures for the workshop book

Copy the curated figures here. Chapters reference them by these names.

| file | used in | source |
|---|---|---|
| `two_alternatives.svg` | ch. 1, what a coefficient does answer | Christian Steglich, workshop page (`slides/MEworkshopEUSN2026/Figure_two_alternatives.svg`) |
| `choice_recip_shift_all_alts.pdf` / `.png` | ch. 2, shared denominator | paper repo `plots/choice_recip_shift_all_alts_sigma1.pdf` |
| `choice_egoX_masses.pdf` / `.png` | ch. 2, ego aside | paper repo `plots/choice_egoX_masses_sigma1.pdf` |
| `interaction_sign_swap.png` | ch. 2, interactions as second differences | Christian Steglich, workshop page (`slides/MEworkshopEUSN2026/interaction.png`) |

Both raster and vector are worth having where a figure is used in both output
formats: HTML cannot embed PDF, PDF output prefers vector. Quarto picks per format
automatically when the extension is omitted in the `![](figs/name)` reference — the
two files above that carry an explicit extension exist in one format only.

The logistic-curve figure in ch. 1 is generated in the chapter itself rather than
copied, so no file is needed for it.
