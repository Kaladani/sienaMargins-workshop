# Data

**The data are not distributed with this repository.** This directory is where
they go once you have them; it ships empty apart from this file.

## What you need

Three files, which the code reads:

| file | contents |
|---|---|
| `AK_friendship.RData` | friendship networks, `friendship[[class]][["W"]]` etc. |
| `AK_gender.RData` | `male[[class]]` |
| `AK_primary.RData` | `primary[[class]]` |

The pupil questionnaire (`QuestionnairePupilsWaveV-1.pdf`) is available from the
same place. Nothing in the code reads it, but it documents the instrument the
data come from.

## Getting them

They are hosted alongside the other workshop material:

<https://steglich.gmw.rug.nl/workshops/>

Put all three `.RData` files in this directory. Chapter 4 and
`scripts/estimation.R` will also fetch them for you on first use — both call
`ensure_ak_data()` from `scripts/download_data.R`, which downloads whatever is
missing and leaves anything already present alone.

> The address in `scripts/download_data.R` is currently a placeholder. Until it
> is set, download the files by hand and put them here.

Chapters 1, 3 and 5 do not need them: chapter 3 uses the `s50` data that ships
with RSiena, and chapters 1 and 5 read the precomputed results in `results/`.

## Source and terms

- **Dataset:** Knecht, A. (2004). *Network and actor attributes in early
  adolescence*. DANS. <https://doi.org/10.17026/dans-z9b-h2bp>
- **Codebook:** Knecht, A. (2006). *Networks and actor attributes in early
  adolescence [2003/04]*. ICS-Codebook no. 61, University of Utrecht.
  Persistent identifier `urn:nbn:nl:ui:13-ehzl-c6`.

The questionnaire (*Vragenlijst Scholierenonderzoek*, Projectgroep
Scholierenonderzoek, Universiteit Utrecht; C. Baerveldt and A. Knecht) is a
separate work by those authors.

The DANS deposit is published under **Restricted access**, governed by the
[DANS Licence](https://doi.org/10.17026/fp39-0x58). Two obligations matter to
anyone using these data here:

- **Citation (Art. 2).** Published research using the data must cite the dataset,
  including DANS and the persistent identifier as a full URL. The chapters cite
  `@Knecht2006`; `@KnechtData` in `references.bib` carries the DOI form.
- **Redistribution (Art. 3).** Distributing the dataset or substantial parts of
  it requires prior permission from the rights holder. This is why the files are
  not included here, and why you should not commit them back into the repository
  or pass them on yourself — `.gitignore` covers the first of those, but not the
  second.

Nothing in this repository's licence applies to these data.
