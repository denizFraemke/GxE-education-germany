# Construction of the East / West migration variables in SHIP

This document describes how the three regional / migration variables
(`east_west`, `east_west_mig`, `east_west_mig_edu`) used in the GxE
analyses are derived from the SHIP data delivery. The
implementation is in §6 of
`02_Phenotype/SHIP/EduGxE_DE_SHIP_DataMining.R`.

## Geographic context

All SHIP participants currently live in **Mecklenburg-Vorpommern**,
a federal state in the former GDR. The reference assumption is
therefore that every SHIP participant is East-German unless there
is explicit evidence that they were West-German pre-1989. SHIP
delivers two sub-samples — SHIP-Start and SHIP-TREND — which differ
in whether such evidence is available.

## Source variables

| Variable | Sub-sample | Type | Meaning |
|---|---|---|---|
| `wohn89` | SHIP-Start only | categorical | Residence in 1989; numeric coding with `wohn89 == 2` indicating West-Germany. |
| `AGE_SHIP0` | SHIP-Start only | continuous | Age at baseline (SHIP-0). Used to date a `wohn89`-derived West-to-East move relative to age 25. |

SHIP-TREND has neither `wohn89` nor an analogous pre-1989 location
item. Every TREND participant therefore falls through to the
East-default branch on all three variables.

## `east_west` — region of origin

For SHIP-Start, the variable is determined by the 1989 residence
question. SHIP-TREND has no 1989-location data and defaults to East.

| Sub-sample | `wohn89` | → `east_west` |
|---|---|---|
| SHIP-Start | 2 (West-resident 1989) | west |
| SHIP-Start | anything else (incl. NA) | east |
| SHIP-TREND | n/a | east |

The "anything else" branch silently absorbs East-1989 residents,
non-responses, and missing-value codes. The East default is
informative because all SHIP participants live in Mecklenburg-
Vorpommern at the time of survey, so an unobserved 1989 residence
is more likely East than West.

## `east_west_mig` — lifetime migration

A SHIP-Start participant who reports `wohn89 == 2` (West-resident
1989) but currently lives in East-Germany is, by definition, a
West-to-East mover. All other SHIP participants are assumed to have
remained in East-Germany throughout.

| Sub-sample | `wohn89` | → `east_west_mig` |
|---|---|---|
| SHIP-Start | 2 | west_to_e_move |
| SHIP-Start | anything else (incl. NA) | east_stay |
| SHIP-TREND | n/a | east_stay |

## `east_west_mig_edu` — migration during education

Same baseline logic as `east_west_mig`, but with a single timing
override: a SHIP-Start respondent classified as `west_to_e_move`
who was **already ≥ 25 at baseline** is reclassified as `west_stay`
on this variable, on the assumption that a move occurring after age
25 falls outside the formative-education window and that the person
was instead West-resident throughout their education.

`AGE_SHIP0` is required for this rule; participants with missing
age receive `NA` on this variable only.

| Sub-sample | `wohn89` | age at baseline | → `east_west_mig_edu` |
|---|---|---|---|
| SHIP-Start | 2 | ≥ 25 | west_stay |
| SHIP-Start | 2 | < 25 | west_to_e_move |
| SHIP-Start | 2 | NA | NA |
| SHIP-Start | anything else | any | east_stay |
| SHIP-TREND | n/a | n/a | east_stay |

## Caveats

- **TREND has no pre-1989 location data.** Every SHIP-TREND
  participant is classified as `east_stay` on all three variables.
  TREND participants who were in fact West-resident pre-1989 cannot
  be distinguished from genuine East-origin participants from the
  available data.
- **SHIP-Start `wohn89` is the only pre-1989 anchor.** A respondent
  whose `wohn89` is missing falls through to the East default;
  this absorbs an unknown share of West-1989 residents whose 1989
  location was not reported.
- **No East-to-West moves are possible in SHIP.** By construction
  (currently resident in Mecklenburg-Vorpommern, East), the
  migration directions reduce to `west_to_e_move` or `east_stay`;
  the categories `east_to_w_move` and `west_stay` cannot occur in
  the cohort.
