# Construction of the East / West migration variables in BASE-II

This document describes how the three regional / migration variables
(`east_west`, `east_west_mig`, `east_west_mig_edu`) used in the GxE
analyses are derived from the BASE-II data delivery. The
implementation is in §6 of
`02_Phenotype/BASE-II/EduGxE_DE_BASE_DataMining.R`.

## Scope: older cohort only

All three regional variables are populated for the **older BASE-II
cohort only** (`agegr == 2`, the 60–80 year-old subsample at
baseline). For the younger cohort (`agegr == 1`, 20–35 year-olds at
baseline) all three columns are left as `NA`. No informative
pre-1989 signal — neither a GDR-residence self-report nor a
schooling-system indicator — is available for participants who were
born after reunification or grew up entirely under the unified
school system.

## Source variables

| Variable | Type | Meaning |
|---|---|---|
| `pddr01` | binary self-report | "Lived ≥ 1 year in the GDR" (1 = yes, 2 = no). Asked of the older cohort only. |
| `lsab3` | categorical | Country of last school visit: 1 = FRG, 2 = GDR. |
| `lsab7x` / `lsab7y` | per-respondent | Final school certificate received in the FRG (`lsab7x`) or GDR (`lsab7y`). Used as presence/absence flags (`!is.na(...)`). |
| `pwohnvb02 … pwohnvb12` | per-wave categorical | Berlin borough of residence in each annual wave 2002–2012. Coded by district. |
| `agegr` | categorical | BASE-II age-group flag: 1 = younger cohort, 2 = older cohort. |

## Derived intermediate signals

Two pairs of derived signals are built before the variables themselves.

**Schooling system.** A respondent is flagged as having attended the
GDR school system if they hold a GDR certificate (`lsab7y` not
missing) OR their last school visit was in the GDR (`lsab3 == 2`).
Symmetrically for the FRG.

```text
east_school = !is.na(lsab7y) | (lsab3 == 2)
west_school = !is.na(lsab7x) | (lsab3 == 1)
```

**Pre-1989 origin (older cohort only).**

```text
east_origin = (agegr == 2) & (pddr01 == 1)
west_origin = (agegr == 2) & (pddr01 == 2)
```

**Berlin borough of first residence (2002–2012).** Each wave's
borough code is classified as East-Berlin, West-Berlin, or neutral.
The first non-neutral, non-missing observation across the eleven
waves is taken as the respondent's first observed Berlin residence.

| Borough code | District | Classification |
|---|---|---|
| 1 | Mitte | Neutral (excluded — split during division) |
| 2 | Friedrichshain-Kreuzberg | Neutral (excluded — split during division) |
| 3 | Pankow | East |
| 9 | Treptow-Köpenick | East |
| 10 | Marzahn-Hellersdorf | East |
| 11 | Lichtenberg | East |
| 4 | Charlottenburg-Wilmersdorf | West |
| 5 | Spandau | West |
| 6 | Steglitz-Zehlendorf | West |
| 7 | Tempelhof-Schöneberg | West |
| 8 | Neukölln | West |
| 12 | Reinickendorf | West |
| −3, −2, −1 | (missing codes) | NA |

Respondents whose only Berlin observations fall in the neutral
districts (Mitte, Friedrichshain-Kreuzberg) receive `FALSE` on both
`first_res_east_ber` and `first_res_west_ber` and are therefore
unobserved on the Berlin axis. The neutral classification reflects
that these districts were administratively split before 1989 and
cannot be cleanly assigned to either side.

## `east_west` — region of origin

Built only for the older cohort. The self-reported GDR-residence
item (`pddr01`) is the primary source; the schooling system is the
secondary source used only when `pddr01` is missing. The schooling
fallback requires the respondent to have an unambiguous one-sided
schooling signal (East-only or West-only); ambiguous or absent
schooling evidence leaves the variable as `NA`.

| `agegr` | `pddr01` | `east_school` | `west_school` | → `east_west` |
|---|---|---|---|---|
| 2 (older) | 1 | any | any | east |
| 2 (older) | 2 | any | any | west |
| 2 (older) | NA | TRUE | FALSE | east |
| 2 (older) | NA | FALSE | TRUE | west |
| 2 (older) | NA | TRUE | TRUE | NA |
| 2 (older) | NA | FALSE | FALSE | NA |
| 1 (younger) | any | any | any | NA |

## `east_west_mig` — lifetime migration

Lifetime migration. A respondent in the older cohort is coded as a
mover if either contrary signal — schooling in the other system, OR
first observed Berlin residence on the other side — is present.
Otherwise they are coded as a stayer in their origin region.

Stayer status does **not** require positive evidence of staying on
the schooling or Berlin axes; absence of contrary signals is
sufficient.

| `agegr` | `east_origin` | `west_origin` | `west_school` OR `first_res_west_ber` | `east_school` OR `first_res_east_ber` | → `east_west_mig` |
|---|---|---|---|---|---|
| 2 | TRUE | — | TRUE | — | east_to_w_move |
| 2 | TRUE | — | FALSE | — | east_stay |
| 2 | — | TRUE | — | TRUE | west_to_e_move |
| 2 | — | TRUE | — | FALSE | west_stay |
| 2 | FALSE | FALSE | — | — | NA |
| 1 | — | — | — | — | NA |

## `east_west_mig_edu` — migration during education

Same conceptual structure as `east_west_mig`, but uses **only the
schooling system** as the contrary-signal axis. Berlin residence
between 2002 and 2012 (i.e. adult residence) is deliberately ignored
because it occurred after the formative-education window.

Stayer status here **does** require positive matching-side schooling
evidence: the variable is informative about migration during
education, so the absence of any schooling information at all leaves
the respondent unclassified.

| `agegr` | `east_origin` | `west_origin` | `east_school` | `west_school` | → `east_west_mig_edu` |
|---|---|---|---|---|---|
| 2 | TRUE | — | — | TRUE | east_to_w_move |
| 2 | TRUE | — | TRUE | — | east_stay |
| 2 | — | TRUE | TRUE | — | west_to_e_move |
| 2 | — | TRUE | — | TRUE | west_stay |
| 2 | TRUE | — | FALSE | FALSE | NA |
| 2 | — | TRUE | FALSE | FALSE | NA |
| 2 | FALSE | FALSE | — | — | NA |
| 1 | — | — | — | — | NA |

## Caveats

- **No regional variables for the younger cohort.** Participants
  with `agegr == 1` (20–35 at baseline, born after reunification or
  growing up entirely under unified schooling) are `NA` on all three
  variables. This is by design — no item in the BASE-II delivery
  provides a pre-1989 regional signal for them.
- **Berlin's central districts (Mitte, Friedrichshain-Kreuzberg)
  are deliberately neutral.** These districts spanned the division
  before 1989 and cannot be cleanly assigned to either side. People
  whose Berlin observations fall only in these districts are
  unobserved on the Berlin axis (but can still be classified via
  `pddr01` / schooling).
- **`pddr01` is the GDR-residence anchor.** It asks whether the
  respondent ever lived at least one year in the GDR. The variable
  is binary; it does not distinguish duration, region within the GDR,
  or whether the residence overlapped the formative-education window.
