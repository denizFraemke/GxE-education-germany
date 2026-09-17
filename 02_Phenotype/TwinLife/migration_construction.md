# Construction of the East / West migration variables in TwinLife

This document describes how the three regional / migration variables
(`east_west`, `east_west_mig`, `east_west_mig_edu`) used in the GxE
analyses are derived from the TwinLife data delivery. The
implementation is in §6 of
`02_Phenotype/TwinLife/EduGxE_DE_TwinLife_DataMining.R`.

## Source variables

| Variable | Type | Meaning |
|---|---|---|
| `mig2001` | binary, time-constant | "Born in the GDR" (1 = yes, 0 = no). Defined only for births before German reunification (1990); from 1990 onward the GDR no longer existed and the question has no defined meaning. |
| `ewi_wid1..wid9` | per-wave categorical | Current federal state at the time of each survey wave: 1 = Eastern Bundesland, 2 = Western Bundesland. |
| `birth` | continuous | Birth year + birth month / 12, derived in §3.1 of the data-mining script. Renamed `birth_year` in the final dataset. |

These regions are normalised to bare `"east"` / `"west"` via the
`norm_region()` helper. `mig2001` becomes `born_reg`; `ewi_wid*` is
collapsed to the first observed wave (`first_ewi`) and becomes
`live_reg`.

**Caveat on Berlin.** The TwinLife residential indicator `ewi`
classifies all of Berlin as part of the Western Bundesländer. This is
a substantive simplification, since Berlin (especially East Berlin)
was a major part — and a major educational hub — of the GDR. Treatment
of Berlin in the residential indicator therefore introduces a known
asymmetry between the birth indicator (`mig2001`, which counts East
Berliners as born in the GDR) and the residential indicator (`ewi`,
which counts them as currently in a Western federal state).

The asymmetry **does not affect `east_west`** for pre-1990 born with
known `mig2001`: `born_reg` wins, so East Berliners are correctly
labelled `"east"`. It does affect the two migration variables — pre-1990
East Berliners are coded as `east_to_w_move` even when they stayed in
Berlin — and it also affects `east_west` for the post-1990 cohort,
which falls back to `live_reg`. The upper bound on the migration
artifact in the current TwinLife analysis sample is ~71 `east_to_w_move`
cases (combined real-mover + East-Berlin-stayer counts); the actual
artifact share is expected to be much smaller. No correction is
applied; the residual is accepted.

## `east_west` — region of formative attachment

`east_west` is a single best-guess assignment of an individual to
"east" or "west". It uses the birth region where that is meaningful,
and falls back to the first observed residential region otherwise.

`mig2001` evidence is restricted to births that pre-date
reunification. `mig2001` is recorded as `0` for nearly all post-1990
births — the question "Born in the GDR" cannot be answered "yes"
after 1990 — so treating it as the primary source regardless of
birth era would label post-1990 persons `"west"` even when they
currently live in Eastern federal states. Residual misclassification
(recall error, the small implausible group of post-1990
`mig2001 = 1` cases, etc.) is accepted; no separate sensitivity
analysis is performed.

### Who ends up in which group, when

| Birth year | `mig2001` | `first_ewi` | → `east_west` |
|---|---|---|---|
| ≤ 1989 | 1 | any | east |
| ≤ 1989 | 0 | any | west |
| ≤ 1989 | NA | 1 | east |
| ≤ 1989 | NA | 2 | west |
| ≤ 1989 | NA | NA | NA |
| ≥ 1990 | any | 1 | east |
| ≥ 1990 | any | 2 | west |
| ≥ 1990 | any | NA | NA |

Rows for `birth_year = NA` (about 5,900 persons in the raw delivery,
mostly outside the analysis sample) are not shown — `mig2001` keeps
precedence for them because the pre-/post-1990 split cannot be
determined.

## `east_west_mig` — lifetime migration relative to the GDR / FRG divide

`east_west_mig` captures whether a person stayed in their birth region
or crossed the GDR / FRG divide at some point in their observed life.
It compares `born_reg` (birth region, from `mig2001`) and `live_reg`
(first observed residential region, from `first_ewi`).

"Stayer" / "mover" labels are assigned only when both `born_reg` and
`live_reg` are known. Anyone whose birth region cannot be determined
is left as `NA` — including every post-1990 birth, whose `born_reg`
is `NA` by construction under the rule above. This is substantively
correct: cross-GDR/FRG migration cannot be defined for people born
after the GDR ceased to exist.

### Construction rule

| `born_reg` | `live_reg` | → `east_west_mig` |
|---|---|---|
| east | east | east_stay |
| west | west | west_stay |
| east | west | east_to_w_move |
| west | east | west_to_e_move |
| NA | any | NA |
| any | NA | NA |

## `east_west_mig_edu` — migration observed during the formative window

`east_west_mig_edu` is the formative-years variant of `east_west_mig`:
a move only counts as a "mover" assignment if it was observed at age
≤ 25 (the educational / developmental window, per the analysis plan).
Where the move was first observed at age > 25, the person is treated
as a stayer in their birth region — they did stay in that region
during their formative years, and the later move is irrelevant for
the variable's intended use as an exposure during education.

Age at first observation is computed as
`first_year + first_month / 12 − birth`, where `first_year` and
`first_month` are `yea_fq` / `mon_fq` at the wave that produced
`first_ewi`. If the wave timing is missing, age at first observation
is `NA` and the move cannot be placed relative to age 25 — in that
case the person is left as `NA` on this variable.

### Construction rule

| `born_reg` | `live_reg` | age at first observation | → `east_west_mig_edu` |
|---|---|---|---|
| east | east | any | east_stay |
| west | west | any | west_stay |
| east | west | ≤ 25 | east_to_w_move |
| west | east | ≤ 25 | west_to_e_move |
| east | west | > 25 | east_stay |
| west | east | > 25 | west_stay |
| east | west | NA | NA |
| west | east | NA | NA |
| NA | any | any | NA |
| any | NA | any | NA |
