# Construction of the East / West migration variables in SOEP

This document describes how the three regional / migration variables
(`east_west`, `east_west_mig`, `east_west_mig_edu`) used in the GxE
analyses are derived from the SOEP data delivery. The
implementation is in §6 of
`02_Phenotype/SOEP/EduGxE_DE_SOEP_DataMining.R`. This is the most
elaborate east / west construction across the four cohorts because
SOEP provides per-year household residency from a long-running panel,
which permits reconstructing residence across the educational window
and identifying within-life region changes.

## Source variables

| Variable | Type | Meaning |
|---|---|---|
| `loc1989` | categorical | Self-reported residence in 1989: 1 = East, 2 = West, 3 = other (treated as NA), all other codes treated as NA. Mapped to `loc1989_ew ∈ {"east", "west", NA}`. |
| `sampreg` (in `hbrutto`) | per-year, per-household | Household region code: 1 = West, 2 = East, every negative code (−1 … −6) treated as NA. Mapped to `ew ∈ {"east", "west", NA}` per household × year. |
| `gebjahr` | year | Birth year. |
| `owneduyears` | continuous | Own years of education completed. Used to derive the schooling window and the education-finish year. |

## Derived intermediate signals

The construction uses three derived quantities per person:

- `finish_yr` — estimated education-finish year, computed as `floor(gebjahr + owneduyears + 6)`. Used by the education-window variable.
- `first_ew` / `first_yr` — the household's first non-missing `sampreg`-derived East/West value and the year it was first observed.
- `first_diff_yr` — the first panel year in which the household's `sampreg`-derived East/West differs from the person's `loc1989_ew`. NA if there is never any difference.
- `first_change_yr` / `first_change_from` / `first_change_to` — the first observed flip in the household's `sampreg` East/West across panel years, with the year, the prior region, and the new region.

## `east_west` — region of origin

`east_west` uses a strict **schooling window**: birth year + 15 up to
birth year + `owneduyears` + 7. For each year in that window, the
respondent's region is taken from `loc1989_ew` if the year is before
1989, and from the household's first observed `sampreg` (`first_ew`)
if the year is 1989 or later. The variable is then assigned only
when **all observed years within the schooling window agree** on
either East or West; any mixed-region trajectory or unknown
trajectory leaves it as NA.

| Years in the schooling window | → `east_west` |
|---|---|
| ≥ 1 observed and all are East | east |
| ≥ 1 observed and all are West | west |
| Mixed regions in the window | NA |
| No observed years (`gebjahr` or `owneduyears` missing, or window starts in adulthood) | NA |

The strict "all-agree" rule is intentional: a respondent who is
observed as East in some schooling years and West in others is not
informative for `east_west` (their formative residence is split),
and they are classified instead via `east_west_mig` and
`east_west_mig_edu` below.

## `east_west_mig` — lifetime migration

Lifetime migration is determined by a four-step priority that
combines the 1989 self-report with later panel changes:

1. **Early mismatch.** If `loc1989_ew` is known and `first_ew` is
   known and they disagree, label as a mover in the direction of the
   disagreement (`east_to_w_move` if `loc1989 = "east"` and
   `first_ew = "west"`, vice versa otherwise).
2. **Later panel change.** Else, if there is a first observed
   panel flip (`first_change_yr`) that occurred strictly after the
   first observed year (`first_yr`), label as a mover in the
   direction of that flip (`first_change_from → first_change_to`).
3. **Agreement.** Else, if `loc1989_ew` and `first_ew` are both
   known and agree, classify as `east_stay` or `west_stay`
   accordingly.
4. **Fallback to schooling residence.** Else, fall back to the
   `east_west` classification from the schooling-window rule
   (`east_stay` if `east_west = "east"`, etc.).
5. NA only if none of the above can be applied (rare; requires
   both panel and schooling-window evidence to be entirely missing).

## `east_west_mig_edu` — migration during education

The education-window variant applies four timing rules, indexed by
`finish_yr` (education-finish year), `first_yr` (first panel year),
and `first_diff_yr` (first panel year of region mismatch).

**A. Education finished before 1989.** If `finish_yr < 1989`,
classify by the schooling-residence fallback (`east_west`). Region
changes that occurred post-1989 cannot have influenced education in
this case.

**B. Education extends past 1989, with early mismatch at first
observation.** If `finish_yr > 1989` and `loc1989_ew ≠ first_ew`
and `first_yr == first_diff_yr` (the very first panel observation
is also the first year of mismatch), classify as a mover in the
`loc1989_ew → first_ew` direction.

**C. Late mismatch (panel flip distinct from first observation).**
If `first_yr ≠ first_diff_yr` (the household was first observed in
the `loc1989` region, then later moved to the other), and
`finish_yr < first_diff_yr`, classify as a mover in the panel
direction (`first_change_from → first_change_to`). Otherwise fall
back to `stayer_from_school`.

**D. Otherwise: stayer fallback.** All cases that fail A–C take the
schooling-residence classification (`east_west`).

The fallback `east_west` may itself be NA, in which case
`east_west_mig_edu` is NA.

## Caveats

- **Schooling window is sensitive to `owneduyears`.** A missing
  `owneduyears` or `gebjahr` leaves the window undefined and
  `east_west` becomes NA. Downstream migration variables fall back
  via rules 4 / D and may also be NA.
- **`loc1989 == 3` ("other") is treated as NA.** Respondents who
  selected the "other" 1989-residence category are unobserved on
  the pre-1989 anchor and rely entirely on panel evidence.
- **Rule C is empirically near-empty.** Rule C classifies as a
  mover when education finished strictly **before** the first
  observed panel change in region (`finish_yr < first_diff_yr`):
  a household region change is taken as evidence of a move that
  may have occurred at any earlier point, including during
  education. Only **2 of 2,298** analysis-sample persons
  (≈0.09 %) fall into the rule-C zone at all, so the rule has
  essentially no effect on the classified sample.
