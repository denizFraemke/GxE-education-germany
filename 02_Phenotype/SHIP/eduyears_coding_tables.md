# SHIP — SOEP-style years-of-education coding tables

This document records the SOEP-style harmonisation of years of
education across SHIP's three sub-deliveries (SHIP-Start at baseline
`t0` and follow-up `t1`, and SHIP-TREND at baseline `t0`). The
implementation is in §3–§4 of
`02_Phenotype/SHIP/EduGxE_DE_SHIP_DataMining.R`.

## Harmonisation rule

Harmonisation across SHIP-Start baseline (`t0`) and follow-up (`t1`)
proceeds by linkage on `pid` with the following preference rule:

- The `t0` value (`eduyrs_soep`) is retained.
- If `t0` is `NA`, it is filled from `t1`.
- **Override (a):** respondents aged < 25 years at `t0` are assigned the
  `t1` value to better reflect ongoing educational accumulation.
- **Override (b):** respondents still in education at `t1`
  (`sozio_8 == 1`) receive `NA` for the final measure.

The resulting variable (`eduyrs_soep_harmonized`) thus prioritises
baseline reports, uses follow-up information only when informative,
and avoids assigning completed schooling to participants still
enrolled.

In each per-source table below, the rule is
`eduyrs_soep = school_base + highest_applicable_vocational_or_tertiary_increment`.

## Table S1 — SHIP-Start `t0` (`schule1`, `ausbild*`)

### School base

| Years | Qualification | Source code | Notes |
|---|---|---|---|
| — | Derzeit Schüler:in | `schule1 == 1` | Base NA (may be filled by fallback) |
| 7 | Ohne Abschluss | `schule1 == 2` | |
| 9 | Volks-/Hauptschule | `schule1 == 3` | |
| 10 | Mittlere Reife / Realschule / Fachschule | `schule1 == 4` | |
| 10 | POS | `schule1 == 5` | |
| 13 | Fachschulreife | `schule1 == 6` | |
| 13 | Abitur / EOS (auch mit Facharbeiter) | `schule1 == 7` | |
| 13 | Fachhochschulreife / "Facharbeiter mit Abitur" | `schule1 == 8` | |
| 10 | Sonstiges / keine Angabe | `schule1 == 9` | SOEP convention |
| — | Fallback base = min(school years, 13) | if base NA and `schule_jahre` exists | Applied only if a years-in-school item exists |

### Vocational / tertiary increments

| Increment | Qualification | Source code | Notes |
|---|---|---|---|
| +1.5 | Lehre mit Abschluss | `ausbild4 == 1` | |
| +2 | Fach-/Berufsfachschule / Fachakademie | `ausbild5 == 1` | |
| +3 | FH / Ingenieurschule / Polytechnikum | `ausbild6 == 1` | |
| +5 | Hochschulabschluss (university) | `ausbild7 == 1` | |
| +1.5 | Anderer beruflicher Abschluss | `ausbild8 == 1` | Other training |
| +0 | Anlernzeit / Teilfacharbeiter / kein Abschluss | `ausbild3 == 1` or `ausbild2 == 1` | |
| — | Unknown → NA | 998 / 999 | |

## Table S2 — SHIP-1 (`sozio_*`)

### School base

| Years | Qualification | Source code | Notes |
|---|---|---|---|
| 7 | No general school degree | `sozio_6 == 1` | |
| 9 | Volks-/Hauptschule / POS 8.–9. | `sozio_6 == 2` | |
| 10 | Realschule / POS / "other" school degree | `sozio_6 == 3` or `6` | SOEP: "other" → 10 |
| 13 | Fachhochschulreife / Fachoberschule | `sozio_6 == 4` | |
| 13 | Abitur / EOS | `sozio_6 == 5` | |
| — | Fallback base = min(school years, 13) | if base NA and `sozio_5` known | Always applied; `sozio_5` excludes vocational/tertiary |

### Vocational / tertiary increments

| Increment | Qualification | Source code | Notes |
|---|---|---|---|
| +1.5 | Lehre (apprenticeship) | `sozio_10 == 1` | |
| +2 | Berufsfachschule / Handelsschule | `sozio_12 == 1` | |
| +3 | Meister / Techniker | `sozio_11 == 1` | |
| +3 | Fachschule / Fachhochschule (higher technical) | `sozio_13 == 1` | |
| +5 | University degree | `sozio_13a == 1` | |
| +0 | No vocational degree | `sozio_9 == 1` | |
| — | Unknown → NA | 998 / 999 | Base/increment derived only from valid codes |

## Table S3 — SHIP-TREND `t0` (`t0_sozio_*`)

Harmonisation between SHIP-TREND `t0` and `t1` uses the same rule as
between SHIP-Start `t0` and `t1`.

### School base

| Years | Qualification | Source code | Notes |
|---|---|---|---|
| — | Still a pupil (no final degree) | `t0_sozio_06 == 1` | Base NA (may be filled by fallback) |
| 7 | School exit without degree | `t0_sozio_06 == 2` | |
| 9 | Volks-/Hauptschule | `t0_sozio_06 == 3` | |
| 10 | Mittlere Reife / Realschule / Fachschulreife | `t0_sozio_06 == 4` | |
| 10 | POS (DDR) | `t0_sozio_06 == 5` | |
| 12 | Fachhochschulreife / FOS | `t0_sozio_06 == 6` | |
| 13 | Abitur / EOS | `t0_sozio_06 == 7` | |
| 10 | Other / keine Angabe | `t0_sozio_06 == 8` | SOEP convention |
| — | Fallback base = min(school years, 13) | if base NA and `t0_sozio_05` exists | Applied only if `t0_sozio_05` present |

### Vocational / tertiary increments

| Increment | Qualification | Source code | Notes |
|---|---|---|---|
| +1.5 | Lehre (apprenticeship) | `t0_sozio_08c == 1` | |
| +2 | Berufsfach-/Handelsschule | `t0_sozio_08d == 1` | |
| +3 | Fachschule / Meister / Techniker / Fachakademie | `t0_sozio_08e == 1` | |
| +3 | Fachhochschulabschluss | `t0_sozio_08f == 1` | |
| +5 | Hochschulabschluss (university) | `t0_sozio_08g == 1` | |
| +1.5 | Other vocational | `t0_sozio_08h == 1` | Other training |
| +0 | No vocational degree | `t0_sozio_08b == 1` | |
| — | Unknown → NA | 998 / 999 | |
