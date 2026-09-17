# =============================================================================
# GxE on Education in Germany — SHIP phenotype data mining
# -----------------------------------------------------------------------------
# Purpose      : Build the harmonised SHIP phenotype file feeding 03_Merge.
#                Combines SHIP-Start (df_0 ... df_4) and SHIP-TREND
#                (df_t0, df_t1) into one per-person table with the
#                cross-cohort harmonised columns.
# Inputs       : - ship_2024_122_d_s0..s4_20250115.sas7bdat   (SHIP-Start)
#                - ship_2024_122_d_t0..t1_20250115.sas7bdat   (SHIP-TREND)
# Output       : EduGxE_DE_SHIP_Phenotype_<YYYYMMDD>.rda
# Dependencies : haven, dplyr, purrr, tidyr, stringr, lubridate
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

rm(list = ls())

library(haven)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(lubridate)

# Paths and date stamp ---------------------------------------------------------
data_dir <- "${SHIP_SHARE}/private/data/006_SHIP_pheno"
out_dir  <- "${SHIP_SHARE}/private/data/007_SHIP_processed_pheno"
date_tag <- format(Sys.Date(), "%Y%m%d")

# Helper functions ------------------------------------------------------------

# Reduce zz_nr to last 5 characters (the rest is metadata) and surface as pid.
reduce_id <- function(df) {
  df %>%
    mutate(
      zz_nr = as.character(zz_nr),
      pid   = substr(zz_nr, pmax(nchar(zz_nr) - 4, 0), nchar(zz_nr))
    ) %>%
    select(-zz_nr) %>%
    select(pid, everything())  # pid first
}

# Harmonize per-wave SEX columns to common 1/2 numeric coding.
sex_map <- function(x) {
  if (is.numeric(x)) return(ifelse(x %in% c(1, 2), x, NA_real_))
  x_chr <- str_trim(str_to_lower(as.character(x)))
  case_when(
    x_chr %in% c("1", "m", "male", "mann", "männlich")              ~ 1,
    x_chr %in% c("2", "f", "female", "frau", "weiblich", "w")       ~ 2,
    TRUE                                                            ~ NA_real_
  )
}


# ---- 1. Load raw data -------------------------------------------------------

df_0 <- read_sas(file.path(data_dir, "ship_2024_122_d_s0_20250115.sas7bdat"))
df_1 <- read_sas(file.path(data_dir, "ship_2024_122_d_s1_20250115.sas7bdat"))
df_2 <- read_sas(file.path(data_dir, "ship_2024_122_d_s2_20250115.sas7bdat"))
df_3 <- read_sas(file.path(data_dir, "ship_2024_122_d_s3_20250115.sas7bdat"))
df_4 <- read_sas(file.path(data_dir, "ship_2024_122_d_s4_20250115.sas7bdat"))

df_t0 <- read_sas(file.path(data_dir, "ship_2024_122_d_t0_20250115.sas7bdat"))
df_t1 <- read_sas(file.path(data_dir, "ship_2024_122_d_t1_20250115.sas7bdat"))

df_0  <- reduce_id(df_0)
df_1  <- reduce_id(df_1)
df_2  <- reduce_id(df_2)
df_3  <- reduce_id(df_3)
df_4  <- reduce_id(df_4)
df_t0 <- reduce_id(df_t0)
df_t1 <- reduce_id(df_t1)


# ---- 2. Person-level scaffold ----------------------------------------------
# Each .sas file is already one row per person; reduce_id() above harmonises
# pid across waves.


# ---- 3. Demographics --------------------------------------------------------
# Builds birth_year for the two carrier dataframes (df_0 for SHIP-Start,
# df_t0 for SHIP-TREND).

# Reduce df_0 to the variables we will use downstream
df_0 <- df_0 %>%
  select(
    pid, SEX, EXDATE_SHIP0, AGE_SHIP0,
    school_s0, edlevel_s0, edyrs_s0,
    schule1, ausbild, ausbild2, ausbild3,
    ausbild4, ausbild5, ausbild6, ausbild7, ausbild8,
    wohn89, wend47, wend48, wend49, agedays,
    som_groe, som_gew
  )

# birth_year for SHIP-Start (fractional year from exam date and age in days)
df_0 <- df_0 %>%
  mutate(
    birth_year = round(
      lubridate::year(EXDATE_SHIP0) +
        (lubridate::yday(EXDATE_SHIP0) - agedays) / 365.25,
      2
    )
  )

# birth_year for SHIP-TREND (fractional year from exam date and age in years)
df_t0 <- df_t0 %>%
  mutate(
    EXDATE_SHIP_T0 = as_date(EXDATE_SHIP_T0),  # ensure date even if POSIXct
    birth_year = round(
      year(EXDATE_SHIP_T0) + yday(EXDATE_SHIP_T0) / 365.25 -
        as.numeric(AGE_SHIP_T0),
      2
    )
  )

# SEX is consistent across all seven waves, so it is read from the carrier
# dataframes df_0 / df_t0 only.


# ---- 4. Education years -----------------------------------------------------
# SOEP-style years of education = school base + one vocational/tertiary
# increment. Computed per wave (df_0, df_1, df_t0, df_t1), then df_0 is
# harmonised with df_1 and df_t0 with df_t1 using the <25-at-baseline rule
# (use t1 if young at t0; otherwise prefer t0 and only fill if t0 missing).

## ---- 4.1 SHIP-Start df_0 ('schule' + 'ausbild1-8') ----
# - Map highest school degree (schule1) to SOEP base years.
# - ALWAYS fall back to reported school years (capped at 13) if base is NA.
# - Add ONE highest vocational increment from ausbild2..ausbild8.
df_0 <- df_0 %>%
  mutate(
    schule1_num = as.numeric(replace(schule1, schule1 %in% c(998, 999), NA)),
    school_base = case_when(
      schule1_num == 1 ~ NA_real_, # derzeit Schüler(in)
      schule1_num == 2 ~ 7,        # ohne Abschluss
      schule1_num == 3 ~ 9,        # Volks-/Hauptschule
      schule1_num == 4 ~ 10,       # Mittlere Reife/Realschule/Fachschule
      schule1_num == 5 ~ 10,       # POS (DDR)
      schule1_num == 6 ~ 10,       # Fachschulreife (mittlerer Abschluss; SOEP convention)
      schule1_num == 8 ~ 13,       # Fachhochschulreife
      schule1_num == 7 ~ 13,       # Abitur/EOS (auch mit Facharbeiter)
      schule1_num == 9 ~ 10,       # Sonstiges/keine Angabe → 10 (SOEP convention)
      TRUE             ~ NA_real_
    )
  ) %>%
  # Optional fallback: if a years-in-school variable is present, use it for
  # missing base values (capped at 13).
  {
    if ("schule_jahre" %in% names(.)) {
      mutate(.,
        schule_jahre_num = as.numeric(replace(schule_jahre,
                                              schule_jahre %in% c(998, 999), NA)),
        school_base = ifelse(is.na(school_base) & !is.na(schule_jahre_num),
                             pmin(schule_jahre_num, 13),
                             school_base)
      ) %>% select(-schule_jahre_num)
    } else .
  } %>%
  mutate(
    # Highest vocational/tertiary increment (choose ONE, prioritise strongest)
    voc_inc = case_when(
      ausbild7 == 1 ~ 5,    # Hochschulabschluss
      ausbild6 == 1 ~ 3,    # FH/Ingenieurschule/Polytechnikum
      ausbild5 == 1 ~ 3,    # Fach-/Berufsfachschule, Handelsschule, Fachakademie
      ausbild4 == 1 ~ 1.5,  # Lehre mit Abschlussprüfung
      ausbild8 == 1 ~ 1.5,  # anderer beruflicher Abschluss
      ausbild3 == 1 ~ 0,    # Anlernzeit/Teilfacharbeiter
      ausbild2 == 1 ~ 0,    # kein beruflicher Abschluss
      TRUE          ~ 0     # ongoing training (ausbild==1) does not add years
    ),
    eduyrs_soep = school_base + voc_inc
  ) %>%
  select(-schule1_num)


## ---- 4.2 SHIP-Start df_1 ('sozio_*') ----
# tolerate "sozio__5" vs "sozio_5"
if ("sozio__5" %in% names(df_1) && !"sozio_5" %in% names(df_1)) {
  df_1 <- dplyr::rename(df_1, sozio_5 = sozio__5)
}

df_1 <- df_1 %>%
  mutate(
    s5 = as.numeric(ifelse(sozio_5 %in% c(998, 999), NA, sozio_5)),
    s6 = as.numeric(ifelse(sozio_6 %in% c(998, 999), NA, sozio_6)),

    school_base_orig = case_when(
      s6 == 1 ~ 7,         # no degree
      s6 == 2 ~ 9,         # Hauptschule/POS 8./9.
      s6 == 3 ~ 10,        # Realschule/POS
      s6 == 4 ~ 13,        # Fachhochschulreife
      s6 == 5 ~ 13,        # Abitur/EOS
      s6 == 6 ~ 10,        # "other" → SOEP uses 10
      s6 == 7 ~ NA_real_,  # still in school (no final degree)
      TRUE    ~ NA_real_
    ),

    # ALWAYS fallback to self-reported school years (clamped to [7, 13])
    school_base = ifelse(is.na(school_base_orig) & !is.na(s5),
                         pmin(pmax(s5, 7), 13),
                         school_base_orig),

    # Highest vocational/tertiary increment
    v9   = as.numeric(sozio_9),
    v10  = as.numeric(sozio_10),
    v11  = as.numeric(sozio_11),
    v12  = as.numeric(sozio_12),
    v13  = as.numeric(sozio_13),
    v13a = as.numeric(sozio_13a),

    voc_inc = case_when(
      v13a == 1 ~ 5,    # University
      v13  == 1 ~ 3,    # Fachschule / Fachhochschule (higher technical)
      v11  == 1 ~ 3,    # Meister / Techniker
      v12  == 1 ~ 2,    # Berufsfachschule / Handelsschule
      v10  == 1 ~ 1.5,  # Apprenticeship (Lehre)
      v9   == 1 ~ 0,    # no vocational degree
      TRUE      ~ 0
    ),

    eduyrs_soep = school_base + voc_inc
  ) %>%
  select(-s5, -s6, -v9, -v10, -v11, -v12, -v13, -v13a)


## ---- 4.3 Harmonise df_0 with df_1 ----
# Combination rules (df_0 is the carrier):
#   1) Prefer eduyrs_soep from df_0; fill from df_1 only if df_0 was NA.
#   2) If age < 25 at wave 0 (agedays < 9125), force df_1's value.
#   3) If still in school at t1 (sozio_8 == 1), set final to NA.
df_0 <- df_0 %>%
  left_join(
    df_1 %>% select(pid,
                    eduyrs_soep_t1 = eduyrs_soep,
                    sozio_8_t1     = sozio_8),
    by = "pid"
  ) %>%
  mutate(
    eduyrs_soep_t1 = as.numeric(eduyrs_soep_t1),
    eduyrs_soep    = as.numeric(eduyrs_soep),

    eduyrs_soep_harmonized = case_when(
      sozio_8_t1 == 1                     ~ NA_real_,
      !is.na(agedays) & agedays < 9125    ~ eduyrs_soep_t1,
      is.na(eduyrs_soep)                  ~ eduyrs_soep_t1,
      TRUE                                ~ eduyrs_soep
    )
  ) %>%
  select(-eduyrs_soep_t1, -sozio_8_t1)


## ---- 4.4 SHIP-TREND df_t0 ('t0_sozio_*') ----
df_t0 <- df_t0 %>%
  mutate(
    t06 = as.numeric(replace(t0_sozio_06, t0_sozio_06 %in% c(998, 999), NA)),
    school_base = case_when(
      t06 == 1 ~ NA_real_, # still pupil
      t06 == 2 ~ 7,        # exit without degree
      t06 == 3 ~ 9,        # Volks-/Hauptschule
      t06 == 4 ~ 10,       # Mittlere Reife/Realschule/Fachschulreife
      t06 == 5 ~ 10,       # POS (DDR) ≈ Realschule
      t06 == 6 ~ 12,       # Fachhochschulreife/FOS
      t06 == 7 ~ 13,       # Abitur/EOS
      t06 == 8 ~ 10,       # other / no answer → 10 (SOEP convention)
      TRUE     ~ NA_real_
    )
  ) %>%
  # Optional fallback: years-in-school if available
  {
    if ("t0_sozio_05" %in% names(.)) {
      mutate(.,
        t05_num = as.numeric(replace(t0_sozio_05,
                                     t0_sozio_05 %in% c(998, 999), NA)),
        school_base = ifelse(is.na(school_base) & !is.na(t05_num),
                             pmin(t05_num, 13),
                             school_base)
      ) %>% select(-t05_num)
    } else .
  } %>%
  mutate(
    voc_inc = case_when(
      t0_sozio_08g == 1 ~ 5,    # Hochschulabschluss
      t0_sozio_08f == 1 ~ 3,    # Fachhochschulabschluss
      t0_sozio_08e == 1 ~ 3,    # Fachschule/Meister/Techniker
      t0_sozio_08d == 1 ~ 2,    # Berufsfach-/Handelsschule
      t0_sozio_08c == 1 ~ 1.5,  # Lehre
      t0_sozio_08h == 1 ~ 1.5,  # other vocational
      t0_sozio_08b == 1 ~ 0,    # explicitly no vocational degree
      TRUE              ~ 0
    ),
    eduyrs_soep = school_base + voc_inc
  ) %>%
  select(-t06)


## ---- 4.5 SHIP-TREND df_t1 (sozio_*) ----
df_t1 <- df_t1 %>%
  mutate(
    s06 = as.numeric(replace(sozio_06, sozio_06 %in% c(998, 999), NA)),
    school_base = case_when(
      s06 == 1 ~ NA_real_,
      s06 == 2 ~ 7,
      s06 == 3 ~ 9,
      s06 == 4 ~ 10,
      s06 == 5 ~ 10,
      s06 == 6 ~ 12,
      s06 == 7 ~ 13,
      s06 == 8 ~ 10,
      TRUE     ~ NA_real_
    ),
    voc_inc = case_when(
      sozio_08g == 1 ~ 5,
      sozio_08f == 1 ~ 3,
      sozio_08e == 1 ~ 3,
      sozio_08d == 1 ~ 2,
      sozio_08c == 1 ~ 1.5,
      sozio_08h == 1 ~ 1.5,
      sozio_08b == 1 ~ 0,
      TRUE           ~ 0
    ),
    eduyrs_soep = school_base + voc_inc
  ) %>%
  select(-s06)


## ---- 4.6 Harmonise df_t0 with df_t1 ----
# Rules (no "still in school" flag at t1):
#   1) Prefer t0; fill from t1 only if t0 missing.
#   2) If age < 25 at t0, take t1 when available; otherwise NA.
df_t0 <- df_t0 %>%
  left_join(
    df_t1 %>% select(pid, eduyrs_soep_t1 = eduyrs_soep),
    by = "pid"
  ) %>%
  mutate(
    eduyrs_soep    = as.numeric(eduyrs_soep),
    eduyrs_soep_t1 = as.numeric(eduyrs_soep_t1),

    under25_t0 = !is.na(AGE_SHIP_T0) & as.numeric(AGE_SHIP_T0) < 25,

    eduyrs_soep_harmonized = case_when(
      under25_t0 & !is.na(eduyrs_soep_t1) ~ eduyrs_soep_t1,
      under25_t0 &  is.na(eduyrs_soep_t1) ~ NA_real_,
      is.na(eduyrs_soep)                  ~ eduyrs_soep_t1,
      TRUE                                ~ eduyrs_soep
    )
  ) %>%
  select(-eduyrs_soep_t1, -under25_t0)


# ---- 5. Parental education --------------------------------------------------
# Not available in the SHIP delivery.


# ---- 6. East / West classification ------------------------------------------
# Unified reasoning across both SHIP cohorts:
#   - All SHIP participants live in Mecklenburg-Vorpommern (East) today.
#   - Therefore the default origin assumption is East.
#   - Only an explicit wohn89 == 2 (West-resident in 1989) overrides this
#     default to mark the participant as a West-to-East mover.
#   - SHIP-Start (df_0) has wohn89, so the override CAN apply.
#   - SHIP-TREND (df_t0) has no wohn89 at all; the override never fires
#     and every t0 participant falls through to the East default.
# Limitation: t0 participants who actually were West-resident pre-1989
# cannot be distinguished from East-origin in the available data.

# SHIP-Start: derive from wohn89 where available, default to east otherwise
df_0 <- df_0 %>%
  mutate(
    wohn89_num = suppressWarnings(as.numeric(wohn89)),

    east_west = case_when(
      wohn89_num == 2 ~ "west",
      TRUE            ~ "east"   # default: lives in MV (East) today
    ),

    east_west_mig = case_when(
      wohn89_num == 2 ~ "west_to_e_move",
      TRUE            ~ "east_stay"   # default: East today, never moved
    )
  ) %>%
  mutate(
    .age_num = suppressWarnings(as.numeric(AGE_SHIP0)),
    # If the respondent was already 25+ at baseline AND has a lifetime
    # west_to_e move, downgrade it to west_stay for the education window
    # (move must have happened post-education). Everyone else keeps the
    # lifetime label. NA only when age is unknown.
    east_west_mig_edu = case_when(
      is.na(.age_num)                                       ~ NA_character_,
      .age_num >= 25 & east_west_mig == "west_to_e_move"    ~ "west_stay",
      TRUE                                                  ~ east_west_mig
    )
  ) %>%
  select(-wohn89_num, -.age_num)

# SHIP-TREND: no wohn89 variable, so the override never fires and the
# East default applies to all participants.
df_t0 <- df_t0 %>%
  mutate(
    east_west         = "east",
    east_west_mig     = "east_stay",
    east_west_mig_edu = "east_stay"
  )


# ---- 7. Anthropometrics -----------------------------------------------------
# All participants gave anthropometric data at the first time point
# (df_0 in SHIP-Start, df_t0 in SHIP-TREND). bmi = kg / m^2.

df_0 <- df_0 %>%
  mutate(
    bmi        = som_gew / ((som_groe / 100) ^ 2),
    bmi_age    = agedays / 365.25,
    height_age = agedays / 365.25
  ) %>%
  rename(height = som_groe) %>%
  select(-som_gew)

df_t0 <- df_t0 %>%
  mutate(
    bmi        = som_gew / ((som_groe / 100) ^ 2),
    bmi_age    = AGE_SHIP_T0,
    height_age = AGE_SHIP_T0
  ) %>%
  rename(height = som_groe) %>%
  select(-som_gew)


# ---- 9. Final harmonised dataset --------------------------------------------

# Standardise sex variable name in df_t0: use "SEX"
if ("SEX_SHIP_T0" %in% names(df_t0) && !"SEX" %in% names(df_t0)) {
  df_t0 <- df_t0 %>% rename(SEX = SEX_SHIP_T0)
} else if ("SEX_SHIP_T0" %in% names(df_t0) && "SEX" %in% names(df_t0)) {
  df_t0 <- df_t0 %>% select(-SEX_SHIP_T0)  # keep existing SEX, drop duplicate
}

keep_cols <- c(
  "pid", "SEX", "birth_year", "eduyrs_soep_harmonized",
  "east_west", "east_west_mig", "east_west_mig_edu",
  "bmi", "bmi_age", "height", "height_age"
)

df_0 <- df_0 %>%
  select(any_of(keep_cols)) %>%
  mutate(pid = paste0("s_", as.character(pid)))  # SHIP-Start prefix

df_t0 <- df_t0 %>%
  select(any_of(keep_cols)) %>%
  mutate(pid = paste0("t_", as.character(pid)))  # SHIP-TREND prefix

df_final <- bind_rows(df_0, df_t0) %>%
  rename(
    gender    = SEX,
    education = eduyrs_soep_harmonized
  ) %>%
  mutate(
    gender = case_when(
      gender == 1 ~ "Male",
      gender == 2 ~ "Female",
      TRUE        ~ NA_character_
    )
  )

# Save final data — .sav keeps column descriptions; .rda is the canonical format
out_sav <- file.path(out_dir,
                     sprintf("EduGxE_DE_SHIP_Phenotype_%s.sav", date_tag))
out_rda <- file.path(out_dir,
                     sprintf("EduGxE_DE_SHIP_Phenotype_%s.rda", date_tag))

# write_sav() can't represent SAS-style character-tagged NAs. Tagged NAs
# from read_sas() are encoded in the IEEE-754 NaN payload, so they survive
# attribute stripping and arithmetic; we have to overwrite every NA cell
# with a plain NA_real_ before writing. The .rda keeps the original.
df_final_for_sav <- df_final %>%
  mutate(across(where(is.numeric),
                ~ replace(as.numeric(.), is.na(.), NA_real_)))
write_sav(df_final_for_sav, out_sav)

save(df_final, file = out_rda, compress = "xz", version = 3)
cat("Saved:", out_sav, "\n")
cat("Saved:", out_rda, "\n")


# ---- 10. Aggregate sanity check ---------------------------------------------
# Quick aggregate-only summary printed on each run.

cat("\n=== SHIP phenotype summary ===\n")
cat("Total N:", nrow(df_final), "\n\n")

cat("gender:\n")
print(table(df_final$gender, useNA = "ifany"))

cat("\neast_west:\n")
print(table(df_final$east_west, useNA = "ifany"))

cat("\neast_west_mig:\n")
print(table(df_final$east_west_mig, useNA = "ifany"))

cat("\n--- By-region education summary ---\n")
df_final %>%
  filter(east_west %in% c("east", "west")) %>%
  group_by(east_west) %>%
  summarise(
    n_nonmiss_edu = sum(!is.na(education)),
    female_pct    = round(100 * sum(gender == "Female" & !is.na(education),
                                    na.rm = TRUE) / sum(!is.na(education)), 1),
    M_education   = round(mean(education, na.rm = TRUE), 2),
    SD_education  = round(sd(education,   na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()
