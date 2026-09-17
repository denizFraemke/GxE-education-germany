# =============================================================================
# GxE on Education in Germany — BASE-II phenotype data mining
# -----------------------------------------------------------------------------
# Purpose      : Build the harmonised BASE-II phenotype file feeding 03_Merge.
# Inputs       : - BASEII_GSOEP12_Fraemke_data2024_3.sav
#                - base2_select_nov25.sas7bdat             (cleaned)
#                - SOEP-BASE Long.dta                      (SOEP-2014 long)
# Output       : Final_data/EduGxE_DE_BASEII_Phenotype_<YYYYMMDD>.{rda,sav}
# Dependencies : dplyr, haven, tidyr
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(haven)
library(tidyr)

# Shared helpers (neg_to_na, earliest_nonmissing, latest_nonmissing, mean_na).
# Run this script from its containing directory (02_Phenotype/BASE-II/).
source(file.path("..", "utils", "edu_helpers.R"))

# Paths and date stamp ---------------------------------------------------------
data_dir <- "${BASE2_SHARE}/private/data/004_BASEII_Socioeceonomic_Pheno"
out_dir  <- file.path(data_dir, "Final_data")
date_tag <- format(Sys.Date(), "%Y%m%d")


# ---- 1. Load raw data -------------------------------------------------------

df_all <- read_sav(file.path(data_dir, "BASEII_GSOEP12_Fraemke_data2024_3.sav"))

# Cleaned data
df_denis <- read_sas(file.path(data_dir, "base2_select_nov25.sas7bdat"))

# Full SOEP-2014 long data
df_long <- read_dta(file.path(data_dir, "SOEP-BASE Long.dta"))

# Pull SOEP-2014 variables of interest into df_soep (larger sample than BASE-II core)
df_soep <- df_long %>%
  select(
    # ------------------------------------------------------------------
    # 1) Core IDs + key variables
    # ------------------------------------------------------------------
    any_of(c(
      # IDs / design
      "nkidpznum", "nkidpz", "pnrfest", "jahr", "hnr","psex", "lgeb", "lgebmo",

      # Anthropometrics
      'pgr', 'pkilo',

      #Parents Education
      "lv05", "lm05", "lm06","lv06",

      # GDR Variables
      'lsab1', #Jahr letzter Schulbesuch
      'lsab3', #Land letzter Schulbesuch incl. GDR & FRG
      "lsab7x", #Schulabschluss in der BRD
      "lsab7y", #Schulabschluss in der DDR

      "pddr01",  # Mind. 1 Jahr in der DDR gelebt
      "pddr02a", # Umzug in die BRD geplant
      "pddr02b", # Ausreiseantrag gestellt
      "pddr02no",# Weder Umzug noch Ausreiseantrag
      "pddr03a", # Selbst inhaftiert aus politischen Gründen
      "pddr03b", # Familienmitglieder inhaftiert aus politischen Gründen
      "pddr03no",# Weder ich selbst noch Familienmitglieder waren inhaftiert
      "pddr04a", # Mitglied einer Oppositionsbewegung
      "pddr04b", # Familienmitglieder in der Oppositionsbewegung
      "pddr04no",# Weder ich selbst noch Familienmitglieder in der Opposition
      "pddr05a", # Familienmitglieder/Freunde als Informanten der Stasi
      "pddr05b", # Arbeitskollegen/Nachbarn als Informanten der Stasi
      "pddr05no",# Keine Stasi-Informanten im Umfeld
      "pddr05wn",# Stasi-Informanten: Weiß nicht
      "pddr06",  # Einsicht in Stasi-Akten genommen, beantragt oder geplant
      "pddr071", # Leben in der DDR war alles in allem gut (Vergleich zu heute)
      "pddr072", # Aufgabe der Staatssicherheit war der Schutz der DDR-Bürger
      "pddr073", # Stasi-Verbrechen sollten härter strafrechtlich verfolgt werden


      # Migration-related
      "lsta2nr",   # birth country excl. GDR
      "lzug01",    # moved to FRG excluding moved from GDR

      # Berlin residence 2002–2012
      "pwohnvb02", "pwohnvb03", "pwohnvb04", "pwohnvb05",
      "pwohnvb06", "pwohnvb07", "pwohnvb08", "pwohnvb09",
      "pwohnvb10", "pwohnvb11", "pwohnvb12"


    )))

# Save variable and value labels from df_long version
var_labels <- lapply(df_soep, function(x) attr(x, "label"))
val_labels <- lapply(df_soep, function(x) attr(x, "labels"))


# ---- 2. Person-level scaffold -----------------------------------------------

# Add birth year and sex from the cleaned dataset to df_all -> df

df <- df_all %>%
  left_join(
    df_denis %>%
      select(nKIDPZnum, sex, birth_yr, birth_mo),
    by = "nKIDPZnum"
  )

df <- df %>%
  mutate(birth_year = birth_yr + birth_mo / 12)

# Clean df_soep: replace negative codes with NA across all columns

df_soep <- df_soep %>%
  mutate(
    across(everything(),neg_to_na)
  )

# Collapse df_soep to one row per person --------------------------------------

vars_special <- c("lgeb",
                  "lgebmo",
                  "jahr", "hnr", "pgr", "pkilo",
                  "lv05", "lm05","lm06","lv06",
                  "lsab1", "lsab3",
                  "lsab7x", "lsab7y",
                  "lzug01")

vars_other <- setdiff(names(df_soep), c("pnrfest", vars_special))

df_soep <- df_soep %>%
  group_by(pnrfest) %>%
  summarise(
    jahr_min = min(as.numeric(jahr), na.rm = TRUE),
    jahr_max = max(as.numeric(jahr), na.rm = TRUE),

    lgeb  = earliest_nonmissing(lgeb,jahr),
    lgebmo  = earliest_nonmissing(lgebmo,jahr),

    hnr   = latest_nonmissing(hnr, jahr),
    pgr   = mean_na(pgr),
    pkilo = mean_na(pkilo),

    lv05   = earliest_nonmissing(lv05, jahr),
    lm05   = earliest_nonmissing(lm05, jahr),

    lv06   = earliest_nonmissing(lv06, jahr),
    lm06   = earliest_nonmissing(lm06, jahr),

    lsab7x = earliest_nonmissing(lsab7x, jahr),
    lsab7y = earliest_nonmissing(lsab7y, jahr),

    lsab1  = latest_nonmissing(lsab1, jahr),
    lsab3  = latest_nonmissing(lsab3, jahr),

    lzug01 = earliest_nonmissing(lzug01, jahr),

    across(all_of(vars_other),
           ~ { v <- .[!is.na(.)]; if (length(v) == 0) NA_real_ else v[1] }),

    .groups = "drop"
  )

# Re-attach variable / value labels stripped by the summarise() above

for (v in names(df_soep)) {
  if (!is.null(var_labels[[v]])) {
    attr(df_soep[[v]], "label") <- var_labels[[v]]
  }
  if (!is.null(val_labels[[v]])) {
    attr(df_soep[[v]], "labels") <- val_labels[[v]]
  }
}

# Merge df_soep extras (anthropometrics + parent training + schooling) into df

df <- df %>%
  left_join(
    df_soep %>%
      mutate(nKIDPZnum = nkidpznum) %>%
      select(nKIDPZnum, pgr, pkilo, lsab3, lsab1, lv06, lm06),
    by = "nKIDPZnum"
  )


# ---- 3. Demographics --------------------------------------------------------
# Already populated above by Section 2 (df$sex, df$birth_year).


# ---- 4. Education years -----------------------------------------------------
# BASE-II's harmonised `education` column is a direct rename of the
# precomputed Educ_final column; the rename happens inside the df_final
# mutate in Section 9.


# ---- 5. Parental education --------------------------------------------------

df <- df %>%
  mutate(

    # SOEP schooling years (pgbilzt logic)
    mother_school_years = case_when(
      lm05 == 1 ~ 7,   # no degree
      lm05 == 2 ~ 9,   # Hauptschule
      lm05 == 3 ~ 10,  # Realschule
      lm05 == 4 ~ 13,  # Abitur/EOS
      lm05 == 5 ~ 10,  # other
      TRUE       ~ NA_real_
    ),
    father_school_years = case_when(
      lv05 == 1 ~ 7,
      lv05 == 2 ~ 9,
      lv05 == 3 ~ 10,
      lv05 == 4 ~ 13,
      lv05 == 5 ~ 10,
      TRUE       ~ NA_real_
    ),

    # Further education / training years (coarse SOEP mapping)
    # 1 = vocational training (apprenticeship) ≈ 3 years
    # 2 = university degree ≈ 5 years
    # 3 = no completed training → 0 years
    # 4 / NA: treated as 0 extra years
    mother_train_years = case_when(
      lm06 == 1 ~ 1.5,
      lm06 == 2 ~ 5,
      lm06 == 3 ~ 0,
      TRUE        ~ 0
    ),
    father_train_years = case_when(
      lv06 == 1 ~ 1.5,
      lv06 == 2 ~ 5,
      lv06 == 3 ~ 0,
      TRUE        ~ 0
    ),

    # Total years of education (7–18 range, SOEP-style)
    mother_education = mother_school_years + mother_train_years,
    father_education = father_school_years + father_train_years
  )

# Row-mean of mother + father
df <- df %>%
  mutate(
    parental_education = rowMeans(
      cbind(mother_education, father_education),
      na.rm = TRUE
    ),
    parental_education = ifelse(
      is.nan(parental_education),
      NA_real_,
      parental_education
    )
    # Not standardized here: standardization happens at analysis time,
    # within parent birth-cohort bins, across pooled data.
  )


# ---- 6. East / West classification ------------------------------------------
# Build three harmonised columns for the older cohort (agegr == 2) only.
# For agegr == 1 (younger cohort) all three columns remain NA — we have no
# informative pre-1989 signal for them.
#
#   east_west          — region of origin / where the person grew up
#   east_west_mig      — lifetime migration (East <-> West, including
#                        post-reunification Berlin residence)
#   east_west_mig_edu  — migration during education only (schooling
#                        switches; ignores post-reunification Berlin)
#
# Inputs used:
#   pddr01           — self-report: lived >= 1 year in the GDR (1=yes, 2=no)
#   lsab3            — country of last school visit (1=FRG, 2=GDR)
#   lsab7x / lsab7y  — final school certificate received in FRG / GDR
#   pwohnvb02..12    — Berlin borough of residence each year 2002-2012
#   agegr            — 1 = younger cohort, 2 = older cohort

## ---- 6.1 Berlin residence signals (2002-2012) ----
# Per respondent, find the first year (2002-2012) they reported a borough
# that is unambiguously East or West Berlin. Used by 6.3 east_west_mig only.
# Borough codes (pwohnvb*):
#   East: Pankow (3), Treptow-Köpenick (9), Marzahn-Hellersdorf (10), Lichtenberg (11)
#   West: Charlottenburg-Wilmersdorf (4), Spandau (5), Steglitz-Zehlendorf (6),
#         Tempelhof-Schöneberg (7), Neukölln (8), Reinickendorf (12)
#   Neutral / not classified: Mitte (1), Friedrichshain-Kreuzberg (2)
#   Missing codes: -3, -2, -1

berlin_vars <- intersect(paste0("pwohnvb", sprintf("%02d", 2:12)), names(df))

if (length(berlin_vars) > 0L) {
  df <- df %>%
    mutate(.rowid_berlin = dplyr::row_number())

  berlin_summary <- df %>%
    select(.rowid_berlin, all_of(berlin_vars)) %>%
    pivot_longer(
      cols      = all_of(berlin_vars),
      names_to  = "wave",
      values_to = "code_raw"
    ) %>%
    mutate(
      code      = ifelse(code_raw %in% c(-3, -2, -1), NA_real_, as.numeric(code_raw)),
      year      = as.integer(sub("pwohnvb", "", wave)) + 2000L,
      east_flag = code %in% c(3, 9, 10, 11),
      west_flag = code %in% c(4, 5, 6, 7, 8, 12)
    ) %>%
    arrange(.rowid_berlin, year) %>%
    group_by(.rowid_berlin) %>%
    summarise(
      first_res_east_ber = {
        idx <- which(east_flag | west_flag)
        if (length(idx) == 0L) FALSE else east_flag[idx[1L]]
      },
      first_res_west_ber = {
        idx <- which(east_flag | west_flag)
        if (length(idx) == 0L) FALSE else west_flag[idx[1L]]
      },
      .groups = "drop"
    )

  df <- df %>%
    left_join(berlin_summary, by = ".rowid_berlin") %>%
    select(-.rowid_berlin)
} else {
  df <- df %>%
    mutate(first_res_east_ber = FALSE, first_res_west_ber = FALSE)
}

## ---- 6.2 Schooling signals + origin signals (older cohort) ----
# - east_school / west_school: TRUE if respondent has a degree from / last
#   attended school in the GDR / FRG
# - east_origin / west_origin: TRUE for older cohort (agegr == 2) by their
#   self-reported region of pre-reunification residence (pddr01)

df <- df %>%
  mutate(
    ddr_degree  = !is.na(lsab7y),
    frg_degree  = !is.na(lsab7x),
    east_school = (ddr_degree | lsab3 == 2),
    west_school = (frg_degree | lsab3 == 1),

    east_origin = (agegr == 2) & (pddr01 == 1),
    west_origin = (agegr == 2) & (pddr01 == 2)
  )

## ---- 6.3 Build east_west, east_west_mig, east_west_mig_edu (older cohort) ----
# All three columns are NA for agegr != 2.

df <- df %>%
  mutate(
    east_west = case_when(
      agegr != 2                                  ~ NA_character_,
      east_origin                                 ~ "east",
      west_origin                                 ~ "west",
      is.na(pddr01) & east_school & !west_school  ~ "east",
      is.na(pddr01) & west_school & !east_school  ~ "west",
      TRUE                                        ~ NA_character_
    ),

    # Lifetime migration: any contrary signal (school OR Berlin) flags the
    # respondent as a mover; otherwise stayer in the origin region.
    # Stayer status requires no positive evidence on either axis — the
    # absence of a contrary signal is sufficient.
    east_west_mig = case_when(
      agegr != 2                                       ~ NA_character_,
      east_origin & (west_school | first_res_west_ber) ~ "east_to_w_move",
      east_origin                                      ~ "east_stay",
      west_origin & (east_school | first_res_east_ber) ~ "west_to_e_move",
      west_origin                                      ~ "west_stay",
      TRUE                                             ~ NA_character_
    ),

    # Migration during education: schooling only (Berlin ignored).
    # Stayer requires positive schooling evidence on the matching side.
    east_west_mig_edu = case_when(
      agegr != 2                ~ NA_character_,
      east_origin & west_school ~ "east_to_w_move",
      east_origin & east_school ~ "east_stay",
      west_origin & east_school ~ "west_to_e_move",
      west_origin & west_school ~ "west_stay",
      TRUE                      ~ NA_character_
    )
  )


# ---- 7. Anthropometrics -----------------------------------------------------

df <- df %>%
  mutate(
    bmi = ifelse(
      !is.na(pkilo) & !is.na(pgr),
      pkilo / ( (pgr / 100)^2 ),
      NA_real_
    ),
    bmi_age   = age_T1_CFGG_U1,
    height    = pgr,
    height_age = age_T1_CFGG_U1
  )


# ---- 9. Final harmonised dataset --------------------------------------------

df_final <- df %>%
  mutate(
    # IDs / structure (to match e.g. TwinLife)
    pid = paste0("b_", nKIDPZnum),

    # gender
    gender = case_when(
      sex == 1 ~ "Male",
      sex == 2 ~ "Female",
      TRUE     ~ NA_character_
    ),

    # Education (years) — rename Educ_final to harmonised name
    education = Educ_final
  ) %>%
  select(
    pid, gender, birth_year,
    education,
    mother_education, father_education,
    parental_education,
    east_west, east_west_mig, east_west_mig_edu,
    bmi, bmi_age, height, height_age
  )

# Save final data — .sav keeps column descriptions; .rda is the canonical format
out_sav <- file.path(out_dir, sprintf("EduGxE_DE_BASEII_Phenotype_%s.sav", date_tag))
out_rda <- file.path(out_dir, sprintf("EduGxE_DE_BASEII_Phenotype_%s.rda", date_tag))

write_sav(df_final, out_sav)
save(df_final, file = out_rda, compress = "xz", version = 3)
cat("Saved:", out_sav, "\n")
cat("Saved:", out_rda, "\n")


# ---- 10. Aggregate sanity check ---------------------------------------------
# Quick aggregate-only summary printed on each run.

cat("\n=== BASE-II phenotype summary ===\n")
cat("Total N:", nrow(df_final), "\n\n")

cat("gender:\n")
print(table(df_final$gender, useNA = "ifany"))

cat("\nEast/West:\n")
print(table(df_final$east_west, useNA = "ifany"))

cat("\nEast/West Migration:\n")
print(table(df_final$east_west_mig, useNA = "ifany"))

cat("\n--- By-region education summary ---\n")
df_final %>%
  filter(east_west %in% c("east", "west")) %>%
  group_by(east_west) %>%
  summarise(
    n_nonmiss_edu = sum(!is.na(education)),
    female_pct    = round(100 * sum(gender == "Female" & !is.na(education), na.rm = TRUE) /
                          sum(!is.na(education)), 1),
    M_education   = round(mean(education, na.rm = TRUE), 2),
    SD_education  = round(sd(education, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()


