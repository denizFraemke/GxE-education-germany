# =============================================================================
# GxE on Education in Germany — TwinLife phenotype data mining
# -----------------------------------------------------------------------------
# Purpose      : Build the harmonised TwinLife phenotype file feeding
#                03_Merge from the long-format delivery.
# Inputs       : - data_request_pheno_2_V1.dta  (TwinLife long-format)
# Output       : EduGxE_DE_TwinLife_Phenotype_<YYYYMMDD>.{rda,sav}
# Dependencies : haven, dplyr, tidyr
# Note on waves:
#   wid: 1=F2F1, 2=CATI1, 3=F2F2, 4=CATI2, 5=F2F3, 6=CATI3,
#        7=F2F4, 8=CATI4, 9=F2F5
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

rm(list = ls())

library(haven)
library(dplyr)
library(tidyr)

# Paths and date stamp ---------------------------------------------------------
data_dir <- "${TWINLIFE_SHARE}/private/data/2025_TwinLife_GermanGenetics/raw_data/Transfer_Bielefeld_20260320"
out_dir  <- "${TWINLIFE_SHARE}/private/data/2025_TwinLife_GermanGenetics/"
date_tag <- format(Sys.Date(), "%Y%m%d")

pheno_path <- file.path(data_dir, "data_request_pheno_2_V1.dta")

# Helper functions ------------------------------------------------------------

# Pivot a single per-wave variable from long → wide with `_wid{N}` suffix.
# Used in Section 2.2.
# names_sort = TRUE forces the per-wave columns into numeric wid order.
# Without it, column order follows first-appearance in the long data —
# which happens to be sorted by wid for the current delivery, but the
# defensive flag makes downstream `tail(..., 1)` "latest non-NA" calls
# robust to delivery-order changes.
pivot_var <- function(data, var_name) {
  data %>%
    filter(!is.na(wid)) %>%
    select(pid, wid, value = all_of(var_name)) %>%
    filter(!is.na(value)) %>%
    group_by(pid, wid) %>%
    summarise(value = first(value), .groups = "drop") %>%
    pivot_wider(names_from = wid, values_from = value,
                names_prefix = paste0(var_name, "_wid"),
                names_sort   = TRUE)
}

# SOEP-style school degree → years of schooling (eca0130).
# Used in Section 3.3 (parent education) and Section 4.1 (own education).
school_years <- function(x) {
  case_when(
    x == 1 ~ 7,    # no degree
    x == 2 ~ 9,    # Hauptschulabschluss
    x == 3 ~ 10,   # Realschulabschluss
    x == 4 ~ 12,   # Fachoberschule
    x == 5 ~ 13,   # Abitur
    TRUE   ~ NA_real_
  )
}

# SOEP-style vocational degree → additional years (eca0230).
# Used in Section 3.3 and Section 4.1.
voc_years <- function(x) {
  case_when(
    x == 1  ~ 0,    # no vocational
    x == 2  ~ 1.5,  # Lehre / vocational
    x == 3  ~ 2,    # health-care schools
    x == 4  ~ 2,    # technical college / Meister
    x == 5  ~ 1.5,  # civil servants
    x == 6  ~ 3,    # FH
    x == 7  ~ 3,    # cooperative education
    x == 8  ~ 3,    # FH / cooperative (older code)
    x == 9  ~ 5,    # university
    x == 10 ~ 5,    # PhD
    TRUE    ~ NA_real_  # 11 = other → NA (unclear)
  )
}

# Normalise free-form region strings to bare "east" / "west".
# Used in Section 6.1.
norm_region <- function(x) case_when(
  is.na(x)                              ~ NA_character_,
  grepl("east", x, ignore.case = TRUE)  ~ "east",
  grepl("west", x, ignore.case = TRUE)  ~ "west",
  TRUE                                  ~ NA_character_
)


# ---- 1. Load raw data -------------------------------------------------------
# The delivery is in LONG format: one row per person × wave (wid).

df_long <- read_dta(pheno_path)

# Convert haven_labelled columns to plain numeric (strip labels)
df_long <- df_long %>% mutate(across(where(is.labelled), zap_labels))

# Replace all negative values with NA (TwinLife missing codes:
# -83, -87, -89, -90, -92, -93, -94, -95, -98, -99)
df_long <- df_long %>% mutate(across(where(is.numeric), ~ ifelse(. < 0, NA, .)))


# ---- 2. Person-level scaffold -----------------------------------------------

## ---- 2.1 Time-constant identifiers (one row per pid) ----
# fid, cgr, ptyp_1, sex, wav0100, mig2001, mig2201, fpr0104 (birth month),
# fpr0105 (birth year). All constant across waves.

# first(.x, default = NA, na_rm = TRUE) returns NA for a pid group that is
# all-NA on a column (e.g. the pid == NA group) instead of a zero-length
# vector, so those groups collapse to one NA row.
id_vars <- df_long %>%
  group_by(pid) %>%
  summarise(
    across(
      c(fid, cgr, ptyp_1, sex, wav0100, mig2001, mig2201),
      ~ dplyr::first(.x, default = NA, na_rm = TRUE)
    ),
    .groups = "drop"
  )

cat("Unique individuals in new data:", nrow(id_vars), "\n")

## ---- 2.2 Pivot long → wide (per-wave columns) ----
# Per-wave values for: age0100, yea_fq, mon_fq, bdy0100_hgt, bdy0300, ewi,
# ptyp, plus the education variables eca0130 / eca0230.
# (pivot_var() is defined in Section 0.)

wide_age     <- pivot_var(df_long, "age0100")
wide_yea_fq  <- pivot_var(df_long, "yea_fq")
wide_mon_fq  <- pivot_var(df_long, "mon_fq")
wide_hgt     <- pivot_var(df_long, "bdy0100_hgt")
wide_bmi     <- pivot_var(df_long, "bdy0300")
wide_ewi     <- pivot_var(df_long, "ewi")
wide_ptyp    <- pivot_var(df_long, "ptyp")
wide_eca0130 <- pivot_var(df_long, "eca0130")
wide_eca0230 <- pivot_var(df_long, "eca0230")

df_twinlife <- id_vars %>%
  left_join(wide_age,     by = "pid") %>%
  left_join(wide_yea_fq,  by = "pid") %>%
  left_join(wide_mon_fq,  by = "pid") %>%
  left_join(wide_hgt,     by = "pid") %>%
  left_join(wide_bmi,     by = "pid") %>%
  left_join(wide_ewi,     by = "pid") %>%
  left_join(wide_ptyp,    by = "pid") %>%
  left_join(wide_eca0130, by = "pid") %>%
  left_join(wide_eca0230, by = "pid")


# ---- 3. Demographics --------------------------------------------------------

## ---- 3.1 birth_year ----
# birth = yea_fq_wid1 + (mon_fq_wid1 / 12) - age0100_wid1
df_twinlife <- df_twinlife %>%
  mutate(birth = yea_fq_wid1 + (mon_fq_wid1 / 12) - age0100_wid1)

## ---- 3.2 gender ----
df_twinlife <- df_twinlife %>%
  mutate(
    gender = case_when(
      sex == 1 ~ "Male",
      sex == 2 ~ "Female",
      TRUE     ~ NA_character_
    )
  )

## ---- 3.3 Parent birth years + parent education from ptyp 300 / 400 ----
# Parents appear in the same long file with ptyp == 300 (mother) or 400
# (father). We extract per-parent education (latest eca0130 + eca0230)
# and birth year, then pivot to one row per family with mother/father
# columns. Joining to df_twinlife by fid spreads those values to every
# child in the family. The parental_education row-mean is built later
# in Section 5.

# (school_years and voc_years are defined in Section 0.)

parent_info <- df_long %>%
  filter(ptyp %in% c(300, 400)) %>%
  group_by(pid) %>%
  summarise(
    fid   = first(na.omit(fid)),
    ptyp  = first(na.omit(ptyp)),
    # Take latest non-NA education values across waves
    eca0130_val = {
      v <- eca0130[!is.na(eca0130)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    },
    eca0230_val = {
      v <- eca0230[!is.na(eca0230)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    },
    # Birth year: from first wave with non-NA timing + age
    birth_year_parent = {
      idx <- which(!is.na(yea_fq) & !is.na(mon_fq) & !is.na(age0100))
      if (length(idx) == 0) NA_real_
      else yea_fq[idx[1]] + mon_fq[idx[1]] / 12 - age0100[idx[1]]
    },
    .groups = "drop"
  ) %>%
  mutate(
    parent_edu_years = school_years(eca0130_val) + voc_years(eca0230_val),
    parent_role      = ifelse(ptyp == 300, "mother", "father")
  )

# Pivot to get mother / father columns per family
parent_wide <- parent_info %>%
  select(fid, parent_role, parent_edu_years, birth_year_parent) %>%
  pivot_wider(
    names_from  = parent_role,
    values_from = c(parent_edu_years, birth_year_parent),
    values_fn   = first   # handle possible duplicates per role
  ) %>%
  rename(
    mother_education  = parent_edu_years_mother,
    father_education  = parent_edu_years_father,
    mother_birth_year = birth_year_parent_mother,
    father_birth_year = birth_year_parent_father
  )

df_twinlife <- df_twinlife %>%
  left_join(parent_wide, by = "fid")


# ---- 4. Education years -----------------------------------------------------

## ---- 4.1 SOEP-style years (eca0130 + eca0230) ----
# Mapping definitions in school_years() / voc_years() (Section 0).

# Use the LATEST available wave for education (prefer higher wid).
eca0130_cols <- grep("^eca0130_wid", names(df_twinlife), value = TRUE)
eca0230_cols <- grep("^eca0230_wid", names(df_twinlife), value = TRUE)

df_twinlife <- df_twinlife %>%
  rowwise() %>%
  mutate(
    latest_eca0130 = {
      vals <- c_across(all_of(eca0130_cols))
      v <- vals[!is.na(vals)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    },
    latest_eca0230 = {
      vals <- c_across(all_of(eca0230_cols))
      v <- vals[!is.na(vals)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    }
  ) %>%
  ungroup()

df_twinlife <- df_twinlife %>%
  mutate(edu_soep = school_years(latest_eca0130) + voc_years(latest_eca0230))


## ---- 4.2 Restrict to age > 25 at education measurement ----
# Education is read from the latest available wave, so the latest wave's
# yea_fq + mon_fq determine the age at education assessment.

# Latest non-NA yea_fq and mon_fq per person
yea_fq_cols <- grep("^yea_fq_wid", names(df_twinlife), value = TRUE)
mon_fq_cols <- grep("^mon_fq_wid", names(df_twinlife), value = TRUE)

df_twinlife <- df_twinlife %>%
  rowwise() %>%
  mutate(
    latest_yea_fq = {
      vals <- c_across(all_of(yea_fq_cols))
      v <- vals[!is.na(vals)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    },
    latest_mon_fq = {
      vals <- c_across(all_of(mon_fq_cols))
      v <- vals[!is.na(vals)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    }
  ) %>%
  ungroup()

# Apply age > 25 filter for SOEP-style education
df_twinlife <- df_twinlife %>%
  mutate(
    education = ifelse(
      !is.na(latest_yea_fq) & !is.na(latest_mon_fq) & !is.na(birth) &
        (latest_yea_fq + latest_mon_fq / 12 - birth) > 25,
      edu_soep,
      NA_real_
    )
  )

# Fallback using wav0100 subsample timing if explicit timing is missing.
# Both subsamples were measured up to ~2024 (F2F5), so the same year is used.
df_twinlife <- df_twinlife %>%
  mutate(
    education = ifelse(
      is.na(education) & !is.na(edu_soep),
      case_when(
        wav0100 == 1 & !is.na(birth) & ((2024 - birth) > 25) ~ edu_soep,
        wav0100 == 2 & !is.na(birth) & ((2024 - birth) > 25) ~ edu_soep,
        TRUE ~ NA_real_
      ),
      education
    )
  )

# ---- 5. Parental education --------------------------------------------------
# Mother / father education years are joined in from Section 3.3. The
# row-mean parental_education is computed here.

df_twinlife <- df_twinlife %>%
  mutate(
    parental_education = rowMeans(
      cbind(mother_education, father_education),
      na.rm = TRUE
    ),
    parental_education = ifelse(
      is.nan(parental_education), NA_real_, parental_education
    )
  )


# ---- 6. East / West classification ------------------------------------------
# Inputs:
#   mig2001 — born in the GDR (0 = no/west, 1 = yes/east); time-constant
#   ewi     — current eastern (1) / western (2) federal state, per wave

## ---- 6.1 born_reg / live_reg ----
# First observed ewi (earliest wave with non-NA ewi)
ewi_cols <- grep("^ewi_wid", names(df_twinlife), value = TRUE)

df_twinlife <- df_twinlife %>%
  rowwise() %>%
  mutate(
    first_ewi = {
      vals <- c_across(all_of(ewi_cols))
      v <- vals[!is.na(vals)]
      if (length(v) == 0) NA_real_ else head(v, 1)
    }
  ) %>%
  ungroup()

df_twinlife <- df_twinlife %>%
  mutate(
    # mig2001 = "born in the GDR (yes/no)" only carries information for
    # births before German reunification. From 1990 onward the GDR no
    # longer existed, so the question has no defined meaning and in the
    # raw TwinLife delivery post-1990 births are overwhelmingly
    # mig2001 = 0 regardless of the federal state the person lived in.
    # mig2001 is therefore suppressed for births in or after 1990 and
    # east_west falls back to live_reg. Where the birth year is unknown,
    # mig2001 keeps precedence (small group, mostly outside the analysis
    # sample). Residual misclassification is accepted; no separate
    # sensitivity analysis is planned.
    # See 02_Phenotype/TwinLife/migration_construction.md.
    born_region = case_when(
      is.na(mig2001)                ~ NA_character_,
      !is.na(birth) & birth >= 1990 ~ NA_character_,
      mig2001 == 1                  ~ "born_east",
      mig2001 == 0                  ~ "born_west",
      TRUE                          ~ NA_character_
    ),
    live_region = case_when(
      is.na(first_ewi) ~ NA_character_,
      first_ewi == 1   ~ "live_east",
      TRUE             ~ "live_west"
    )
  )

# Normalise region strings to bare "east" / "west" via norm_region() (Section 0)
df_twinlife <- df_twinlife %>%
  mutate(
    born_reg = norm_region(born_region),
    live_reg = norm_region(live_region)
  )

## ---- 6.2 east_west (origin if known, else first observed living region) ----
df_twinlife <- df_twinlife %>%
  mutate(
    east_west = case_when(
      !is.na(born_reg)                   ~ born_reg,
      is.na(born_reg) & !is.na(live_reg) ~ live_reg,
      TRUE                               ~ NA_character_
    )
  )

## ---- 6.3 east_west_mig (lifetime: stayer vs mover from born vs live) ----
# Only classify as stayer/mover when BOTH born_reg and live_reg are known.
# If either side is missing, leave as NA rather than fall back to a one-
# sided guess at staying. Every post-1990 birth is therefore NA on this
# variable (their born region is unknown by construction), which is
# correct: cross-GDR/FRG migration cannot be defined for people born
# after the GDR ceased to exist.
# See 02_Phenotype/TwinLife/migration_construction.md.
df_twinlife <- df_twinlife %>%
  mutate(
    east_west_mig = case_when(
      !is.na(born_reg) & !is.na(live_reg) & born_reg == "east" & live_reg == "east" ~ "east_stay",
      !is.na(born_reg) & !is.na(live_reg) & born_reg == "west" & live_reg == "west" ~ "west_stay",
      !is.na(born_reg) & !is.na(live_reg) & born_reg == "east" & live_reg == "west" ~ "east_to_w_move",
      !is.na(born_reg) & !is.na(live_reg) & born_reg == "west" & live_reg == "east" ~ "west_to_e_move",
      TRUE                                                                          ~ NA_character_
    )
  )

## ---- 6.4 east_west_mig_edu (move only counts if observed by age <= 25) ----
# Derive the wid of the first non-missing ewi observation, then look up
# the timing (yea_fq, mon_fq) at that wave to compute age at first ewi.
df_twinlife <- df_twinlife %>%
  rowwise() %>%
  mutate(
    first_wid_ewi = {
      vals   <- c_across(all_of(ewi_cols))
      wids   <- as.numeric(gsub("ewi_wid", "", ewi_cols))
      non_na <- which(!is.na(vals))
      if (length(non_na) == 0) NA_real_ else wids[non_na[1]]
    }
  ) %>%
  ungroup()

# Look up yea_fq and mon_fq at first_wid_ewi via a join to the long data
first_ewi_timing <- df_long %>%
  select(pid, wid, yea_fq, mon_fq) %>%
  filter(!is.na(wid))

df_twinlife <- df_twinlife %>%
  left_join(
    first_ewi_timing %>%
      rename(first_wid_ewi_check = wid,
             first_year          = yea_fq,
             first_month         = mon_fq),
    by = c("pid", "first_wid_ewi" = "first_wid_ewi_check")
  )

df_twinlife <- df_twinlife %>%
  mutate(
    age_at_first_ewi = if_else(
      !is.na(first_year) & !is.na(first_month) & !is.na(birth),
      (first_year + first_month / 12) - birth,
      NA_real_
    )
  )

# Only classify as stayer/mover when we have positive evidence:
#   - stayer requires born_reg == live_reg (both known)  OR  born_reg != live_reg
#     with a known age at first observation > 25 (move happened post-formative)
#   - mover  requires born_reg != live_reg with a known age at first
#     observation ≤ 25 (move observed during the formative window)
# Everything else (either side unknown, or timing unknown for an
# observed move) is NA. See 02_Phenotype/TwinLife/migration_construction.md.
df_twinlife <- df_twinlife %>%
  mutate(
    east_west_mig_edu = case_when(
      # same region in birth and current → stayer in that region
      !is.na(born_reg) & !is.na(live_reg) & born_reg == live_reg ~
        if_else(born_reg == "east", "east_stay", "west_stay"),

      # different region AND move observed at age ≤ 25 → mover during formative years
      !is.na(born_reg) & !is.na(live_reg) & born_reg != live_reg &
        !is.na(age_at_first_ewi) & age_at_first_ewi <= 25 ~
        case_when(
          born_reg == "east" & live_reg == "west" ~ "east_to_w_move",
          born_reg == "west" & live_reg == "east" ~ "west_to_e_move",
          TRUE                                    ~ NA_character_
        ),

      # different region AND move first observed at age > 25 → stayer in birth region
      # (the move happened after the formative-years window)
      !is.na(born_reg) & !is.na(live_reg) & born_reg != live_reg &
        !is.na(age_at_first_ewi) & age_at_first_ewi > 25 ~
        if_else(born_reg == "east", "east_stay", "west_stay"),

      # different region with unknown timing, or either side missing → cannot classify
      TRUE ~ NA_character_
    )
  )


# ---- 7. Anthropometrics -----------------------------------------------------
# Physical measurements are at F2F waves only:
#   bmi (bdy0300):        wids 3 (F2F2), 5 (F2F3, very few), 7 (F2F4), 9 (F2F5)
#   height (bdy0100_hgt): wids 3 (F2F2), 7 (F2F4), 9 (F2F5)
#
# The delivery does not include corrected bmi/height at F2F1, so wids
# 3, 7, 9 are used to maximise coverage.
# Logic:
#   bmi    = average across available F2F waves (3, 7, 9)
#   height = if age at the latest measured wave (wid9) < 20 → take latest
#            available; otherwise average across waves (3, 7, 9)

# Ensure columns exist (create as NA if missing in delivery)
for (col in c("bdy0300_wid3","bdy0300_wid7","bdy0300_wid9",
              "bdy0100_hgt_wid3","bdy0100_hgt_wid7","bdy0100_hgt_wid9",
              "age0100_wid3","age0100_wid7","age0100_wid9")) {
  if (!col %in% names(df_twinlife)) df_twinlife[[col]] <- NA_real_
}

df_twinlife <- df_twinlife %>%
  mutate(
    bmi     = rowMeans(pick(bdy0300_wid3, bdy0300_wid7, bdy0300_wid9), na.rm = TRUE),
    bmi_age = rowMeans(pick(age0100_wid3, age0100_wid7, age0100_wid9), na.rm = TRUE)
  ) %>%
  mutate(
    bmi     = ifelse(is.nan(bmi),     NA_real_, bmi),
    bmi_age = ifelse(is.nan(bmi_age), NA_real_, bmi_age)
  )

df_twinlife <- df_twinlife %>%
  rowwise() %>%
  mutate(
    height = case_when(
      !is.na(age0100_wid9) & age0100_wid9 < 20 ~
        coalesce(bdy0100_hgt_wid9, bdy0100_hgt_wid7, bdy0100_hgt_wid3),
      TRUE ~ {
        m <- mean(c_across(c(bdy0100_hgt_wid3, bdy0100_hgt_wid7, bdy0100_hgt_wid9)),
                  na.rm = TRUE)
        if (is.nan(m)) NA_real_ else m
      }
    ),
    height_age = case_when(
      !is.na(age0100_wid9) & age0100_wid9 < 20 ~
        coalesce(age0100_wid9, age0100_wid7, age0100_wid3),
      TRUE ~ {
        m <- mean(c_across(c(age0100_wid3, age0100_wid7, age0100_wid9)),
                  na.rm = TRUE)
        if (is.nan(m)) NA_real_ else m
      }
    )
  ) %>%
  ungroup()


# ---- 9. Final harmonised dataset --------------------------------------------

names(df_twinlife)[names(df_twinlife) == "birth"]  <- "birth_year"
names(df_twinlife)[names(df_twinlife) == "ptyp_1"] <- "ptyp"

df_final <- df_twinlife[, c(
  "fid", "pid", "cgr", "ptyp", "gender", "birth_year",
  "education",
  "mother_education", "father_education",
  "mother_birth_year", "father_birth_year",
  "parental_education",
  "east_west", "east_west_mig", "east_west_mig_edu",
  "bmi", "bmi_age", "height", "height_age"
)]

# Save final data — .sav keeps column descriptions; .rda is canonical
out_sav <- file.path(out_dir, sprintf("EduGxE_DE_TwinLife_Phenotype_%s.sav", date_tag))
out_rda <- file.path(out_dir, sprintf("EduGxE_DE_TwinLife_Phenotype_%s.rda", date_tag))

write_sav(df_final, out_sav)
save(df_final, file = out_rda, compress = "xz", version = 3)

cat("\nSaved:", out_sav, "\n")
cat("Saved:", out_rda, "\n")


# ---- 10. Aggregate sanity check ---------------------------------------------

cat("\n======================================================\n")
cat("=== TwinLife phenotype summary ===\n")
cat("======================================================\n")
cat("Total N:", nrow(df_final), "\n\n")

numeric_vars <- c("birth_year", "education",
                  "mother_education", "father_education",
                  "mother_birth_year", "father_birth_year",
                  "parental_education",
                  "bmi", "bmi_age", "height", "height_age")

for (v in numeric_vars) {
  vals <- df_final[[v]]
  cat(sprintf("  %-25s  N=%5d  M=%8.3f  SD=%7.3f\n",
              v,
              sum(!is.na(vals)),
              mean(vals, na.rm = TRUE),
              sd(vals,   na.rm = TRUE)))
}

cat("\nGender:\n")
print(table(df_final$gender, useNA = "ifany"))

cat("\nEast/West:\n")
print(table(df_final$east_west, useNA = "ifany"))

# east_west split by the GDR-existence cutoff (born <1990 vs >=1990):
# post-1990 born are classified via live_reg (first_ewi) only, since
# mig2001 has no defined meaning after reunification.
cat("\nEast/West x birth era (pre/post-1990 GDR-existence cutoff):\n")
.era_1990 <- cut(df_final$birth_year,
                 breaks = c(-Inf, 1989.99, Inf),
                 labels = c("born <1990", "born >=1990"))
print(addmargins(table(east_west = df_final$east_west,
                       era      = .era_1990,
                       useNA    = "ifany")))

cat("\nEast/West Migration:\n")
print(table(df_final$east_west_mig, useNA = "ifany"))

cat("\nEast/West Migration (edu):\n")
print(table(df_final$east_west_mig_edu, useNA = "ifany"))

# By-region education summary
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

cat("\n=== Done ===\n")
