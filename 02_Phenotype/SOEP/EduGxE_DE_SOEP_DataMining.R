# =============================================================================
# GxE on Education in Germany — SOEP phenotype data mining
# -----------------------------------------------------------------------------
# Purpose      : Build the harmonised SOEP-IS phenotype file feeding 03_Merge.
# Inputs       : - GSOEP_merge_20230314.csv  (processed SOEP-IS sample)
#                - bioparen.dta              (parent SES, education codes, birth)
#                - p.dta                     (panel-respondent education item)
#                - bioage.dta                (school-age trajectories)
#                - ppfad.dta                 (birth year/month, 1989 location)
#                - hbrutto.dta               (household sampreg per year)
# Output       : EduGxE_DE_SOEP_Phenotype_<YYYYMMDD>.rda
# Dependencies : dplyr, data.table, haven, tidyr, purrr
# Note on structure:
#   The canonical 10-section template is non-monotonic for SOEP because the
#   data flow chains: demographics (Section 3) need columns that only land
#   in `df` after the ppfad merge that closes Section 2. Section dividers
#   therefore mark the canonical block currently being built, not strict
#   left-to-right order.
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(data.table)
library(haven)
library(tidyr)
library(purrr)

# Paths and date stamp ---------------------------------------------------------
soep_processed <- "${SOEP_SHARE}/private/data/001_SOEP_processed_data/007_GermanGenetics"
soep_raw       <- "${SOEP_SHARE}/private/data/001_SOEP_rawdata/002_SOEP_Pheno/001_SOEP-IS/soep-is.2020_stata_en"
out_dir        <- soep_processed
date_tag       <- format(Sys.Date(), "%Y%m%d")


# ---- 1. Load raw data -------------------------------------------------------

    soepis <- fread(file.path(soep_processed, "GSOEP_merge_20230314.csv"))

    # Drop rows 2495–2507: they carry only IDs (everything else NA) but
    # the methylation pipeline still has data for them. Keep IDs for
    # cross-reference if needed, but exclude from downstream analysis.
    soepis <- soepis[-c(2495:2507), ]

    # Recode and rename 'male' to 'gender'
    soepis$gender <- ifelse(soepis$male == 1, "Male", "Female")

    # Remove the source 'male' variable
    soepis$male <- NULL

    # Parental SES and education codes
    bioparen <- read_dta(file.path(soep_raw, "bioparen.dta"))

    # Birth month and 1989 East/West location
    ppfad <- read_dta(file.path(soep_raw, "ppfad.dta"))

    # Household sampreg per year (for east/west panel reconstruction)
    hbrutto <- read_dta(file.path(soep_raw, "hbrutto.dta"))
    
      
# Variable selections from the auxiliary .dta files ----

# Parental SES + education codes from bioparen
parents <- select(
  bioparen, pid,
  mpid, fpid,           # parent IDs
  valter, malter,       # parent birth years (or ages-at-child-birth, see Section 3)
  vsbil, msbil,         # level of education (general schooling)
  vbbil, mbbil          # type of education (vocational / tertiary)
)

# Birth year/month + 1989 location from ppfad
ppfad <- select(
  ppfad, pid,
  gebjahr, gebmonat, loc1989
)




# ---- 2. Person-level scaffold (build parent_genes; soepis stays the carrier) ----

dim(soepis)
parent_genes <- left_join(soepis, parents, by = "pid")
dim(parent_genes)

# ---- 5. Parental education --------------------------------------------------
        # Generate years-of-education score for both parents.
        #
        # SOEP/bioparen convention: m* = mother (Mutter), v* = father (Vater).
        parent_edu <- mutate(parent_genes,
                                   m.edu.years =  case_when(msbil == 6 ~ 7,
                                                            msbil == 1 ~ 9,
                                                            msbil == 2 ~ 10,
                                                            msbil == 3 ~ 12,
                                                            msbil == 4 ~ 13,
                                                            msbil == 5 ~ 10,

                                                            msbil == 0 ~ 0),
                                   v.edu.years =  case_when(vsbil == 6 ~ 7,
                                                            vsbil == 1 ~ 9,
                                                            vsbil == 2 ~ 10,
                                                            vsbil == 3 ~ 12,
                                                            vsbil == 4 ~ 13,
                                                            vsbil == 5 ~ 10,

                                                            vsbil == 0 ~ 0))
        
        #add years of VOCATIONAL or UNIVERSITAL education score for both parents
        
        #add 1.5 years for apprenticeship and civil servants
        parent_edu$m.edu.years <- ifelse(parent_edu$mbbil %in% c(20:25,28,40) & parent_edu$m.edu.years != 0, 
                                          parent_edu$m.edu.years +1.5, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$vbbil %in% c(20:25,28,40) & parent_edu$v.edu.years != 0,
                                          parent_edu$v.edu.years +1.5, parent_edu$v.edu.years)
        
        #add 2 years for technical and health schools
        parent_edu$m.edu.years <- ifelse(parent_edu$mbbil %in% c(26:27) & parent_edu$m.edu.years != 0, 
                                          parent_edu$m.edu.years +2, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$vbbil %in% c(26:27) & parent_edu$v.edu.years != 0,
                                          parent_edu$v.edu.years +2, parent_edu$v.edu.years)
        
        # Add 3 years for higher technical college (FHS / Ingenieurschule).
        # SOEP/bioparen code 30 = Fachhochschule / Ingenieurschule
        # (3-year qualification beyond apprenticeship). Code 28
        # (Beamtenausbildung / civil-service training) takes +1.5 in the
        # apprenticeship block above, per SOEP's pgbilzeit/pgbilzt rule.
        parent_edu$m.edu.years <- ifelse(parent_edu$mbbil %in% c(30) & parent_edu$m.edu.years != 0,
                                          parent_edu$m.edu.years + 3, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$vbbil %in% c(30) & parent_edu$v.edu.years != 0,
                                          parent_edu$v.edu.years + 3, parent_edu$v.edu.years)
        
        #add 5 years for university and college
        parent_edu$m.edu.years <- ifelse(parent_edu$mbbil %in% c(31:32) & parent_edu$m.edu.years != 0, 
                                          parent_edu$m.edu.years +5, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$vbbil %in% c(31:32) & parent_edu$v.edu.years != 0,
                                          parent_edu$v.edu.years +5, parent_edu$v.edu.years)
        
        
        #infer "do not know" responses from mbbil/ vbbil
        #add minimum years of education required for occupational training; 
        # eg. 13+5 for university or 9 + 1.5 for apprentice ship
        
        #university 18 years
        parent_edu$m.edu.years <- ifelse(parent_edu$mbbil == 32 & parent_edu$m.edu.years == 0, 18, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$vbbil == 32 & parent_edu$v.edu.years == 0, 18, parent_edu$v.edu.years)
        
        #apprenticeship / vocational degree
        parent_edu$m.edu.years <- ifelse(parent_edu$mbbil %in% c(20:25) & parent_edu$m.edu.years == 0, 10.5, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$vbbil %in% c(20:25) & parent_edu$v.edu.years == 0, 10.5, parent_edu$v.edu.years)
        
        #Replace remaining 0 (no vocational training, or do not know) with NA
        parent_edu$m.edu.years <- ifelse(parent_edu$m.edu.years == 0, NA, parent_edu$m.edu.years)
        parent_edu$v.edu.years <- ifelse(parent_edu$v.edu.years == 0, NA, parent_edu$v.edu.years)
        
        
        
        parent_edu <-
        parent_edu %>%
          rowwise() %>%
        mutate( eduparents = mean(c(m.edu.years , v.edu.years), na.rm = T))
        parent_edu$eduparents[is.na(parent_edu$eduparents) ] <- NA #weird format of NA (NaN). This code makes it congurent with other variables
        
        
        parent_edu$mpid[parent_edu$mpid < 0] <- NA
        parent_edu$fpid[parent_edu$fpid < 0] <- NA
        
        #Add parental education to main dataset
        soepis  <- left_join(soepis,parent_edu[,c("pid","m.edu.years","v.edu.years", "eduparents", "mpid","fpid", "malter", "valter")], by = "pid")
        

# ---- 2. (cont.) close the scaffold: ppfad -> df -----------------------------
# After this, `df` is the working dataframe with everything Section 6 below
# needs (loc1989, gebjahr, gebmonat).

df <- left_join(soepis, ppfad, by = "pid")


# ---- 6. East / West classification ------------------------------------------
# The most complex east/west logic of any cohort: per-person residency
# reconstructed from year-by-year household sampreg observations between
# age 15 and end of schooling (owneduyears + 7). Sub-sections below.

## ---- 6.1 Helpers + hbrutto sampreg → ew ----

loc1989_to_ew <- function(code) {
  dplyr::case_when(
    code == 1 ~ "east",
    code == 2 ~ "west",
    code == 3 ~ NA_character_,
    TRUE      ~ NA_character_
  )
}

# Restrict hbrutto; map sampreg -> ew (household × year)
hbrutto_sub <- hbrutto %>%
  semi_join(df %>% distinct(hid), by = "hid") %>%
  select(hid, syear, sampreg) %>%
  mutate(
    sampreg = na_if(sampreg, -1L),
    sampreg = na_if(sampreg, -2L),
    sampreg = na_if(sampreg, -3L),
    sampreg = na_if(sampreg, -4L),
    sampreg = na_if(sampreg, -5L),
    sampreg = na_if(sampreg, -6L),
    ew = case_when(
      sampreg == 1 ~ "west",
      sampreg == 2 ~ "east",
      TRUE         ~ NA_character_
    )
  ) %>%
  arrange(hid, syear)

# First available sampreg EW and its year per household
hh_first <- hbrutto_sub %>%
  group_by(hid) %>%
  summarise(
    first_ew  = dplyr::first(ew[!is.na(ew)], default = NA_character_),
    first_yr  = dplyr::first(syear[!is.na(ew)], default = NA_integer_),
    .groups = "drop"
  )

## ---- 6.2 Person-level scaffold for east/west: df_person ----
df_person <- df %>%
  transmute(
    pid, hid,
    loc1989     = suppressWarnings(as.numeric(loc1989)),
    loc1989_ew  = loc1989_to_ew(loc1989),
    gebjahr     = suppressWarnings(as.numeric(gebjahr)),
    owneduyears = suppressWarnings(as.numeric(owneduyears)),
    # finish year for the “education migration” rule (owneduyears + 6)
    finish_yr   = if_else(!is.na(gebjahr) & !is.na(owneduyears),
                          as.integer(gebjahr + floor(owneduyears + 6)),
                          NA_integer_)
  ) %>%
  left_join(hh_first, by = "hid")

# Distinct person–household key
pid_hh <- df %>% distinct(pid, hid, loc1989)

## ---- 6.3 First year where sampreg EW ≠ loc1989 EW (per person) ----
first_diff <- pid_hh %>%
  left_join(hbrutto_sub, by = "hid", relationship = "many-to-many") %>%  # expected m:m
  mutate(loc1989_ew = loc1989_to_ew(as.numeric(loc1989))) %>%
  arrange(pid, syear) %>%
  group_by(pid) %>%
  summarise(
    first_diff_yr = suppressWarnings(min(syear[!is.na(ew) & ew != loc1989_ew], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(first_diff_yr = ifelse(is.infinite(first_diff_yr), NA_integer_, first_diff_yr))

df_person <- df_person %>%
  left_join(first_diff, by = "pid")

## ---- 6.4 east_west (strict schooling window: age 15 .. owneduyears+7) ----
# For pre-1989 years use loc1989; for 1989+ proxy with household first_ew

pers_years <- df_person %>%
  filter(!is.na(gebjahr), !is.na(owneduyears)) %>%
  mutate(end_age_ew = as.integer(floor(owneduyears + 7))) %>%
  filter(end_age_ew >= 15) %>%
  mutate(syear_seq = map2(gebjahr, end_age_ew, ~ seq(from = .x + 15L, to = .x + .y, by = 1L))) %>%
  select(pid, loc1989, loc1989_ew, hh_first_ew = first_ew, syear_seq) %>%
  unnest_longer(syear_seq, values_to = "syear") %>%
  mutate(
    ew = if_else(syear < 1989L, loc1989_ew, hh_first_ew)
  )

pers_ew_class <- pers_years %>%
  group_by(pid) %>%
  summarise(
    n_obs_school = sum(!is.na(ew)),
    all_east     = n_obs_school > 0 && all(ew == "east", na.rm = TRUE),
    all_west     = n_obs_school > 0 && all(ew == "west", na.rm = TRUE),
    east_west    = case_when(
      all_east ~ "east",
      all_west ~ "west",
      TRUE     ~ NA_character_     # any mixture/unknown -> exclude
    ),
    .groups = "drop"
  )

excluded_n <- sum(is.na(pers_ew_class$east_west))

## ---- 6.5 east_west_mig (loc1989 vs first sampreg EW + later moves) ----
# Harmonized NA handling with Step 5:
# - Use early mismatch if both loc1989_ew and first_ew observed.
# - Else, if a later panel change after first_yr exists, classify as mover.
# - Else, fallback to schooling-residence (east_west) to classify as stayer.
# - Only remain NA if no rule can classify (e.g., no panel change and no east_west).

# --- Panel changes without pbrutto (household-based) --------------------------
# 1) First EW flip per household
hh_change <- hbrutto_sub %>%
  arrange(hid, syear) %>%
  group_by(hid) %>%
  mutate(prev_ew = dplyr::lag(ew)) %>%
  filter(!is.na(prev_ew), !is.na(ew), ew != prev_ew) %>%
  slice_head(n = 1) %>%
  transmute(
    hid,
    first_change_yr   = as.integer(syear),
    first_change_from = prev_ew,
    first_change_to   = ew
  ) %>%
  ungroup()

# 2) Map to persons (earliest flip across any household they appear in)
panel_changes <- df %>%
  distinct(pid, hid) %>%
  left_join(hh_change, by = "hid") %>%
  arrange(pid, first_change_yr) %>%
  group_by(pid) %>%
  slice_head(n = 1) %>%                 # earliest observed change per person
  ungroup()

# 3) Ensure all pids are present (NA when no change observed)
panel_changes <- df %>% distinct(pid) %>%
  left_join(panel_changes, by = "pid")

mig_any <- df_person %>%
  left_join(
    panel_changes %>%
      select(pid, first_change_yr, first_change_from, first_change_to),
    by = "pid"
  ) %>%
  # fallback stayer from schooling-residence (as in Step 5)
  left_join(pers_ew_class %>% select(pid, east_west), by = "pid") %>%
  mutate(
    move_dir_diff = case_when(
      loc1989_ew == "east" & first_ew == "west" ~ "east_to_w_move",
      loc1989_ew == "west" & first_ew == "east" ~ "west_to_e_move",
      TRUE                                      ~ NA_character_
    ),
    move_dir_panel = case_when(
      first_change_from == "east" & first_change_to == "west" ~ "east_to_w_move",
      first_change_from == "west" & first_change_to == "east" ~ "west_to_e_move",
      TRUE                                                    ~ NA_character_
    ),
    moved_after_first = !is.na(first_change_yr) & !is.na(first_yr) & (first_change_yr > first_yr),
    stayer_from_school = case_when(
      east_west == "east" ~ "east_stay",
      east_west == "west" ~ "west_stay",
      TRUE                ~ NA_character_
    )
  ) %>%
  transmute(
    pid,
    east_west_mig = case_when(
      # (1) Early mismatch available -> mover by loc1989_ew -> first_ew
      !is.na(move_dir_diff) ~ move_dir_diff,
      
      # (2) Otherwise, if a later panel change AFTER first_yr -> mover by panel direction
      moved_after_first & !is.na(move_dir_panel) ~ move_dir_panel,
      
      # (3) Otherwise, if both loc1989_ew and first_ew agree -> stayer by that agreement
      !is.na(loc1989_ew) & !is.na(first_ew) & loc1989_ew == first_ew ~
        if_else(first_ew == "east", "east_stay", "west_stay"),
      
      # (4) Fallback: stayer from schooling-residence classification
      !is.na(stayer_from_school) ~ stayer_from_school,
      
      # If nothing applies, remain NA
      TRUE ~ NA_character_
    )
  )


## ---- 6.6 east_west_mig_edu (timed rules per spec) ----
# Rules:
# A) finish_yr < 1989       -> stayer
# B) finish_yr > 1989 & loc1989_ew != first_ew & first_yr == first_diff_yr
#                           -> mover (dir: loc1989_ew -> first_ew)
# C) first_yr != first_diff_yr:
#       mover only if finish_yr < first_diff_yr (dir: panel from->to)
#       else stayer
# Fallback: stayer from schooling-residence (east_west)

mig_edu <- df_person %>%
  # bring panel change info for direction in rule C
  left_join(
    panel_changes %>% select(pid, first_change_yr, first_change_from, first_change_to),
    by = "pid"
  ) %>%
  # schooling-residence fallback
  left_join(pers_ew_class %>% select(pid, east_west), by = "pid") %>%
  mutate(
    finished_pre1989 = !is.na(finish_yr) & finish_yr < 1989L,
    early_mismatch   = !is.na(loc1989_ew) & !is.na(first_ew) & (loc1989_ew != first_ew),
    same_time        = !is.na(first_yr) & !is.na(first_diff_yr) & (first_yr == first_diff_yr),
    diff_time        = !is.na(first_yr) & !is.na(first_diff_yr) & (first_yr != first_diff_yr),
    
    move_dir_diff = case_when(
      loc1989_ew == "east" & first_ew == "west" ~ "east_to_w_move",
      loc1989_ew == "west" & first_ew == "east" ~ "west_to_e_move",
      TRUE                                      ~ NA_character_
    ),
    move_dir_panel = case_when(
      first_change_from == "east" & first_change_to == "west" ~ "east_to_w_move",
      first_change_from == "west" & first_change_to == "east" ~ "west_to_e_move",
      TRUE                                                    ~ NA_character_
    ),
    stayer_from_school = case_when(
      east_west == "east" ~ "east_stay",
      east_west == "west" ~ "west_stay",
      TRUE                ~ NA_character_
    ),
    
    east_west_mig_edu = case_when(
      # A) finished education before 1989 -> stayer
      finished_pre1989 ~ stayer_from_school,
      
      # B) finish after 1989 AND immediate early mismatch at first_yr -> mover (diff direction)
      !is.na(finish_yr) & finish_yr > 1989L & early_mismatch & same_time ~ move_dir_diff,
      
      # C) first_yr != first_diff_yr -> mover only if finish_yr < first_diff_yr (panel direction)
      diff_time & !is.na(finish_yr) & finish_yr < first_diff_yr ~ move_dir_panel,
      
      # Otherwise -> stayer (fallback to schooling-residence)
      TRUE ~ stayer_from_school
    )
  ) %>%
  select(pid, east_west_mig_edu)


## ---- 6.7 Merge east/west columns back into df ----
df <- df %>%
  left_join(pers_ew_class %>% select(pid, east_west), by = "pid") %>%
  left_join(mig_any, by = "pid") %>%
  left_join(mig_edu, by = "pid")

# Quick QA counts
table(df$east_west,         useNA = "ifany")
table(df$east_west_mig,     useNA = "ifany")
table(df$east_west_mig_edu, useNA = "ifany")

# Excluded from east_west (any mixture/unknown during schooling window):
excluded_n


# ---- 3. Demographics --------------------------------------------------------
# birth_year derived once df has gebjahr / gebmonat from the ppfad merge.
# (gender already set in Section 1; mother/father birth years computed in
# Section 9 below from malter/valter.)

df <- df %>%
  mutate(birth_yr = gebjahr + gebmonat / 12)


# ---- 4. Education years (rename) + 9. (renames toward harmonised columns) ---
# Rename internal column names to the cross-cohort harmonised names.
names(df)[names(df) == 'maxedu']   <- 'education'        # Section 4
names(df)[names(df) == 'birth_yr'] <- 'birth_year'
names(df)[names(df) == 'cmheight'] <- 'height'

# Rename individual parent education for harmonization across cohorts
names(df)[names(df) == 'm.edu.years'] <- 'mother_education'
names(df)[names(df) == 'v.edu.years'] <- 'father_education'

# ---- 3. Demographics (cont.): mother / father birth years ------------------
# Compute parent birth years from bioparen's malter/valter
      # bioparen: malter = Geburtsjahr der Mutter, valter = Geburtsjahr des Vaters
      # Values >100 are birth years directly; values <100 would be ages at child's birth
      df <- df %>%
        mutate(
          malter_clean = ifelse(!is.na(malter) & malter > 0, malter, NA_real_),
          valter_clean = ifelse(!is.na(valter) & valter > 0, valter, NA_real_),
          mother_birth_year = case_when(
            is.na(malter_clean) ~ NA_real_,
            malter_clean > 100  ~ malter_clean,
            TRUE                ~ gebjahr - malter_clean
          ),
          father_birth_year = case_when(
            is.na(valter_clean) ~ NA_real_,
            valter_clean > 100  ~ valter_clean,
            TRUE                ~ gebjahr - valter_clean
          )
        ) %>%
        select(-malter_clean, -valter_clean)

# ---- 5. (cont.) parental_education = row-wise mean -------------------------
df <- df %>%
  mutate(
    # row-wise mean of available parental education (uses one if the other is missing)
    parental_education = rowMeans(select(., mother_education, father_education), na.rm = TRUE)
    # Not standardized here: standardization happens at analysis time, within
    # parent birth-cohort bins, across pooled datasets
  )


# ---- 7. Anthropometrics -----------------------------------------------------
# bmi and height are inherited from the upstream SOEP-IS file (renamed
# above). Only the age-at-measurement is set here; SOEP records BMI/height
# centrally so we adopt ageJan12019 as the reference age.
df <- df %>%
  mutate(
    bmi_age    = as.numeric(ageJan12019),
    height_age = as.numeric(ageJan12019)
  )
      
      
      
# ---- 9. Final harmonised dataset --------------------------------------------

soepis_final <- df[, c(
  "pid", "ID", "hid",
  "gender", "birth_year", "education",
  "mother_education", "father_education",
  "mother_birth_year", "father_birth_year",
  "parental_education",
  "height", "bmi", "height_age", "bmi_age",
  "east_west", "east_west_mig", "east_west_mig_edu"
)]

# Save final data — .sav keeps column descriptions; .rda is the canonical format
out_sav <- file.path(out_dir,
                     sprintf("EduGxE_DE_SOEP_Phenotype_%s.sav", date_tag))
out_rda <- file.path(out_dir,
                     sprintf("EduGxE_DE_SOEP_Phenotype_%s.rda", date_tag))

# write_sav() can't represent SAS-style character-tagged NAs that may have
# travelled in from the Stata files. Strip them by replacing every NA cell
# with a plain NA_real_; the .rda keeps the original.
soepis_final_for_sav <- soepis_final %>%
  mutate(across(where(is.numeric),
                ~ replace(as.numeric(.), is.na(.), NA_real_)))
write_sav(soepis_final_for_sav, out_sav)

save(soepis_final, file = out_rda, compress = "xz", version = 3)
cat("Saved:", out_sav, "\n")
cat("Saved:", out_rda, "\n")


# ---- 10. Aggregate sanity check --------------------------------------------
# Quick aggregate-only summary printed on each run.

cat("\n=== SOEP phenotype summary ===\n")
cat("Total N:", nrow(soepis_final), "\n\n")

cat("gender:\n")
print(table(soepis_final$gender, useNA = "ifany"))

cat("\neast_west:\n")
print(table(soepis_final$east_west, useNA = "ifany"))

cat("\neast_west_mig:\n")
print(table(soepis_final$east_west_mig, useNA = "ifany"))

cat("\n--- By-region education summary ---\n")
soepis_final %>%
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
