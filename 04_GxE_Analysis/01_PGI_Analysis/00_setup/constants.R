# constants.R
# Shared analytic constants for the 01_PGI_Analysis workflow.
# Sourced by every section's analysis and report script.

# ---- Standardization settings ----
KERNEL_BW             <- 5       # Gaussian kernel bandwidth (years) for z-scoring
PARENT_GENERATION_GAP <- 28      # approximate years between parent and child cohorts
NP_BANDWIDTH          <- 8       # bandwidth for non-parametric kernel regression overlays
MAX_ABS_PGI_EDU       <- 3       # exclude extreme global-PGI outliers from analysis

# ---- Descriptives settings ----
# Minimum finite N per group for a descriptive two-group test (Welch t /
# Hedges' g) to be computed in 00_Descriptives/. Cells below this floor
# report "n/a" rather than an unstable estimate.
DESC_MIN_CELL_N <- 30

# ---- GAM settings (analysis plan: k = 8, tp default, cr sensitivity) ----
K_DEFAULT  <- 8
BS_DEFAULT <- "tp"   # thin-plate spline
K_SENS     <- 8
BS_SENS    <- "cr"   # cubic regression spline (sensitivity)
N_SIM      <- 5000   # draws for simultaneous CIs

# ---- Plot colour / linetype constants ----
COL_EAST <- "#D62728"   # red  for East (slope plots)
COL_WEST <- "#1F77B4"   # blue for West (slope plots)
# Kernel-diagnostic plots use a separate palette so they don't read as
# East/West slope plots at a glance (orange/purple, distinct from red/blue).
COL_KERNEL_EAST <- "#E66101"   # orange
COL_KERNEL_WEST <- "#5E3C99"   # purple
LT_MALE   <- "solid"
LT_FEMALE <- "dashed"
