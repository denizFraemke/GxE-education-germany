# paths.R — single source of truth for DATA_ROOT in the 02_SNP_Heritability
# stage.
#
# Every script here resolves the analysis-share root through data_root():
# set the DATA_ROOT environment variable to override it (e.g. point at a
# replication folder); otherwise the placeholder default below is used.
#
# Change the default in THIS one place only; never inline it in individual
# scripts.
data_root <- function() {
  Sys.getenv(
    "DATA_ROOT",
    unset = "/path/to/data_root"
  )
}
