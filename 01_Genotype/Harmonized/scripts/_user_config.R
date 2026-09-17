# =============================================================================
# _user_config.R — read per-user settings (ANALYSIS_ROOT, ...) for R scripts
# =============================================================================
# Resolves a user_config key (e.g. ANALYSIS_ROOT) by checking:
#   1. The process environment (Sys.getenv).
#   2. user_config.sh next to this file's parent (../user_config.sh from
#      scripts/, since this file lives in scripts/).
# Errors with a clear instruction if the key is missing in both places.
#
# Usage from another R script in scripts/:
#   source(file.path(dirname(sys.frame(1)$ofile %||% commandArgs(trailingOnly = FALSE)[
#     grep("^--file=", commandArgs(trailingOnly = FALSE))]), "_user_config.R"))
#   analysis_dir <- get_user_config_value("ANALYSIS_ROOT")
# =============================================================================

get_user_config_value <- function(key, required = TRUE) {
  v <- Sys.getenv(key, unset = "")
  if (nzchar(v)) return(v)

  # Locate this script's directory; user_config.sh lives one level above.
  this_file <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      sub("^--file=", "", file_arg[[1]])
    } else if (!is.null(sys.frame(1)$ofile)) {
      sys.frame(1)$ofile
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)

  candidates <- character(0)
  if (!is.na(this_file)) {
    script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
    candidates <- c(candidates,
                    file.path(script_dir, "user_config.sh"),
                    file.path(script_dir, "..", "user_config.sh"),
                    file.path(script_dir, "..", "..", "user_config.sh"))
  }
  candidates <- c(candidates,
                  "user_config.sh",
                  "01_Genotype/Harmonized/user_config.sh")

  for (path in candidates) {
    if (file.exists(path)) {
      lines <- tryCatch(readLines(path, warn = FALSE),
                        error = function(e) character(0))
      pattern <- sprintf('^[[:space:]]*export[[:space:]]+%s=("?)([^"]*)\\1[[:space:]]*$', key)
      hits <- regmatches(lines, regexec(pattern, lines))
      vals <- vapply(hits, function(m) if (length(m) >= 3) m[[3]] else NA_character_,
                     character(1))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) > 0) {
        # Expand $HOME / ${HOME} so the returned path is usable.
        v <- tail(vals, 1)
        v <- gsub("\\$\\{?HOME\\}?", Sys.getenv("HOME"), v)
        v <- gsub("\\$\\{?USER\\}?", Sys.getenv("USER"), v)
        return(v)
      }
    }
  }

  if (required) {
    stop(sprintf(
      paste0(
        "%s is not configured. To run this script, either:\n",
        "  1. Export it in your shell: `export %s=\"/path/to/...\"`\n",
        "  2. Or set %s in 01_Genotype/Harmonized/user_config.sh\n",
        "     (written by _bootstrap_user_config.sh on first run).\n",
        "  3. Or `source` the existing user_config.sh in your shell first:\n",
        "       source 01_Genotype/Harmonized/user_config.sh\n",
        "     then re-run this R script."
      ),
      key, key, key
    ), call. = FALSE)
  }
  ""
}
