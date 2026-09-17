#!/usr/bin/env bash
# =============================================================================
# setup_tardis_data.sh — Set up data directories on Tardis HPC
# =============================================================================
# This script sets up the data directory structure and checks for required
# files. It provides two strategies for getting data onto the cluster:
#
#   Strategy A (recommended): Mount the project shares via CIFS
#   Strategy B (fallback):    Copy data from the RStudio server via scp
#
# Run this ONCE before starting the pipeline. If files are missing, follow
# the instructions printed at the end.
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# Load TARDIS_USER, TARDIS_BIN_DIR, optional MPIB_SERVER, etc. from
# user_config.sh (or run the interactive bootstrap on first invocation).
# shellcheck source=../_bootstrap_user_config.sh
source "${PIPELINE_DIR}/_bootstrap_user_config.sh"

source "${PIPELINE_DIR}/config.sh"

# =============================================================================
# RSTUDIO SERVER PATHS (where the data lives on the RStudio server)
# =============================================================================
# These are used to generate scp commands if CIFS mounts aren't available.
# MPIB_SERVER comes from user_config.sh; the default below is a placeholder.
MPIB_SERVER="${MPIB_SERVER:-${TARDIS_USER}@<institution-rstudio-host>}"

# Formatted summary statistics on the RStudio server (scp fallback only)
# Set this to wherever your formatted sumstats live on the RStudio server
# if you want to scp them to the cluster. Leave empty to disable the fallback.
#
# CANONICAL SOURCE (the pre-formatted per-chr SBayesR files for every trait):
#   data share -> Projects/04_data_analysis/013-PGS/data/formatted/
# Its subdirs match NEEDED_SUMSTATS_DIRS below exactly (EA4_excl_SHIP,
# EA4_excl_BASEII, Cog_Malanchini, NCog_Malanchini, Height_GIANT). The raw
# Malanchini Cog/NCog GWAS that feed Cog_/NCog_Malanchini live alongside them
# in ../sumstats/ (Malanchini_Cog_ext.txt, Malanchini_NCog_ext.txt). Easiest
# staging: rsync those subdirs from the mounted share straight to the cluster.
#
# EXCEPTION — Height_GIANT: the pre-formatted copy under 013-PGS/data/formatted
# is MALFORMED (extra leading CHR column + unfiltered NA rows). Do NOT stage it
# from there. Instead regenerate it from the RAW GIANT file with
# scripts/format_height_sumstats.sh (raw lives on the share at
# Projects/03_data/009_SUMSTATS/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz).
MPIB_FORMATTED=""

# Genotype data on the RStudio server (project shares), which typically mounts
# the shares under ${SHARE_MOUNT_ROOT}/. Override in the environment; the
# default below is a placeholder.
_MPIB_SRT="${SHARE_MOUNT_ROOT:-/path/to/share_mount_root}"
MPIB_GENO_BASEII="${_MPIB_SRT}/<base2-share>/private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY"
MPIB_GENO_SHIP0="${_MPIB_SRT}/<ship-share>/private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose"
MPIB_GENO_SHIPTD="${_MPIB_SRT}/<ship-share>/private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall"
MPIB_GENO_SOEP="${_MPIB_SRT}/<soep-share>/private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr1.dose.vcf"  # encrypted; see setup_soep_geno.sh
MPIB_GENO_TWINLIFE="${_MPIB_SRT}/<twinlife-share>/private/data/Gendata/001_processed/TwinLife.b37.chrall"
unset _MPIB_SRT

# LD reference on the RStudio server (scp fallback only). Set this if you want the
# setup to propose an scp command; leave empty otherwise.
MPIB_LD_REF=""

# =============================================================================
# CIFS MOUNT PATHS ON TARDIS
# =============================================================================
# If you've requested CIFS mounts (via IT), they'll appear here.
# Defaults to ${HOME} on Tardis — override if your mounts live elsewhere.
TARDIS_MOUNT_BASE="${TARDIS_MOUNT_BASE:-${HOME}}"

# Possible locations for formatted sumstats via CIFS. Add extra candidates
# here if your site mounts the sumstats somewhere else.
CIFS_FORMATTED=(
    "${TARDIS_MOUNT_BASE}/formatted_sumstats"
)

# Possible locations for genotype data via CIFS
CIFS_SRT_BASE_CANDIDATES=(
    "${TARDIS_MOUNT_BASE}/shares"
    "${TARDIS_MOUNT_BASE}/SRT"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
found_count=0
error_count=0
scp_commands=()

create_symlink() {
    local src="$1"
    local dst_prefix="$2"
    local label="$3"

    local linked=0
    for ext in .bed .bim .fam; do
        if [[ -f "${src}${ext}" ]]; then
            ln -sf "${src}${ext}" "${dst_prefix}${ext}" 2>/dev/null || \
                cp "${src}${ext}" "${dst_prefix}${ext}"
            linked=1
        fi
    done
    # Also link .log and .nosex if they exist
    for ext in .log .nosex; do
        [[ -f "${src}${ext}" ]] && ln -sf "${src}${ext}" "${dst_prefix}${ext}" 2>/dev/null || true
    done
    if [[ $linked -eq 1 ]]; then
        log_info "  OK: $label"
        ((found_count++))
    else
        log_error "  MISSING: $label — no .bed/.bim/.fam at ${src}.*"
        ((error_count++))
    fi
}

# =============================================================================
# START SETUP
# =============================================================================

log_info "============================================================"
log_info "Setting up Tardis data directories for PGI Computation"
log_info "============================================================"
log_info ""

# --- Create target directories -----------------------------------------------
GENO_DIR="${PROJ_ROOT}/data/geno"
SUMSTATS_PARENT="${PROJ_ROOT}/data/formatted_sumstats"
BIN_DIR="${TARDIS_BIN_DIR}"
REF_PARENT="$(dirname "$LD_REF_DIR")"

mkdir -p "$GENO_DIR" "$SUMSTATS_PARENT" "$BIN_DIR" "$REF_PARENT"
ensure_dirs

# =============================================================================
# 1. TOOL BINARIES
# =============================================================================
log_info "--- Checking tool binaries ---"

# GCTB
if [[ -x "$GCTB" ]]; then
    log_info "  GCTB OK: $GCTB"
    ((found_count++))
else
    log_error "  GCTB not found at $GCTB"
    log_info "  To install GCTB 2.5.2:"
    log_info "    cd ${BIN_DIR}"
    log_info "    wget https://cnsgenomics.com/software/gctb/download/gctb_2.5.2_Linux.zip"
    log_info "    unzip gctb_2.5.2_Linux.zip"
    log_info "    cp gctb_2.5.2_Linux/gctb ."
    log_info "    chmod +x gctb"
    scp_commands+=("# GCTB: download on Tardis (see above)")
    ((error_count++))
fi

# PLINK 1.9
if [[ -x "$PLINK19" ]]; then
    log_info "  PLINK 1.9 OK: $PLINK19"
    ((found_count++))
else
    # Check PATH
    if command -v plink &>/dev/null; then
        PLINK19_PATH=$(command -v plink)
        ln -sf "$PLINK19_PATH" "$PLINK19" 2>/dev/null || cp "$PLINK19_PATH" "$PLINK19"
        log_info "  PLINK 1.9 linked from PATH: $PLINK19_PATH"
        ((found_count++))
    else
        log_error "  PLINK 1.9 not found at $PLINK19"
        log_info "  To install PLINK 1.9:"
        log_info "    cd ${BIN_DIR}"
        log_info "    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip"
        log_info "    unzip plink_linux_x86_64_20231018.zip plink"
        log_info "    chmod +x plink"
        scp_commands+=("# PLINK 1.9: download on Tardis (see above)")
        ((error_count++))
    fi
fi

# PLINK 2
if [[ -x "$PLINK2" ]]; then
    log_info "  PLINK 2 OK: $PLINK2"
    ((found_count++))
else
    log_error "  PLINK 2 not found at $PLINK2"
    ((error_count++))
fi

# =============================================================================
# 2. LD REFERENCE MATRICES (UKB 50k shrunk sparse)
# =============================================================================
log_info ""
log_info "--- Checking LD reference matrices ---"

LD_TEST_FILE="${LD_REF_DIR}/ukb50k_shrunk_chr1_mafpt01.ldm.sparse"
if [[ -f "$LD_TEST_FILE" ]] || [[ -f "${LD_TEST_FILE}.bin" ]]; then
    log_info "  LD reference OK: $LD_REF_DIR"
    ((found_count++))
else
    log_error "  LD reference NOT found at $LD_REF_DIR"
    log_info ""
    log_info "  Either — download directly on the cluster (large, ~50GB):"
    log_info "    mkdir -p $LD_REF_DIR && cd $LD_REF_DIR"
    log_info "    wget https://zenodo.org/record/3350914/files/ukb50k_2.8M_shrunk_sparse.zip"
    log_info "    unzip ukb50k_2.8M_shrunk_sparse.zip"
    log_info ""
    log_info "  Or — copy from the RStudio server:"
    scp_commands+=("scp -r ${MPIB_SERVER}:${MPIB_LD_REF} ${REF_PARENT}/")
    log_info "    scp -r ${MPIB_SERVER}:${MPIB_LD_REF} ${REF_PARENT}/"
    ((error_count++))
fi

# =============================================================================
# 3. FORMATTED SUMMARY STATISTICS
# =============================================================================
log_info ""
log_info "--- Checking formatted summary statistics ---"

NEEDED_SUMSTATS_DIRS=(
    "EA4_excl_SHIP"
    "EA4_excl_BASEII"
    "Cog_Malanchini"
    "NCog_Malanchini"
    "Height_GIANT"
)

# Try to find a CIFS mount with the formatted sumstats
FMT_ROOT=""
for candidate in "${CIFS_FORMATTED[@]}"; do
    if [[ -d "$candidate" ]]; then
        for subdir in "${NEEDED_SUMSTATS_DIRS[@]}"; do
            if [[ -d "${candidate}/${subdir}" ]]; then
                FMT_ROOT="$candidate"
                log_info "  Found formatted sumstats via CIFS: $FMT_ROOT"
                break 2
            fi
        done
    fi
done

if [[ -n "$FMT_ROOT" ]]; then
    # Link from CIFS mount
    for subdir in "${NEEDED_SUMSTATS_DIRS[@]}"; do
        src="${FMT_ROOT}/${subdir}"
        dst="${SUMSTATS_PARENT}/${subdir}"
        if [[ -d "$src" ]]; then
            [[ -L "$dst" ]] && rm -f "$dst"
            ln -sf "$src" "$dst" 2>/dev/null || cp -r "$src" "$dst"
            nfiles=$(ls "$dst"/*SBayesR* 2>/dev/null | wc -l || echo "0")
            log_info "  OK: $subdir ($nfiles chr files)"
            ((found_count++))
        else
            log_error "  MISSING: $subdir not at CIFS source"
            ((error_count++))
        fi
    done
else
    # Check if already copied
    all_present=true
    for subdir in "${NEEDED_SUMSTATS_DIRS[@]}"; do
        dst="${SUMSTATS_PARENT}/${subdir}"
        if [[ -d "$dst" ]]; then
            nfiles=$(ls "$dst"/*SBayesR* 2>/dev/null | wc -l || echo "0")
            if [[ "$nfiles" -ge 22 ]]; then
                log_info "  OK: $subdir (already present, $nfiles files)"
                ((found_count++))
            else
                log_warn "  INCOMPLETE: $subdir ($nfiles/22 chr files)"
                all_present=false
                ((error_count++))
            fi
        else
            log_error "  MISSING: $subdir"
            all_present=false
            ((error_count++))
        fi
    done

    if [[ "$all_present" == false ]]; then
        _subdir_csv=$(IFS=,; echo "${NEEDED_SUMSTATS_DIRS[*]}")
        log_info ""
        log_info "  Stage the missing formatted sumstats."
        log_info "  CANONICAL SOURCE = the data share:"
        log_info "    Projects/04_data_analysis/013-PGS/data/formatted/<SUBDIR>/"
        log_info "  From any machine where the share is mounted, rsync straight to Tardis:"
        log_info "    rsync -av \"<SHARE>/Projects/04_data_analysis/013-PGS/data/formatted/\"{${_subdir_csv}} \\"
        log_info "      ${TARDIS_USER}@<tardis>:${SUMSTATS_PARENT}/"
        if [[ -n "$MPIB_FORMATTED" ]]; then
            log_info "  Or scp from the configured MPIB server (\$MPIB_FORMATTED):"
            for subdir in "${NEEDED_SUMSTATS_DIRS[@]}"; do
                dst="${SUMSTATS_PARENT}/${subdir}"
                if [[ ! -d "$dst" ]] || [[ $(ls "$dst"/*SBayesR* 2>/dev/null | wc -l || echo 0) -lt 22 ]]; then
                    cmd="scp -r ${MPIB_SERVER}:${MPIB_FORMATTED}/${subdir} ${SUMSTATS_PARENT}/"
                    log_info "    $cmd"
                    scp_commands+=("$cmd")
                fi
            done
        fi
    fi
fi

# =============================================================================
# 4. GENOTYPE DATA
# =============================================================================
log_info ""
log_info "--- Checking genotype data ---"

# Try to find SRT CIFS mount
SRT_ROOT=""
for candidate in "${CIFS_SRT_BASE_CANDIDATES[@]}"; do
    if [[ -d "$candidate" ]]; then
        SRT_ROOT="$candidate"
        log_info "  Found SRT CIFS mount: $SRT_ROOT"
        break
    fi
done

# Define genotype entries: LABEL|CONFIG_PREFIX|SRT_SUBPATH|MPIB_FULL_PATH
GENO_ENTRIES=(
    "BASE-II|${GENO_BASEII}|<base2-share>/private/data/003_BASEII_geno/LIFEBRAIN_BASEII_2015_AFFY|${MPIB_GENO_BASEII}"
    "SHIP-0|${GENO_SHIP0}|<ship-share>/private/data/008_SHIP_processed_geno/Final-SHIP-0_R4a.chrall.dose|${MPIB_GENO_SHIP0}"
    "SHIP-Td|${GENO_SHIPTD}|<ship-share>/private/data/008_SHIP_processed_geno/SHIP-Td-SHIP-Td_B2_merged.chrall|${MPIB_GENO_SHIPTD}"
    "SOEP-G|${GENO_SOEP}|<soep-share>/private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/bed_bim_fam/SOEP-G.b37.mildQC.hrc1-1_imp.chr1.dose.vcf|${MPIB_GENO_SOEP}"
    "TwinLife|${GENO_TWINLIFE}|<twinlife-share>/private/data/Gendata/001_processed/TwinLife.b37.chrall|${MPIB_GENO_TWINLIFE}"
)

for entry in "${GENO_ENTRIES[@]}"; do
    IFS='|' read -r label target srt_rel mpib_path <<< "$entry"

    # Already present?
    if [[ -f "${target}.bed" ]] && [[ -f "${target}.bim" ]] && [[ -f "${target}.fam" ]]; then
        log_info "  OK: $label (already present)"
        ((found_count++))
        continue
    fi

    linked=0
    # Try SRT CIFS mount
    if [[ -n "$SRT_ROOT" ]]; then
        src="${SRT_ROOT}/${srt_rel}"
        if [[ -f "${src}.bed" ]]; then
            create_symlink "$src" "$target" "$label (CIFS)"
            linked=1
        fi
    fi

    if [[ $linked -eq 0 ]]; then
        bname=$(basename "$target")
        log_error "  MISSING: $label"
        # Build scp commands for .bed/.bim/.fam
        for ext in .bed .bim .fam; do
            cmd="scp ${MPIB_SERVER}:${mpib_path}${ext} ${GENO_DIR}/"
            scp_commands+=("$cmd")
        done
        log_info "    scp ${MPIB_SERVER}:${mpib_path}.{bed,bim,fam} ${GENO_DIR}/"
        ((error_count++))
    fi
done

# =============================================================================
# SUMMARY
# =============================================================================
log_info ""
log_info "============================================================"
log_info "SETUP SUMMARY"
log_info "============================================================"
log_info "  Items OK:      $found_count"

if [[ $error_count -gt 0 ]]; then
    log_error "  Items MISSING: $error_count"
    log_info ""
    log_info "============================================================"
    log_info "TO FIX: Run these commands on Tardis"
    log_info "============================================================"
    log_info ""

    # Print all scp commands as a copy-paste block
    if [[ ${#scp_commands[@]} -gt 0 ]]; then
        log_info "# --- Copy data from the RStudio server ---"
        for cmd in "${scp_commands[@]}"; do
            log_info "$cmd"
        done
    fi

    log_info ""
    log_info "After fixing, re-run: bash scripts/setup_tardis_data.sh"
    log_info "Then verify:          bash run_pipeline.sh --step 00"
else
    log_info ""
    log_info "  All items found!"
    log_info "  Next: bash run_pipeline.sh --step 00"
fi
