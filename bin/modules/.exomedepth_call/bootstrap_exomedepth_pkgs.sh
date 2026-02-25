#!/usr/bin/env bash
set -euo pipefail

msg(){ echo "[bootstrap_exomedepth] $*"; }
die(){ msg "ERROR: $*"; exit 1; }

R_MODULES="${R_MODULES:-"gcc/11.4.0 openmpi/4.1.4 R/4.4.1"}"
R_AUTOINSTALL="${R_AUTOINSTALL:-1}"

maybe_enable_modules() {
  if command -v module >/dev/null 2>&1; then return 0; fi
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh && command -v module >/dev/null 2>&1 && return 0
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash && command -v module >/dev/null 2>&1 && return 0
  return 1
}

pick_r_user_lib_root() {
  if [[ -n "${MAPI_R_LIBS_USER_ROOT:-}" ]]; then echo "$MAPI_R_LIBS_USER_ROOT"; return 0; fi
  if [[ -d "/standard" && -n "${HPC_ALLOCATION:-}" ]]; then
    local cand="/standard/${HPC_ALLOCATION}/${USER}/R/goolf"
    if mkdir -p "$cand" >/dev/null 2>&1; then echo "$cand"; return 0; fi
  fi
  echo "$HOME/R/goolf"
}

msg "Enabling modules + loading: $R_MODULES"
maybe_enable_modules || die "module system not available"
module purge >/dev/null 2>&1 || true
# shellcheck disable=SC2086
module load $R_MODULES || die "failed to load: $R_MODULES"

command -v Rscript >/dev/null 2>&1 || die "Rscript not found after module load"
r_mm="$(Rscript --vanilla -e 'cat(paste0(R.version$major,".",strsplit(R.version$minor,"\\.")[[1]][1]))')"
root="$(pick_r_user_lib_root)"
export R_LIBS_USER="${root}/${r_mm}"
mkdir -p "$R_LIBS_USER" || die "Failed to create R_LIBS_USER: $R_LIBS_USER"
msg "R_LIBS_USER=$R_LIBS_USER"
msg "Rscript=$(command -v Rscript)"

# Install + verify
Rscript --vanilla - <<'RS'
repos <- c(CRAN="https://cloud.r-project.org")

# CRAN bootstrap
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos=repos)
}

# Install Bioconductor deps for ExomeDepth
pkgs <- c("ExomeDepth","GenomicRanges","rtracklayer","Rsamtools")
BiocManager::install(pkgs, ask=FALSE, update=FALSE)

miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  cat("MISSING:", paste(miss, collapse=" "), "\n")
  print(.libPaths())
  quit(status=2)
} else {
  cat("OK: installed + loadable\n")
  print(.libPaths())
}
RS

msg "done"
