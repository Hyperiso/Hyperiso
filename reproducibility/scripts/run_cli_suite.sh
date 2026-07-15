#!/usr/bin/env bash
set -euo pipefail

UPDATE_EXPECTED=0
if [[ "${1:-}" == "--update-expected" ]]; then
  UPDATE_EXPECTED=1
  shift
fi
if [[ "$#" -ne 0 ]]; then
  echo "Usage: $0 [--update-expected]" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPRO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${REPRO_DIR}/.." && pwd)"
OUT_DIR="${REPRO_DIR}/outputs"
EXP_DIR="${REPRO_DIR}/expected_outputs"
INPUT="${REPRO_DIR}/inputs/sm_reference.flha"
mkdir -p "${OUT_DIR}" "${EXP_DIR}"

find_hyperiso_bin() {
  if [[ -n "${HYPERISO_BIN:-}" ]]; then
    printf '%s\n' "${HYPERISO_BIN}"
    return
  fi
  local candidates=(
    "${ROOT_DIR}/Hyperiso/Hyperiso/core/Test/build/UserInterfaceLib/hyperiso-ui"
    "${ROOT_DIR}/Hyperiso/Hyperiso/core/Test/build/hyperiso-ui"
    "${ROOT_DIR}/build/hyperiso-ui"
    "${ROOT_DIR}/build/UserInterfaceLib/hyperiso-ui"
    "${ROOT_DIR}/build/Hyperiso/Hyperiso/core/src/UserInterface/hyperiso-ui"
    "${ROOT_DIR}/build/Hyperiso/Hyperiso/core/src/UserInterfaceLib/hyperiso-ui"
  )
  for candidate in "${candidates[@]}"; do
    if [[ -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return
    fi
  done
  if command -v hyperiso-ui >/dev/null 2>&1; then
    command -v hyperiso-ui
    return
  fi
  echo "Cannot find hyperiso-ui. Set HYPERISO_BIN=/path/to/hyperiso-ui." >&2
  return 1
}

BIN="$(find_hyperiso_bin)"
if [[ ! -x "${BIN}" ]]; then
  echo "hyperiso-ui is not executable: ${BIN}" >&2
  exit 1
fi

echo "Using hyperiso-ui: $(realpath "${BIN}")"
EXPECTED_VERSION="$(python3 - "${REPRO_DIR}/manifest.json" <<'PY_VERSION'
import json
import sys
from pathlib import Path

print(json.loads(Path(sys.argv[1]).read_text())["hyperiso_version"])
PY_VERSION
)"
CLI_VERSION="$("${BIN}" --version)"
if [[ "${CLI_VERSION}" != "${EXPECTED_VERSION}" ]]; then
  echo "Unexpected hyperiso-ui version: ${CLI_VERSION} (expected ${EXPECTED_VERSION})" >&2
  exit 1
fi

run_case() {
  local id="$1"
  local outfile="$2"
  shift 2
  local raw_stdout raw_stderr
  raw_stdout="$(mktemp)"
  raw_stderr="$(mktemp)"
  echo "[${id}] $*"
  if ! "$@" >"${raw_stdout}" 2>"${raw_stderr}"; then
    cat "${raw_stdout}" >&2
    cat "${raw_stderr}" >&2
    rm -f "${raw_stdout}" "${raw_stderr}"
    return 1
  fi
  if [[ -s "${raw_stderr}" ]]; then
    echo "[${id}] unexpected stderr from reference command:" >&2
    cat "${raw_stderr}" >&2
    rm -f "${raw_stdout}" "${raw_stderr}"
    return 1
  fi
  python3 "${SCRIPT_DIR}/normalize_cli_output.py" "${raw_stdout}" "${OUT_DIR}/${outfile}"
  rm -f "${raw_stdout}" "${raw_stderr}"
  cat "${OUT_DIR}/${outfile}"
}

run_case R1 wilson_b_sm_nnlo.txt \
  "${BIN}" wilson summary \
  --model SM --lha "${INPUT}" --groups BCoefficients \
  --coeffs C7,C8,C9,C10 --qmatch 81 --q 4.8 --order NNLO

run_case R2 wilson_bscalar_sm_lo.txt \
  "${BIN}" wilson summary \
  --model SM --lha "${INPUT}" --groups BScalarCoefficients \
  --coeffs CQ1_MU,CQ2_MU --qmatch 81 --q 4.8 --order LO

run_case R3 observables_sm_nnlo.txt \
  "${BIN}" observable summary \
  --model SM --lha "${INPUT}" \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma --order NNLO

run_case R4 observable_binned_fl_kstar_sm_nnlo.txt \
  "${BIN}" observable summary \
  --model SM --lha "${INPUT}" --observables BR_Bs__mu_mu \
  --bins 'F_L_B0__K*0_mu_mu:1.1:6.0' --order NNLO

run_case R5 statistics_sm_seed123456.txt \
  "${BIN}" statistic summary \
  --model SM --lha "${INPUT}" \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma \
  --uncertainties --draws 200 --seed 123456 \
  --samples-csv "${OUT_DIR}/statistics_samples.csv" --order NNLO

python3 "${SCRIPT_DIR}/check_expected_outputs.py" \
  --manifest "${REPRO_DIR}/manifest.json" \
  --outputs "${OUT_DIR}" \
  --validate-only

if [[ "${UPDATE_EXPECTED}" -eq 1 ]]; then
  python3 "${SCRIPT_DIR}/freeze_reference_outputs.py" --root "${ROOT_DIR}"
  echo "References updated. Review the numerical diff before committing."
else
  python3 "${SCRIPT_DIR}/check_expected_outputs.py" \
    --manifest "${REPRO_DIR}/manifest.json" \
    --outputs "${OUT_DIR}" \
    --expected "${EXP_DIR}"
fi
