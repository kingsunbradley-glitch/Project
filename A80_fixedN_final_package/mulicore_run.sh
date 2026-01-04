#!/usr/bin/env bash
set -euo pipefail

# 用法：
#   ./mulicore_run.sh NPROCS INDIR RUN_FIRST RUN_LAST [OUTDIR]
#
# 例子：
#   ./mulicore_run.sh 8 ../ 1 385
#   ./mulicore_run.sh 8 ../ 1 385 out
#
# 环境变量：
#   CLEAN=1  运行前清理旧输出（默认 1）
#   CLEAN=0  不清理

NPROCS="${1:-8}"
INDIR="${2:-$(pwd)}"
RUN_FIRST="${3:-1}"
RUN_LAST="${4:-385}"
OUTDIR="${5:-out}"

# 脚本/源码所在目录（不要用 pwd，避免从别的目录调用时路径错）
BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 规范化 INDIR（避免相对路径 + 子进程 cd 导致找不到数据）
INDIR="$(cd "$INDIR" && pwd)"

# 规范化 OUTDIR（默认放在脚本目录下 out/）
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

JOBS_DIR="${OUTDIR}/jobs"
PDF_DIR="${OUTDIR}/pdf"
BIN="${OUTDIR}/process_runs_A80"

if [[ "${CLEAN:-1}" == "1" ]]; then
  rm -rf "$JOBS_DIR" "$PDF_DIR" \
         "${OUTDIR}/summary.root" "${OUTDIR}/qa_plots.root" \
         "${OUTDIR}/hadd_summary.log" "${OUTDIR}/hadd_qa.log" \
         "${OUTDIR}"/summary_*.root "${OUTDIR}"/qa_plots_*.root
fi

mkdir -p "$JOBS_DIR" "$PDF_DIR"

echo "[1/3] Build binary..."
# 小工程：直接编译源码目录下的 .cpp，输出到 out/
g++ -O2 -std=c++17 \
  "${BASEDIR}/main.cpp" \
  "${BASEDIR}/RunProcessor.cpp" \
  "${BASEDIR}/PeakFinder.cpp" \
  "${BASEDIR}/Metrics.cpp" \
  "${BASEDIR}/QaIO.cpp" \
  $(root-config --cflags --libs) \
  -o "$BIN"

TOTAL=$((RUN_LAST - RUN_FIRST + 1))
CHUNK=$(( (TOTAL + NPROCS - 1) / NPROCS ))

echo "[2/3] Launch jobs: NPROCS=$NPROCS runs=[$RUN_FIRST,$RUN_LAST] chunk=$CHUNK"
echo "  INDIR : $INDIR"
echo "  OUTDIR: $OUTDIR"

pids=()
out_sums=()
out_qas=()

for ((i=0; i<NPROCS; i++)); do
  s=$((RUN_FIRST + i*CHUNK))
  e=$((s + CHUNK - 1))
  if (( s > RUN_LAST )); then break; fi
  if (( e > RUN_LAST )); then e=$RUN_LAST; fi

  s5=$(printf "%05d" "$s")
  e5=$(printf "%05d" "$e")

  jobdir="${JOBS_DIR}/job_${i}_${s5}_${e5}"
  mkdir -p "$jobdir"

  outsum="${OUTDIR}/summary_${s5}_${e5}.root"
  outqa="${OUTDIR}/qa_plots_${s5}_${e5}.root"
  out_sums+=("$outsum")
  out_qas+=("$outqa")

  (
    cd "$jobdir"
    echo "JOB $i: runs [$s,$e]" > job.info
    "$BIN" "$INDIR" "$s" "$e" \
      "$outsum" "$outqa" "$PDF_DIR" \
      > run.log 2>&1
    echo "DONE" > DONE
  ) &

  pids+=($!)
  echo "  started job $i runs [$s,$e] pid=${pids[-1]}"
done

fail=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then fail=1; fi
done
if (( fail != 0 )); then
  echo "ERROR: some jobs failed. Check ${JOBS_DIR}/job_*/run.log"
  exit 1
fi

echo "[3/3] hadd merge..."
hadd -f "${OUTDIR}/summary.root" "${out_sums[@]}" > "${OUTDIR}/hadd_summary.log" 2>&1
hadd -f "${OUTDIR}/qa_plots.root" "${out_qas[@]}" > "${OUTDIR}/hadd_qa.log" 2>&1

echo "Cleanup chunk outputs..."
rm -f "${out_sums[@]}" "${out_qas[@]}"

echo "DONE."
echo "Outputs in: ${OUTDIR}"
echo "  summary.root"
echo "  qa_plots.root"
echo "  pdf/ (all per-run PDFs)"
echo "Logs:"
echo "  jobs/job_*/run.log"
echo "  hadd_summary.log  hadd_qa.log"
