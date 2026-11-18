#!/bin/bash
#
# submit.sh -- submit the matrix-element calculation to LSF as a job array.
#
# One array element per diagram template: array index N computes template N
# (1-based, read from LSB_JOBINDEX in src/main.cpp), sweeping all 25 xi values
# with 8 threads, and appending its results to $LOG_DIR.
#
# USAGE
#   ./submit.sh                 all 82 matrix elements
#   ./submit.sh 17              only matrix element 17
#   ./submit.sh 17,23,40        only those three (resubmitting stragglers)
#   ./submit.sh 40-82           a range
#
#   EASY_MODE=1 ./submit.sh     quick pass: 1000x less statistics, for testing
#
# AFTERWARDS
#   bjobs -A                    progress of the array as a whole
#   bjobs -a -J numerics        per-element state; look for EXIT
#   bpeek <jobid>[<index>]      stdout of a running element
#   ls $LOG_DIR | wc -l         how many result files exist so far
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BINARY="$SCRIPT_DIR/build/Viktor_MC_integral"
LOG_DIR="$HOME/logs"      # MUST match log_dir in src/main.cpp
LSF_DIR="$HOME/lsf"       # batch system's own stdout/stderr, created below

N_TASKS=130               # number of entries in vector_diagram_template
N_THREADS=8               # MUST match N_threads in src/main.cpp
WALLTIME="23:30"          # hh:mm; a killed job keeps the xi values it finished
MEM_MB=16                 # per slot
QUEUE="physics-short"          # empty would use your default queue
JOB_NAME="numerics"

EASY_MODE="${EASY_MODE:-0}"

# ---------------------------------------------------------------------------
# Preflight -- every one of these failures is silent if you skip the check
# ---------------------------------------------------------------------------
if [[ ! -x "$BINARY" ]]; then
    echo "No executable at $BINARY" >&2
    echo "Build it first:" >&2
    echo "  cmake -S '$SCRIPT_DIR' -B '$SCRIPT_DIR/build' && cmake --build '$SCRIPT_DIR/build' -j" >&2
    exit 1
fi

# src/main.cpp aborts with exit 1 if this is missing, but that costs you a
# scheduling round trip per element to find out.
if [[ ! -d "$LOG_DIR" ]]; then
    echo "Results directory $LOG_DIR does not exist." >&2
    echo "Create it, or change log_dir in src/main.cpp to match." >&2
    exit 1
fi

# LSF does not create the -o/-e directory; it silently mails the output instead.
mkdir -p "$LSF_DIR"

INDICES="${1:-1-$N_TASKS}"
if [[ ! "$INDICES" =~ ^[0-9]+([,-][0-9]+)*$ ]]; then
    echo "Bad index specification '$INDICES' -- expected e.g. 1-82, 17, or 17,23,40" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
bsub_args=(
    -J "${JOB_NAME}[${INDICES}]"
    -n "$N_THREADS"
    -R "span[hosts=1]"          # all 8 slots on one machine: the threads share memory
    -R "rusage[mem=${MEM_MB}]"
    -W "$WALLTIME"
    -o "$LSF_DIR/%J_%I.out"     # %J = job id, %I = array index
    -e "$LSF_DIR/%J_%I.err"
)
[[ -n "$QUEUE" ]] && bsub_args+=(-q "$QUEUE")

echo "binary    : $BINARY"
echo "elements  : $INDICES  ($N_THREADS threads each, span[hosts=1])"
echo "walltime  : $WALLTIME"
echo "results   : $LOG_DIR"
echo "batch logs: $LSF_DIR"
echo "EASY_MODE : $EASY_MODE$([[ "$EASY_MODE" != 0 ]] && echo '   <-- reduced statistics, not physics results')"
echo

# bsub snapshots the current environment into the job, so exporting is enough.
# If your site disables that, add:  -env "all, EASY_MODE=$EASY_MODE"
export EASY_MODE
bsub "${bsub_args[@]}" "$BINARY"
