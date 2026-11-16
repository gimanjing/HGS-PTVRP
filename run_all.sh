#!/usr/bin/env bash
# Runs HGS-PTVRP on every instance of the three sets and collects one line per run in Results/results.csv.
#
#   ./run_all.sh [TIME] [SEEDS] [JOBS]
#
#   TIME   time limit per run, in seconds of CPU time   (default 60)
#   SEEDS  number of seeds per instance, 1..SEEDS        (default 1)
#   JOBS   runs executed in parallel, one core each      (default 1)
#
# Optional environment variables:
#   FILTER   only the instances whose file name matches this pattern, e.g. FILTER='X-n1*'
#   HGS_BIN  path of the solver                          (default build/hgs)
#   OUT      output folder                               (default Results)
#
# A run that is already in results.csv is skipped, so an interrupted batch can simply be started again.

set -euo pipefail

TIME="${1:-60}"
SEEDS="${2:-1}"
JOBS="${3:-1}"

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="${HGS_BIN:-$ROOT/build/hgs}"
OUT="${OUT:-$ROOT/Results}"
FILTER="${FILTER:-*}"
CSV="$OUT/results.csv"

if [ ! -x "$BIN" ]; then
	echo "Solver not found at $BIN. Build it first:  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build" >&2
	exit 1
fi

mkdir -p "$OUT"
export BIN OUT CSV TIME ROOT

run_one() {
	local set="$1" name="$2" seed="$3"
	mkdir -p "$OUT/$set"
	local base="$OUT/$set/${name}_seed${seed}"
	"$BIN" "$ROOT/Instances/$set/$name.vrp" "$base.sol" -seed "$seed" -t "$TIME" -csv "$CSV" > "$base.log" 2>&1
	if grep -q "^EXCEPTION" "$base.log"; then
		echo "FAILED  $set / $name / seed $seed : $(grep -m1 '^EXCEPTION' "$base.log")"
	else
		echo "done    $set / $name / seed $seed"
	fi
}
export -f run_one

already_done() {
	[ -f "$CSV" ] && awk -F, -v s="$1" -v i="$2" -v d="$3" '$1==s && $2==i && $6==d {found=1} END {exit !found}' "$CSV"
}

for set in "Instances 1" "Instances 2" "Instances 3"; do
	for file in "$ROOT/Instances/$set"/$FILTER.vrp; do
		[ -e "$file" ] || continue
		name="$(basename "$file" .vrp)"
		for seed in $(seq 1 "$SEEDS"); do
			if ! already_done "$set" "$name" "$seed"; then
				printf '%s\0%s\0%s\0' "$set" "$name" "$seed"
			fi
		done
	done
done | xargs -0 -r -n 3 -P "$JOBS" bash -c 'run_one "$@"' _

echo "Results in $CSV"
