#!/bin/bash
## one (tissue, chr) of stage S3.2: run_conditional_v3m.R with its own log (env RS, CODE, NC, LOGDIR; cwd = root)
T=$1; C=$2; L=$LOGDIR/cond/$T.$C.log
if "$RS" "$CODE/run_conditional_v3m.R" freeze "$T" real "$C" "$NC" > "$L" 2>&1; then
  echo "[S3.2] $T $C PASS: $(grep -E 'units, signals' "$L" | tail -n 1)"
else echo "[S3.2] $T $C FAIL (see $L)"; exit 1; fi
