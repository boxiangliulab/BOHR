#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR="/home/users/nus/e1327989/scratch/gtex_v8/split_junctions"
SJ_DIR="/home/users/nus/e1327989/scratch/gtex_v8/star_junctions"
NUM_THREADS=48
mkdir -p "$SJ_DIR"

convert_file() {
  local INPUT_FILE="$1"
  local SAMPLE OUTPUT_FILE TMP LOG

  SAMPLE="$(basename "$INPUT_FILE" .junction.gct.gz)"
  OUTPUT_FILE="$SJ_DIR/${SAMPLE}.SJ.out.tab"
  TMP="${OUTPUT_FILE}.tmp"
  LOG="$SJ_DIR/${SAMPLE}.log"

  # ✅ Skip automatically: the result file already exists and is non-empty
  if [[ -s "$OUTPUT_FILE" ]]; then
    echo "[SKIP] $SAMPLE already done: $OUTPUT_FILE"
    return 0
  fi

  # Remove any leftover temporary file (so that an earlier interruption does not affect the check)
  rm -f "$TMP"

  # Use `if ! pipeline` to swallow a single-file error and keep running the other files
  if ! gzip -cd -- "$INPUT_FILE" 2> "$LOG" | awk 'BEGIN{OFS="\t"}
    NR<=3 {next}
    $1=="Name" {next}
    NF<2 {next}
    {
      id=$1
      n=split(id,a,"_")
      if(n<3) next

      chr=a[1]
      start=a[2]
      end=a[3]

      if(start !~ /^[0-9]+$/) next
      if(end   !~ /^[0-9]+$/) next

      strand=0
      motif=0
      annotated=0
      unique=$2
      multi=0
      maxOverhang=0

      if(unique !~ /^-?[0-9]+$/) unique=0

      print chr,start,end,strand,motif,annotated,unique,multi,maxOverhang
    }' > "$TMP"
  then
    echo "[WARN] failed for $SAMPLE (input=$INPUT_FILE). See log: $LOG" >&2
    rm -f "$TMP"
    : > "$OUTPUT_FILE"   # Keep an empty file as a placeholder for the failure (the original logic)
    return 0
  fi

  if [[ ! -s "$TMP" ]]; then
    echo "[WARN] empty output for $SAMPLE (input=$INPUT_FILE). See log: $LOG" >&2
    : > "$OUTPUT_FILE"
    rm -f "$TMP"
    return 0
  fi

  mv -f "$TMP" "$OUTPUT_FILE"
}

export -f convert_file
export OUTPUT_DIR SJ_DIR

find "$OUTPUT_DIR" -type f -name "*.junction.gct.gz" -print0 \
  | xargs -0 -P "$NUM_THREADS" -I{} bash -c 'convert_file "$@"' _ {}

echo "All files converted."