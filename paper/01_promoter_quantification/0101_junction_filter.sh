#!/bin/bash

# Configuration parameters
SJ_DIR="star_junctions"          # Directory containing the SJ.out.tab files
FILTERED_DIR="filtered_junctions" # Output directory for the filtered files
THRESHOLD=5                       # Filtering threshold
NUM_THREADS=48                    # Number of parallel threads

# Create the output directory
mkdir -p $FILTERED_DIR

# Filter function: remove records with count < 5
filter_file() {
    local INPUT_FILE=$1
    local SAMPLE=$(basename "$INPUT_FILE")
    local OUTPUT_FILE="$FILTERED_DIR/$SAMPLE"
    
    echo "Filtering file: $INPUT_FILE -> $OUTPUT_FILE"
    
    # Filter records: those with uniqueReads < 5 are removed (column 7)
    awk -v threshold=$THRESHOLD '$7 >= threshold {print $0}' $INPUT_FILE > $OUTPUT_FILE
    
    # Count the records before and after filtering
    BEFORE=$(wc -l < $INPUT_FILE)
    AFTER=$(wc -l < $OUTPUT_FILE)
    REMOVED=$((BEFORE - AFTER))
    
    echo "Sample $SAMPLE: $BEFORE records before filtering, $AFTER after filtering, $REMOVED removed"
}

export -f filter_file
export FILTERED_DIR THRESHOLD

# Find all SJ.out.tab files and process them in parallel
find "$SJ_DIR" -name "*.SJ.out.tab" | xargs -P $NUM_THREADS -I{} bash -c 'filter_file "{}"'

# Summary statistics
TOTAL_FILES=$(find "$SJ_DIR" -name "*.SJ.out.tab" | wc -l)
echo "Done: $TOTAL_FILES files filtered in total"

# Optional: average filtering ratio over all files
if [ $TOTAL_FILES -gt 0 ]; then
    echo "Computing the average filtering rate..."
    TOTAL_BEFORE=0
    TOTAL_AFTER=0
    
    for file in $SJ_DIR/*.SJ.out.tab; do
        base_file=$(basename "$file")
        before=$(wc -l < "$file")
        after=$(wc -l < "$FILTERED_DIR/$base_file")
        
        TOTAL_BEFORE=$((TOTAL_BEFORE + before))
        TOTAL_AFTER=$((TOTAL_AFTER + after))
    done
    
    REMOVED=$((TOTAL_BEFORE - TOTAL_AFTER))
    PERCENT=$(echo "scale=2; $REMOVED * 100 / $TOTAL_BEFORE" | bc)
    
    echo "Total: $TOTAL_BEFORE records before filtering, $TOTAL_AFTER after filtering"
    echo "Removed $REMOVED records, $PERCENT% of the total"
fi