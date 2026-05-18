#!/bin/bash

# Make sure the environment is properly set
source /lustre/r3b/mjena/G249/hpc_stuff/config.sh
source ./check_environment.sh

# Create log directories if they don't exist
mkdir -p /lustre/r3b/mjena/G249/rootfiles/logs/output
mkdir -p /lustre/r3b/mjena/G249/rootfiles/logs/errors

# Count all root files starting with main*
FILES=($(find "$OUTROOTPATH/../5cmshift" -type f -name 'main*.root' -printf "%f\n" | sort))
FILE_COUNT=${#FILES[@]}

if [ "$FILE_COUNT" -eq 0 ]; then
    log_error "No root files found in $OUTROOTPATH"
    exit 1
else 
    log_success "Found $FILE_COUNT root files in $OUTROOTPATH"
fi

# Save the sorted list of files to a temporary file
FILE_LIST="$PWD/file_list.txt"
rm -f "$FILE_LIST"
for file in "${FILES[@]}"; do
    echo "$file" >> "$FILE_LIST"
done

log_success "File list created at $FILE_LIST"

# Check how many files already successfully processed
OUTPUT_DIR="/lustre/r3b/mjena/G249/rootfiles/analyzed/ts_4mus/with_randomization"
ALREADY_DONE=0
if [ -d "$OUTPUT_DIR" ]; then
    for file in "${FILES[@]}"; do
        output_file="$OUTPUT_DIR/${file/.root/_analysed.root}"
        if [ -f "$output_file" ]; then
            # Check if it has content (not just empty)
            size=$(stat -c%s "$output_file" 2>/dev/null || echo 0)
            if [ $size -gt 1000 ]; then
                ((ALREADY_DONE++))
            fi
        fi
    done
fi

log_message "Already processed successfully: $ALREADY_DONE / $FILE_COUNT"
log_message "Will process: $((FILE_COUNT - ALREADY_DONE)) files (including re-processing failed ones)"

# IMPORTANT: Limit concurrent jobs to avoid overwhelming the filesystem
MAX_CONCURRENT=100

log_message "Submitting array job for $FILE_COUNT files (max $MAX_CONCURRENT concurrent)..."

# Submit the SLURM array job with throttling
JOB_ID=$(sbatch --export=ALL,FILE_LIST="$FILE_LIST" \
    --array=0-$((FILE_COUNT - 1))%${MAX_CONCURRENT} \
    --parsable \
    analyse_files.sh)

if [ $? -eq 0 ]; then
    log_success "Job submitted successfully! Job ID: $JOB_ID"
    log_message "Monitor with: squeue -j $JOB_ID"
    log_message "Or: watch -n 5 'squeue -j $JOB_ID | tail -20'"
    log_message "Logs in: /lustre/r3b/mjena/G249/rootfiles/logs/"
else
    log_error "Failed to submit job!"
    exit 1
fi
