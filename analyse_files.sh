#!/bin/bash
#SBATCH --job-name=analyse_WR
#SBATCH --time=2:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --partition=highmem
#SBATCH --output=/lustre/r3b/mjena/G249/rootfiles/logs/output/temp_%A_%a.out
#SBATCH --error=/lustre/r3b/mjena/G249/rootfiles/logs/errors/temp_%A_%a.err

# Load environment
source /lustre/r3b/mjena/G249/hpc_stuff/config.sh
source ./check_environment.sh

# Get the current file from the list
FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FILE_LIST")

if [ -z "$FILE" ]; then
    echo "ERROR: Could not read file for task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract basename for log naming
FILEBASE=$(basename "$FILE" .root)

# Define final log paths
FINAL_OUT="/lustre/r3b/mjena/G249/rootfiles/logs/output/${FILEBASE}.out"
FINAL_ERR="/lustre/r3b/mjena/G249/rootfiles/logs/errors/${FILEBASE}.err"

# Redirect output to final locations
exec 1>"$FINAL_OUT"
exec 2>"$FINAL_ERR"

echo "=========================================="
echo "Processing: $FILE"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Started: $(date)"
echo "Node: $SLURM_NODELIST"
echo "=========================================="

# Check if analysis executable exists
if [ ! -f "/lustre/r3b/mjena/G249/analyse_all/analysis" ]; then
    echo "ERROR: Analysis executable not found at /lustre/r3b/mjena/G249/analyse_all/analysis"
    exit 1
fi

# Run analysis
echo "Running: /lustre/r3b/mjena/G249/analyse_all/analysis $FILE"
/lustre/r3b/mjena/G249/analyse_all/analysis "$FILE"

EXIT_CODE=$?
echo "Analysis finished with exit code: $EXIT_CODE at $(date)"

# Clean up temp files (they'll be empty since we redirected immediately)
rm -f "/lustre/r3b/mjena/G249/rootfiles/logs/output/temp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
rm -f "/lustre/r3b/mjena/G249/rootfiles/logs/errors/temp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

exit $EXIT_CODE
