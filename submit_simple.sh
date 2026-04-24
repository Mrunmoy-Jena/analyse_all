#!/bin/bash
# submit_simple.sh

source /u/mjena/config.sh
mkdir -p /u/mjena/rootfiles/logs/output
mkdir -p /u/mjena/rootfiles/logs/errors

for file in $OUTROOTPATH/main*.root; do
    basename=$(basename "$file" .root)
    
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=analyze_${basename}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --partition=long
#SBATCH --time=8:00:00
#SBATCH --output=/u/mjena/rootfiles/logs/output/${basename}.out
#SBATCH --error=/u/mjena/rootfiles/logs/errors/${basename}.err

source /u/mjena/config.sh
/u/mjena/analyse_all/analysis "$file"
EOF
    
    echo "Submitted job for $basename"
done
