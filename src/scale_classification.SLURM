#!/bin/bash
#SBATCH --job-name=tif_to_df    # Job name
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --cpus-per-task=1            # Number of cores per MPI rank
#SBATCH --mem=14gb

#SBATCH --time=24:00:00       #Time limit hrs:min:sec
#SBATCH --output=scale.out   # Standard output and error log
#SBATCH --error=scale.err


pwd; hostname; date
ml conda
source activate quetzal3

python3 /blue/ewhite/s.marconi/NeonSpeciesClassification/src/prediction_at_tile.py $1 $2 $3 $4 $5
echo "done"

date
