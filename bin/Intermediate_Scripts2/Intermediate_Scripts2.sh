#!/bin/bash
##SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=IS2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mnaidu@gradcenter.cuny.edu#  Created by Praveena  on 2/25/24.
#

source ~/miniconda3/etc/profile.d/conda.sh

conda activate VenomFlowAnalysis

#There are 6 intermediate scripts
current_dir=$(pwd)

echo "GENOME_ID is: $GENOME_ID"
echo "SPECIES is: $SPECIES"
cd Intermediate_Scripts2
Rscript Figure_generation_SignalP.R
Rscript Figure_generation_Trinity.R
Rscript Generating_TopTables_Trinity.R $GENOME_ID $SPECIES
Rscript Figure_generation_Transdecoder.R
Rscript Generating_Tables_Transdecoder_SignalP.R $GENOME_ID $SPECIES


