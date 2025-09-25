#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=Rmarkdown
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mnaidu@gradcenter.cuny.edu
#  Created by Praveena Naidu on 2/2/25.
#  
current_dir=$(pwd)
source ~/miniconda3/etc/profile.d/conda.sh

conda activate r_env


#Menupages
bash "$current_dir/Rmarkdown_scripts/A.sh"

bash "$current_dir/Rmarkdown_scripts/B.sh"

#bash "$current_dir/Rmarkdown_scripts/H.sh"

#bash "$current_dir/Rmarkdown_scripts/J.sh"

bash "$current_dir/Rmarkdown_scripts/L.sh"

bash "$current_dir/Rmarkdown_scripts/M.sh"

bash "$current_dir/Rmarkdown_scripts/N.sh"

bash "$current_dir/Rmarkdown_scripts/O.sh"

#bash "$current_dir/Rmarkdown_scripts/P.sh"

bash "$current_dir/Rmarkdown_scripts/Q.sh"

bash "$current_dir/Rmarkdown_scripts/R.sh"

bash "$current_dir/Rmarkdown_scripts/S.sh"

bash "$current_dir/Rmarkdown_scripts/V.sh"

bash "$current_dir/Rmarkdown_scripts/W.sh"

bash "$current_dir/Rmarkdown_scripts/X.sh"

