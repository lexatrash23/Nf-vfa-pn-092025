#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mnaidu@gradcenter.cuny.edu
#  Created by Praveena Naidu on 2/2/25.
#  
current_dir=$(pwd)
Rscript -e "rmarkdown::render('$current_dir/Rmarkdown_scripts/R.Rmd')"
