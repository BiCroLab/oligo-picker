#!/bin/bash -l
#SBATCH -A b2011210
#SBATCH -p devel -n 1
#SBATCH -t 00:01:00
#SBATCH -J compress
#SBATCH --open-mode=append
#SBATCH -o compressDir.%j.log -e compressDir.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mihaela.martis@bils.se



echo "yes, it worked"
