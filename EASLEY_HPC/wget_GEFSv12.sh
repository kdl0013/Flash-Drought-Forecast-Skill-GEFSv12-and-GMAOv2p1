#!/bin/bash
#SBATCH --job-name=gefs_preprocess             # job name
#SBATCH --nodes=1                    #number of nodes
#SBATCH --ntasks=20
#SBATCH --partition=davisj8_std          # name of partition to submit job
#SBATCH --time=100:00:00              # Run time (D-HH:MM:SS)
#SBATCH --output=gefs.out          # Output file. %j is replaced with job ID
#SBATCH --error=gefs.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=kdl0013@auburn.edu
#SBATCH --mem=20gb

conda init bash #for easley to start new shell
source ~/.bashrc #to open a new session
conda activate cfgrib #to activate the environment for GEFS pre-processing

#required modules
module load cdo
module load ncl
module load python/anaconda/3.8.6

#create bash wget scripts
python3 wget_create_GEFS.py

#Download wget scripts
for variable in wget_GEFS*.sh;do
bash $variable;done








