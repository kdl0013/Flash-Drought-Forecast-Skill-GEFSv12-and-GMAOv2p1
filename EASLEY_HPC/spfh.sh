#!/bin/bash
#SBATCH --job-name=spfh             # job name
#SBATCH --nodes=1                    #number of nodes
#SBATCH --ntasks=20
#SBATCH --partition=davisj8_std          # name of partition to submit job
#SBATCH --time=100:00:00              # Run time (D-HH:MM:SS)
#SBATCH --output=spfh.out          # Output file. %j is replaced with job ID
#SBATCH --error=spfh.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=kdl0013@auburn.edu
#SBATCH --mem=2gb

conda init bash #for easley to start new shell
source ~/.bashrc #to open a new session

#required modules
module load cdo
module load ncl
module load python/anaconda/3.8.6

conda activate cfgrib #to activate the environment for GEFS pre-processing

#export path for cfgrib
newpath='~/.conda/envs/cfgrib/bin/cfgrib'
export $PATH:$newpath

cd ~/process_GEFS
#create bash wget scripts
python3 ~/process_GEFS/spfh_test_compression.py









