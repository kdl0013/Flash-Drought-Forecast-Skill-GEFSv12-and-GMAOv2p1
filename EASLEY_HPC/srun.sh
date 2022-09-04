
srun -N1 -n1 --time=24:00:00 --pty bash

conda init bash #for easley to start new shell
source ~/.bashrc #to open a new session

#required modules
module load cdo
module load ncl
module load python/anaconda/3.8.6

conda activate cfgrib #to activate the environment for GEFS pre-processing










