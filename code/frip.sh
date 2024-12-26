#!/bin/bash
#FRiP score

#SBATCH --time 20:00:00
#SBATCH --job-name=frip  
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=frip.%J.out
#SBATCH --error=frip.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload python


#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>frip.log 2>&1


python frip.py

