#!/bin/bash

#SBATCH --time 20:00:00
#SBATCH --job-name=reg_diff_exp_peaks 
#SBATCH -p knl_usr_prod
#SBATCH -N 1
#SBATCH -n 20 
#SBATCH --mem=86000
#SBATCH --output=reg_diff_exp_peaks.%J.out
#SBATCH --error=reg_diff_exp_peaks.%J.err
#SBATCH --account=<ACCOUNT_NAME>


module load profile/global
module load profile/bioinf

module load autoload macs

 
MACS(){
	/cineca/prod/opt/applications/macs/2.1.0/python--2.7.12/bin/macs2 "$@"
}
export MACS


workDir=$WORK/drip_project/
tmpDir=$WORK/drip_project/tmp
aligDir=$WORK/drip_project/alignment
logDir=$WORK/drip_project/logs
qualDir=$WORK/drip_project/qual
peakDir=$WORK/drip_project/peak_calling
rawDir=$WORK/drip_raw

cd ${workDir}

#create a logfile 
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>${logDir}/diff_exp_peaks_reg.log 2>&1


MACS bdgdiff --t1 ${peakDir}/regular_Model_rep1_treat_pileup.bdg --c1 ${peakDir}/regular_Model_rep1_control_lambda.bdg --t2 ${peakDir}/regular_Model_rep2_treat_pileup.bdg\
   --c2 ${peakDir}/regular_Model_rep2_control_lambda.bdg --d1 26592730 --d2 26592730 -g 76 -l 200 --o-prefix ${peakDir}/reg_Model_Diff_fragCtrl
  

MACS bdgdiff --t1 ${peakDir}/regular_noModel_rep1_treat_pileup.bdg --c1 ${peakDir}/regular_noModel_rep1_control_lambda.bdg --t2 ${peakDir}/regular_noModel_rep2_treat_pileup.bdg\
   --c2 ${peakDir}/regular_noModel_rep2_control_lambda.bdg --d1 26592730 --d2 26592730 -g 76 -l 200 --o-prefix ${peakDir}/reg_noModel_Diff_fragCtrl


