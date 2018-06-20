#!/bin/bash
#!
#! Name of the job:
#SBATCH -J exp-analyser
#SBATCH --mail-type=FAIL
#SBATCH -m cyclic:block
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem-per-cpu=1500
#SBATCH --cpus-per-task=1
#SBATCH -p CIAO
#SBATCH --mail-user=fp305@cam.ac.uk

CMD="hostname"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
echo -e "\nExecuting command:\n==================\n$SBATCH_CMD_CIAO\n"

echo -e "\n$SBATCH_CMD_CIAO\n"

eval $SBATCH_CMD_CIAO

#hostname
#sleep 100
