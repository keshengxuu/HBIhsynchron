#!/bin/bash
#SBATCH --output=output_%A_%a.out
#SBATCH --error=error_%A_%a.out
#SBATCH --cpus-per-task=3        # Number of CPU cores per task
#SBATCH --time=16:40:00            # Time limit hrs:min:sec
#SBATCH --array=0-23%12

export OMP_NUM_THREADS=3

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

srun -u python $1  $2 


