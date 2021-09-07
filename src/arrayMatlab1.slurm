#!/bin/bash
#SBATCH --mail-user=benjamin.k.kesler@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=0-1500
#SBATCH --output=Out_%A_%a.out

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load MATLAB
matlab -nodisplay -nosplash -r  "A0_segment_define_variables_streamlined_non_GUI(${SLURM_ARRAY_TASK_ID});"
