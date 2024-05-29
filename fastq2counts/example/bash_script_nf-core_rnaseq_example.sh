#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=10GB
#SBATCH --partition=shortterm

# Make conda available within a slurm job within which Nextflow is installed
PATH=$WORK/.omics/miniforge3/bin:$PATH
source $WORK/.omics/miniforge3/etc/profile.d/conda.sh
conda activate YOUR_NEXTFLOW_ENV     #you have to activate your environment with a Nextflow version 22.10.1 or higher

module load singularity

cd PATH_TO_YOUR_LAUNCHDIR
export NXF_SINGULARITY_HOME_MOUNT=true
nextflow run nf-core/rnaseq \
    -profile uzl_omics \
    -params-file PATH_TO_YOUR_PARAMS_FILE