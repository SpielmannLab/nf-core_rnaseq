#!/usr/bin/bash

#SBATCH -c 1
#SBATCH --mem=10GB
#SBATCH --partition=shortterm

PATH=$WORK/.omics/miniforge3/bin:$PATH 
source $WORK/.omics/miniforge3/etc/profile.d/conda.sh 
conda activate deseq2
module load singularity
cd /data/humangen_external/test_area/sygo/project37/for_GitHub

nextflow deseq2.nf -params-file /data/humangen_external/test_area/sygo/project37/for_GitHub/other_files/deseq2_params.yaml