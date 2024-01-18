# This repo contains several workflows to perform several tasks:

1. [fastq2counts](fastq2counts), which is our configuration of [nf-core/rnaseq](https://nf-co.re/rnaseq/3.13.2) to generate count matrices from fastq files
2. [deseq2](deseq2) to perform differential expression analysis using the count matrices

## fastq2counts

### nf-core/configs: UzL OMICS Cluster Configuration

The rnaseq nf-core pipeline has been successfully configured for use on the UzL OMICS cluster at the University of Luebeck.

To use, run any nf-core pipeline with `-profile uzl_omics`. This will download and launch the `uzl_omics.config` which has been pre-configured with a setup suitable for the UzL OMICS cluster. Using this profile, docker images containing the required softwares will be downloaded, and converted to singularity images before execution of the pipeline.

### Software dependencies

Nextflow and Singularity are required to run nf-core scripts. Singularity is available using the environment module system on UzL OMICS cluster and can be made available by issuing the command below:

```bash
module load singularity
```

But nf-core requires Nextflow version 22.10.1 or higher, so you have to intall a more recent version first (as of Nov 2023, higher versions are not installed in the cluster). This can be done within a conda environment for example. Note: for Nextflow versions newer than 23.07.0-edge, it is necessary to mount the home directory using the command:

```bash
NXF_SINGULARITY_HOME_MOUNT=true
```

### Running the pipeline

You can use the following batch script as an example to run the pipeline (the example shown for the nf-core/rnaseq pipeline can also be extended to other nf-core pipelines):

```bash
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
    --outdir PATH_TO_YOUR_OUTPUT \
    -profile uzl_omics \
    -params-file PATH_TO_YOUR_PARAMS_FILE
```

## Below are non-mandatory information

> note:
> You will need access to the UzL OMICS cluster in order to run the pipeline. In doubt contact IT.

## deseq2

Fill out the [parameters file](./deseq2/deseq2_params.yaml) and start the job as follows:

```bash
sbatch deseq2_sbatch.sh
```
