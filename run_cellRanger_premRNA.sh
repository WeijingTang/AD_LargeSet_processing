#!/bin/bash
## No need to change this unless you could also run with a different PI account. 
## Replace MY_PI_SUNetID_or_Project_ID with the PI/project to be charged.
#SBATCH --account=icobos

## Sherlock partitions are normal, sfgf and ... 
#SBATCH --partition=sfgf,owners

## Set job time to 24 hours
#SBATCH --time=48:00:00

## Set a name for the job, visible in `squeue`
#SBATCH --job-name="run_cellranger_with_intron"

## The following settings are optimal for *most* software, we want one task 
## to have one or more cores available for that task to fork or use threads.
## One node, -N also works
#SBATCH --nodes=1
## Number of tasks, -n also works. For 10x local mode use only 1 taks per node.
#SBATCH --ntasks=1
## To take advantage of multi-threading, use 16 CPU/core per task, -c also works.
#SBATCH --cpus-per-task=16

## There are to ways to specify memory, --mem= and --mem-per-cpu=
## --mem is usually enough since total job memory is easy to specify 
## this way.
## 256GB of RAM is recommended for 10x runs.
#SBATCH --mem=256G

## Open an job array for multiple 10x Channel
## Always remember to check the number here!!!
## The range of numbers to use for the array, equal to number of 10x channels

## Specify log file location to help with better logging of errors/outputs.
#SBATCH -o CRcount-%A-%a.out
#SBATCH -e CRcount-%A-%a.err

## Put any module here, anaconda environment example shown here:
## Modules needed for the pipeline to run
#SBATCH --array=1-8
## Optionally activate environment
# source  activate py3.6

## First specify the reference and base path
#base_path=/oak/stanford/groups/icobos/Single_Cell_data/Raw_FASTQ

## enter base path so make sure 10x sample list is stored here!
#cd ${base_path}

## Use bash command to get input sample/channel name, fastqpath, expected cell number from a SAMPLELIST, i.e., 10x sample list
## Sample list should be comma-separated (csv) file that have no header
## Sample list should contain 3 columns for ID, fastq_path(absolute or relative to base_path), expected_cell_number (no comma please) for each channel

#SAMPLELIST=sample_sheet.txt
#SEED=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SAMPLELIST)
#name=$(echo "$SEED" | cut -d/ -f9)
#| cut -d/ -f8
#bowtie2-build $name".fa" $name
#mkdir "$name"
#mv *.bt2 "$name"
#bowtie2 -p 8 -x $name"/"$name $SEED > $name".sam"
FOLDER='/scratch/groups/icobos/Marazzi_lab/Quantification/Exon_and_Intron/'
cd ${FOLDER}
SAMPLELIST=/scratch/groups/icobos/Marazzi_lab/Quantification/sample_sheet.txt
SEED=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SAMPLELIST)
name=$(echo "$SEED" | cut -d/ -f9)
#mkdir ${name}
#cd ${name}
for SP in $SEED
do
SAMPLE=$(echo "$SP" | cut -d/ -f7)
mkdir ${SAMPLE}_concat
cd ${SAMPLE_concat}
cat ${SEED}/${SAMPLE}_S*_L001_I1_001.fastq.gz ${SEED}/${SAMPLE}_S*_L002_I1_001.fastq.gz ${SEED}/${SAMPLE}_S*_L003_I1_001.fastq.gz ${SEED}/${SAMPLE}_S*_L004_I1_001.fastq.gz \
  > ${FOLDER}/${SAMPLE}_concat/${SAMPLE}_S0_L000_I1_001.fastq.gz
cat ${SEED}/${SAMPLE}_S*_L001_R1_001.fastq.gz ${SEED}/${SAMPLE}_S*_L002_R1_001.fastq.gz ${SEED}/${SAMPLE}_S*_L003_R1_001.fastq.gz ${SEED}/${SAMPLE}_S*_L004_R1_001.fastq.gz \
  > ${FOLDER}/${SAMPLE}_concat/${SAMPLE}_S0_L000_R1_001.fastq.gz
cat ${SEED}/${SAMPLE}_S*_L001_R2_001.fastq.gz ${SEED}/${SAMPLE}_S*_L002_R2_001.fastq.gz ${SEED}/${SAMPLE}_S*_L003_R2_001.fastq.gz ${SEED}/${SAMPLE}_S*_L004_R2_001.fastq.gz \
  > ${FOLDER}/${SAMPLE}_concat/${SAMPLE}_S0_L000_R2_001.fastq.gz
cd ${FOLDER}/${SAMPLE}_concat
/scratch/groups/icobos/Marazzi_lab/Cell_Ranger_5.0/cellranger-5.0.0/cellranger count --id=$SAMPLE --transcriptome=/scratch/groups/icobos/Marazzi_lab/Mouse_reference/mm10 --fastqs=${FOLDER}/${SAMPLE}_concat --sample=$SAMPLE --expect-cells=10000 --include-introns

echo "Completed successfully."

done

## deactivating python virtual environment
# source deactivate

# What variables does SLURM set?
env | grep SLURM | sort -V


