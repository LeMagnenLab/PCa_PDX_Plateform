#!/bin/bash
#SBATCH --job-name=final_export
#SBATCH --error=/scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/Projects/Prostate/Xenograft_Models/exp/All_PDOXs_PDOXOs/10_Final_Export/final_export_SLURM.err
#SBATCH --output=/scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/Projects/Prostate/Xenograft_Models/exp/All_PDOXs_PDOXOs/10_Final_Export/final_export_SLURM.out
#SBATCH --cpus-per-task=4      # Number of threads reserved to perform parallel operations (max = 20)
#SBATCH --mem-per-cpu=32GB     #This is the memory reserved per core. (4x52 = 208GB) (max = 256GB total for one node, but -20% for the Operating System)
#SBATCH --time=1:00:00       #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=romuald.parmentier@usb.ch

echo ${SLURM_JOBID}
echo ${TMPDIR}
echo ${id}

ml purge					            # Purge before loading modules

echo "Load LoupeR env"
START=$(date +%s)

export R_LIBS="/scicore/home/wykopa75/parmen0000/miniconda3/envs/LoupeR_ubuntu/lib/R/library"

# Initialize conda for the shell
eval "$(conda shell.bash hook)"

conda activate LoupeR_ubuntu

# Verify that the environment is activated and infercnv is installed
echo "Current Conda Environment: $(conda env list | grep \* | awk '{print $1}')"
if [[ $(conda env list | grep \* | awk '{print $1}') != "LoupeR_ubuntu" ]]; then
echo "Conda environment 'LoupeR_ubuntu' not activated!"
exit 1
fi

cd /scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/Projects/Prostate/Xenograft_Models/bin/All_PDOXs_PDOXOs

Rscript 10_Final_Export.R

END=$(date +%s)
DIFF=$(( $END - $START ))
UHR=$(( $DIFF / 3600 )); MIN=$(( $(( $DIFF % 3600 )) / 60 )); SEC=$(( $(( $DIFF % 3600 )) % 60 ))  # Register script duration
echo " ---->  it took  ${UHR}h ${MIN}m ${SEC}s for this run to complete <---- "

ml purge					     # Purge at the end

echo "End LoupeR final export"