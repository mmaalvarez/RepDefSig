#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=10G

conda activate singularity

# Running Variational Autoencoder with dataset split (90 % training and 10% test) to find optimal parameters. Recommended to be run in Singularity environment
singularity exec tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python VAE_tensorflow1.py \
	--learning_rate 0.0005 \
	--batch_size 200 \
	--epochs 200 \
	--num_components 14 \
	--dataset_training 'TCGA_Hartwig_PCAWG_least45of56_14832samples_90split_balanced.txt' \
	--dataset_test 'TCGA_Hartwig_PCAWG_least45of56_14832samples_10split_balanced.txt' \
	--kappa 0.5 \
	--depth 1

# Running on complete dataset after finding optimal parameters
singularity exec tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python VAE_tensorflow1_nosplit_alldata.py \
	--learning_rate 0.0005 \
	--batch_size 200 \
	--epochs 200 \
	--num_components 14 \
	--dataset_training 'TCGA_Hartwig_PCAWG_least45of56_14832samples_nosplit_balanced.txt' \
	--kappa 0.5 \
	--depth 1
