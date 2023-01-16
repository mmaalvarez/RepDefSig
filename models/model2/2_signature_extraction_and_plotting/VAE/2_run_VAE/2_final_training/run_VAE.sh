#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=10G

conda activate singularity


# Running Variational Autoencoder without dataset split, and optimal parameters
singularity exec ../singularity_container/tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python VAE_tensorflow1.py \
	--learning_rate  \
	--batch_size  \
	--epochs  \
	--num_components  \
	--dataset_training '../../1_coefficient_scaling/VAE_input_1000iters.tsv' \
	--dataset_real '../../1_coefficient_scaling/original_data_scaled.tsv' \
	--kappa  \
	--depth 
