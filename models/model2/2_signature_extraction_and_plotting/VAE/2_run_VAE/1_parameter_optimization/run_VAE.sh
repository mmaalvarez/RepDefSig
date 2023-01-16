#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=10G

conda activate singularity


# Running Variational Autoencoder with dataset split (90 % training and 10% validation) to find optimal parameters
singularity exec ../singularity_container/tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg python VAE_tensorflow1.py \
	--learning_rate 0.0005 \
	--batch_size 200 \
	--epochs 200 \
	--num_components 6 \
	--dataset_training '../../1_coefficient_scaling/VAE_input_1000iters.tsv' \
	--validation 0.1 \
	--kappa 0.5 \
	--depth 1
