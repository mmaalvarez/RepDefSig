# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# Modified from Mischan Vali Pour's Variational Autoencoder: github.com/lehner-lab/RDGVassociation/tree/main/somatic_component_extraction
#
#     adapted/modified from github.com/greenelab/tybalt/blob/master/tybalt_vae.ipynb
#
#     adapted/inspired from "Hand-On Maschine Learning with Scikit-Learn, Keras & Tensorflow" by Aurélien Géron December 2020
#
# Outputs different performance parameters
#
# Runs with tensorflow 1.15.5

# +
####import all important stuff

## important modules
# Python ≥3.5 is required
import sys
import argparse #for parsing
print(sys.version) #print version

# Scikit-Learn ≥0.20 is required
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

##keras stuff
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Input, Dense, Lambda, Layer, Activation, BatchNormalization
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras import backend as K
from tensorflow.keras import metrics, optimizers, losses
from tensorflow.keras.callbacks import Callback

print(keras.__version__)
print(tf.__version__)

# Common imports
import numpy as np
import os

#pandas
import pandas as pd

# for pearson
from scipy.stats import pearsonr

# +
####start getting info

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--learning_rate',
                    default=0.0005,
                    help='learning rate of the Adam optimizer')
parser.add_argument('-b', '--batch_size',
                    default=200,
                    help='number of samples to include in each learning batch')
parser.add_argument('-e', '--epochs',         
                    default=200,
                    help='how many times to cycle through the full dataset')
parser.add_argument('-n', '--num_components',
                    default=14,
                    help='latent space dimensionality (size)')
parser.add_argument('-t', '--dataset_training', 
                    default='../../1_coefficient_scaling/VAE_input_1000iters.tsv',
                    help='training dataset, put in full name + direc')
parser.add_argument('-r', '--dataset_real', 
                    default='../../1_coefficient_scaling/original_data_scaled.tsv',
                    help='real dataset (i.e. no permutations, just [-1,1]-scaled), put in full name + direc')
parser.add_argument('-k', '--kappa',
                    default=0.5,
                    help='kappa, how strongly to linearly ramp up the KL loss after each epoch')
parser.add_argument('-d', '--depth',
                    default=1,
                    help='define whether there should be a layer between latent and input/ouput layer, if yes depth=2, else depth=1')
args = parser.parse_args(args=[])

# +
######## Set hyper parameters

learning_rate = float(args.learning_rate)
batch_size = int(args.batch_size)
epochs = int(args.epochs)
latent_dim = int(args.num_components)
dataset_training = args.dataset_training
dataset_real = args.dataset_real
kappa = float(args.kappa)
beta = K.variable(0) #KL loss weighting at first epoch

#decide on whether there should be a layer before input and latent layer
#make this layer
depth = int(args.depth)
hidden_dim =  latent_dim*2 #douple size of latent dimensions

# Random seed
seed = int(np.random.randint(low=0, high=10000, size=1))
np.random.seed(seed)


# +
######## upload and parse input data

## upload [-1,1]-scaled full PERMUTED data (for training)
full_df = pd.read_csv(dataset_training, sep='\t')

## upload [-1,1]-scaled ORIGINAL data (for final signature extraction)
real_df = pd.read_csv(dataset_real, sep='\t')
# store sample names column
sample_id = real_df.drop(real_df.columns[1:], axis=1)

# convert to numpy array
train_df_input = np.array(full_df.drop(full_df.columns[0], axis=1))
real_df_input = np.array(real_df.drop(real_df.columns[0], axis=1))

print(train_df_input.shape)
print(real_df_input.shape)

# Set architecture dimensions
original_dim = train_df_input.shape[1]


# +
######## def functions and classes

# Function for reparameterization trick to make model differentiable
# custom layer to sample the codings given mean and log_var
# samples from a normal distribution
def sampling(args):
    # Function with args required for Keras Lambda function
    z_mean, z_log_var = args
    # sample epsilon for a normal distribution with same shape as input columns aka phenotypes aka somatic features
    epsilon = K.random_normal(shape=tf.shape(z_mean), mean=0.,
                              stddev=1)
    # estimate latent codings by adding epsilon to mu and standard deviation which was learned
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z

## custom layer with loss
class CustomVariationalLayer(Layer):
    def __init__(self, **kwargs):
        self.is_placeholder = True
        super(CustomVariationalLayer, self).__init__(**kwargs)
    def vae_loss(self, x_input, x_decoded):
        reconstruction_loss = original_dim * \
                              metrics.mse(x_input, x_decoded) #using here mean squared error, multiplying with original dim since tensor already uses mean
        kl_loss = -0.5 * K.sum(1 + latent_log_encoded -
                               K.square(latent_mean_encoded) -
                               K.exp(latent_log_encoded), axis=-1)
        return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss)) #combining reconstruction and KL loss and taking the mean
    def call(self, inputs):
        x = inputs[0]
        x_decoded = inputs[1]
        loss = self.vae_loss(x, x_decoded)
        self.add_loss(loss, inputs=inputs)
        return x

# implement warmup to slowly ramp up the KL loss (beta=0 is vanilla autoencoder and beta=1 full VAE)
# KL will be weighted by KL*beta
# modified code from https://github.com/keras-team/keras/issues/2595 from https://github.com/greenelab/tybalt/blob/master/tybalt_vae.ipynb
# idea of ladder VAE from https://arxiv.org/abs/1602.02282class WarmUpCallback(Callback):
class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        self.beta = beta
        self.kappa = kappa 
    def on_epoch_end(self, epoch, logs={}): #on_epoch_begin, # Behavior on each epoch
        if K.get_value(self.beta) <= 1:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)

# +
#### ENCODER ####
#first dense layer to get mean and log var
#batch normalization
#relu as activation function

pheno_input = Input(shape=(original_dim, )) # shape == number of neurons in input layer (==DNA repair mark coefficients)
z_shape = latent_dim # number of signatures (==neurons) in latent layer

if depth == 1:
    latent_mean = Dense(latent_dim,
                        kernel_initializer='glorot_uniform')(pheno_input)
    latent_log_var = Dense(latent_dim,
                           kernel_initializer='glorot_uniform')(pheno_input)
elif depth == 2:
    hidden_dense = Dense(hidden_dim,
                         kernel_initializer='glorot_uniform')(pheno_input)
    hidden_dense_batchnorm = BatchNormalization()(hidden_dense)
    hidden_enc = Activation('relu')(hidden_dense_batchnorm)
    latent_mean = Dense(latent_dim,
                         kernel_initializer='glorot_uniform')(hidden_enc)
    latent_log_var = Dense(latent_dim,
                            kernel_initializer='glorot_uniform')(hidden_enc)

# batch normalization and activation for mu and log variance
latent_mean_batchnorm = BatchNormalization()(latent_mean)
latent_mean_encoded = Activation('relu')(latent_mean_batchnorm)

latent_log_batchnorm = BatchNormalization()(latent_log_var)
latent_log_encoded = Activation('relu')(latent_log_batchnorm)

# reparameterization
latent_codings = Lambda(sampling,
           output_shape=(z_shape, ))([latent_mean_encoded, latent_log_encoded])

# +
#### DECODER ####
#need tanh as activation function since our output values have to be between -1 and 1 (as the input)
if depth == 1:
    decoder_to_reconstruct = Dense(original_dim,
                                   kernel_initializer='glorot_uniform',
                                   activation='tanh')
elif depth == 2:
    decoder_to_reconstruct = Sequential()
    decoder_to_reconstruct.add(Dense(hidden_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='relu',
                                     input_dim=latent_dim))
    decoder_to_reconstruct.add(Dense(original_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='tanh'))

phenos_reconstruct = decoder_to_reconstruct(latent_codings)

# +
############### build up autoencoder ###############

adam = optimizers.Adam(lr=learning_rate)
vae_layer = CustomVariationalLayer()([pheno_input, phenos_reconstruct])
vae = Model(pheno_input, vae_layer)
vae.compile(optimizer=adam, loss=None, loss_weights=[beta])

# +
############### fit Model ##########################

history = vae.fit(train_df_input,
                  shuffle=True,
                  epochs=epochs,
                  batch_size=batch_size,
                  # validation_data=(validation_df_input, None),
                  callbacks=[WarmUpCallback(beta, kappa)])

# evaluate final loss
training_loss= vae.evaluate(train_df_input)
# validation_loss= vae.evaluate(validation_df_input)

print(training_loss)
# print(validation_loss)
# +
######################################################################
## QC
######################################################################


############### output "signatures" from real (i.e. not permuted) data #############

# Mischan: "latent_mean_encoded layer of the VAE is equivalent to the exposures from the NMF"

# extract 'signatures' (latent layer means) from original (-1,1 scaled) data
encoder = Model(pheno_input, latent_mean_encoded)
encoded_real_df = encoder.predict_on_batch(real_df_input)
# to pandas, and rename columns ('signatures')
encoded_real_df = pd.DataFrame(encoded_real_df, columns= range(1, latent_dim+1)).add_prefix('sig')
# normalize every row (sample) to sum up to 1 for plotting
encoded_real_df = encoded_real_df.div(encoded_real_df.sum(axis=1), axis=0)
# append sample names column
encoded_real_df = pd.concat([sample_id, encoded_real_df], axis=1)

output_df_direc = os.path.join("signatures__" + str(latent_dim) + "_components_" + 
                               str(learning_rate) + "_lr_" + str(batch_size) + "_batchsize_" + 
                               str(epochs) + "_epochs_" + str(kappa) + "_kappa_" +
                               str(1) + "_beta_" + str(hidden_dim) + "_hiddenDim_" + str(depth) + "_depth" + ".tsv")
#save output
encoded_real_df.to_csv(output_df_direc, sep='\t', index= False)




# ############# pearson reconstruction vs input #################################################################

# # encoding again
# val_encoded = encoder.predict_on_batch(validation_df_input)

# ####decoder generative model
# decoder_input = Input(shape=(latent_dim, ))  # can generate from any sampled z vector
# _x_decoded_mean = decoder_to_reconstruct(decoder_input)
# decoder = Model(decoder_input, _x_decoded_mean)

# # reconstruction
# val_reconstructed = decoder.predict(val_encoded) 
# val_reconstructed = pd.DataFrame(val_reconstructed, columns=train_df.drop(['sample_short'], axis=1).columns)
# print(val_reconstructed)

# validation_df= pd.DataFrame(validation_df_input, columns=train_df.drop(['sample_short'], axis=1).columns)

# ## check the mean pearson between input and reconstructed, pearson for each sample
# r = [pearsonr(val_reconstructed.iloc[x, :],
#               validation_df.iloc[x, :])[0] for x in range(val_reconstructed.shape[0])]
# r_mean = np.mean(np.array(r))
# print(r_mean)



# ############### check Correlation with NMF signatures #############

# ##upload file
# file_NMF_sig = '../../../NMF/good_model_score_table.tsv'
# NMF_sig = pd.read_csv(file_NMF_sig,sep='\t')

# ##take the latent encoding matrix and this information to do pearson correlations
# ###### upload latent encodings of the mean
# encoded_mean_train = encoder.predict_on_batch(train_df_input)
# encoded_mean_validation = encoder.predict_on_batch(validation_df_input)

# encoded_mean_train_df = pd.DataFrame(encoded_mean_train, columns= range(1, latent_dim+1))
# encoded_mean_validation_df = pd.DataFrame(encoded_mean_validation, columns= range(1, latent_dim+1))

# #add sample short info
# train_df_newindex = train_df.copy()
# validation_df_newindex = validation_df.copy()

# train_df_newindex.index =  range(train_df_newindex.shape[0])
# encoded_mean_train_df = encoded_mean_train_df.assign(sample_short= train_df_newindex.sample_short)

# validation_df_newindex.index =  range(validation_df_newindex.shape[0])
# encoded_mean_validation_df = encoded_mean_validation_df.assign(sample_short=validation_df_newindex.sample_short)

# #combine training and test mean encodings
# encoded_mean_dataset_df = pd.concat([encoded_mean_train_df,encoded_mean_validation_df])

# ##add NMF_sig to encodings
# encoded_encodings_NMF_sig_ID = pd.merge(encoded_mean_dataset_df, NMF_sig, on='sample_short')

# ###### do pearson correlation between all columns except for sample_short ####
# pearson_encodings_vs_NMF_sigs = encoded_encodings_NMF_sig_ID.corr(method='pearson')

# ###estimate maximum correlation with the VAE_mean_encodings with each of the NMF_sig
# sig1 = max(pearson_encodings_vs_NMF_sigs.sig1.iloc[0:latent_dim])
# sig2 = max(pearson_encodings_vs_NMF_sigs.sig2.iloc[0:latent_dim])


# ##create df from it
# max_correlation_NMF_sigs = pd.DataFrame(np.array([[sig1,sig2]]),
#                                           columns=['sig1', 'sig2'])
# max_correlation_NMF_sigs = np.absolute(max_correlation_NMF_sigs) ##take absolute values/no minus

# ##get mean value
# mean_pearson_NMF_sigs = max_correlation_NMF_sigs.mean(axis=1)
