## conda activate myenv



# Import statements


import matplotlib.pyplot as plt
import numpy as np

from keras.datasets import mnist
from keras.layers import (
        Activation, BatchNormalization, concatenate, Dense,
        Embedding, Flatten, Input, Multiply, Reshape, UpSampling3D, Concatenate)
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import Conv3D, Conv3DTranspose
from keras.layers.pooling import MaxPool3D
from keras.models import Model, Sequential
from keras.optimizers import Adam


import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(15)
run_opts = tf.compat.v1.RunOptions(report_tensor_allocations_upon_oom = True)

from keras.utils.vis_utils import plot_model
from random import randint
from numpy import (zeros, ones)


import nibabel as nib
img = nib.load("../data/brainweb_all_6.nii")
a = np.array(img.dataobj)

import glob as glob
from pandas import read_csv

# VegList = sorted(glob.glob('../brainweb_new/result/generated_r_5_*.csv'))
VegList = sorted(glob.glob('../brainweb_new/result/generated_r_5_1[12].csv'))
data = []

data = []
for f in VegList:
    print(f)
    dat_iter = np.asarray(read_csv(f, header=None))  ## np.asarray gives different format
    data.append(dat_iter.reshape(181, 217, 181, 12))

# Maybe add pd.concat if you wish

data_array = np.asarray(data)
# print(data_array)

data_orig = data_array[:,:,:,:,7]
data_array = data_array[:,:,:,:,0:5]
print(data_array.shape)

## data contains data with 2 people, 181 n_x, 217 n_y, 181 n_z, 6 TE/TR settings




## Model input dimensions:
img_rows = 128
img_cols = 128
img_depth = 128
channels = 1
settings_no = 6

img_shape = (img_rows, img_cols, img_depth, channels)
img_shape_8 = (img_rows, img_cols, img_depth, settings_no)
pix_shape = (1, 1, 1, settings_no)

z_dim = 100				## Size of noise vector z, used as input to the generator


TE_seq = (0.1, 0.2, 0.6, 0.7, 0.1, 0.2)
TR_seq = (0.2, .2, 0.2, 0.2, 0.03, .03)



## cGAN Discriminator:
def build_discriminator(img_shape):
    print("\n\n\nDiscriminator:")
    visible = Input(img_shape)
    
    print("\nImage shape is:")
    print(visible.shape)
    
    x1_n = Conv3D(filters=1, kernel_size=64, strides=32, activation="elu",padding='same',name="x1_n")(visible)    
    x2_n = Conv3D(filters=1, kernel_size=64, strides=1, activation="elu",padding='same',name="x2_n")(visible)
    x3_n = Conv3D(filters=1, kernel_size=64, strides=2, activation="elu",padding='same',name="x3_n")(visible)
    x4_n = Conv3D(filters=1, kernel_size=64, strides=4, activation="elu",padding='same',name="x4_n")(visible)
    x5_n = Conv3D(filters=1, kernel_size=64, strides=8, activation="elu",padding='same',name="x5_n")(visible)
    x6_n = Conv3D(filters=1, kernel_size=64, strides=16, activation="elu",padding='same',name="x6_n")(visible)
    
    ## My fix: faltten and then connect using - https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8662660
    x1_flat = Flatten()(x1_n)
    x2_flat = Flatten()(x2_n)
    x3_flat = Flatten()(x3_n)
    x4_flat = Flatten()(x4_n)
    x5_flat = Flatten()(x5_n)
    x6_flat = Flatten()(x6_n)

    x1_flat = Dense(units=1, activation="sigmoid")(x1_flat)
    x2_flat = Dense(units=1, activation="sigmoid")(x2_flat)
    x3_flat = Dense(units=1, activation="sigmoid")(x3_flat)
    x4_flat = Dense(units=1, activation="sigmoid")(x4_flat)
    x5_flat = Dense(units=1, activation="sigmoid")(x5_flat)
    x6_flat = Dense(units=1, activation="sigmoid")(x6_flat)
    
    final_1 = concatenate([x1_flat, x2_flat, x3_flat, x4_flat, x5_flat, x6_flat])
    final_2 = Dense(units=1, activation="sigmoid")(final_1)
    
    model = Model(inputs=visible, outputs=final_2)
    
    # compile model
    opt = Adam(lr=0.0002, beta_1=0.5)
    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
    plot_model(model, to_file='gan_dis_test.pdf', show_shapes=True, show_layer_names=True) # , expand_nested=True
    #print(model.summary())
    return model


#build_discriminator(img_shape)

# select real samples
def generate_real_samples(dataset, n_batch):
    # choose random instances
    ix = randint(0, dataset.shape[1]-128); iy = randint(0, dataset.shape[2]-128); iz = randint(0, dataset.shape[3]-128)
    # select images and labels
    X = dataset[:,ix:(ix+128), iy:(iy+128), iz:(iz+128)]	## (no_sample, X, Y, Z, settings_no)
    return X, ix, iy, iz








## cGAN generator: 
def build_generator(img_shape_gen):
    print("\n\n\nGenerator:")
    print("img_shape_gen")
    print(img_shape_gen)
    		
    visible = Input(img_shape_gen, name="img")
    TE_seq_1 = Input(img_shape_gen, name="TE")
    TR_seq_1 = Input(img_shape_gen, name="TR")
    
    vis_final = concatenate([visible, TE_seq_1, TR_seq_1])
    print(vis_final.shape)
    print("Check 1")
    
    vis_final = Reshape((128,128,128,settings_no*3))(vis_final)								## Check the formation
    
    print(vis_final.shape)
    print("\n\n\n\n\n")
    print(vis_final)

    x1 = Dense(units=128, activation="sigmoid")(vis_final)
    x1 = Dense(units=1, activation="sigmoid")(x1)
    final = Dense(units=1, activation="elu")(x1)

    model = Model(inputs=[visible, TE_seq_1, TR_seq_1], outputs=final, name="try_gen")  ## No, this would have same dim as ??

    plot_model(model, to_file='gan_gen_test.pdf', show_shapes=True, show_layer_names=True)
    # print(model.summary())				# Added by Subrata
    return model


# create the discriminator
d_model = build_discriminator(img_shape)
# create the generator
g_model = build_generator(img_shape_8)



dataset = data_array

X_real, ix, iy, iz = generate_real_samples(data_orig, 128)
print("X_real[:,:,:,:].shape: "); print(X_real[1:3,:,:,:].shape)
y_real = ones(((X_real[1:3,:,:,:]).shape[0], 1))


#print("Train: ")
# d_model.train_on_batch(X_real[1:3,:,:,:], y_real)

z_input = dataset[:,ix:(ix+128), iy:(iy+128), iz:(iz+128),:]
print(z_input[1,1,1,:].shape)
images = g_model.predict([z_input[1,1,1], TE_seq, TR_seq])						# predict outputs




#print("See\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")


#def define_gan(g_model, d_model):
#    d_model.trainable = False
#    gen_raw, gen_TE, gen_TR = g_model.input			    # get noise and label inputs from generator model
#    gen_output = g_model.output							# get image output from the generator model
#    print("gen_output.shape:"); print(gen_output.shape)  
#    ## (None, 1, 1, 1) -- not matching with input of d_model in the next line wants img_shape (128, 128, 1)
#    
#    gan_output = d_model(gen_output)					# connect image output and label input from generator as inputs to discriminator
#    model = Model([gen_raw, gen_TE, gen_TR], gan_output)# define gan model as taking noise and label and outputting a classification
#    model.compile(loss='binary_crossentropy', optimizer=Adam(lr=0.0002, beta_1=0.5))  	# compile model
#    return model



#def generate_fake_samples(generator, dataset, ix, iy, iz):
#	# z_input, labels_input = generate_latent_points(latent_dim, n_samples)		# generate points in latent space
#	z_input = dataset[:,ix:(ix+128), iy:(iy+128), iz:(iz+128),:]
#	images = generator.predict([z_input, TE_seq, TR_seq])						# predict outputs
#	y = zeros((n_samples, 1))													# create class labels
#	return images, y






#dataset = data_array

#X_real, ix, iy, iz = generate_real_samples(dataset, 128)
#print("X_real[:,:,:,:,0].shape: "); print(X_real[1:3,:,:,:,0].shape)
#y_real = ones(((X_real[1:3,:,:,:,0]).shape[0], 1))


##print("Train: ")
##d_model.train_on_batch(X_real[1:3,:,:,:,0], y_real)

## generate 'fake' examples
#[X_fake, labels], y_fake = generate_fake_samples(g_model, dataset, ix, iy, iz)
#print("Check 3")

## update discriminator model weights
#d_loss2, _ = d_model.train_on_batch([X_fake, labels], y_fake)
#print("Check 4")

## prepare points in latent space as input for the generator
#[z_input, labels_input] = generate_latent_points(latent_dim, n_batch)
#print("Check 5")

## create inverted labels for the fake samples
#y_gan = ones((n_batch, 1))
#print("Check 6")

## update the generator via the discriminator's error
#g_loss = gan_model.train_on_batch([z_input, labels_input], y_gan)







#def train1(g_model, d_model, gan_model, dataset, n_batch=128, n_iter=10):
#    # manually enumerate epochs
#    for i in range(n_iter):
#        # enumerate batches over the training set
#        for j in range(n_batch):
#            # get a 128x128x128 patch from the main data
#            # No, only take the flair image from the data. 
#            X_real = generate_real_samples(dataset, n_batch)
#            print("X_real[:,:,:,:,0].shape: ")
#            print(X_real[:,:,:,:,0].shape)
#            
#            # update discriminator model weights
#            y_real = ones((X_real.shape[0], 1))
#            print("y_real shape: ")
#            print(y_real.shape)
#            
#            
#            ## Problem here: Changed the format of output. 
#            ## See this also: https://github.com/fengwang/MCNN/blob/master/tutorial/tutorial.ipynb            
#            d_loss1, _ = d_model.train_on_batch(X_real[:,:,:,:,0], y_real)
#            # discriminator model should get the pd image - check which one is the rho density
#            print("Check 2")
#            
#            # generate 'fake' examples
#            [X_fake, labels], y_fake = generate_fake_samples(g_model, latent_dim, half_batch)
#            print("Check 3")
#            
#            # update discriminator model weights
#            d_loss2, _ = d_model.train_on_batch([X_fake, labels], y_fake)
#            print("Check 4")
#            
#            # prepare points in latent space as input for the generator
#            [z_input, labels_input] = generate_latent_points(latent_dim, n_batch)
#            print("Check 5")
#            
#            # create inverted labels for the fake samples
#            y_gan = ones((n_batch, 1))
#            print("Check 6")
#            
#            # update the generator via the discriminator's error
#            g_loss = gan_model.train_on_batch([z_input, labels_input], y_gan)





### size of the latent space
## latent_dim = 100
### create the discriminator
##d_model = build_discriminator(img_shape)
### create the generator
##g_model = build_generator(pix_shape)
### create the gan
##gan_model = define_gan(g_model, d_model)
### load image data
## dataset = load_real_samples()
### train model
## train1(g_model, d_model, gan_model, data_array, 128, 10)








