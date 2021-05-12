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



#arr = np.random.random([4,3,2,3])
#for idx, value in np.ndenumerate(arr):
#    print(idx, value)





#import nibabel as nib
#img = nib.load("../data/brainweb_all_6.nii")
#a = np.array(img.dataobj)

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
data_array = data_array[:,:,:,:,0:6]
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

TE_seq = np.array([0.1, 0.2, 0.6, 0.7, 0.1, 0.2])
TR_seq = np.array([0.2, .2, 0.2, 0.2, 0.03, .03])



## cGAN Discriminator:
def build_discriminator(img_shape):
    print("\n\n\nDiscriminator:")
    visible = Input(img_shape)
    
    print("\nImage shape is:")
    print(visible.shape)
    
    x1_n = Conv3D(filters=1, kernel_size=4, strides=2, activation="elu",padding='same',name="x1_n")(visible)    
    x2_n = Conv3D(filters=1, kernel_size=4, strides=1, activation="elu",padding='same',name="x2_n")(visible)
    
    ## My fix: faltten and then connect using - https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8662660
    x1_flat = Flatten()(x1_n)
    x2_flat = Flatten()(x2_n)
    
    x1_flat = Dense(units=1, activation="sigmoid")(x1_flat)
    x2_flat = Dense(units=1, activation="sigmoid")(x2_flat)
    
    final_1 = concatenate([x1_flat, x2_flat])
    final_2 = Dense(units=1, activation="sigmoid")(final_1)
    
    model = Model(inputs=visible, outputs=final_2)
    
    # compile model
    opt = Adam(lr=0.02, beta_1=0.5)
    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
    plot_model(model, to_file='gan_dis_test.pdf', show_shapes=True, show_layer_names=True) # , expand_nested=True
    #print(model.summary())
    return model


#build_discriminator(img_shape)



## cGAN generator: 
def build_generator(img_shape_gen):
    print("\n\n\nGenerator:")
    print("img_shape_gen")
    print(img_shape_gen)
    
    visible = Input(img_shape_gen, name="img")
    TE_seq_1 = Input(img_shape_gen, name="TE")
    TR_seq_1 = Input(img_shape_gen, name="TR")
    
    vis_final = concatenate([visible, TE_seq_1, TR_seq_1])
    
    # vis_final = Reshape((128,128,128,settings_no*3))(vis_final)								## Check the formation
    
    x1 = Dense(units=1, activation="sigmoid")(vis_final)
    final = Dense(units=1, activation="elu")(x1)

    model = Model(inputs=[visible, TE_seq_1, TR_seq_1], outputs=final, name="try_gen")  ## No, this would have same dim as ??

    plot_model(model, to_file='gan_gen_test.pdf', show_shapes=True, show_layer_names=True)
    # print(model.summary())				# Added by Subrata
    return model





# create the discriminator
d_model = build_discriminator(img_shape)
# create the generator
g_model = build_generator(settings_no)




dataset = data_array



# select samples
def generate_real_samples(dataset, i, n_batch):
    # choose random instances
    ix = randint(0, dataset.shape[1]-128); iy = randint(0, dataset.shape[2]-128); iz = randint(0, dataset.shape[3]-128)
    # select images and labels
    X = dataset[i, ix:(ix+128), iy:(iy+128), iz:(iz+128)]	## (no_sample, X, Y, Z, (settings_no))
    return X, np.array(ix), np.array(iy), np.array(iz)





X_real, ix, iy, iz = generate_real_samples(data_orig, 0, 128) 
X_real = np.expand_dims(X_real, axis=0)

for i in range(3):			## should be range(15-1)
	X_real_tmp, ix_tmp, iy_tmp, iz_tmp = generate_real_samples(data_orig, 0, 128)
	X_real_tmp = np.expand_dims(X_real_tmp, axis=0)
	X_real = np.concatenate([X_real, X_real_tmp], axis=0)
	ix = np.append(ix, ix_tmp); iy = np.append(iy, iy_tmp); iz = np.append(iz, iz_tmp);


print("X_real.shape: "); print(X_real.shape)

y_real = ones(X_real.shape[0], dtype=int)   ### original images - 15x128x128x128 sample


print("Train: ")
d_model.train_on_batch(X_real, y_real)
print("Train ended\n\n\n\n\n")

print("data_array.shape")
print(data_array.shape)
print(ix[0]); print(iy[0]); print(iz[0])
X = data_array[0, ix[0]:(ix[0]+128), iy[0]:(iy[0]+128), iz[0]:(iz[0]+128),:]
print("X.shape"); print(X.shape)

## Feed the whole X by reshaping: 



print("X[0,0,0,:].shape"); print((X[0,0,0,:]).shape); print(TE_seq.shape); 

tmp_X = X[0,0,0,:]

tmp_X = np.expand_dims(tmp_X, axis=0)
TE_seq_expand = np.expand_dims(TE_seq, axis=0)
TR_seq_expand = np.expand_dims(TR_seq, axis=0)
print("tmp_X.shape"); print(tmp_X.shape)

g_model.predict([tmp_X, TE_seq_expand, TR_seq_expand])
print("Check seperate\n\n\n\n\n")




tmp_X_all = np.reshape(X, (128*128*128, settings_no))
TE_seq_all = np.repeat(TE_seq_expand, 128*128*128, axis=0)
TR_seq_all = np.repeat(TR_seq_expand, 128*128*128, axis=0)
print("tmp_X_all.shape"); print(tmp_X_all.shape)
print("TE_seq_all.shape"); print(TE_seq_all.shape)

#fake_data = g_model.predict([tmp_X_all, TE_seq_all, TR_seq_all])
#print("Check seperate\n\n\n\n\n")

#fake_data = np.reshape(fake_data, X.shape[:-1])




# Generated sample: 
def generate_fake_samples(generator, dataset, i, ix, iy, iz, n_batch):
    X = dataset[i, ix:(ix+128), iy:(iy+128), iz:(iz+128),:]## (no_sample, X, Y, Z, (settings_no))
    shape_x = X.shape
    for j in range(shape_x[0]):
        for k in range(shape_x[1]):
            for l in range(shape_x[2]):
                tmp_X = X[j,k,l,:]
                tmp_X = np.expand_dims(tmp_X, axis=0)
                fake_data[j,k,l] = generator.predict(X[j,k,l,:], TE_seq_expand, TR_seq_expand)
    print("Check5")
    return fake_data


# Generated sample: 
def generate_fake_samples_2(generator, dataset, i, ix, iy, iz, n_batch):
    X = dataset[i, ix:(ix+128), iy:(iy+128), iz:(iz+128),:]## (no_sample, X, Y, Z, (settings_no))
    shape_x = X.shape
    tmp_X_all = np.reshape(X, (128*128*128, settings_no))
    TE_seq_all = np.repeat(TE_seq, 128*128*128, axis=0)
    TR_seq_all = np.repeat(TR_seq, 128*128*128, axis=0)

    fake_data = generator.predict([tmp_X_all, TE_seq_all, TR_seq_all])

    fake_data = np.reshape(fake_data, shape_x[:-1])
    return fake_data






def define_gan_2(g_model, d_model):
    d_model.trainable = False
    gen_raw, gen_TE, gen_TR = g_model.input				# get noise and label inputs from generator model
    print("gen_raw.shape:"); print(gen_raw.shape)
    gen_output = g_model.output							# get image output from the generator model
    print("gen_output.shape:"); print(gen_output.shape)
    
    gan_output = d_model(gen_output)					# connect image output and label input from generator as inputs to discriminator
    model = Model([gen_raw, gen_TE, gen_TR], gan_output)# define gan model as taking noise and label and outputting a classification
    model.compile(loss='binary_crossentropy', optimizer=Adam(lr=0.0002, beta_1=0.5))  	# compile model
    return model



gan_model = define_gan(g_model, d_model)		# create the gan




	#for i in range(3):			## should be range(15-1)
	#	X_fake_tmp = generate_fake_samples(g_model, data_array, 0, ix[i], iy[i], iz[i], 128)
	#	X_fake_tmp = np.expand_dims(X_fake_tmp, axis=0)
	#	print("X_fake_tmp.shape"); print(X_fake_tmp.shape)
	#	X_real = np.concatenate([X_real, X_fake_tmp], axis=0)
	#	

	#print("X_real.shape: "); print(X_real.shape)

	#y_real = zeros(X_real.shape[0], dtype=int)   ### fake images


	#print("Train: ")
	#d_model.train_on_batch(X_real, y_real)
	##### 




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








