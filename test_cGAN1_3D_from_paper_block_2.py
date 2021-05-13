## conda activate myenv

"""
Block version with generator also. 
So there would be a huge array in generator
"""

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


import glob as glob
from pandas import read_csv

VegList = sorted(glob.glob('../brainweb_new/result/generated_r_5_1[12].csv'))
data = []
for f in VegList:
    print(f)
    dat_iter = np.asarray(read_csv(f, header=None))  ## np.asarray gives different format
    data.append(dat_iter.reshape(181, 217, 181, 12))

data_array = np.asarray(data)

data_orig = data_array[:,:,:,:,7]			## Suppose this is rho-weighted settings image
data_array = data_array[:,:,:,:,0:6]		## Suppose these are all other train settings image
print(data_array.shape)

## data contains data with 2 people, 181 n_x, 217 n_y, 181 n_z, 6 TE/TR settings




## Model input dimensions:
img_rows = 128; img_cols = 128; img_depth = 128
channels = 1
settings_no = 6
no_of_sample = data_orig.shape[0]
no_of_patch = 3  		## should be 15 according to the paper

img_shape = (img_rows, img_cols, img_depth, channels)
img_shape_8 = (img_rows, img_cols, img_depth, settings_no)
pix_shape = (1, 1, 1, settings_no)

TE_seq = np.array([0.1, 0.2, 0.6, 0.7, 0.1, 0.2])
TR_seq = np.array([0.2, .2, 0.2, 0.2, 0.03, .03])






## cGAN Discriminator:
def build_discriminator(img_shape_dis):
    print("\n\n\nDiscriminator:"); 
    visible = Input(img_shape_dis)
    print("\nImage shape is:"); print(visible.shape)		# (None, 128, 128, 128, 1)
    
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


subject_no = 0; 		## i-th patient's data

d_model = build_discriminator(img_shape)



"""
Generate real patch with corresponding indices
i-th person data
"""
def generate_real_samples(dataset, subject_no, n_batch):
    
    ix = randint(0, dataset.shape[1]-128); iy = randint(0, dataset.shape[2]-128); iz = randint(0, dataset.shape[3]-128)
    X_real = dataset[subject_no, ix:(ix+128), iy:(iy+128), iz:(iz+128)]	## (no_sample, X, Y, Z)
    X_real = np.expand_dims(X_real, axis=0)
    ix = np.array(ix); iy = np.array(iy); iz = np.array(iz)
    
    
    for i in range(1,no_of_patch):
        ix_tmp = randint(0, dataset.shape[1]-128); iy_tmp = randint(0, dataset.shape[2]-128); iz_tmp = randint(0, dataset.shape[3]-128)
        X_real_tmp = dataset[subject_no, ix_tmp:(ix_tmp+128), iy_tmp:(iy_tmp+128), iz_tmp:(iz_tmp+128)]	## (no_sample, X, Y, Z)
        X_real_tmp = np.expand_dims(X_real_tmp, axis=0)
    
        X_real = np.concatenate([X_real, X_real_tmp], axis=0)
        ix = np.append(ix, ix_tmp); iy = np.append(iy, iy_tmp); iz = np.append(iz, iz_tmp);

    return X_real, np.array(ix), np.array(iy), np.array(iz)


X_real, ix, iy, iz = generate_real_samples(data_orig, 0, 128)



print("X_real.shape: "); print(X_real.shape)

y_real = ones(X_real.shape[0], dtype=int)   ### original images - 15x128x128x128 sample
print("y_real.shape: "); print(y_real.shape)

print("Train: ")
d_model.train_on_batch(X_real, y_real)					## no_of_patch is the sample size
print("Real Train ended\n\n\n\n\n")








## cGAN generator: 
def build_generator(img_shape_gen):
    print("\n\n\nGenerator:")
    print("img_shape_gen"); print(img_shape_gen)
    
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


g_model = build_generator(img_shape_8)

print("Check 1")






"""
Generated fake sample: 
input: 	whole data: no_of_sample  x  n_x  x  n_y  x  n_z  x  settings_no
		subject_no: the number of subject over which it is taken.
		ix, iy, iz: numpy arrays, each of length no_of_patch
		
inter:	tmp_X_all: (128x128x128) x settings_no

output: fake_data:  128 x 128 x 128 
"""

TE_seq_expand = np.expand_dims(TE_seq, axis=0)
TR_seq_expand = np.expand_dims(TR_seq, axis=0)
TE_seq_all = np.tile(TE_seq_expand, (no_of_patch, 128, 128, 128, 1))
TR_seq_all = np.tile(TR_seq_expand, (no_of_patch, 128, 128, 128, 1))
#TE_seq_all = np.expand_dims(TE_seq_all, axis=0)
#TR_seq_all = np.expand_dims(TR_seq_all, axis=0)



def generate_inputs(generator, dataset, subject_no, ix, iy, iz, n_batch):
    for i in range(len(ix)):
        ix_tmp = ix[i]; iy_tmp = iy[i]; iz_tmp = iz[i];
        X_fake_tmp = dataset[subject_no, ix_tmp:(ix_tmp+128), iy_tmp:(iy_tmp+128), iz_tmp:(iz_tmp+128),:]  
        X_fake_tmp = np.expand_dims(X_fake_tmp, axis=0)
        
        if i==0:
            tmp_X_all = X_fake_tmp
        else:
            tmp_X_all = np.concatenate([tmp_X_all, X_fake_tmp], axis=0)
    
    print("tmp_X_all"); print(tmp_X_all.shape)
    return tmp_X_all



def generate_fake_samples_2(generator, dataset, subject_no, ix, iy, iz, n_batch):
    tmp_X_all = generate_inputs(generator, dataset, subject_no, ix, iy, iz, n_batch)
    fake_data = generator.predict([tmp_X_all, TE_seq_all, TR_seq_all])
    return fake_data






abc = generate_fake_samples_2(g_model, data_array, 0, ix, iy, iz, 128)
print("abc.shape"); print(abc.shape)

y_fake = ones(abc.shape[0], dtype=int)   ### original images - 15x128x128x128 sample

print("Train: ")
d_model.train_on_batch(abc, y_fake)					## no_of_patch is the sample size
print("Fake Train ended\n\n\n\n\n")










def define_gan_2(g_model, d_model):
    d_model.trainable = False
    gen_raw, gen_TE, gen_TR = g_model.input				# get noise and label inputs from generator model
    print("gen_raw.shape:"); print(gen_raw.shape)
    gen_output = g_model.output							# get image output from the generator model
    print("gen_output.shape:"); print(gen_output.shape)
    
    gan_output = d_model(gen_output)					# connect image output and label input from generator as inputs to discriminator
    model = Model([gen_raw, gen_TE, gen_TR], gan_output)# define gan model as taking noise and label and outputting a classification
    model.compile(loss='binary_crossentropy', optimizer=Adam(lr=0.0002, beta_1=0.5))  	# compile model
    plot_model(model, to_file='gan_whole_test.pdf', show_shapes=True, show_layer_names=True)
    return model



gan_model = define_gan_2(g_model, d_model)		# create the gan







X_real, ix, iy, iz = generate_real_samples(data_orig, 0, 128)
y_real = ones(X_real.shape[0], dtype=int)   ### original images - 15x128x128x128 sample

# update discriminator model weights
d_loss1, _ = d_model.train_on_batch(X_real, y_real)					## no_of_patch is the sample size

# generate 'fake' examples
X_fake = generate_fake_samples_2(g_model, data_array, 0, ix, iy, iz, 128)
y_fake = ones(X_fake.shape[0], dtype=int)   ### original images - 15x128x128x128 sample

# update discriminator model weights
d_loss2, _ = d_model.train_on_batch(X_fake, y_fake)					## no_of_patch is the sample size


# create inverted labels for the fake samples
y_gan = ones(X_fake.shape[0], dtype=int)


# generate 'fake' examples
X_inputs = generate_inputs(g_model, data_array, 0, ix, iy, iz, 128)

print("X_inputs"); print(X_inputs.shape)
print("TE_seq_all"); print(TE_seq_all.shape)
print("y_gan"); print(y_gan.shape)

# update the generator via the discriminator's error
g_loss, _ = gan_model.train_on_batch([X_inputs, TE_seq_all, TR_seq_all], y_gan)






