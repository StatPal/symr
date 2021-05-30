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
from datetime import datetime



import glob as glob
from pandas import read_csv

VegList = sorted(glob.glob('../brainweb_new/result/generated_r_5_[123]?.csv'))
data = []
for f in VegList:
    print(f)
    dat_iter = np.asarray(read_csv(f, header=None))  ## np.asarray gives different format
    data.append(dat_iter.reshape(181, 217, 181, 12))

data_array = np.asarray(data)


idx = [0,8,9]
data_orig = data_array[:,:,:,:,1]			## Suppose this is rho-weighted settings image
data_array = data_array[:,:,:,:,idx]		## Suppose these are all other train settings image
print("data_array.shape"); print(data_array.shape)

## data contains data with 2 people, 181 n_x, 217 n_y, 181 n_z, 6 TE/TR settings






## Model input dimensions:
img_rows = 64; img_cols = 64; img_depth = 64
channels = 1
settings_no = data_array.shape[4]
no_of_sample = data_orig.shape[0]
no_of_patch = 3  		## should be 15 according to the paper

img_shape = (img_rows, img_cols, img_depth, channels)
img_shape_8 = (img_rows, img_cols, img_depth, settings_no)
pix_shape = (1, 1, 1, settings_no)

TE_seq = np.array([0.01, 0.08, 0.01])
TR_seq = np.array([0.6, 2.0, 3.0])





## cGAN Discriminator:

def build_discriminator(img_shape_dis):
    print("\n\n\nDiscriminator:")
    visible = Input(img_shape_dis)
    
    u2 = Conv3D(filters=16, kernel_size=3, activation="elu",padding='same',name="conv_21")(visible)
    u2 = BatchNormalization(name="BN2")(u2)
    u2 = Conv3D(filters=16, kernel_size=3, activation="elu",padding='same',name="conv_22")(u2)
    u2 = BatchNormalization(name="u2")(u2)
    
    v2 = Conv3D(filters=32, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_2")(u2)
    v2 = BatchNormalization(name="v2")(v2)
    y32 = u2
    
    u3 = Conv3D(filters=32, kernel_size=3, activation="elu",padding='same',name="conv_31")(v2)
    u3 = BatchNormalization(name="BN3")(u3)
    u3 = Conv3D(filters=32, kernel_size=3, activation="elu",padding='same',name="conv_32")(u3)
    u3 = BatchNormalization(name="u3")(u3)
    
    v3 = Conv3D(filters=64, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_3")(u3)
    v3 = BatchNormalization(name="v3")(v3)
    y33 = UpSampling3D(size = 2)(u3)
    
    u4 = Conv3D(filters=64, kernel_size=3, activation="elu",padding='same',name="conv_41")(v3)
    u4 = BatchNormalization(name="BN4")(u4)
    u4 = Conv3D(filters=64, kernel_size=3, activation="elu",padding='same',name="conv_42")(u4)
    u4 = BatchNormalization(name="u4")(u4)
    
    v4 = Conv3D(filters=128, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_4")(u4)
    v4 = BatchNormalization(name="v4")(v4)
    y34 = UpSampling3D(size = 4)(u4)
    
    u5 = Conv3D(filters=128, kernel_size=3, activation="elu",padding='same',name="conv_51")(v4)
    u5 = BatchNormalization(name="BN5")(u5)
    u5 = Conv3D(filters=128, kernel_size=3, activation="elu",padding='same',name="conv_52")(u5)
    u5 = BatchNormalization(name="u5")(u5)
    
    v5 = Conv3D(filters=256, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_5")(u5)
    v5 = BatchNormalization(name="v5")(v5)
    y35 = UpSampling3D(size = 8)(u5)
    y36 = UpSampling3D(size = 16)(v5)
    
    x1 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x1")(v5)
    
        
    z3 = Concatenate(axis=-1)([y32, y33, y34, y35, y36])
    x3 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x3")(z3)
    
    y43 = MaxPool3D(pool_size = 2,name="y43")(y33)
    y44 = MaxPool3D(pool_size = 2,name="y44")(y34)
    y45 = MaxPool3D(pool_size = 2,name="y45")(y35)
    y46 = MaxPool3D(pool_size = 2,name="y46")(y36)

    z4 = Concatenate(axis=-1)([y43, y44, y45, y46])
    x4 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x4")(z4)
    
    
    y54 = MaxPool3D(pool_size = 2,name="y54")(y44)
    y55 = MaxPool3D(pool_size = 2,name="y55")(y45)
    y56 = MaxPool3D(pool_size = 2,name="y56")(y46)
    z5 = Concatenate(axis=-1)([y54, y55, y56])
    x5 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x5")(z5)
    
    
    y65 = MaxPool3D(pool_size = 2,name="y65")(y55)
    y66 = MaxPool3D(pool_size = 2,name="y66")(y56)
        
    z6 = Concatenate(axis=-1)([y65, y66])
    x6 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x6")(z6)
    
    x1_flat = Flatten()(x1)
    x3_flat = Flatten()(x3)
    x4_flat = Flatten()(x4)
    x5_flat = Flatten()(x5)
    x6_flat = Flatten()(x6)

    x1_flat = Dense(units=1, activation="sigmoid")(x1_flat)
    x3_flat = Dense(units=1, activation="sigmoid")(x3_flat)
    x4_flat = Dense(units=1, activation="sigmoid")(x4_flat)
    x5_flat = Dense(units=1, activation="sigmoid")(x5_flat)
    x6_flat = Dense(units=1, activation="sigmoid")(x6_flat)
    
    final_1 = concatenate([x1_flat, x3_flat, x4_flat, x5_flat, x6_flat])
    final_2 = Dense(units=1, activation="sigmoid")(final_1)
    
    model = Model(inputs=visible, outputs=final_2)
    
    # compile model
    opt = Adam(lr=0.000001)
    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
    
    plot_model(model, to_file='gan_dis.pdf', show_shapes=True, show_layer_names=True) # , expand_nested=True
    #print(model.summary())
    return model



d_model = build_discriminator(img_shape)



"""
Generate real patch with corresponding indices
i-th person data
"""
def generate_real_samples(dataset, subject_no, n_batch):
    
    ix = randint(0, dataset.shape[1]-64); iy = randint(0, dataset.shape[2]-64); iz = randint(0, dataset.shape[3]-64)
    X_real = dataset[subject_no, ix:(ix+64), iy:(iy+64), iz:(iz+64)]	## (no_sample, X, Y, Z)
    X_real = np.expand_dims(X_real, axis=0)
    ix = np.array(ix); iy = np.array(iy); iz = np.array(iz)
    
    
    for i in range(1,no_of_patch):
        ix_tmp = randint(0, dataset.shape[1]-64); iy_tmp = randint(0, dataset.shape[2]-64); iz_tmp = randint(0, dataset.shape[3]-64)
        X_real_tmp = dataset[subject_no, ix_tmp:(ix_tmp+64), iy_tmp:(iy_tmp+64), iz_tmp:(iz_tmp+64)]	## (no_sample, X, Y, Z)
        X_real_tmp = np.expand_dims(X_real_tmp, axis=0)
    
        X_real = np.concatenate([X_real, X_real_tmp], axis=0)
        ix = np.append(ix, ix_tmp); iy = np.append(iy, iy_tmp); iz = np.append(iz, iz_tmp);

    return X_real, np.array(ix), np.array(iy), np.array(iz)


X_real, ix, iy, iz = generate_real_samples(data_orig, 0, 64)
y_real = ones(X_real.shape[0], dtype=int)   ### original images - 15x64x64x64 sample

print("X_real.shape: ", X_real.shape)

d_model.train_on_batch(X_real, y_real)					## no_of_patch is the sample size
print("Real Train ended\n\n\n\n\n")









### https://stackoverflow.com/questions/43915482/how-do-you-create-a-custom-activation-function-with-keras
### See also, https://stackoverflow.com/questions/57023350/implementing-the-square-non-linearity-sqnl-activation-function-in-keras
### and this: https://stackoverflow.com/questions/49982438/how-to-restrict-the-output-of-neural-network-to-be-positive-in-python-keras

## Custom activation function
#from keras.layers import Activation
#from keras import backend as K
#from keras.utils.generic_utils import get_custom_objects


#def custom_activation(x):
#    return (x * x)

#get_custom_objects().update({'sq_activation': Activation(custom_activation)})



## cGAN generator: 


def build_generator(img_shape_gen):
    print("\n\n\nGenerator:")
    
    visible = Input(img_shape_gen, name="img")
    TE_seq_1 = Input(img_shape_gen, name="TE")
    TR_seq_1 = Input(img_shape_gen, name="TR")
    
    vis_final = concatenate([visible, TE_seq_1, TR_seq_1])
    
    # vis_final = Reshape((64,64,64,settings_no*3))(vis_final)								## Check the formation
    # vis_final = Concatenate(axis=-1)([visible, TE_seq_1, TR_seq_1])
    
    
    print("vis_final.shape"); print(vis_final.shape)
    
    x1 = Dense(units=64, activation="sigmoid")(vis_final)
    x1 = Dense(units=1, activation="sigmoid")(x1)

    x2 = Dense(units=64, activation="sigmoid")(vis_final)
    x2 = Dense(units=1, activation="sigmoid")(x2)

    joined_representation = Multiply()([x1, x2])
    final = Dense(units=1, activation="elu")(joined_representation)

    model = Model(inputs=[visible, TE_seq_1, TR_seq_1], outputs=final, name="try_gen")

    plot_model(model, to_file='gan_gen.pdf', show_shapes=True, show_layer_names=True)
    # print(model.summary())
    return model



g_model = build_generator(img_shape_8)




"""
Generated fake sample: 
input: 	whole data: no_of_sample  x  n_x  x  n_y  x  n_z  x  settings_no
		subject_no: the number of subject over which it is taken.
		ix, iy, iz: numpy arrays, each of length no_of_patch
		
inter:	tmp_X_all: (64x64x64) x settings_no

output: fake_data:  64 x 64 x 64 
"""

TE_seq_expand = np.expand_dims(TE_seq, axis=0)
TR_seq_expand = np.expand_dims(TR_seq, axis=0)
TE_seq_all = np.tile(TE_seq_expand, (no_of_patch, 64, 64, 64, 1))
TR_seq_all = np.tile(TR_seq_expand, (no_of_patch, 64, 64, 64, 1))



def generate_inputs(generator, dataset, subject_no, ix, iy, iz, n_batch):
    for i in range(len(ix)):
        ix_tmp = ix[i]; iy_tmp = iy[i]; iz_tmp = iz[i];
        X_fake_tmp = dataset[subject_no, ix_tmp:(ix_tmp+64), iy_tmp:(iy_tmp+64), iz_tmp:(iz_tmp+64),:]  
        X_fake_tmp = np.expand_dims(X_fake_tmp, axis=0)
        
        if i==0:
            tmp_X_all = X_fake_tmp
        else:
            tmp_X_all = np.concatenate([tmp_X_all, X_fake_tmp], axis=0)
    
    return tmp_X_all


def generate_fake_samples_2(generator, dataset, subject_no, ix, iy, iz, n_batch):
    tmp_X_all = generate_inputs(generator, dataset, subject_no, ix, iy, iz, n_batch)
    fake_data = generator.predict([tmp_X_all, TE_seq_all, TR_seq_all])
    return fake_data[:,:,:,:,0]




X_fake = generate_fake_samples_2(g_model, data_array, 0, ix, iy, iz, 64)
y_fake = ones(X_fake.shape[0], dtype=int)   ### original images - 15x64x64x64 sample

print("X_fake.shape: ", X_fake.shape)

d_model.train_on_batch(X_fake, y_fake)					## no_of_patch is the sample size
print("Fake Train ended\n\n\n\n\n")








def define_gan_2(g_model, d_model):
    d_model.trainable = False
    gen_raw, gen_TE, gen_TR = g_model.input				# get noise and label inputs from generator model
    gen_output = g_model.output							# get image output from the generator model
    
    gan_output = d_model(gen_output)					# connect image output and label input from generator as inputs to discriminator
    model = Model([gen_raw, gen_TE, gen_TR], gan_output)# define gan model as taking noise and label and outputting a classification
    model.compile(loss='binary_crossentropy', optimizer=Adam(lr=0.0005))  	# compile model
    plot_model(model, to_file='gan_whole_test.pdf', show_shapes=True, show_layer_names=True)
    return model


gan_model = define_gan_2(g_model, d_model)		# create the gan

print("GAN model defined")




# train the generator and discriminator
def train(g_model, d_model, gan_model, dataset_total, dataset_orig, n_epochs=100, n_batch=64):
	for i in range(n_epochs):
		# enumerate batches over the training set
		for j in range(no_of_sample):
			X_real, ix, iy, iz = generate_real_samples(dataset_orig, j, 64)
			y_real = ones(X_real.shape[0], dtype=int)   ### original images - 15x64x64x64 sample
			# update discriminator model weights
			d_loss1, _ = d_model.train_on_batch(X_real, y_real)					## no_of_patch is the sample size
			# generate 'fake' examples
			X_fake = generate_fake_samples_2(g_model, dataset_total, j, ix, iy, iz, 64)
			y_fake = ones(X_fake.shape[0], dtype=int)   ### original images - 15x64x64x64 sample
			# update discriminator model weights
			d_loss2, _ = d_model.train_on_batch(X_fake, y_fake)					## no_of_patch is the sample size

			# create inverted labels for the fake samples
			y_gan = ones(X_fake.shape[0], dtype=int)

			# generate 'fake' examples
			X_inputs = generate_inputs(g_model, dataset_total, j, ix, iy, iz, 64)

			# update the generator via the discriminator's error
			g_loss = gan_model.train_on_batch([X_inputs, TE_seq_all, TR_seq_all], y_gan)

			# summarize loss on this batch
			print('>%d, %d, d1=%.3f, d2=%.3f g=%.3f' %
				(i+1, j+1, d_loss1, d_loss2, g_loss)); print("now=", datetime.now())
	# save the generator model
	g_model.save('cgan_generator_5_1.h5')



# train model
iter=100000/8
epoch_no=int(iter/no_of_sample)
train(g_model, d_model, gan_model, data_array, data_orig, epoch_no, 5)




