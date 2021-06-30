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


#from keras import backend as K
#K.set_session(K.tf.Session(config=K.tf.ConfigProto(intra_op_parallelism_threads=15 , inter_op_parallelism_threads=15)))


import tensorflow as tf
#config = tf.ConfigProto(intra_op_parallelism_threads=15, 
#                        inter_op_parallelism_threads=15, 
#                        allow_soft_placement=True,
#                        device_count = {'CPU': 15})
#session=tf.Session(config=config) 

tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(15)

run_opts = tf.compat.v1.RunOptions(report_tensor_allocations_upon_oom = True)


from keras.utils.vis_utils import plot_model
# or use from keras_visualizer import visualizer 
# visualizer(model, format='png', view=True)
## from keras.utils import model_to_dot
from random import randint


from numpy import (zeros, ones)


import nibabel as nib

img = nib.load("../data/brainweb_all_6.nii")

a = np.array(img.dataobj)




import glob as glob
from numpy import genfromtxt
import csv

from pandas import read_csv


# create list of vegetation files to be opened
# VegList = sorted(glob.glob('../brainweb_new/result/generated_r_5_*.csv'))
VegList = sorted(glob.glob('../brainweb_new/result/generated_r_5_1[1234].csv'))
data = []

#f = VegList[0]
#with open(f, 'r') as file1:
#    dat_iter = csv.reader(file1, skipinitialspace=True)
#    data = [data for data in dat_iter]
##    for row in dat_iter:
##        print(row)

#df = read_csv(f, header=None)
#df[:10]
#df2 = np.asarray(df)
#df3 = df2.reshape(181,217,181,12)


data = []
for f in VegList:
    print(f)
    # reader = csv.reader(f)
    # dat = genfromtxt(f, delimiter=',')
    #with open(f, 'r') as file:
    #    dat_iter = csv.reader(file)
    #all_dat = [all_dat, dat]
    dat_iter = np.asarray(read_csv(f, header=None))  ## np.asarray gives different format
    data.append(dat_iter.reshape(181, 217, 181, 12))
    # VegListCoords = np.vstack((Veg.x, Veg.y, Veg.z)).transpose()

# Maybe add pd.concat if you wish
# data


data_array = np.asarray(data)
# print(data_array)
print(data_array.shape)






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


def build_generator(img_shape_8):
    print("\n\n\nGenerator:")
    print("img_shape_8")
    print(img_shape_8)
    
    visible = Input(img_shape_8, name="img")
    TE_seq_1 = Input(img_shape_8, name="TE")
    TR_seq_1 = Input(img_shape_8, name="TR")
    
    vis_final = concatenate([visible, TE_seq_1, TR_seq_1])
    print(vis_final.shape)
    print("Check 1")
    
    vis_final = Reshape((128,128,128,settings_no*3))(vis_final)								## Check the formation
    # vis_final = Concatenate(axis=-1)([visible, TE_seq_1, TR_seq_1])
    
    print(vis_final.shape)
    print("\n\n\n\n\n")
    # vis_final = visible
    print(vis_final)

    x1 = Dense(units=128, activation="sigmoid")(vis_final)
    x1 = Dense(units=1, activation="sigmoid")(x1)

    x2 = Dense(units=128, activation="sigmoid")(vis_final)
    x2 = Dense(units=1, activation="sigmoid")(x2)

    joined_representation = Multiply()([x1, x2])
    final = Dense(units=1, activation="elu")(joined_representation)

    model = Model(inputs=[visible, TE_seq_1, TR_seq_1], outputs=final, name="try_gen")  ## No, this would have same dim as ??

    plot_model(model, to_file='gan_gen.pdf', show_shapes=True, show_layer_names=True)
    # print(model.summary())				# Added by Subrata
    return model



build_generator(img_shape_8)







## cGAN Discriminator:

def build_discriminator(img_shape):
    print("\n\n\nDiscriminator:")
    visible = Input(img_shape)
    
    print("\nImage shape is:")
    print(visible.shape)

    
    u1 = Conv3D(filters=16, kernel_size=3, activation="elu",padding='same',name="conv_11")(visible)
    u1 = BatchNormalization(name="BN1")(u1)
    u1 = Conv3D(filters=16, kernel_size=3, activation="elu",padding='same',name="conv_12")(u1)
    u1 = BatchNormalization(name="u1")(u1)
    
    v1 = Conv3D(filters=32, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_1")(u1)
    v1 = BatchNormalization(name="v1")(v1)
    y21 = u1
    
        
    u2 = Conv3D(filters=32, kernel_size=3, activation="elu",padding='same',name="conv_21")(v1)
    u2 = BatchNormalization(name="BN2")(u2)
    u2 = Conv3D(filters=32, kernel_size=3, activation="elu",padding='same',name="conv_22")(u2)
    u2 = BatchNormalization(name="u2")(u2)
    
    v2 = Conv3D(filters=64, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_2")(u2)
    v2 = BatchNormalization(name="v2")(v2)
    y22 = UpSampling3D(size = 2)(u2)
    
    
    u3 = Conv3D(filters=64, kernel_size=3, activation="elu",padding='same',name="conv_31")(v2)
    u3 = BatchNormalization(name="BN3")(u3)
    u3 = Conv3D(filters=64, kernel_size=3, activation="elu",padding='same',name="conv_32")(u3)
    u3 = BatchNormalization(name="u3")(u3)
    
    v3 = Conv3D(filters=128, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_3")(u3)
    v3 = BatchNormalization(name="v3")(v3)
    y23 = UpSampling3D(size = 4)(u3)
    
    
    u4 = Conv3D(filters=128, kernel_size=3, activation="elu",padding='same',name="conv_41")(v3)
    u4 = BatchNormalization(name="BN4")(u4)
    u4 = Conv3D(filters=128, kernel_size=3, activation="elu",padding='same',name="conv_42")(u4)
    u4 = BatchNormalization(name="u4")(u4)
    
    v4 = Conv3D(filters=256, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_4")(u4)
    v4 = BatchNormalization(name="v4")(v4)
    y24 = UpSampling3D(size = 8)(u4)
    
    
    u5 = Conv3D(filters=256, kernel_size=3, activation="elu",padding='same',name="conv_51")(v4)
    u5 = BatchNormalization(name="BN5")(u5)
    u5 = Conv3D(filters=256, kernel_size=3, activation="elu",padding='same',name="conv_52")(u5)
    u5 = BatchNormalization(name="u5")(u5)
    
    v5 = Conv3D(filters=512, kernel_size=3, strides=2, activation="elu",padding='same',name="conv_5")(u5)
    v5 = BatchNormalization(name="v5")(v5)
    y25 = UpSampling3D(size = 16)(u5)
    y26 = UpSampling3D(size = 32)(v5)
    
    x1 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x1")(v5)
    
    z2 = Concatenate(axis=-1)([y21, y22, y23, y24, y25, y26]) # Concatenates images with their label embeddings
    x2 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x2")(z2)
    
    
    y32 = MaxPool3D(pool_size = 2,name="y32")(y22)
    y33 = MaxPool3D(pool_size = 2,name="y33")(y23)
    y34 = MaxPool3D(pool_size = 2,name="y34")(y24)
    y35 = MaxPool3D(pool_size = 2,name="y35")(y25)
    y36 = MaxPool3D(pool_size = 2,name="y36")(y26)
        
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
    
    print("y65: "); print(y65.shape)
    print("y66: "); print(y66.shape)
    
    z6 = Concatenate(axis=-1)([y65, y66])
    print("z6: "); print(z6.shape)
    x6 = Conv3D(filters=1, kernel_size=3, strides=1, padding='same',name="x6")(z6)
    
    
    
    # model = Model(inputs=visible, outputs=[x1, x2, x3, x4, x5, x6])
    ### Problem is in this line. 
    
    
    ## My fix: faltten and then connect using - https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8662660
    
#    x1_flat = Flatten()(x1)
#    x2_flat = Flatten()(x2)
#    x3_flat = Flatten()(x3)
#    x4_flat = Flatten()(x4)
#    x5_flat = Flatten()(x5)
#    x6_flat = Flatten()(x6)
#    
#    
#    final_1 = concatenate([x1_flat, x2_flat, x3_flat, x4_flat, x5_flat, x6_flat])
#    final_2 = Dense(units=1, activation="sigmoid")(final_1)
#    
#    model = Model(inputs=visible, outputs=final_2)


    x1_flat = Flatten()(x1)
    x2_flat = Flatten()(x2)
    x3_flat = Flatten()(x3)
    x4_flat = Flatten()(x4)
    x5_flat = Flatten()(x5)
    x6_flat = Flatten()(x6)

    
    
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
    
    plot_model(model, to_file='gan_dis.pdf', show_shapes=True, show_layer_names=True) # , expand_nested=True
    # model_to_dot(model, show_shapes=True, show_layer_names=True, expand_nested=True)
    #print(model.summary())
    return model



build_discriminator(img_shape)





def define_gan(g_model, d_model):
    # make weights in the discriminator not trainable
    d_model.trainable = False
    # get noise and label inputs from generator model
    gen_raw, gen_TE, gen_TR = g_model.input
    # get image output from the generator model
    gen_output = g_model.output
    print(gen_output.shape)  ## (None, 1, 1, 1) -- not matching with input of d_model in the next line
    # wants img_shape (128, 128, 1)
    
    # connect image output and label input from generator as inputs to discriminator
    gan_output = d_model(gen_output)
    # define gan model as taking noise and label and outputting a classification
    model = Model([gen_raw, gen_TE, gen_TR], gan_output)
    # compile model
    opt = Adam(lr=0.0002, beta_1=0.5)
    model.compile(loss='binary_crossentropy', optimizer=opt)
    plot_model(model, to_file='gan_all.pdf', show_shapes=True, show_layer_names=True)
    return model




# create the discriminator
d_model = build_discriminator(img_shape)
# create the generator
g_model = build_generator(img_shape_8)

d_model

g_model

# create the gan
gan_model = define_gan(g_model, d_model)
# load image data



# select real samples
def generate_real_samples(dataset, n_batch):
    # choose random instances
    ix = randint(0, dataset.shape[1]-128); iy = randint(0, dataset.shape[2]-128); iz = randint(0, dataset.shape[3]-128)
    # select images and labels
    X = dataset[:,ix:(ix+128), iy:(iy+128), iz:(iz+128),:]
    return X



X_real = generate_real_samples(data_array, 128)
print("X_real[0].shape")
print(X_real[0].shape)
print("X_real[1:3,:,:,:,0].shape: ")
print(X_real[1:3,:,:,:,0].shape)

# update discriminator model weights
y_real = ones(((X_real[1:3,:,:,:,0]).shape[0], 1))
print("y_real shape: ")
print(y_real.shape)

d_model = build_discriminator(img_shape)


print("Train: ")
d_model.train_on_batch(X_real[1:3,:,:,:,0], y_real)


print("Check 2")

# generate 'fake' examples
[X_fake, labels], y_fake = generate_fake_samples(g_model, latent_dim, half_batch)
print("Check 3")

# update discriminator model weights
d_loss2, _ = d_model.train_on_batch([X_fake, labels], y_fake)
print("Check 4")

# prepare points in latent space as input for the generator
[z_input, labels_input] = generate_latent_points(latent_dim, n_batch)
print("Check 5")

# create inverted labels for the fake samples
y_gan = ones((n_batch, 1))
print("Check 6")

# update the generator via the discriminator's error
g_loss = gan_model.train_on_batch([z_input, labels_input], y_gan)









def train1(g_model, d_model, gan_model, dataset, n_batch=128, n_iter=10):
    # manually enumerate epochs
    for i in range(n_iter):
        # enumerate batches over the training set
        for j in range(n_batch):
            # get a 128x128x128 patch from the main data
            # No, only take the flair image from the data. 
            X_real = generate_real_samples(dataset, n_batch)
            print("X_real.shape: ")
            print(X_real.shape)
            
            # update discriminator model weights
            y_real = ones((X_real.shape[0], 1))
            print("y_real shape: ")
            print(y_real.shape)
            d_loss1, _ = d_model.train_on_batch(X_real[:,:,:,:,0], y_real)
            # discriminator model should get the pd image - check which one is the rho density
            print("Check 2")
            
            # generate 'fake' examples
            [X_fake, labels], y_fake = generate_fake_samples(g_model, latent_dim, half_batch)
            print("Check 3")
            
            # update discriminator model weights
            d_loss2, _ = d_model.train_on_batch([X_fake, labels], y_fake)
            print("Check 4")
            
            # prepare points in latent space as input for the generator
            [z_input, labels_input] = generate_latent_points(latent_dim, n_batch)
            print("Check 5")
            
            # create inverted labels for the fake samples
            y_gan = ones((n_batch, 1))
            print("Check 6")
            
            # update the generator via the discriminator's error
            g_loss = gan_model.train_on_batch([z_input, labels_input], y_gan)





## size of the latent space
# latent_dim = 100
## create the discriminator
#d_model = build_discriminator(img_shape)
## create the generator
#g_model = build_generator(pix_shape)
## create the gan
#gan_model = define_gan(g_model, d_model)
## load image data
# dataset = load_real_samples()
## train model
# train1(g_model, d_model, gan_model, data_array, 128, 10)





