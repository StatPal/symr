############################ Import libs: ##########################
from __future__ import print_function
import warnings
warnings.filterwarnings('ignore')





import matplotlib.pyplot as plt
import os

import numpy as np
from models import *

import torch
import torch.optim

#from skimage.measure import compare_psnr		## changed to
from skimage.metrics import peak_signal_noise_ratio
from utils.denoising_utils import *

torch.backends.cudnn.enabled = True
torch.backends.cudnn.benchmark =True
dtype = torch.cuda.FloatTensor

dtype = torch.cuda.FloatTensor

imsize =-1
PLOT = True
sigma = 25
sigma_ = sigma/255.




############## LOAD IMAGE: ############################
import pandas as pd
f = '../LS_with_deep/result/r_LS_brainweb_next_stage_1_pred.csv'
dat_iter = np.asarray(pd.read_csv(f, header=None))  ## np.asarray gives different format
dat_2 = dat_iter.reshape(181, 217, 181, 9)


import PIL
from PIL import Image
from matplotlib import cm


target = 0
dat_target = dat_2[:,:,:,target]

scaling_factor = np.amax(dat_target)
print("Scaling factor:", scaling_factor)
dat_target = dat_target*255/scaling_factor   ## Check this
print("Scaling max value:", np.amax(dat_target), "\n\n\n")

print("Check, before np. cip. dat_target.shape", dat_target.shape)


ar = np.clip(dat_target,0,255).astype(np.uint8)

print("ar done")

if dat_target.shape[0] == 1:
	ar = ar[0]
else:
	#ar = ar.transpose(1, 2, 0)
	ar = ar.transpose(0, 1, 2)
#act_image = Image.fromarray(ar)

act_image = ar

print("Check, before img_noisy_pil .cip")
## img_noisy_pil = crop_image(act_image, d=32)   # Don't crop, pad



## Pad images:
new_shape = (act_image.shape[0] - act_image.shape[0] % 32, act_image.shape[1] - act_image.shape[1] % 32, 
			 act_image.shape[2] - act_image.shape[2] % 32)
if act_image.shape[0] % 32 != 0:
	tmp_1 = new_shape[0]+32
if act_image.shape[1] % 32 != 0:
	tmp_2 = new_shape[1]+32
if act_image.shape[2] % 32 != 0:
	tmp_3 = new_shape[2]+32
new_shape = (tmp_1, tmp_2, tmp_3)

img_noisy_pil = np.zeros(new_shape)
img_noisy_pil[0:181,0:217,0:181] = act_image


img_noisy_np = pil_to_np(img_noisy_pil)

# As we don't have ground truth
img_pil = img_noisy_pil
img_np = img_noisy_np

plot_image_grid([img_np], 4, 5, name="original_paper_images/first_snail.pdf");
        


print("Check, before setup")




############### SETUP: #######################

INPUT = 'noise' # 'meshgrid'
pad = 'reflection'					## Subrata, 'ReflectionPad3d' is joined just a few days ago - should be in models/common.py 
## https://github.com/pytorch/pytorch/issues/27655
## https://github.com/cyyever/pytorch/commit/86ee83211949d188ab93b85e6647920232a48b72
## https://github.com/pytorch/pytorch/pull/59791
## That means it's possibly in the nightly-build
## https://anaconda.org/pytorch-nightly/repo
pad = 'zero'


OPT_OVER = 'net' # 'net,input'

reg_noise_std = 1./30. # set to 1./20. for sigma=50
LR = 0.01

OPTIMIZER='adam' # 'LBFGS'
show_every = 200
exp_weight=0.99


num_iter = 2400
num_iter = 150  # Subrata
input_depth = 3
figsize = 5 

net = skip(
        input_depth, 3, 
        num_channels_down = [8, 16, 32, 64, 128], 
        num_channels_up   = [8, 16, 32, 64, 128],
        num_channels_skip = [0, 0, 0, 4, 4], 
        upsample_mode='bilinear',
        need_sigmoid=True, need_bias=True, pad=pad, act_fun='LeakyReLU')

net = net.type(dtype)



net_input = get_noise(input_depth, INPUT, (img_pil.size[1], img_pil.size[0])).type(dtype).detach()

# Compute number of parameters
s  = sum([np.prod(list(p.size())) for p in net.parameters()]);
print ('Number of params: %d' % s)

# Loss
mse = torch.nn.MSELoss().type(dtype)
img_noisy_torch = np_to_torch(img_noisy_np).type(dtype)







################# OPTIMIZE ######################


net_input_saved = net_input.detach().clone()
noise = net_input.detach().clone()
out_avg = None
last_net = None
psrn_noisy_last = 0

i = 0
def closure():
    
    global i, out_avg, psrn_noisy_last, last_net, net_input
    
    if reg_noise_std > 0:
        net_input = net_input_saved + (noise.normal_() * reg_noise_std)
    
    out = net(net_input)
    
    # Smoothing
    if out_avg is None:
        out_avg = out.detach()
    else:
        out_avg = out_avg * exp_weight + out.detach() * (1 - exp_weight)
            
    total_loss = mse(out, img_noisy_torch)
    total_loss.backward()
        
    
    psrn_noisy = peak_signal_noise_ratio(img_noisy_np, out.detach().cpu().numpy()[0]) 
    psrn_gt    = peak_signal_noise_ratio(img_np, out.detach().cpu().numpy()[0]) 
    psrn_gt_sm = peak_signal_noise_ratio(img_np, out_avg.detach().cpu().numpy()[0]) 
    
    # Note that we do not have GT for the "snail" example
    # So 'PSRN_gt', 'PSNR_gt_sm' make no sense
    print ('Iteration %05d    Loss %f   PSNR_noisy: %f   PSRN_gt: %f PSNR_gt_sm: %f' % (i, total_loss.item(), psrn_noisy, psrn_gt, psrn_gt_sm), '\r', end='')
    if  PLOT and i % show_every == 0:
        out_np = torch_to_np(out)
        plot_image_grid([np.clip(out_np, 0, 1), 
                         np.clip(torch_to_np(out_avg), 0, 1)], factor=figsize, nrow=1, name="original_paper_images/steps"+str(i)+".pdf")
        
        
    
    # Backtracking
    if i % show_every:
        if psrn_noisy - psrn_noisy_last < -5: 
            print('Falling back to previous checkpoint.')

            for new_param, net_param in zip(last_net, net.parameters()):
                net_param.data.copy_(new_param.cuda())

            return total_loss*0
        else:
            last_net = [x.detach().cpu() for x in net.parameters()]
            psrn_noisy_last = psrn_noisy
            
    i += 1

    return total_loss

p = get_params(OPT_OVER, net, net_input)
optimize(OPTIMIZER, p, closure, LR, num_iter)







###############################################


out_np = torch_to_np(net(net_input))
q = plot_image_grid([np.clip(out_np, 0, 1), img_np], factor=13, name="original_paper_images/last.pdf");
