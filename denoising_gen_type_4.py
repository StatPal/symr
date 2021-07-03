## Used the denoise example and simplifying the corresponding settings
## Number 4, trying to do without the whole PIL business

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
from skimage.metrics import peak_signal_noise_ratio
from utils.denoising_utils import *

import time

torch.backends.cudnn.enabled = True
torch.backends.cudnn.benchmark =True
dtype = torch.cuda.FloatTensor

imsize =-1
PLOT = True #--original True
sigma = 25
sigma_ = sigma/255.


############## LOAD IMAGE: ############################
import pandas as pd
f = 'result/r_LS_brainweb_next_stage_1_pred.csv'
dat_iter = np.asarray(pd.read_csv(f, header=None))  ## np.asarray gives different format
dat_2 = dat_iter.reshape(181, 217, 181, 9)



target = 0

print("\n\ntarget image/settings no:", target)

dat_target = dat_2[:,:,:,target]
dat_target_refined = dat_target

slice_no = 91

print("slice_no:", slice_no)
slice_2D = dat_target[:,:,slice_no]
print("slice_2D.shape", slice_2D.shape)
scaling_factor = np.amax(slice_2D)


print("Scaling factor:", scaling_factor)
slice_2D = slice_2D*255/scaling_factor   ## Check this

ar = np.clip(slice_2D,0,255).astype(np.uint8)
if slice_2D.shape[0] == 1:
	ar = ar[0]
else:
	#ar = ar.transpose(1, 2, 0)
	ar = ar.transpose(0, 1)
act_image = ar
#act_image = Image.fromarray(ar)

## Pad images:
new_shape = (act_image.shape[0] - act_image.shape[0] % 32, act_image.shape[1] - act_image.shape[1] % 32)
if act_image.shape[0] % 32 != 0:
	tmp_1 = new_shape[0]+32
if act_image.shape[1] % 32 != 0:
	tmp_2 = new_shape[1]+32
new_shape = (tmp_1, tmp_2)

img_noisy_pil = np.zeros(new_shape)
img_noisy_pil[0:181,0:217] = act_image

#img_noisy_pil = crop_image(act_image, d=32)			## This changes the dimension - don't crop - pad
img_noisy_np = pil_to_np(img_noisy_pil)

# As we don't have ground truth
img_pil = img_noisy_pil
img_np = img_noisy_np

if PLOT:
	plot_image_grid([img_np], 4, 5, name="test/init_img_type_4.pdf");


############### SETUP: #######################


INPUT = 'noise' # 'meshgrid'
pad = 'reflection'
OPT_OVER = 'net' # 'net,input'

reg_noise_std = 1./30. # set to 1./20. for sigma=50
reg_noise_std = 0. 	   #  Subrata try
LR = 0.01

OPTIMIZER='adam' # 'LBFGS' -original 'adam'
show_every = 100  ## - original? 
exp_weight=0.99

num_iter = 240		# 2400
input_depth = 32
figsize = 4 


net = skip(input_depth, 1, 
				num_channels_down = [128, 128, 128, 128, 128],
                num_channels_up =   [128, 128, 128, 128, 128],
                num_channels_skip = [4, 4, 4, 4, 4], 
                upsample_mode='bilinear', downsample_mode='stride',
                need_sigmoid=True, need_bias=True, pad=pad, act_fun='LeakyReLU')

#~7 min






net = net.type(dtype)
net_input = get_noise(input_depth, INPUT, (img_pil.shape[0], img_pil.shape[1])).type(dtype).detach()

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
	print ('Iteration %05d    Loss %f   PSNR_noisy: %f   PSRN_gt: %f PSNR_gt_sm: %f' 
			% (i, total_loss.item(), psrn_noisy, psrn_gt, psrn_gt_sm), '\r', end='')
	if  PLOT and i % show_every == 0:
	    out_np = torch_to_np(out)
	    plot_image_grid([np.clip(out_np, 0, 1), 
	                     np.clip(torch_to_np(out_avg), 0, 1)], factor=figsize, nrow=1, name="test/steps"+str(i)+"_type_4.pdf")
	
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
q = plot_image_grid([np.clip(out_np, 0, 1), img_np], factor=13, name="test/err_1_"+str(target)+"_last_type_4.pdf");
print("out_np.shape:", (out_np).shape)

df = pd.DataFrame(out_np[0,0:181,0:217])

out_np = out_np*scaling_factor/255.0
dat_target[:,:,slice_no] = out_np[0,0:181,0:217]
df.to_csv('test/final_'+str(slice_no)+'_type_4.csv', header=None, index=None)



