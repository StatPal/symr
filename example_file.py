import os
os.chdir("./src")

import sympy as sim

import numpy as np





## Sample row of parameters, W
W_row = np.array([50, 0.01, 0.003])
W = np.array([[50, 0.01, 0.003], [36, 0.02, 0.04]])	## Two sample rows

## Design parameters
TE = np.array([0.01, 0.03, 0.04, 0.01])
TR = np.array([0.6, 0.6, 1, 0.8])
sig = np.array([1.2, 1.1, 1.4, 1.2])

## Forward transformed values: 
sim.Bloch_eqn(W_row, TE, TR)
sim.v_mat(W, TE, TR)
sim.Generate_r(W, TE, TR, sig)

sim.dee_v_ij_dee_W_ik(W_row, TE, TR, 1, 1)					## dv_i1/dW_i1
sim.dee_2_v_ij_dee_W_ik_dee_W_ik1(W_row, TE, TR, 1, 1, 2)	## d^2v_i1/(dW_i1 dW_i2)
