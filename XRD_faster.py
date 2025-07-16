#!/usr/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import time
#from deco import *
from multiprocessing import Pool

# get the current time in seconds since the epoch
start = time.time()

# Cromer-Mann coefficients for Ga3+ and O2- (a_i, b_i, c):
CM_Ga = [12.692, 6.69883, 6.06692, 1.0066, 2.81262, 0.22789, 6.36441, 14.4122, 1.53545]
CM_O = [3.7504, 2.84294, 1.54298, 1.62091, 16.5151, 6.59203, 0.319201, 43.3486, 0.242060]

# Constant
const = 16. * np.pi * np.pi

# Atomic scattering factor
def fat(Q, atom):
	CM = CM_Ga if atom == 1 else CM_O

	Qabs = np.linalg.norm(Q)
	ans = CM[8]

	for i in range(4):
		ans += CM[i]*np.exp(-CM[i+4]*Qabs*Qabs/const)

	return ans

# Atomic scattering factor Ga
def fGa(Q):
	return fat(Q, 1)

# Atomic scattering factor O
def fO(Q):
	return fat(Q, 0)

# X-ray amplitude as calculated by MD
def A_MD(Q, dfGa, dfO):
	ans = fGa(Q)*np.sum(np.exp(1.j*np.dot(dfGa,Q))) + fO(Q)*np.sum(np.exp(1.j*np.dot(dfO,Q)))

	return ans

def compute_intensity(args):
    """Compute intensity for a given Qx and Qz."""
    Qx, Qy, Qz, dfGa, dfO = args
    Q = [Qx, Qy, Qz]
    amp = A_MD(Q, dfGa, dfO)
    ans = np.abs(amp)**2

    return Qx * 10., Qy * 10., Qz * 10., ans

def QSpace(Qxs, Qys, Qzs, fileName, df):
    dfGa = df[df['Type'] == 1][['X', 'Y', 'Z']].to_numpy()
    dfO = df[df['Type'] == 2][['X', 'Y', 'Z']].to_numpy()

    # Prepare arguments for parallel computation
    tasks = [(Qx, Qy, Qz, dfGa, dfO) for Qx in Qxs for Qy in Qys for Qz in Qzs]
    Ntotal = len(tasks)

    with open(fileName, 'w') as file:
        with Pool() as pool:
            for n, result in enumerate(pool.map(compute_intensity, tasks), 1):
                Qx, Qy, Qz, intensity = result
                file.write(f'{Qx}\t{Qy}\t{Qz}\t{intensity}\n')
                print(f'Computed Q-points: {n}/{Ntotal} ({n/Ntotal*100:.2f} %)', end='\r')


if __name__ == '__main__':
    # Lê os dados da simulação
    for i in range(0,101,1):
        col_names = ['Id', 'Type', 'X', 'Y', 'Z']
        df = pd.read_csv(f'./{i}-frame.data', names=col_names, sep='\s+', skiprows=12, nrows=320000) 
        print(df)
        Qxs = np.linspace(-8.2, 8.2, num=51)
        Qys = np.linspace(-8.2, 8.2, num=51)
        Qzs = np.linspace(-8.2, 8.2, num=51)

        fileName= f'RSM_{i}-frame.dat'


        QSpace(Qxs, Qys, Qzs, fileName, df)
        # Print running time
        print(f'Elapsed time: {time.time() - start} s')
