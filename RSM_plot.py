import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
import matplotlib.colors as colors
from matplotlib import cm, ticker

def pol_to_cart(r,theta):
	x = r*np.cos(theta)
	y = r*np.sin(theta)
	return x, y


wavelength = .15406 #nm

#sns.set_context("poster")

fig,axs = plt.subplots(2,4)

names = [[0,25,51,76],[101,126,151,201]]
overlaps = [[0,25,50,75],[100,125,150,200]]

for i in range(2):
	for j in range(4):
		file = f'RSM_{names[i][j]}-frame.dat'


		df = pd.read_csv(file, sep='\t', header=None, names=['Qx','Qy','Qz','Int'], dtype='float')
		df['Q'] = np.sqrt(df['Qx']**2 + df['Qy']**2 + df['Qz']**2)

		df = df[(df['Qy']==0.) & (df['Q']> -1.)] # edit the RSM slice and remove the "direct" beam component

		max_int = df['Int'].max()
		max_ind = df['Int'].idxmax()

		# span theta
		theta = np.linspace(0, np.pi, 1000)

		# Ewald sphere limits for a Cu Kalpha1 source
		r1 = 4.*np.pi/wavelength #* np.sin(120./2. * np.pi/180.) # gives the cutoff for a given max 2theta
		r2 = 4.*np.pi/wavelength*np.cos(theta)*(theta < np.pi/2) -4.*np.pi/wavelength*np.cos(theta)*(theta >= np.pi/2)


		# compute x1 and x2
		x1, y1 = pol_to_cart(r1, theta)
		x2, y2 = pol_to_cart(r2, theta)


		axs[i,j].plot(x1, y1, linewidth=2, color='white')
		axs[i,j].plot(x2, y2, linewidth=2, color='white')



		plot = axs[i,j].tricontourf(df['Qx'], df['Qz'], df['Int']/max_int, vmax=1., cmap='nipy_spectral', levels=np.logspace(-6., 0., 13), norm=colors.LogNorm())


		cbar = plt.colorbar(plot,ax=axs[i,j])
		cbar.set_label('Normalised Intensity', rotation=270., labelpad=25.)

		axs[i,j].set_title(f'{overlaps[i][j]}')


		axs[i,j].set_xlabel(r'$Q_{\parallel}$ (nm$^{-1}$)')
		axs[i,j].set_ylabel(r'$Q_{\perp}$ (nm$^{-1}$)')
		axs[i,j].set_facecolor("black")
		#axs.set_xlim(-60.,60.)
		#axs.set_ylim(-60.,60.)
		#ax.set_ylim(29.,32.)

		max_Qx = df['Qx'][max_ind]
		max_Qy = df['Qy'][max_ind]
		max_Qz = df['Qz'][max_ind]


		print(f'Left:\t Max intensity of {max_int}\t Qx = {max_Qx} /nm\t Qy = {max_Qy} /nm\t Qz = {max_Qz} /nm')

		axs[i,j].set_aspect('equal')

plt.show()
