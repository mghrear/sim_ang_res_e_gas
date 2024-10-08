## This script computes the dN/dx for each of the selected energies / gases

import sys
sys.path.append('../')


from ROOT import TFile, TVector3, TMath
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from math import sqrt, fabs
import SimTools
from matplotlib.colors import LogNorm


###############################################################
###############################################################

################   Variables to be specified  ################

###############################################################
###############################################################

#Select the gas mixture
#Options: 'he_co2', 'he_cf4'
gases = ['he_co2', 'he_cf4']

# Energies in keV
Energies = [30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]


# distances on which to compute dN/dx
dists = [ 0.1 , 0.2 , 0.3 , 0.4 , 0.5 ]


###############################################################
###############################################################


df_results = pd.DataFrame(columns = ['gas', 'energy', 'dNdx', 'dNdx_err', 'energy_loss'])


for gas in gases:
	for Energy in Energies:

		# There are two possible W-values
		if gas == 'he_co2':
			W = 0.0352
		else:
			W = 0.0304

		# List of electron file(s)
		file_path = '~/data/e_ang_res/MS_test_sims/'+gas+'/'+str(Energy)+'_keV/'+gas+'_'+str(Energy)+'keV_' + str(0) + '/'+gas+'_'+str(Energy)+'keV_' + str(0) + '.root'

		#read tracks
		tracks,times = SimTools.read_degrad(file_path)

		N_vals = []
		x_vals = []

		for track,time in zip(tracks,times):

			#Get linear direction
			charge_counts = SimTools.get_N_vs_x(track,time,dists)

			N_vals += charge_counts
			x_vals += dists

		# Fit line to get dNdx
		popt, pcov = curve_fit(SimTools.slope_fit, x_vals, N_vals, sigma=np.sqrt(N_vals), absolute_sigma=True)
		perr = np.sqrt(np.diag(pcov))

		# Plot results
		plt.figure()
		x_bins = np.arange(0,np.max(N_vals),10)

		h = plt.hist2d(x_vals, N_vals, bins=[len(dists),40], norm=LogNorm())
		cbar1 = plt.colorbar(h[3])
		cbar1.set_label("Count") 

		xs = np.arange(0,max(x_vals),max(x_vals)/1000.0)
		plt.plot(xs,SimTools.slope_fit(xs,*popt),'r-',label='dN/dX =  = %5.3f electon/cm' % tuple([popt[0]]))
		plt.legend()

		plt.xlabel("Fit Length")
		plt.ylabel("No. Electrons")

		plt.savefig("./plots/dNdx_"+gas+"_"+str(Energy)+"keV.pdf")

		# Save data
		df_results = df_results.append({'gas' : gas, 'energy' : Energy, 'dNdx' : popt[0], 'dNdx_err' : perr[0], 'energy_loss' :   SimTools.slope_fit(max(dists),*popt)*W}, ignore_index = True)

df_results.to_pickle('./dNdx_data.pk')

