# This script is used to process the raw degrad simulation files into angular distributions which are fit to a guassian to determine the angular resolution
# Results are stored in MS_fitting_data.pk

import sys
sys.path.append('../')
import os
import matplotlib.pyplot as plt
import numpy as np
import root_pandas as rp
import pandas as pd
from math import sqrt, fabs
import SimTools
from scipy.optimize import curve_fit

###############################################################
###############################################################

################   Variables to be specified  ################

###############################################################
###############################################################

#Select the gas
gases = ['cf4','co2','ch4','c2h6','ne','ar','xe']

# Energies to calculate
Energies = [50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]

# Select range in cm
dists = [0.5,1.0,2.0]

#Number of electron file to analyze
efiles = 100

# Percent of distribution to trim
Trim = 2.0


###############################################################
###############################################################


# dataframe to store all results
df_results = pd.DataFrame(columns = ['gas', 'energy', 'fit_length', 'ang_dist', 'sigma', 'sigma_err'])

# Loop through all gases, energies, and fit lengths
for gas in gases:
    for Energy in Energies:

        # Only proceed if directory exists for this Gas / Energy
        if not os.path.isdir('/Users/majdghrear/data/e_ang_res/MS_fit_sims/'+gas+'/'+str(Energy)+'_keV/'):
            continue

        for dist in dists:

            # List to store angle for each recoil track
            angles = []

            # Loop through different simulation files (100 for each case)
            for i in range(efiles):

                # Path to degrad dataframe
                file_path = '~/data/e_ang_res/MS_fit_sims/'+gas+'/'+str(Energy)+'_keV/'+gas+'_'+str(Energy)+'keV_' + str(i) + '/'+gas+'_'+str(Energy)+'keV_' + str(i) + '.root'

                # Read the recoil tracks and corresponding time ordering
                tracks,times = SimTools.read_degrad(file_path)

                for track,time in zip(tracks,times):

                    #Get a linear direction of the recoil up to the distance specified
                    dir1,x,y,z = SimTools.get_dir(track,time,dist)

                    #The true direction is the z-direction
                    dir0 = np.array([0.0,0.0,1.0])

                    #Compute angle between true direction and fit direction
                    ang = np.arccos(np.dot(dir0,dir1))

                    #Assign - sign to directions along - y direction
                    if dir1[1] < 0.0:
                        ang = -ang

                    angles += [ang]

            # Fit to a Gaussian to determine sigma
            sigma, sigma_err = SimTools.fit_gauss(angles, gas, Energy, dist)

            # Add results to dataframe
            df_results = df_results.append({'gas' : gas, 'energy' : Energy, 'fit_length' : dist, 'ang_dist' : angles, 'sigma' :  sigma, 'sigma_err' : sigma_err }, ignore_index = True)

# Svae the dataframe
df_results.to_pickle('./MS_fitting_data.pk')
