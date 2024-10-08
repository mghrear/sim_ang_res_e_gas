# This script reads the raw electron simulation root files (from degrad)
# and outputs processed files where the recoils are diffused
# Must select 1 recoil energy at a time 

import sys
sys.path.append('../')

from ROOT import TFile, TVector3, TMath
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, fabs
import SimTools
from SimTools import WriteNtuple

he_cf4 = SimTools.he_cf4
he_co2 = SimTools.he_co2

###############################################################
###############################################################

################   Variables to be specified  ################

###############################################################
###############################################################

# Select electron recoil energy
Energy = 200

# Drift [cm]
drift = 25.0

#Select gas mixture he_cf4 or he_co2 as defined in SimTools
mixture = he_cf4

###############################################################
###############################################################

sigma_T = mixture['sigma_T']
sigma_r = mixture['sigma_r']
gas = mixture['name']

# Total sigma in cm
if 'overide_sigma' in list(mixture):
	sigma = mixture['overide_sigma'] * float(1e-4)
else:
	sigma = sqrt( sigma_r**2 + (sqrt(drift) * sigma_T )**2 ) * float(1e-4)

# Get file path
file_path = '~/data/e_ang_res/MS_test_sims/'+gas+'/'+str(Energy)+'_keV/'+gas+'_'+str(Energy)+'keV_' + str(0) + '/'+gas+'_'+str(Energy)+'keV_' + str(0) + '.root'

# Read raw degrad file
tracks,times = SimTools.read_degrad(file_path)

# Loop through the degrad file and diffuse the tracks 
tracks2 = []
for track in tracks:

	track2 = []

	for point in track:

		charge = TVector3(point[0],point[1],point[2]) + TVector3(sigma* np.random.normal(),sigma* np.random.normal(),sigma* np.random.normal())

		track2+= [charge]

	tracks2 += [track2]


# Save the diffused file
tree_name = 'recoils'
root_name = '~/data/e_ang_res/MS_Diff_test_sims/'+gas+'/'+str( int(drift) )+'cm_drift_'+str(Energy)+'keV'
WriteNtuple(tracks2, times, tree_name, root_name)


