# SYS IMPORTS
import sys, os, time
import re
# PACKAGE IMPORTS
import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt

def analysis():
	print('Hello Analysis!')

# plots a galaxy given filtered MassModel data
def plot_single_galaxy(mm_data, gal_id):
	gal_data = mm_data[mm_data['ID'] == gal_id]
	x = gal_data['R']
	y = gal_data['Vobs']
	y_err = gal_data['e_Vobs']

	plt.errorbar(x,y,yerr=y_err,xerr=None,fmt='b.')
	plt.title(gal_id)
	plt.xlabel('Radius (kpc)')
	plt.ylabel('V_obs (km/s)')
	plt.show()