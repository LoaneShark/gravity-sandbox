# PYTHON IMPORTS
import sys, os, time
import re
# DEPENDENCY IMPORTS
import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
#from mrtparse import *

def main():
	massmodels = import_massmodel_data()

	print(massmodels)
	plot_single_galaxy(massmodels[massmodels['ID'] == 'UGCA442'])

# plots a galaxy given filtered MassModel data
def plot_single_galaxy(gal_data):
	gal_id = gal_data['ID'][0]
	# if there is more than one galaxy in the dataset, filter out all but the first
	gal_data = gal_data[gal_data['ID'] == gal_id]
	x = gal_data['R']
	y = gal_data['Vobs']
	y_err = gal_data['e_Vobs']

	plt.errorbar(x,y,yerr=y_err,xerr=None,fmt='b.')
	plt.title(gal_id)
	plt.xlabel('Radius (kpc)')
	plt.ylabel('V_obs (km/s)')
	plt.show()

# imports the MassModel data from the locally saved copy of the SPARC dataset
def import_massmodel_data(filepath='../data/MassModels_Lelli2016c.mrt'):
	mm_dt = np.dtype([('ID','U255'),				#	    ID 	: Galaxy Identifier
					  ('D',np.float64),			#	     D	: Assumed distance 					(Mpc)
					  ('R',np.float64),			#	     R	: Galactocentric radius 			(kpc)
					  ('Vobs',np.float64),		#	  Vobs	: Observed circular velocity 		(km/s)
					  ('e_Vobs',np.float64),	#	e_Vobs	: Error in Vobs 					(km/s)
					  ('Vgas',np.float64),		#	  Vgas	: Gas velocity contribution 		(km/s)
					  ('Vdisk',np.float64),		#	 Vdisk	: Disk velocity contribution 		(km/s)
					  ('Vbul',np.float64),		#	  Vbul	: Bulge velocity contribution 		(km/s)
					  ('SBdisk',np.float64),	#	SBdisk	: Disk surface brightness 			(solLum/pc2)
					  ('SBbul',np.float64)])	#	 SBbul	: Bulge surface brightness 			(solLum/pc2)

	mm_data = np.empty((3391,),dtype=mm_dt)
	i = 0
	h_len = 25	# number of header lines to skip in the file
	for line in open(filepath):
		# ignore first few lines of header and file description
		if i >= h_len:
			line_data = re.split('\\s+',line)[:10]

			mm_data[i-h_len]['ID'] 		= line_data[0]
			mm_data[i-h_len]['D'] 		= line_data[1]
			mm_data[i-h_len]['R'] 		= line_data[2]
			mm_data[i-h_len]['Vobs'] 	= line_data[3]
			mm_data[i-h_len]['e_Vobs'] 	= line_data[4]
			mm_data[i-h_len]['Vgas'] 	= line_data[5]
			mm_data[i-h_len]['Vdisk'] 	= line_data[6]
			mm_data[i-h_len]['Vbul'] 	= line_data[7]
			mm_data[i-h_len]['SBdisk'] 	= line_data[8]
			mm_data[i-h_len]['SBbul'] 	= line_data[9]
		i += 1

	return mm_data

if __name__ == '__main__':
	main()