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
# - mm_data : numpy array of imported mass model data
# - gal_id  : string containing the galaxy's identifier
def plot_single_galaxy_mm(mm_data, gal_id):
	gal_data = mm_data[mm_data['ID'] == gal_id]
	x = gal_data['R']
	y = gal_data['Vobs']
	y_err = gal_data['e_Vobs']

	plt.errorbar(x,y,yerr=y_err,xerr=None,fmt='r.')
	plt.title('Newtonian Mass Model ['+gal_id+']')
	plt.xlabel('Radius (kpc)')
	plt.ylabel('V_obs (km/s)')
	plt.show()

# plots a galaxy given its photometric profile
# - photo_data : numpy array of imported photometric profile data
# - gal_id     : string containing the galaxy's identifier
def plot_single_galaxy_photo_profile(photo_data, gal_id):
	x = photo_data['Radius']
	y = photo_data['mu']
	y_err = photo_data['Error']

	plt.errorbar(x,y,yerr=y_err,xerr=None,fmt='g.')
	plt.title('Photometric Profile ['+gal_id+']')
	plt.xlabel('Radius (???)')
	plt.ylabel('mu (mag/arcsec^2)')
	plt.show()

# plots a galaxy given its bulge/disk decomposition
# - dec_data : numpy array of imported decomposition data
# - gal_id   : string containing the galaxy's identifier
def plot_single_galaxy_decomp(dec_data, gal_id):
	x = dec_data['Radius']
	y = dec_data['SBdisk']

	plt.plot(x,y,'b.')
	plt.title('Bulge/Disk Decomposition ['+gal_id+']')
	plt.xlabel('Radius (kpc)')
	plt.ylabel('SBdisk (solLum/pc^2)')
	plt.show()

# imports the surface photometric profile data for a given galaxy identifer
# optimized for (Lelli et al, 2016) and ID = UGCA442
# - datapath : string containing the path to the photometric profile data directory
# - gal_ID   : string containing the galaxy identifier
# - num_header_rows : the number of rows containing header and comments to ignore. Typically 1 in this case.
# - num_rows : the number of total rows that contain data in the document, including the header.
# TODO: Verify data units and descriptions. Radius in particular seems to be in a different value than most.
def import_photometric_profile(datapath='../sparc/sfb_LTG/', gal_ID='UGCA442', num_header_rows=1, num_rows=27):
	filepath = datapath + str(gal_ID) + '.sfb'
	pp_dt    = [('Radius',  		np.float64),	#	Radius	:  ???						
				('mu',				np.float64),	#	mu	    :  Absolute Sky Value (?)			(mag/arcsec^2)			
				('Kill',			bool),      	#	Kill	:  ???								bool
				('Error',			np.float64)]	#	Error	:  Error in mu (?)					(mag/arcsec^2)
	pp_data  = np.empty((num_rows-num_header_rows,),dtype=pp_dt)

	i = 0
	h_len = num_header_rows	# number of header lines to skip in the file
	for line in open(filepath):
		# ignore first line(s) as header
		if i >= num_header_rows and i <= num_rows-1:
			line_data = line.split()

			pp_data[i-h_len]['Radius']	= line_data[0]
			pp_data[i-h_len]['mu']    	= line_data[1]
			pp_data[i-h_len]['Kill']  	= np.nan if int(line_data[2]) not in [0,1] else (int(line_data[2]) == 1)
			pp_data[i-h_len]['Error'] 	= line_data[3]
		i += 1

	return pp_data

# imports the Bulge/Disk Decompositions for a given galaxy identifier
# optimized for (Lelli et al, 2016) and ID = UGCA442
# - datapath : string containing the path to the photometric profile data directory
# - gal_ID   : string containing the galaxy identifier
# - num_header_rows : the number of rows containing header and comments to ignore. Typically 1 in this case.
# - num_rows : the number of total rows that contain data in the document, including the header.
def import_bulge_disk_decomps(datapath='../sparc/BulgeDiskDec_LTG/', gal_ID='UGCA442', num_header_rows=1, num_rows=27):
	filepath = datapath + str(gal_ID) + '.dens'
	bdd_dt   = [('Radius',  		np.float64),	#	Radius	: Measurement(?) radius 			(kpc)
				('SBdisk',			np.float64),	#	SBdisk	: Disk surface brightness 			(solLum/pc^2)
				('SBbulge',			np.float64)]    #	SBbulge	: Bulge surface brightness			(solLum/pc^2)
	bdd_data = np.empty((num_rows-num_header_rows,),dtype=bdd_dt)

	i = 0
	h_len = num_header_rows	# number of header lines to skip in the file
	for line in open(filepath):
		# ignore first line(s) as header
		if i >= num_header_rows and i <= num_rows-1:
			line_data = line.split()

			bdd_data[i-h_len]['Radius']	 = line_data[0]
			bdd_data[i-h_len]['SBdisk']  = line_data[1]
			bdd_data[i-h_len]['SBbulge'] = line_data[2]
		i += 1

	return bdd_data
