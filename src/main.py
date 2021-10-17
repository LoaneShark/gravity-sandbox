# SYS IMPORTS
import sys, os, time
import re
# PACKAGE IMPORTS
import numpy as np 
import scipy as sp
#from mrtparse import *
# LOCAL IMPORTS
from src.tools import *

hubble_types = {0: 'S0', 1: 'Sa', 2: 'Sab', 3: 'Sb', 4: 'Sbc', 5: 'Sc', 6: 'Scd', 7: 'Sd', 8: 'Sdm', 9: 'Sm', 10: 'Im', 11: 'BCD'}
dist_methods = {1: 'Hubble-Flow', 2: 'RGB Tip Magnitude', 3: 'Cepheid Magnitude-Period Relation', 4: 'Ursa Major Cluster', 5: 'SN Light Curve'}
quality_flag = {1: 'High', 2: 'Medium', 3: 'Low'}

def main():
	massmodels = import_massmodel_data()
	metadata, refs = import_galaxy_sample_data()
	mm_dict = {mm['ID']: massmodels[massmodels['ID'] == mm['ID']] for mm in massmodels}

	print(True)

def debug(gal_ID):
	massmodels = import_massmodel_data()
	metadata, refs = import_galaxy_sample_data()
	mm_dict = {mm['ID']: massmodels[massmodels['ID'] == mm['ID']] for mm in massmodels}

	# sample data
	print(mm_dict[gal_ID])
	print(metadata[gal_ID])
	print(refs[metadata[gal_ID]['Ref']])

	# sample plot
	plot_single_galaxy(massmodels,gal_ID)

# imports the Newtonian Mass Model data from the locally saved copy of the SPARC dataset
# optimized for (Lelli et al, 2016)
def import_massmodel_data(filepath='../sparc/MassModels_Lelli2016c.mrt', num_header_rows=25, num_rows=3391):
	mm_dt = np.dtype([('ID',		    'U255'),	#	    ID 	: Galaxy Identifier
					  ('D',			np.float64),	#	     D	: Assumed distance 					(Mpc)
					  ('R',			np.float64),	#	     R	: Galactocentric radius 			(kpc)
					  ('Vobs',		np.float64),	#	  Vobs	: Observed circular velocity 		(km/s)
					  ('e_Vobs',	np.float64),	#	e_Vobs	: Error in Vobs 					(km/s)
					  ('Vgas',		np.float64),	#	  Vgas	: Gas velocity contribution 		(km/s)
					  ('Vdisk',		np.float64),	#	 Vdisk	: Disk velocity contribution 		(km/s)
					  ('Vbul',		np.float64),	#	  Vbul	: Bulge velocity contribution 		(km/s)
					  ('SBdisk',	np.float64),	#	SBdisk	: Disk surface brightness 			(solLum/pc^2)
					  ('SBbul',		np.float64)])	#	 SBbul	: Bulge surface brightness 			(solLum/pc^2)
	mm_data = np.empty((num_rows,), dtype=mm_dt)

	i = 0
	h_len = num_header_rows	# number of header lines to skip in the file
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

# imports the galaxy sample data, providing attributes and classifications for each sample
# optimized for SPARC (Lelli et al, 2016c)
def import_galaxy_sample_data(filepath='../sparc/SPARC_Lelli2016c.mrt', num_header_rows=41, data_start_row=98):
	citations = {}
	sample_data = {}
	meta_dt = np.dtype([('ID',			'U255'),		#	    ID 	: Galaxy Identifier
						('T',			np.int32),		#	     T	: Hubble Type 						hubble_types[T]
						('D',			np.float64),	#	     D	: Assumed distance 					(Mpc)
						('e_D',			np.float64),	#	   e_D	: Mean error in D 					(Mpc)
						('f_D',			np.int32),		#	   f_D	: Distance Method 					dist_methods[f_D]
						('Inc',			np.float64),	#	   Inc	: Inclination 						(degrees)
						('e_Inc',		np.float64),	#	 e_Inc	: Mean error in Inclination 		(degrees)
						('L',			np.float64),	#	     L	: Total Luminosity at [3.6]			(10+9solLum)
						('e_L',			np.float64),	#	   e_L	: Mean error in Luminosity			(10+9solLum)
						('Reff',		np.float64),	#	  Reff	: Effective Radius at [3.6] 		(kpc)
						('SBeff',		np.float64),	#	 SBeff	: Eff. Surface Brightness at [3.6] 	(solLum/pc2)
						('Rdisk',		np.float64),	#	 Rdisk	: Disk scale length at [3.6] 		(kpc)
						('SBdisk',		np.float64),	#	SBdisk	: Disk Surface Brightness at [3.6] 	(solLum/pc2)
						('MHI',			np.float64),	#	   MHI	: Total HI Mass 					(10+9solMass)
						('RHI',			np.float64),	#	   RHI	: HI Radius at 1 solMass/pc2 		(kpc)
						('Vflat',		np.float64),	#	 Vflat	: Asymptotic flat rotation vel.		(km/s)
						('e_Vflat',		np.float64),	#  e_Vflat	: Mean error in Vflat 				(km/s)
						('Q',			np.int32),		#	     Q	: Quality flag 						quality_flag[T]
						('Ref',			'U255')])		#	   Ref	: Reference citation				citations[Ref]

	i = 0
	for line in open(filepath):
		# begin citation reference data
		if i >= num_header_rows and i <= (data_start_row-2):
			ref, citation = line.split(' = ')
			ref = ref.split(' ')[1]
			citation = citation.split('\n')[0]
			citations[ref] = citation
		# begin galaxy sample data
		if i >= data_start_row:
			line_data = re.split('\\s+',line)
			# ignore leading/trailing whitespace
			if line_data[0] == '':
				line_data = line_data[1:]
			if line_data[-1] == '':
				line_data = line_data[:-1]
			# create a numpy array with labeled and type-checked fields
			line_arr = np.empty((1,), dtype=meta_dt)
			for dtype_field,val in zip(meta_dt.fields.items(), line_data):
				dtype_label = dtype_field[0]
				line_arr[dtype_label] = val

			sample_data[line_arr[0]['ID']] = line_arr[0]
		i += 1

	return sample_data,citations

if __name__ == '__main__':
	main()

if __name__ == '__debug__':
	debug('UGCA442')