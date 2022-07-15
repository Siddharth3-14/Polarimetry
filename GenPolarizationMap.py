from astropy.io import fits
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from math import *
import pandas as pd
from astropy.wcs import WCS

plank_I_fits = fits.open('FITS_files/skyview/Stokes_I.fits')[0]
plank_Q_fits = fits.open('FITS_files/skyview/Stokes_Q.fits')[0]
plank_U_fits = fits.open('FITS_files/skyview/Stokes_U.fits')[0]

# wcs = WCS(plank_I_fits.header)
# print(wcs)

RA_delt = plank_I_fits.header['CDELT1']
DEC_delt = plank_I_fits.header['CDELT2']
RA_ref = int(plank_I_fits.header['CRPIX1'])
DEC_ref = int(plank_I_fits.header['CRPIX2'])
RA_ref_value = plank_I_fits.header['CRVAL1']
DEC_ref_value = plank_I_fits.header['CRVAL2']
RA_axis_len = plank_I_fits.header['NAXIS1']
DEC_axis_len = plank_I_fits.header['NAXIS2']

RA_axis = np.arange(0,RA_axis_len,1)
DEC_axis = np.arange(0,DEC_axis_len,1)
RA_axis = RA_ref_value - RA_delt*(RA_ref - RA_axis)
DEC_axis = DEC_ref_value + DEC_delt*(DEC_ref - DEC_axis)

#making a meshgrid from the arrays
DEC_grid,RA_grid = np.meshgrid(DEC_axis,RA_axis , sparse=False, indexing='ij')

plank_I_data = plank_I_fits.data
plank_Q_data = plank_Q_fits.data
plank_U_data = plank_U_fits.data


# #calculating polarization and polarization angle
plank_polarization = np.sqrt(plank_Q_data*plank_Q_data + plank_U_data*plank_U_data)/plank_I_data
ratio = plank_U_data/plank_Q_data
plank_theta_mod = (180/np.pi)*0.5*np.arctan(plank_U_data/plank_Q_data) + 90 - 63
plt.imshow(plank_theta_mod,origin='lower')
plt.show()
# #averages out te matrix
# def bining_avg(matrix,avg):
#      shape_y = matrix.shape[0]
#      shape_x = matrix.shape[1]
#      shape_x_mod = int(shape_x/avg) 
#      shape_y_mod = int(shape_y/avg)
#      modded_matrix  = np.zeros((shape_y_mod,shape_x_mod))
#      error_matrix = np.zeros((shape_y_mod,shape_x_mod))
#      for i in range(shape_y_mod):
#           for j in range(shape_x_mod):
#                modded_matrix[i,j] = np.nanmean(matrix[i*avg:(i+1)*avg,j*avg:(j+1)*avg])
#                error_matrix[i,j] = np.nanstd(matrix[i*avg:(i+1)*avg,j*avg:(j+1)*avg])
#      return modded_matrix,error_matrix

# DEC_grid_avg,_ = bining_avg(DEC_grid,5)
# RA_grid_avg,_ = bining_avg(RA_grid,5)
# polarization_avg , polarization_error = bining_avg(plank_polarization,5)
# theta_avg,theta_error = bining_avg(plank_theta_mod,5)

# # #makes the matrix into one column
# RA_array = RA_grid_avg.flatten()
# DEC_array = DEC_grid_avg.flatten()
# polarization_array = polarization_avg.flatten()
# polarization_error_array = polarization_avg.flatten()
# theta_array = theta_avg.flatten()
# theta_error_array = theta_avg.flatten()

# # #makes the panda dataframe and saves it into a csv file
# output_df = pd.DataFrame({'ra':RA_array,'dec':DEC_array,'P':polarization_array,'eP':polarization_error_array,'PA':theta_array,'ePA':theta_error_array})
# print(output_df.head)
# output_df = output_df.dropna()
# print(output_df.head)
# output_df.to_csv('data_files_generated/whole_region_plank_polarizationv5.csv')

