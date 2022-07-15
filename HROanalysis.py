import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image
import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PolarSi import *
import matplotlib
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 18})

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    

cloud = 'L328'

path_RA_grid = 'HRO_images/'+cloud + '_RA_low_res_grid.npy'
path_DEC_grid = 'HRO_images/'+cloud + '_DEC_low_res_grid.npy'
path_image = 'HRO_images/' + cloud + '_low_res.png'

polarization_file = 'data_gaia_added/' + cloud + '_pol_mod2.csv'
RA_grid = np.load(path_RA_grid)
DEC_grid = np.load(path_DEC_grid)

Polarization_data = pd.read_csv(polarization_file,delimiter=',')
polarization_RA = Polarization_data['ra']
polarization_DEC = Polarization_data['dec']
optical_pa_mod = Polarization_data['PA']

a = mapping_func(RA_grid,DEC_grid,polarization_RA,polarization_DEC,optical_pa_mod,0.005)
matplotlib.image.imsave(cloud +'.png', a)
# a,b,c,d=Gaussian_derivative(path_image,25,1)
# plt.imshow(d)
# # plt.show()
# relative_orientation_mod = HRO_analysis(RA_grid, DEC_grid,polarization_file,path_image)
# bins_array = np.arange(-180,220,10)
# gauss_array = np.arange(-180,220,0.5)
# histogram, bin_edges = np.histogram(relative_orientation_mod,bins_array)

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.hist(relative_orientation_mod, bins_array ,color='k',histtype = 'step')

# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
# p0 = [(np.amax(histogram)+10), np.mean(relative_orientation_mod), np.std(relative_orientation_mod)]
# coeff, var_matrix = curve_fit(gauss, bin_centres, histogram, p0=p0)
# hist_fit = gauss(gauss_array, *coeff)

# top_textbox = cloud + '{linebreak} $\mu$: {mean:.0f}{linebreak}$\sigma$: {std:.0f}'.format(mean=coeff[1],std=coeff[2],linebreak='\n')
# props = dict(boxstyle='round', facecolor='white', alpha=0.5)
# ax1.text(0.05, 0.95, top_textbox, transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)
# ax1.plot(gauss_array, hist_fit,'k-')
# ax1.set_xlabel('$\Phi$',fontsize= 20)
# ax1.set_ylabel('Number',fontsize= 20)
# plt.tight_layout()
# plt.show()

