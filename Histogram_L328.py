from numpy.lib.npyio import load
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})



#The input file
cloud = 'L328'
input_data = '../data_gaia_added/' +cloud +'_pol_mod2.csv'
df = pd.read_csv(input_data, sep=",", header=0)
Polarization_angle = df['PA']
Polarization = df['P']
bins_array = np.arange(-180,220,10)
gauss_array = np.arange(-180,220,0.5)
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

hist, bin_edges = np.histogram(Polarization_angle,bins_array)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
# plt.hist(Polarization_angle,bins_array)


local_minima = [23,30]
local_PA_minima = 90


fig = plt.figure()


#plotting the multgaussian fit plot
ax1 = fig.add_subplot(121)
ax1.hist(Polarization_angle, bins_array ,color='k',histtype = 'step')
ax3 = ax1.twinx()
ax3.scatter(Polarization_angle,Polarization,s=15,color = 'grey')
#fit the left gaussian
left_hist = hist[0:24]
left_bin_centres = bin_centres[0:24]
left_p0 = [np.amax(left_hist), np.mean(left_bin_centres), np.std(left_bin_centres)]
left_coeff, var_matrix_left = curve_fit(gauss, left_bin_centres, left_hist, p0=left_p0)
left_hist_fit = gauss(gauss_array, *left_coeff)

#fit the centre gaussian
centre_bin_centres = bin_centres[local_minima[0]:local_minima[1]]
centre_hist = hist[local_minima[0]:local_minima[1]]
centre_p0 = [np.amax(centre_hist), np.mean(centre_bin_centres), np.std(centre_bin_centres)]
centre_coeff, var_matrix_left = curve_fit(gauss, centre_bin_centres, centre_hist, p0=centre_p0)
centre_hist_fit = gauss(gauss_array, *centre_coeff)


#fit the right gaussian
right_bin_centres = bin_centres[local_minima[1]:]
right_hist = hist[local_minima[1]:]
right_p0 = [np.amax(right_hist), np.mean(right_bin_centres), np.std(right_bin_centres)]
right_coeff, var_matrix_left = curve_fit(gauss, right_bin_centres, right_hist, p0=right_p0)
right_hist_fit = gauss(gauss_array, *right_coeff)

#plotting fitted gauss curves
ax1.plot(gauss_array, left_hist_fit, 'k-',label='$\mu$ :{left_mean:.0f}{linebreak}$\sigma$ :{left_std:.0f}'.format(linebreak = '\n', left_mean =left_coeff[1],left_std=left_coeff[2]))
ax1.plot(gauss_array,centre_hist_fit,'k-.',label='$\mu$ : {centre_mean:.0f}{linebreak}$\sigma$ :{centre_std:.0f}'.format(linebreak = '\n',centre_mean =centre_coeff[1],centre_std=centre_coeff[2]))
ax1.plot(gauss_array, right_hist_fit,'k--', label='$\mu$ :{right_mean:.0f}{linebreak}$\sigma$ :{right_std:.0f}'.format(linebreak='\n',right_mean =right_coeff[1],right_std=right_coeff[2]))


# ax1.title.set_text('L328 polarization angle {linebreak} Left fitted gauss mean and std :{left_mean:.0f}, {left_std:.0f} {linebreak}Centre fitted gauss mean and std :{centre_mean:.0f}, {centre_std:.0f}{linebreak} Right fitted gauss mean and std :{right_mean:.0f}, {right_std:.0f}'.format(linebreak = '\n', left_mean =left_coeff[1],left_std=left_coeff[2], right_mean =right_coeff[1],right_std=right_coeff[2],centre_mean =centre_coeff[1],centre_std=centre_coeff[2] ))
ax1.set_xlabel('Polarization Angle (deg)')
ax3.set_ylabel('Polarization')
ax1.set_ylabel('Number')
ax1.legend()


#plotting single gauss fit
ax2 = fig.add_subplot(122)
# modifying the polarization angles
Polarization_angle_mod = Polarization_angle -180*(Polarization_angle>local_PA_minima)*np.ones_like(Polarization_angle)

hist, bin_edges = np.histogram(Polarization_angle_mod,bins_array)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
ax2.hist(Polarization_angle_mod,bins_array,color='k',histtype = 'step')
ax4 = ax2.twinx()
ax4.scatter(Polarization_angle_mod,Polarization,s=15,color = 'grey')
#fit the left gaussian
left_hist = hist[10:21]
left_bin_centres = bin_centres[10:21]
left_p0 = [np.amax(left_hist), np.mean(left_bin_centres), np.std(left_bin_centres)]
left_coeff, var_matrix_left = curve_fit(gauss, left_bin_centres, left_hist, p0=left_p0)
left_hist_fit = gauss(gauss_array, *left_coeff)

#fit the centre gaussian
# centre_bin_centres = bin_centres[local_minima[0]:local_minima[1]]
# centre_hist = hist[local_minima[0]:local_minima[1]]
# centre_p0 = [np.amax(centre_hist), np.mean(centre_bin_centres), np.std(centre_bin_centres)]
# centre_coeff, var_matrix_left = curve_fit(gauss, centre_bin_centres, centre_hist, p0=centre_p0)
# centre_hist_fit = gauss(gauss_array, *centre_coeff)

#plotting the fitted curve
ax2.plot(gauss_array, left_hist_fit, 'k-',label='$\mu$ :{left_mean:.0f}{linebreak}$\sigma$ :{left_std:.0f}'.format(linebreak = '\n', left_mean =left_coeff[1],left_std=left_coeff[2]))
# ax2.plot(gauss_array, centre_hist_fit, 'k--',label='$\mu$ :{centre_mean:.0f}{linebreak}$\sigma$ :{centre_std:.0f}'.format(linebreak = '\n',centre_mean =right_coeff[1],centre_std=right_coeff[2]))
# ax2.title.set_text('L328 polarization angle modded {linebreak} Left fitted gauss mean and std :{left_mean:.0f}, {left_std:.0f} {linebreak} Centre fitted gauss mean and std :{centre_mean:.0f}, {centre_std:.0f}'.format(linebreak = '\n', left_mean =left_coeff[1],left_std=left_coeff[2], centre_mean =right_coeff[1],centre_std=right_coeff[2] ))
ax2.set_xlabel('Polarization Angle (deg)')
ax2.set_ylabel('Number')
ax4.set_ylabel('Polarization')
ax2.legend()
plt.tight_layout()
plt.show()









