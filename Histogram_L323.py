import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

#The input file
cloud = 'L323'
input_data = '../data_gaia_added/' + cloud + '_pol_mod2.csv'
df = pd.read_csv(input_data, sep=",", header=0)
Polarization_angle = df['PA']
Polarization = df['P']
bins_array = np.arange(-180,220,10)
gauss_array = np.arange(-180,220,0.5)


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# plt.hist(Polarization_angle, bins_array ,alpha=0.5, color='y')

hist, bin_edges = np.histogram(Polarization_angle,bins_array)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
local_minima = 26
local_PA_minima = 90
fig = plt.figure()


#plotting the multgaussian fit plot
ax1 = fig.add_subplot(121)

ax1.hist(Polarization_angle, bins_array ,alpha=0.5, color='k',histtype = 'step')
ax3 = ax1.twinx()
ax3.scatter(Polarization_angle,Polarization,s=20,color = 'grey')

#fit the left gaussian
left_hist = hist[18:local_minima]
left_bin_centres = bin_centres[18:local_minima]
left_p0 = [np.amax(left_hist), np.mean(left_bin_centres), np.std(left_bin_centres)]
# print(left_p0)
left_coeff, var_matrix_left = curve_fit(gauss, left_bin_centres, left_hist, p0=left_p0)
left_hist_fit = gauss(gauss_array, *left_coeff)

#fit the right gaussian
right_bin_centres = bin_centres[local_minima:]
right_hist = hist[local_minima:]
right_p0 = [np.amax(right_hist), np.mean(right_bin_centres), np.std(right_bin_centres)]
right_coeff, var_matrix_left = curve_fit(gauss, right_bin_centres, right_hist, p0=right_p0)
right_hist_fit = gauss(gauss_array, *right_coeff)

#plotting fitted gauss curves
ax1.plot(gauss_array, left_hist_fit,'k--', label='$\mu$:{left_mean:.0f}{linebreak} $\sigma$:{left_std:.0f}'.format(left_mean =left_coeff[1],left_std=left_coeff[2],linebreak='\n'))
ax1.plot(gauss_array, right_hist_fit,'k-',label='$\mu$:{right_mean:.0f}{linebreak} $\sigma$:{right_std:.0f}'.format(right_mean =right_coeff[1],right_std=right_coeff[2],linebreak='\n'))
# ax1.title.set_text('L323 polarization angle {linebreak} Left fitted gauss mean and std :{left_mean:.0f}, {left_std:.0f} {linebreak} Right fitted gauss mean and std :{right_mean:.0f}, {right_std:.0f}'.format(linebreak = '\n', left_mean =left_coeff[1],left_std=left_coeff[2], right_mean =right_coeff[1],right_std=right_coeff[2] ))
ax1.set_xlabel('Polarization Angle (deg)')
ax3.set_ylabel('Polarization')
ax1.set_ylabel('Number')
ax1.legend()


#plotting single gauss fit
ax2 = fig.add_subplot(122)

# modifying the polarization angles
Polarization_angle_mod = Polarization_angle -180*(Polarization_angle>local_PA_minima)*np.ones_like(Polarization_angle)

# #plotting the histogram
hist, bin_edges, patches = ax2.hist(Polarization_angle_mod,bins_array ,color='k',histtype = 'step')
ax4 = ax2.twinx()
ax4.scatter( Polarization_angle_mod,Polarization,s=20,color = 'grey')
#fitting the gauss function
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
p0 = [np.amax(hist), np.mean(Polarization_angle_mod), np.std(Polarization_angle_mod)]
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
hist_fit = gauss(gauss_array, *coeff)

#plotting the fitted curve
ax2.plot(gauss_array, hist_fit,'k-', label='$\mu$:{mean:.0f}{linebreak}$\sigma$:{std:.0f}'.format(mean=coeff[1],std=coeff[2],linebreak='\n'))
# ax2.title.set_text('L323 polarization angle modded{linebreak} Fitted mean = {mean:.0f} {linebreak} Fitted standard deviation = {std:.0f}'.format(mean=coeff[1],std=coeff[2],linebreak = '\n'))
ax2.set_xlabel('Polarization Angle (deg)')
ax2.set_ylabel('Number')
ax4.set_ylabel('Polarization')
ax2.legend()
plt.tight_layout()
plt.show()














