import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 18})


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    

L328_divider_bottom_coord_x = 274.3584501
L328_divider_bottom_coord_y = -18.1788278
L328_divider_top_coord_x = 274.0540910
L328_divider_top_coord_y = -18.0579074


def L328_divider(x, y):
    slope = (L328_divider_top_coord_y - L328_divider_bottom_coord_y) / \
        (L328_divider_top_coord_x-L328_divider_bottom_coord_x)
    return slope*(x-L328_divider_top_coord_x) - (y - L328_divider_top_coord_y)

datafile= "data_gaia_added/L328_pol_full_mod.csv"
L328_df = pd.read_csv(datafile, delimiter=',')


######################### L328 top part ################################################
################## filterring L328 top data #######################################

L328_top_df = L328_df[L328_divider(L328_df['ra'],L328_df['dec'])<0]
# ############ calculating the shape angle ####################
# local_PA_minima = 125
# L328_top_axis_angle =40
# Pol_Angle_top = np.array(L328_top_df['PA'])
# Pol_Angle_top_mod = Pol_Angle_top -180*(Pol_Angle_top>local_PA_minima)*np.ones_like(Pol_Angle_top)
# top_shape_angle = abs(Pol_Angle_top_mod - L328_top_axis_angle)
# Polarization = L328_top_df['P']


# ############### Generate bins for modelling ###########

# bins_array = np.arange(-180,220,10)
# gauss_array = np.arange(-180,220,0.5)
# top_hist, bin_edges = np.histogram(top_shape_angle,bins_array)
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# ################## plotting histogram #################
# ################## plotting histogram #################
# fig = plt.figure()
# ax1 = fig.add_subplot(121)
# ax1.hist(top_shape_angle, bins_array ,alpha=0.5, color='k',histtype = 'step')
# # ax2 = ax1.twinx()
# # ax2.scatter(top_shape_angle,Polarization,s=5,color = 'grey')

# ##################### Modelling the histogram ######################
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
# p0 = [np.amax(top_hist), np.mean(top_shape_angle), np.std(top_shape_angle)]
# top_coeff, var_matrix = curve_fit(gauss, bin_centres, top_hist, p0=p0)
# top_hist_fit = gauss(gauss_array, *top_coeff)


# #################### plotting the modelled gaussian 
# top_textbox = 'L328 Top {linebreak} $\mu$: {mean:.0f}{linebreak}$\sigma$: {std:.0f}'.format(mean=top_coeff[1],std=top_coeff[2],linebreak='\n')
# props = dict(boxstyle='round', facecolor='white', alpha=0.5)
# ax1.text(0.05, 0.95, top_textbox, transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)
# ax1.plot(gauss_array, top_hist_fit,'k-')
# ax1.set_xlabel('$\\theta_{pol} - \\theta_{cloud}$ (deg)',fontsize= 20)
# ax1.set_ylabel('Number',fontsize= 20)
# plt.tight_layout()

# # ########################### L328 bottom part #########################################
# L328_bottom_df = L328_df[L328_divider(L328_df['ra'],L328_df['dec'])>0]
# L328_bottom_axis_angle =3
# Pol_Angle_bottom = np.array(L328_bottom_df['PA'])
# Pol_Angle_bottom_mod = Pol_Angle_bottom -180*(Pol_Angle_bottom>local_PA_minima)*np.ones_like(Pol_Angle_bottom)
# bottom_shape_angle = abs(Pol_Angle_bottom_mod - L328_bottom_axis_angle)
# Polarization = L328_bottom_df['P']

# ############### Generate bins for modelling ###########

# bins_array = np.arange(-180,220,10)
# gauss_array = np.arange(-180,220,0.5)
# bottom_hist, bin_edges = np.histogram(bottom_shape_angle,bins_array)
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# ################## plotting histogram #################
# ################## plotting histogram #################
# ax3 = fig.add_subplot(122)
# ax3.hist(bottom_shape_angle, bins_array ,alpha=0.5, color='k',histtype = 'step')

# ##################### Modelling the histogram ######################
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
# p0 = [np.amax(bottom_hist), np.mean(bottom_shape_angle), np.std(bottom_shape_angle)]
# bottom_coeff, var_matrix = curve_fit(gauss, bin_centres, bottom_hist, p0=p0)
# bottom_hist_fit = gauss(gauss_array, *bottom_coeff)


# # #################### plotting the modelled gaussian 
# top_textbox = 'L328 bottom {linebreak}$\mu$: {mean:.0f}{linebreak}$\sigma$: {std:.0f}'.format(mean=bottom_coeff[1],std=bottom_coeff[2],linebreak='\n')
# props = dict(boxstyle='round', facecolor='white', alpha=0.5)
# ax3.text(0.05, 0.95, top_textbox, transform=ax3.transAxes, fontsize=14,verticalalignment='top', bbox=props)

# ax3.plot(gauss_array, bottom_hist_fit,'k-')
# ax3.set_xlabel('$\\theta_{pol} - \\theta_{cloud}$ (deg)',fontsize= 20)
# ax3.set_ylabel('Number',fontsize= 20)

# plt.tight_layout()
# # plt.show()
# ################################# L323 ############################################
# ########### importing L323 data ##################
# datafile= "data_gaia_added/L323_Pol_ra_dec_mod2.csv"
# L323_df = pd.read_csv(datafile, delimiter=',')

# ############ calculating the shape angle ####################
# L323_main_axis_angle =61
# Shape_angle = abs(np.array(L323_df['PA']) - L323_main_axis_angle)
# Polarization = L323_df['P']

# ############### Generate bins for modelling ###########
# bins_array = np.arange(-180,220,10)
# gauss_array = np.arange(-180,220,0.5)
# hist, bin_edges = np.histogram(Shape_angle,bins_array)
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2


# ################## plotting histogram #################
# fig = plt.figure()
# ax7 = fig.add_subplot(111)
# ax7.hist(Shape_angle, bins_array ,alpha=0.5, color='k',histtype = 'step')
# # ax8 = ax7.twinx()
# # ax8.scatter(Shape_angle,Polarization,s=5,color = 'grey')

# ##################### Modelling the histogram ######################
# bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
# p0 = [np.amax(hist), np.mean(Shape_angle), np.std(Shape_angle)]
# coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
# hist_fit = gauss(gauss_array, *coeff)

# #################### plotting the modelled gaussian 
# top_textbox = 'L323 Top {linebreak} $\mu$: {mean:.0f}{linebreak}$\sigma$: {std:.0f}'.format(mean=coeff[1],std=coeff[2],linebreak='\n')
# props = dict(boxstyle='round', facecolor='white', alpha=0.5)
# ax7.text(0.05, 0.95, top_textbox, transform=ax7.transAxes, fontsize=14,verticalalignment='top', bbox=props)
# ax7.plot(gauss_array, hist_fit,'k-')
# # ax2.title.set_text('L323 polarization angle modded{linebreak} Fitted mean = {mean:.3f} {linebreak} Fitted standard deviation = {std:.3f}'.format(mean=coeff[1],std=coeff[2],linebreak = '\n'))
# ax7.set_xlabel('$\\theta_{pol} - \\theta_{cloud}$ (deg)',fontsize= 20)
# ax7.set_ylabel('Number',fontsize= 20)
# # ax8.set_ylabel('Polarization')
# # ax7.legend()
# plt.tight_layout()
# plt.show()

