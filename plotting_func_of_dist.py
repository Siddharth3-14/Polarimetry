import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})


# plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1.5



cloud = 'L323'

input_data = 'data_gaia_added/' + cloud + '_pol_mod2.csv'
data = pd.read_csv(input_data, sep=",", header=0)
data = data[(data['P']/data['eP']) > 2]
data = data[data['distance_parallax']/1000 < 6]
data = data[0<data['distance_parallax']/1000]
dist_array = np.array(data['distance_parallax'])/1000
P_array = data['P']
P_error_array = data['eP']
PA_array = data['PA']
PA_error_array = data['ePA']

PA_array = PA_array-180*(PA_array>90)



# plt.figure()
# plt.scatter(dist_array,P_array,s=30)
# plt.plot(distance_center,polarization)
# plt.ylabel('Polarization (%)',fontsize= 20)
# plt.xlabel('distance (Kpc)',fontsize= 20)
# plt.errorbar(dist_array, P_array,color='dimgray',yerr=P_error_array, fmt="o")
# plt.text(0.3,6,cloud,fontsize = 18,weight="bold")
# # plt.title('Polarization vs distance Hpol Data',fontsize= 18)
# plt.tight_layout
# plt.show()
# # plt.savefig('figures/function_of_dist/P_dist_L323.png')

# plt.figure()
# plt.scatter(dist_array,PA_array,s=30)
# plt.ylabel('Polarization angle ($\degree$)',fontsize= 20)
# plt.xlabel('distance (Kpc)',fontsize= 20)
# plt.errorbar(dist_array, PA_array,yerr=PA_error_array, fmt="o",color = 'dimgray')
# plt.text(.3,178,cloud,fontsize = 18,weight="bold")

# # plt.title('Polarization Angle vs distance Hpol Data')
# # plt.savefig('figures/function_of_dist/PA_dist_L323.png')
# plt.tight_layout()
# plt.show()



dist_min = 0
dist_max = 6

P_max = np.max(P_array)
P_min = np.min(P_array)

PA_min = np.min(PA_array)
PA_max = np.max(PA_array)
  
dist_bins = np.arange(dist_min, dist_max, 0.25)
P_bins = np.arange(P_min, P_max, 0.25)
PA_bins = np.arange(PA_min, PA_max, 10)


# Creating plot
fig = plt.subplots(figsize =(20, 10))
plt.hist2d(dist_array, P_array,bins =[dist_bins, P_bins])
plt.title("P X D 2D histogram")
plt.ylabel('P [%]')
plt.xlabel('distance [kpc]')
plt.text(.3,5.5,cloud,fontsize = 18,weight="bold")

fig = plt.subplots(figsize =(20, 10))
plt.hist2d(dist_array, PA_array,bins =[dist_bins, PA_bins])
plt.title("PA X D 2D histogram")
plt.ylabel('PA [deg]')
plt.xlabel('distance [kpc]')
plt.text(.3,85,cloud,fontsize = 18,weight="bold")
plt.show()


# plt.hist2d(dist_array, PA_array,bins =[x_bins, y_bins])
# plt.title("PA X D 2D histogram")
# plt.ylabel('P [%]')
# plt.xlabel('distance [kpc]')
  
# show plot
plt.show()
