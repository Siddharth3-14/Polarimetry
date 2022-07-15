from math import degrees
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 18})

# plt.rcParams["figure.autolayout"] = True
# plt.rcParams["axes.edgecolor"] = "black"
# plt.rcParams["axes.linewidth"] = 1.5

# l323 = np.array([22.48, 1.31972924])*np.pi/180
# l328 = np.array([26.90, 0.3437104])*np.pi/180
# l331 = np.array([13.97, 2.10144931])*np.pi/180
# print(l331[0]/(np.sqrt(2 - l331[0]**2)) - (l331[0] - l331[1])/(np.sqrt(2 - (l331[0] - l331[1])**2)))

   # path to datafile
cloud ='L328'
datafile = 'data_gaia_added/' + cloud + '_pol_mod2.csv'
data = pd.read_csv(datafile,delimiter=',')   #import data      
data = data[(data['P']/data['eP']) > 2]   #only taking data with P/eP more than 2
data = data[data['distance_parallax']/1000 < 4]   #only taking data with distance less than 4 Kpc
data = data[0<data['distance_parallax']/1000]   #only taking data withh positive distance
ra = np.array(data['ra'])   #Right Ascension
dec = np.array(data['dec'])   # Declination
PA = np.array(data['PA'])   #Polarization Angle
ePA = np.array(data['ePA'])   # Error in Polarization Angle

def Calc_l(ra1,dec1,ra2,dec2):

    c1 = SkyCoord(ra1,dec1,unit = 'deg')
    c2 = SkyCoord(ra2,dec2,unit = 'deg')
    sep = c1.separation(c2)
    return sep.arcminute

distance_mat = np.zeros((ra.shape[0],ra.shape[0]))
PolAngleSqr_diff = np.zeros((ra.shape[0],ra.shape[0]))
ePolAngleSqr_mat = np.zeros((ra.shape[0],ra.shape[0]))

for i in range(ra.shape[0]):
    for j in range(i):
        distance_mat[i,j] = Calc_l(ra[i],dec[i],ra[j],dec[j])
        if (abs(PA[i] - PA[j])<90):
            PolAngleSqr_diff[i,j] = (PA[i] - PA[j])**2
            ePolAngleSqr_mat[i,j] = ePA[i]**2+ePA[j]**2
            distance_mat[i,j] = Calc_l(ra[i],dec[i],ra[j],dec[j])



l_bins_length = 0.75
l_bins = np.arange(0,np.max(distance_mat)+l_bins_length,l_bins_length)
l_bins_centre = (l_bins[:-1] + l_bins[1:])/2

struc_func = []
error_func = []
error1_func = []
count_array = []

for k in range(l_bins.shape[0]-1):
    count = 0
    temp_PolAngle = []
    temp_ePolAngle = []
    for i in range(PolAngleSqr_diff.shape[0]):
        for j in range(i):
            if l_bins[k]< distance_mat[i,j]<=l_bins[k+1]:
                count += 1
                temp_PolAngle.append(PolAngleSqr_diff[i,j])
                temp_ePolAngle.append(ePolAngleSqr_mat[i,j])
            else:
                pass
    if count != 0:
        temp_PolAngle = np.array(temp_PolAngle)
        temp_ePolAngle = np.array(temp_ePolAngle)
        std_PA = np.mean(temp_PolAngle) - np.mean(temp_ePolAngle)
        mean_PA = np.sqrt(std_PA)
        a = np.sqrt(temp_ePolAngle)
        err_PA1 = np.std(a)
        err_PA = err_PA1/np.sqrt(count)
        struc_func.append(mean_PA)
        error1_func.append(err_PA1)
        error_func.append(err_PA)
        count_array.append(count)
    else:
        struc_func.append(0)
        error1_func.append(0)
        error_func.append(0)

struc_func = np.array(struc_func)
error_func = np.array(error_func)
struc_func = np.sqrt(struc_func**2 - error_func**2)

l_bins_centre_mod = l_bins_centre-l_bins_centre*(l_bins_centre>10)
# struc_func_mod = struc_func-struc_func*(l_bins_centre>10)
# error_func_mod = error_func - error_func*(l_bins_centre>10)

l_bins_centre_mod = l_bins_centre_mod[l_bins_centre_mod != 0]
struc_func_mod = struc_func[:l_bins_centre_mod.shape[0]]
error_func_mod = error_func[:l_bins_centre_mod.shape[0]]


##################### calculating b paramter and magnitude of magnetic field 
def fitfun(l, m, b):
  y = np.sqrt(b**2 + (m**2) * (l**2))
  return y

popt, pcov = curve_fit(fitfun, l_bins_centre_mod[0:3], struc_func_mod[0:3])
line = np.arange(0,np.max(l_bins),0.5)
fitted_line = fitfun(line,popt[0],popt[1])

print(l_bins_centre_mod.shape,struc_func_mod.shape)
b = popt[1]
b1 = popt[1]*np.pi/180
perr = np.sqrt(np.diag(pcov))
Bmag = b1/(np.sqrt(2-(b1**2)))
column_density = 500
velocity_dispersion_FWHM = 2.35*0.5
Bmag_modCF = 9.3*(np.sqrt(2*column_density))*velocity_dispersion_FWHM*(1/b)
print('perr',perr)
print('Bmag',Bmag_modCF)
print('b',b)

#n = 100/cc or 500
################### plotting 
print( cloud + ' {linebreak}b: {b:.2f} deg{linebreak}b1: {b1:.2f} radian{linebreak}B_ratio: {Bmag:.2f}{linebreak}B_modCF: {Bmag_modCF:.2f}$\mu Gauss$'.format(b=b,b1=b1,Bmag=Bmag,Bmag_modCF=Bmag_modCF,linebreak='\n'))
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.step(l_bins_centre ,struc_func,'k',where='mid',linewidth=2.0)
ax1.scatter(l_bins_centre_mod ,struc_func_mod,c = 'k',s=40)
plt.errorbar(l_bins_centre_mod, struc_func_mod,color='dimgray',yerr=error_func_mod, fmt="o")
ax1.plot(line,fitted_line,'gray',linewidth=2.0)
plt.text(0.5,35,cloud,fontsize = 18,weight="bold")

textbox = cloud + ' {linebreak}b: {b:.2f} deg{linebreak}b1: {b1:.2f} radian{linebreak}B_ratio: {Bmag:.2f}{linebreak}B_modCF: {Bmag_modCF:.2f}$\mu Gauss$'.format(b=b,b1=b1,Bmag=Bmag,Bmag_modCF=Bmag_modCF,linebreak='\n')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax1.text(0.05, 0.95, textbox, transform=ax1.transAxes, fontsize=10,verticalalignment='top', bbox=props)
ax1.set_ylim(20,40)
ax1.set_xlim(0,5)
ax1.set_xlabel(r'$\ell$ (arcmin)',fontsize = 20)
ax1.set_ylabel(r'Angular Dispersion $\langle (\Delta \phi)^2 \rangle^{1/2}$' + '(\N{DEGREE SIGN})',fontsize =20)
plt.show()

