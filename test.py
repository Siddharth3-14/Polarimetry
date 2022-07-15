from astropy.io import fits
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib 
import numpy as np
from astropy.table import Table
from math import *
import pandas as pd
import aplpy
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import cv2


cloud  = 'L328'
datafile = 'data_gaia_added/' + cloud + '_pol_mod2.csv'
FITs_file = 'FITS_files/Optical/Halp_rsubt_L328_323_331_new_header.fits'

df = pd.read_csv(datafile,delimiter=',')
print(df.head)
ra_data = np.array(df['ra'])
dec_data = np.array(df['dec'])
p_data = np.array(df['P'])
ep_data = np.array(df['eP'])
pa_data = np.array(df['PA'])
epa_data = np.array(df['ePA'])

ra_max = np.amax(ra_data)
ra_min = np.amin(ra_data)

dec_max = np.amax(dec_data)
dec_min = np.amin(dec_data)

ra_range = np.arange(ra_min,ra_max,0.1)
dec_range = np.arange(dec_min,dec_max,0.1)

ra_array_mod = []
dec_array_mod = []
p_array_mod = []
pa_array_mod = []

for i in range(ra_range.shape[0]-1):
    for j in range(dec_range.shape[0]-1):
        temp_df = df[df['ra']<ra_range[i+1]]
        temp_df = temp_df[ra_range[i]<temp_df['ra']]
        temp_df = temp_df[temp_df['dec']<dec_range[j+1]]
        temp_df = temp_df[dec_range[j]<temp_df['dec']]

        ra_array_mod.append(np.mean(temp_df['ra']))
        dec_array_mod.append(np.mean(temp_df['dec']))
        p_array_mod.append(np.mean(temp_df['P']))
        pa_array_mod.append(np.mean(temp_df['PA']))

output_df = pd.DataFrame({'ra':ra_array_mod,'dec':dec_array_mod,'P':p_array_mod,'eP':p_array_mod,'PA':pa_array_mod,'ePA':pa_array_mod})
output_df = output_df.dropna()
print(output_df.head)
# output_df.to_csv('data_files_generated/whole_region_plank_polarizationv5.csv')
ra_array_mod = np.array(ra_array_mod)
dec_array_mod = np.array(dec_array_mod)
p_array_mod = np.array(p_array_mod)
pa_array_mod = np.array(pa_array_mod)

scuba = aplpy.FITSFigure(FITs_file)
scuba.ticks.set_xspacing(0.1) # number in degrees
scuba.show_grayscale()

scuba.add_scalebar(0.0057, color='k', corner='top left',lw=2)
scuba.scalebar.set_font(size='xx-large')
# scuba.add_colorbar()
# scuba.colorbar.set_width(0.25)
# scuba.colorbar.set_location('right')
# scuba.colorbar.set_axis_label_font(size=14, weight='medium')
# scuba.colorbar.set_axis_label_text('$\\rm Intensity_{850\mu m}$(mJy/beam)')
scuba.ticks.show()
scuba.ticks.show_x()
scuba.ticks.show_y()


rat = ra_array_mod
det = dec_array_mod 
prt = p_array_mod
eprt = p_array_mod
pa = pa_array_mod
epart = pa_array_mod
prt1 = (prt/prt)*100
prt = prt1
part=pa
meanp = 75 * np.mean(prt)

data= fits.open(FITs_file)
p_lo=-data[0].header['CDELT1']
p_la=data[0].header['CDELT2']
ra_ref=data[0].header['CRVAL1']
de_ref=data[0].header['CRVAL2']
xcen=data[0].header['NAXIS1']
ycen=data[0].header['NAXIS2']
xc      = xcen/2
yc      = ycen/2

# #------POL-2data-----------
p_loo1=p_lo/np.cos(det*0.01745)
x1      = (ra_ref-rat)/p_loo1+xc
y1      = (det-de_ref)/p_la+yc
x2      = (ra_ref-(rat+(prt/meanp)*np.sin(0.0175*part)))/p_loo1+xc+10
y2      = ((det+(prt/meanp)*np.cos(0.0175*part))-de_ref)/p_la+yc-10
x3      = (ra_ref-(rat-(prt/meanp)*np.sin(0.0175*part)))/p_loo1+xc+10
y3      = ((det-(prt/meanp)*np.cos(0.0175*part))-de_ref)/p_la+yc-10

#------------plotting-------------

xx2, yy2=scuba.pixel2world(x2, y2)
xx3, yy3=scuba.pixel2world(x3, y3)
scuba.show_arrows(xx3,yy3,xx2-xx3,yy2-yy3, color='cyan', head_length=0,head_width=0, length_includes_head=False, width = 0.5)  # Polarization vector
plt.show()