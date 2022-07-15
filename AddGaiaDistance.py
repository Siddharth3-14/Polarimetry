from astropy.coordinates.jparser import search
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.utils import data
import pandas as pd
import numpy as np
from PolarSi import *

'''
This Program uses the astroquery.gaia to find the gaia parallaxes to the stars with given RA and DEC
'''

#The input file
data_file = 'HpolPA_180added'
input_data = 'csv/' + data_file +'.csv'
df = pd.read_csv(input_data, sep=",", header=0)
ra = df['ra']
dec = df['dec']

#Initiating the columns to be added
parralax_array = np.zeros_like(ra) # parallax column
dist_array = np.zeros_like(ra) # Distance column
search_width_array = np.zeros_like(ra) #Search width column

#Executes the searh for RA and Dec
for i in range(ra.shape[0]):
    print(i)
    temp_parallax,temp_dist,temp_search_width = paramters_search(ra[i],dec[i])
    parralax_array[i] = temp_parallax
    dist_array[i] = temp_dist
    search_width_array[i] = temp_search_width

#Adds the new column
new_df = df.assign(parralax = parralax_array)
new_df = new_df.assign(dist=dist_array)
new_df = new_df.assign(search_width=search_width_array)

#Saves the modded csv file
new_data = 'gaia_added/' + data_file + '_mod.csv'
new_df.to_csv(new_data)



