
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import pandas as pd
from astropy.io import fits

FITS1 = 'FITS_files/for_subtraction/Halpha.fits'
FITS2 = 'FITS_files/for_subtraction/DSSR2.fits'
FITS3 = 'FITS_files/for_subtraction/DSS.fits'
FITS1_image = fits.open(FITS1)[0].data
FITS2_image = fits.open(FITS2)[0].data
FITS3_image = fits.open(FITS3)[0].data

data_subtraction = (FITS2_image-FITS3_image)
plt.imshow(data_subtraction)
plt.show()