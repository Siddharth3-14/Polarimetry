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


#best way to find the planck polarization would be to get the galactic map
#change to equitoriol and find the angle between the galactic and equitorial using edge detection
#then use this to change the polarizatin angles
#ez

def Gal2EQ(GLON,GLAT):
    GAL_coord =  SkyCoord(l=GLON*u.degree, b=GLAT*u.degree, frame='galactic')
    EQ_coord = GAL_coord.transform_to('fk5')  
    return EQ_coord.ra.degree, EQ_coord.dec.degree 

def Gaussian_derivative(path_to_image,size,sigma):
    image = cv2.imread(path_to_image) # load image
    image = image[:,:,0]
    def convolution(derivative,gaussian,image):
        gaussian_derivative = cv2.filter2D(derivative, -1,gaussian)
        image_convolution = cv2.filter2D(image,-1,gaussian_derivative)
        return image_convolution

    #defining all the derivative kernels
    derivativex_right = (1/6)*np.array([[-1,0,1],[-1, 0, 1],[-1, 0, 1]])
    derivativex_left = (1/6)*np.array([[1,0,-1],[1, 0, -1],[1, 0, -1]])
    derivativey_bottom = (1/6)*np.array([[-1,-1,-1],[0, 0, 0],[1, 1, 1]])
    derivativey_top = (1/6)*np.array([[1,1,1],[0, 0, 0],[-1,-1,-1]])

    #initialising the gaussian kernel
    Gaussian1D = cv2.getGaussianKernel(size, sigma) 
    Gaussian2D = np.matmul(Gaussian1D,np.transpose(Gaussian1D))

    image_gradient_top = convolution(derivativey_top,Gaussian2D,image)
    image_gradient_bottom = convolution(derivativey_bottom,Gaussian2D,image)
    image_gradient_left = convolution(derivativex_left,Gaussian2D,image)
    image_gradient_right = convolution(derivativex_right,Gaussian2D,image)
    
    image_gradient_x = image_gradient_left+image_gradient_right
    image_gradient_y = image_gradient_top+image_gradient_bottom
   
    return image_gradient_x,image_gradient_y

def mod_arctan(x,y):
    if x==0 and y==0:
        return np.nan
    else:
        return (180/np.pi)*np.arctan(x/y)

plank_I_fits = fits.open('FITS_files/PLCK/PLCKI.fits')[0]
plank_Q_fits = fits.open('FITS_files/PLCK/PLCKQ.fits')[0]
plank_U_fits = fits.open('FITS_files/PLCK/PLCKU.fits')[0]
# print(plank_I_fits.header)
FITs_file = 'FITS_files/Optical/Halp_rsubt_L328_323_331_new_header.fits'

GLON_delt = plank_I_fits.header['CDELT1']
GLAT_delt = plank_I_fits.header['CDELT2']
GLON_ref = (plank_I_fits.header['CRPIX1'])
GLAT_ref = (plank_I_fits.header['CRPIX2'])
GLON_ref_value = plank_I_fits.header['CRVAL1']
GLAT_ref_value = plank_I_fits.header['CRVAL2']
GLON_axis_len = plank_I_fits.header['NAXIS1']
GLAT_axis_len = plank_I_fits.header['NAXIS2']
GLON_axis = np.arange(1,GLON_axis_len+1)
GLAT_axis = np.arange(1,GLAT_axis_len+1)
GLON_ref_value = GLON_ref_value - 0.5*GLON_delt
GLAT_ref_value = GLAT_ref_value - 0.5*GLAT_delt


GLAT_array = (GLAT_axis - GLAT_axis_len/2)*GLAT_delt + GLAT_ref_value
GLON_array = GLON_ref_value-(GLON_axis - GLON_axis_len/2)*(GLON_delt*(-1)/np.cos(GLAT_array*0.01745))

GLAT_grid,GLON_grid = np.meshgrid(GLAT_array,GLON_array , sparse=False, indexing='ij')


plt.figure()
ax1 = plt.subplot(121)
ax1.imshow(GLON_grid,origin='lower')
ax2 = plt.subplot(122)
ax2.imshow(GLAT_grid,origin='lower')
plt.show()
RA_grid,DEC_grid = Gal2EQ(GLON_grid,GLAT_grid)


# plt.figure()
# ax1 = plt.subplot(121)
# ax1.imshow(RA_grid,origin='lower')
# ax2 = plt.subplot(122)
# ax2.imshow(DEC_grid,origin='lower')
# plt.show()

# matplotlib.image.imsave('DEC.png', DEC_grid)
# gradient_x,gradient_y = Gaussian_derivative('DEC.png',4,0.001)
# gradient_x = gradient_x*(gradient_y>0)
# gradient_y = gradient_y*(gradient_x>0)

# plt.figure()
# ax1 = plt.subplot(121)
# ax1.imshow(gradient_x,origin='lower')
# ax2 = plt.subplot(122)
# ax2.imshow(gradient_y,origin='lower')
# plt.show()

# temp_angle = []
# for i in range(gradient_x.shape[0]):
#     for j in range(gradient_x.shape[1]):
#         temp_angle.append(mod_arctan(gradient_x[i,j],gradient_y[i,j]))

# angle_offset = np.nanmean(np.array(temp_angle))
print(angle_offset)
plank_I_data = plank_I_fits.data
plank_Q_data = plank_Q_fits.data
plank_U_data = plank_U_fits.data


# #calculating polarization and polarization angle
plank_polarization = np.sqrt(plank_Q_data*plank_Q_data + plank_U_data*plank_U_data)/plank_I_data
plank_theta_mod = (180/np.pi)*0.5*np.arctan(plank_U_data/plank_Q_data) + 90 - angle_offset


plt.figure()
ax1 = plt.subplot(131)
ax1.imshow(RA_grid,origin='lower')
ax2 = plt.subplot(132)
ax2.imshow(DEC_grid,origin='lower')
ax3 = plt.subplot(133)
ax3.imshow(plank_theta_mod,origin='lower')
plt.show()

#averages out te matrix
def bining_avg(matrix,avg):
     shape_y = matrix.shape[0]
     shape_x = matrix.shape[1]
     shape_x_mod = int(shape_x/avg) 
     shape_y_mod = int(shape_y/avg)
     modded_matrix  = np.zeros((shape_y_mod,shape_x_mod))
     error_matrix = np.zeros((shape_y_mod,shape_x_mod))
     for i in range(shape_y_mod):
          for j in range(shape_x_mod):
               modded_matrix[i,j] = np.nanmean(matrix[i*avg:(i+1)*avg,j*avg:(j+1)*avg])
               error_matrix[i,j] = np.nanstd(matrix[i*avg:(i+1)*avg,j*avg:(j+1)*avg])
     return modded_matrix,error_matrix

DEC_grid_avg,_ = bining_avg(DEC_grid,5)
RA_grid_avg,_ = bining_avg(RA_grid,5)
polarization_avg , polarization_error = bining_avg(plank_polarization,5)
theta_avg,theta_error = bining_avg(plank_theta_mod,5)

# #makes the matrix into one column
RA_array = RA_grid_avg.flatten()
DEC_array = DEC_grid_avg.flatten()
polarization_array = polarization_avg.flatten()
polarization_error_array = polarization_avg.flatten()
theta_array = theta_avg.flatten()
theta_error_array = theta_avg.flatten()



scuba = aplpy.FITSFigure(FITs_file)
scuba.ticks.set_xspacing(0.1) # number in degrees
scuba.show_grayscale()

scuba.add_scalebar(0.0057, color='k', corner='top left',lw=2)
scuba.scalebar.set_label("10 %")
scuba.scalebar.set_font(size='xx-large')
scuba.add_colorbar()
scuba.colorbar.set_width(0.25)
scuba.colorbar.set_location('right')
scuba.colorbar.set_axis_label_font(size=14, weight='medium')
scuba.colorbar.set_axis_label_text('$\\rm Intensity_{850\mu m}$(mJy/beam)')
scuba.ticks.show()
scuba.ticks.show_x()
scuba.ticks.show_y()



rat = RA_array
det = DEC_array 
prt = polarization_array
eprt = polarization_error_array
pa = theta_array
epart = theta_error_array
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