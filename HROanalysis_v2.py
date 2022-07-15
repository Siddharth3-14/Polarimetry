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


def mapping_func(grid_xx,grid_yy,x,y,pol,pol_angle,tol = 0.001):

    mappped_pol = np.zeros((grid_xx.shape[0],grid_xx.shape[1]))
    mappped_pol_angle = np.zeros((grid_xx.shape[0],grid_xx.shape[1]))

    for i in range(pol.shape[0]):
        temp1 = (abs(grid_xx - x[i])<tol)*(abs(grid_yy - y[i])<tol)*pol[i]
        temp2 = (abs(grid_xx - x[i])<tol)*(abs(grid_yy - y[i])<tol)*pol_angle[i].real

        mappped_pol += temp1*(mappped_pol==0)
        mappped_pol_angle += temp2*(mappped_pol_angle == 0)

    return mappped_pol,mappped_pol_angle


def L328_divider(x, y):
    L328_divider_bottom_coord_x = 274.2166155
    L328_divider_bottom_coord_y = -18.1189
    L328_divider_top_coord_x = 274.17644
    L328_divider_top_coord_y = -18.08743
    slope = (L328_divider_top_coord_y - L328_divider_bottom_coord_y) / \
        (L328_divider_top_coord_x-L328_divider_bottom_coord_x)
    return slope*(x-L328_divider_top_coord_x) - (y - L328_divider_top_coord_y)


cloud = 'L328'
intensity = 2000
polarization_file = 'data_gaia_added/' + cloud + '_pol_mod2.csv'      

if cloud =='L328':
    RA_bottom_left=274.3282075
    DEC_bottom_left=-18.2947917
    RA_top_right= 274.1062916
    DEC_top_right=-17.9863183
elif cloud == 'L323':
    RA_bottom_left=273.9858242
    DEC_bottom_left=-18.3837631
    RA_top_right= 273.7027558
    DEC_top_right=-18.0594761
elif cloud == 'L331':
    RA_bottom_left=274.3091877
    DEC_bottom_left=-17.9927452
    RA_top_right= 274.1308738
    DEC_top_right=-17.8184623


L328_df = pd.read_csv(polarization_file, delimiter=',')
L328_top_df = L328_df[L328_divider(L328_df['ra'],L328_df['dec'])<0]
L328_bottom_df = L328_df[L328_divider(L328_df['ra'],L328_df['dec'])>0]




region_image = FITS('FITS_files/Optical/Halp_rsubt_mod.fits')
region_image.generate_RA_DEC_mesh()

cloud_image,cloud_RA_mesh,cloud_DEC_mesh = region_image.slicing_image(RA_bottom_left,DEC_bottom_left,RA_top_right, DEC_top_right)
# print((cloud_DEC_mesh[1,0] - cloud_DEC_mesh[0,0])*60*60)


############## optical polarization

# Polarization_data = pd.read_csv(polarization_file,delimiter=',')
# Polarization_data = L328_top_df
Polarization_data = L328_bottom_df

polarization_RA = np.array(Polarization_data['ra'])
polarization_DEC = np.array(Polarization_data['dec'])
optical_pol = np.array(Polarization_data['P'])
optical_pa_mod = np.array(Polarization_data['PA'] + 1j)
optical_mapped_pol,optical_mapped_pol_angle = mapping_func(cloud_RA_mesh,cloud_DEC_mesh,polarization_RA,polarization_DEC,optical_pol,optical_pa_mod,1*0.001661841331264088)


optical_polarization_x = optical_mapped_pol*np.cos(np.radians(optical_mapped_pol_angle))
optical_polarization_y = optical_mapped_pol*np.sin(np.radians(optical_mapped_pol_angle))

cloud_image_mod = cloud_image*(cloud_image < intensity)
matplotlib.image.imsave('HRO_images/cloud_mod.png', cloud_image_mod)
gradient_x,gradient_y = Gaussian_derivative('HRO_images/cloud_mod.png',35,1)

plt.figure()
ax1 = plt.subplot(121)
ax1.imshow(gradient_x)
ax2 = plt.subplot(122)
ax2.imshow(gradient_y)
plt.show()

def gen_relative_orientation(x1,y1,x2,y2):
    def mod_arctan(x,y):
        if x==0 and y==0:
            return np.nan
        else:
            return (180/np.pi)*np.arctan(x/y)
    scalar_product = x1*x2 + y1*y2
    cross_product = x1*y2 - x2*y1
    temp = np.zeros_like(cross_product)
    for i in range(cross_product.shape[0]):
        for j in range(cross_product.shape[1]):
            temp[i,j] = mod_arctan(cross_product[i,j],scalar_product[i,j])
    plt.imshow(temp)
    plt.show()
    temp2 = temp.flatten()
    temp2 = temp2[temp2 != np.nan]
    return temp2   

# plt.figure()
# ax1 = plt.subplot(121)
# ax1.imshow(gradient_x)
# ax2 = plt.subplot(122)
# ax2.imshow(gradient_y)
# plt.show()

relative_orientation = gen_relative_orientation(optical_polarization_x,optical_polarization_y,gradient_x,gradient_y)
bins_array = np.arange(-90,90,0.5)
gauss_array = np.arange(-180,220,10)
histogram, bin_edges = np.histogram(relative_orientation,bins_array)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
width = bin_edges[1] - bin_edges[0]
Ac = 0
Ae = 0
h_c = 0
h_e = 0
h_tot = 0
for i in range(bin_centres.shape[0]):
    if -22.5<bin_centres[i]<22.5:
        Ac += width*histogram[i]
        h_c += histogram[i]
        h_tot += histogram[i]
    elif (-90 < bin_centres[i]<-67.5) or (67.5<bin_centres[i]<90):
        Ae += width*histogram[i]
        h_e += histogram[i]
        h_tot += histogram[i]

E = (Ac - Ae)/(Ac + Ae)
sigma_c_sqr = h_c*(1 - h_c/h_tot)
sigma_e_sqr = h_e*(1 - h_e/h_tot)
sigma_E = np.sqrt((4*(Ae*Ae*sigma_e_sqr + Ac*Ac*sigma_c_sqr))/((Ac + Ae)**4))

print('Ac' , Ac)
print('Ae', Ae)
print('E', E)
print('h_c',h_c)
print('h_e',h_e)
print('h_tot',h_tot)
print('sigma_c_sqr',sigma_c_sqr)
print('sigma_e_sqr',sigma_e_sqr)
print('sigma_E',sigma_E)
print('max_intensity', intensity)
# plt.figure()
# ax1 = plt.subplot(121)
# ax1.imshow(optical_polarization_x)
# ax2 = plt.subplot(122)
# ax2.imshow(optical_polarization_y)
# plt.show()

# x = 600
# E_array = []
# E_error_array = []
# x_array = []
# for i in range(4):
#     print(x)
#     cloud_image_mod = cloud_image*(x<cloud_image)*(cloud_image < (x+300))
#     matplotlib.image.imsave('HRO_images/cloud_mod.png', cloud_image_mod)
#     gradient_x,gradient_y = Gaussian_derivative('HRO_images/cloud_mod.png',35,1)

#     # plt.figure()
#     # ax1 = plt.subplot(121)
#     # ax1.imshow(gradient_x)
#     # ax2 = plt.subplot(122)
#     # ax2.imshow(gradient_y)
#     # plt.show()

#     def gen_relative_orientation(x1,y1,x2,y2):
#         def mod_arctan(x,y):
#             if x==0 and y==0:
#                 return np.nan
#             else:
#                 return (180/np.pi)*np.arctan(x/y)
#         scalar_product = x1*x2 + y1*y2
#         cross_product = x1*y2 - x2*y1
#         temp = np.zeros_like(cross_product)
#         for i in range(cross_product.shape[0]):
#             for j in range(cross_product.shape[1]):
#                 temp[i,j] = mod_arctan(cross_product[i,j],scalar_product[i,j])
#         # plt.imshow(temp)
#         # plt.show()
#         temp2 = temp.flatten()
#         temp2 = temp2[temp2 != np.nan]
#         return temp2   

#     # plt.figure()
#     # ax1 = plt.subplot(121)
#     # ax1.imshow(gradient_x)
#     # ax2 = plt.subplot(122)
#     # ax2.imshow(gradient_y)
#     # plt.show()

#     relative_orientation = gen_relative_orientation(optical_polarization_x,optical_polarization_y,gradient_x,gradient_y)
#     bins_array = np.arange(-90,90,0.5)
#     gauss_array = np.arange(-180,220,10)
#     histogram, bin_edges = np.histogram(relative_orientation,bins_array)
#     bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
#     width = bin_edges[1] - bin_edges[0]
#     Ac = 0
#     Ae = 0
#     h_c = 0
#     h_e = 0
#     h_tot = 0
#     for i in range(bin_centres.shape[0]):
#         if -22.5<bin_centres[i]<22.5:
#             Ac += width*histogram[i]
#             h_c += histogram[i]
#             h_tot += histogram[i]
#         elif (-90 < bin_centres[i]<-67.5) or (67.5<bin_centres[i]<90):
#             Ae += width*histogram[i]
#             h_e += histogram[i]
#             h_tot += histogram[i]

#     E = (Ac - Ae)/(Ac + Ae)
#     sigma_c_sqr = h_c*(1 - h_c/h_tot)
#     sigma_e_sqr = h_e*(1 - h_e/h_tot)
#     sigma_E = np.sqrt((4*(Ae*Ae*sigma_e_sqr + Ac*Ac*sigma_c_sqr))/((Ac + Ae)**4))

#     # print('Ac' , Ac)
#     # print('Ae', Ae)
#     # print('E', E)
#     # print('h_c',h_c)
#     # print('h_e',h_e)
#     # print('h_tot',h_tot)
#     # print('sigma_c_sqr',sigma_c_sqr)
#     # print('sigma_e_sqr',sigma_e_sqr)
#     # print('sigma_E',sigma_E)
#     E_array.append(E)
#     E_error_array.append(sigma_E)
#     x_array.append((x + 300)/2)
#     x+= 300

# plt.plot(x_array,np.zeros_like(x_array),'k')
# plt.plot(x_array ,E_array,'k')
# plt.errorbar(x_array, E_array,color='dimgray',yerr=E_error_array, fmt="o")
# plt.xlabel('Intensity')
# plt.ylabel(r'$ \xi$')
# plt.title(cloud + ' HRO analysis')
# plt.show()
######## Planck polarizaiton


# planck_I_fits = FITS('FITS_files/Stokes/Stokes_I_mod.fits')
# planck_Q_fits = FITS('FITS_files/Stokes/Stokes_Q_mod.fits')
# planck_U_fits = FITS('FITS_files/Stokes/Stokes_U_mod.fits')

# planck_I_fits.generate_RA_DEC_mesh()
# planck_Q_fits.generate_RA_DEC_mesh()
# planck_U_fits.generate_RA_DEC_mesh()

# planck_I_cloud,planckI_RAmesh,planckI_DECmesh = planck_I_fits.slicing_image(274.3282075,-18.2947917,274.1062916,-17.9863183)
# planck_Q_cloud,planckQ_RAmesh,planckQ_DECmesh = planck_Q_fits.slicing_image(274.3282075,-18.2947917,274.1062916,-17.9863183)
# planck_U_cloud,planckU_RAmesh,planckU_DECmesh = planck_U_fits.slicing_image(274.3282075,-18.2947917,274.1062916,-17.9863183)

# plank_polarization = np.sqrt(planck_Q_cloud*planck_Q_cloud + planck_U_cloud*planck_U_cloud)/planck_I_cloud
# plank_theta_mod = (180/np.pi)*0.5*np.arctan(planck_U_cloud/planck_Q_cloud) + 90 + 62

# planck_RA_array = planckI_RAmesh.flatten()
# planck_DEC_array = planckI_DECmesh.flatten()
# planck_pol_array = plank_polarization.flatten()
# # polarization_error_array = plank_polarization.flatten()
# planck_pol_angle_array = plank_theta_mod.flatten()
# # theta_error_array = plank_theta_mod.flatten()



# planck_mapped_pol,planck_mapped_pol_angle = mapping_func(cloud_RA_mesh,cloud_DEC_mesh,planck_RA_array,planck_DEC_array,planck_pol_array,planck_pol_angle_array,1*0.001661841331264088)

# planck_polarization_x = planck_mapped_pol*np.cos(np.radians(planck_mapped_pol_angle))
# planck_polarization_y = planck_mapped_pol*np.sin(np.radians(planck_mapped_pol))

# plt.figure()
# ax1 = plt.subplot(121)
# ax1.imshow(planck_polarization_x)
# ax2 = plt.subplot(122)
# ax2.imshow(planck_polarization_y)
# plt.show()





















