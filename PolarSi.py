import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import cv2
from astroquery.gaia import Gaia 
import matplotlib
# import aplpy


def paramters_search(ra_coord,dec_coord):
    """parameters_search

    Given RA and DEC coordinates this program finds the parallax and distance(this is from hipparcos)
    
    Args:
        ra_coord (float): RA coordinate at which to search 
        dec_coord (float): DEC coordinate at which to search 
    Returns:
        parallax: (float): Gaia parallax of the object at the given RA and DEC
        dist: (float): Hipparcos distance to the object at the given RA and DEC
        search_width (float): Search area in which object was found
    """
    Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # select Gaia release 

    search_width = 0 #initialising the search width
    flag = 0 
    while flag == 0: # initialize the while loop till we find an object within search width
        search_width += 0.25
        coord = SkyCoord(ra=ra_coord, dec=dec_coord, unit=(u.degree, u.degree),frame='icrs')
        width = u.Quantity(search_width, u.arcsecond)
        height = u.Quantity(search_width, u.arcsecond)
        results = Gaia.query_object(coordinate=coord, width=width, height=height)
        flag  = results['parallax'].shape[0]
        if abs(search_width - 30*0.25)< 0.2:  # condition so that the search width doesn't increase 0.25 degrees
            return -1,-1,-1 # returns -1 values as parallax and distances can't be negative 
    print(search_width) 
    return results['parallax'][0],results['dist'][0],search_width

def test_paramters_search(ra_coord,dec_coord,search_width):
    """parameters_search

    This is a test for Gaia parameter search function to search within the given search area 
    
    Args:
        ra_coord (float): RA coordinate at which to search 
        dec_coord (float): DEC coordinate at which to search 
        search_width (float): Search area in which to search for the object
    Returns:
        parallax: (float): Gaia parallax of the object at the given RA and DEC
        dist: (float): Hipparcos distance to the object at the given RA and DEC

    """
    coord = SkyCoord(ra=ra_coord, dec=dec_coord, unit=(u.degree, u.degree),frame='icrs')
    width = u.Quantity(search_width, u.arcsecond)
    height = u.Quantity(search_width, u.arcsecond)
    results = Gaia.query_object(coordinate=coord, width=width, height=height)
    print(results['parallax'][0],results['dist'][0])

def window_avg(matrix,avg):
    """window_avg

    This function averages the image in a window of specified length
    
    Args:
        matrix (2D array): matrix which is to be averaged 
        avg (integer): size of the window to average 
    Returns:
        modded_matrix: (2D array): averaged matrix
        error_matrix: (2D array): Error in the matrix
    """
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

class FITS():
    """FITS

    The class for accessing and manipulating fits file
        
    Attributes:
        hdu (astropy.io PrimaryHDU): Primary HDU file in FITS
        header (astropy.io header): Header of the fits file.
        data (numpy array): Intesity data stored in the fits file.
        DEC_grid (2D array): DEC grid generated from the FITS file header
        RA_grid (2D array): RA grid generated from the FITS file header
    """
    def __init__(self,path_to_datafile):
        self.hdu = fits.open(path_to_datafile)[0]
        self.header = self.hdu.header
        self.data = self.hdu.data

    def print_info(self):
        """print_info

        Prints the header and plots the data stored in the fits file
    
        """
        print('Header info')
        print(self.header)
        print("data info")
        print('shape',self.data.shape)
        plt.imshow(self.data,origin='lower')
        plt.show()
    
    def slicing_image(self,RA_bottom_left,DEC_bottom_left,RA_top_right, DEC_top_right):
        delta = self.DEC_grid[1,0] - self.DEC_grid[0,0]
        x_bottom_left = np.where(abs(RA_bottom_left-self.RA_grid)<delta/2)[1][0]
        x_top_right = np.where(abs(RA_top_right-self.RA_grid)<delta/2)[1][0]
        y_bottom_left = np.where(abs(DEC_bottom_left-self.DEC_grid)<delta/2)[0][0]
        y_top_right = np.where(abs(DEC_top_right-self.DEC_grid)<delta/2)[0][0]
        sliced_data = self.data[y_bottom_left:y_top_right,x_bottom_left:x_top_right]
        sliced_RA_mesh = self.RA_grid[y_bottom_left:y_top_right,x_bottom_left:x_top_right]
        sliced_DEC_mesh = self.DEC_grid[y_bottom_left:y_top_right,x_bottom_left:x_top_right]
        return sliced_data,sliced_RA_mesh,sliced_DEC_mesh

    def generate_RA_DEC_mesh(self):
        """generate_RA_DEC_mesh

        Generates the RA and DEC grid for the intensity map
    
        """
        if 'CDELT1' in self.header:
            RA_delt = self.header['CDELT1']
            DEC_delt = self.header['CDELT2']

        if 'CD1_1' in self.header:
            RA_delt = self.header['CD1_1']
            DEC_delt = self.header['CD2_2']

        RA_ref = int(self.header['CRPIX1'])
        DEC_ref = int(self.header['CRPIX2'])
        RA_ref_value = self.header['CRVAL1']
        DEC_ref_value = self.header['CRVAL2']
        RA_axis_len = self.header['NAXIS1']
        DEC_axis_len = self.header['NAXIS2']

        RA_axis = np.arange(0,RA_axis_len,1)
        DEC_axis = np.arange(0,DEC_axis_len,1)
        RA_axis = RA_ref_value - RA_delt*(RA_ref - RA_axis)
        DEC_axis = DEC_ref_value + DEC_delt*(DEC_axis - DEC_ref)

        #making a meshgrid from the arrays
        self.DEC_grid,self.RA_grid = np.meshgrid(DEC_axis,RA_axis , sparse=False, indexing='ij')

def PolarizationMapCSV(path_to_stokes_I,path_to_stokes_Q,path_to_Stokes_U,window_avg_length):
    """PolarizationMapCSV

    This function creates polarization maps from I,Q and U stokes files averages it and then saves it as CSV
    
    Args:
        path_to_stokes_I (string): path to I stokes files 
        path_to_stokes_Q (string): path to Q stokes files 
        path_to_Stokes_U (string): path to U stokes files
        window_avg_length (integer): length of the binning average 
    Returns:
        None: None
    """
    plank_I = FITS(path_to_stokes_I)
    plank_Q = FITS(path_to_stokes_Q)
    plank_U = FITS(path_to_Stokes_U)

    plank_I.generate_RA_DEC_mesh()
    plank_polarization = np.sqrt(plank_Q.data*plank_Q.data + plank_U.data*plank_U.data)/plank_I.data
    plank_theta_mod = (180/np.pi)*0.5*np.arctan(plank_U.data/plank_Q.data) + 90 + 62
 
    def window_avg(matrix,avg):

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

    DEC_grid_avg,_ = window_avg(plank_I.DEC_grid,window_avg_length)
    RA_grid_avg,_ = window_avg(plank_I.RA_grid,window_avg_length)
    polarization_avg , polarization_error = window_avg(plank_polarization,window_avg_length)
    theta_avg,theta_error = window_avg(plank_theta_mod,window_avg_length)
    
    RA_array = RA_grid_avg.flatten()
    DEC_array = DEC_grid_avg.flatten()
    polarization_array = polarization_avg.flatten()
    polarization_error_array = polarization_error.flatten()
    polarization_angle_array = theta_avg.flatten()
    theta_error_array = theta_error.flatten()

    output_df = pd.DataFrame({'ra':RA_array,'dec':DEC_array,'P':polarization_array,'eP':polarization_error_array,'PA':polarization_angle_array,'ePA':theta_error_array})
    output_df = output_df.dropna()
    print(output_df.head)
    output_df.to_csv('data_files_generated/whole_region_plank_polarizationv2.csv')



def PlankPolarizationMaps(path_to_stokes_I,path_to_stokes_Q,path_to_Stokes_U,path_to_image_for_ploting,window_avg_length):
    """PolarizationMapCSV

    This function creates plank polarization maps from I,Q and U stokes files averages it and then plots it on a specified fits file.
    Note: This requires an old version of astropy comaptible with aplpy
    
    Args:
        path_to_stokes_I (string): path to I stokes files 
        path_to_stokes_Q (string): path to Q stokes files 
        path_to_Stokes_U (string): path to U stokes files
        window_avg_length (integer): length of the binning average 
        path_to_image_for_ploting (string): path to fits file on which plank polarization is to be plotted
    Returns:
        None: None
    """
    plank_I = FITS(path_to_stokes_I)
    plank_Q = FITS(path_to_stokes_Q)
    plank_U = FITS(path_to_Stokes_U)

    plank_I.generate_RA_DEC_mesh()
    plank_polarization = np.sqrt(plank_Q.data*plank_Q.data + plank_U.data*plank_U.data)/plank_I.data
    plank_theta_mod = (180/np.pi)*0.5*np.arctan(plank_U.data/plank_Q.data) + 90 + 62
    def window_avg(matrix,avg):
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

    DEC_grid_avg,_ = window_avg(plank_I.DEC_grid,window_avg_length)
    RA_grid_avg,_ = window_avg(plank_I.RA_grid,window_avg_length)
    polarization_avg , polarization_error = window_avg(plank_polarization,window_avg_length)
    theta_avg,theta_error = window_avg(plank_theta_mod,window_avg_length)
    
    RA_array = RA_grid_avg.flatten()
    DEC_array = DEC_grid_avg.flatten()
    polarization_array = polarization_avg.flatten()
    polarization_error_array = polarization_avg.flatten()
    polarization_angle_array = theta_avg.flatten()
    theta_error_array = theta_error.flatten()
    
    #plotting
    scuba = aplpy.FITSFigure(path_to_image_for_ploting)

    meanp = 125 * np.mean( polarization_array)
  
    p_lo= -plank_I.header['CDELT1'] #why is there a negatic sign
    p_la= plank_I.header['CDELT2']
    ra_ref= plank_I.header['CRVAL1']
    de_ref=plank_I.header['CRVAL2']
    xcen= plank_I.header['NAXIS1']
    ycen= plank_I.header['NAXIS2']
    xc      = xcen/2
    yc      = ycen/2

    p_loo1=p_lo/np.cos(DEC_array*0.01745)

    x2      = (ra_ref-(RA_array+( polarization_array/meanp)*np.sin(0.0175*polarization_angle_array )))/p_loo1+xc+10
    y2      = ((DEC_array+( polarization_array/meanp)*np.cos(0.0175*polarization_angle_array ))-de_ref)/p_la+yc-10
    x3      = (ra_ref-(RA_array-( polarization_array/meanp)*np.sin(0.0175*polarization_angle_array )))/p_loo1+xc+10
    y3      = ((DEC_array-( polarization_array/meanp)*np.cos(0.0175*polarization_angle_array ))-de_ref)/p_la+yc-10

    xx2, yy2=scuba.pixel2world(x2, y2)
    xx3, yy3=scuba.pixel2world(x3, y3)

    scuba.show_arrows(xx3,yy3,xx2-xx3,yy2-yy3, color='y', head_length=0,head_width=0, length_includes_head=False, width = 0.5)  # Polarization vector
    plt.show()

def mapping_func(grid_xx,grid_yy,x,y,f_x_y,tol = 0.001):
    """mapping_func

    This function maps f(x,y) onto (x,y) grid from a data array of x , y and f(x,y)
    
    Args:
        grid_xx (2D array): x coordinates mesh 
        grid-yy (2D array): y coordinates mesh 
        x (1D array): x coordinate value at which f is evaluated
        y (1D array): y coordinate value at which f is evaluated
        f_x_y (1D array): value of f evaluated at x and y
    Returns:
        mappped_func (2D array): f(x,y) mapped into (x,y)
    """
    mappped_func = np.zeros((grid_xx.shape[0],grid_xx.shape[1]))
    for i in range(f_x_y.shape[0]):
        temp = (abs(grid_xx - x[i])<tol)*(abs(grid_yy - y[i])<tol)*f_x_y[i]
        mappped_func += temp*(mappped_func==0)
    return mappped_func

def slicing(bottom_left_xy,top_right_xy,grid):
    """Gaussian_derivative

    This function slices a 2D array using the bottom left and top right coordinates
    
    Args:
        bottom_left_xy ([x,y] x,y integers): x,y coordinates of bottom left point of sliced 2D array
        top_right_xy ([x,y] x,y integers): x,y coordinates of top right point of sliced 2D array
        grid (2D array): 2D array to be sliced
    Returns:
        grid (2D array): 2D array to be sliced
    """
    return grid[bottom_left_xy[1]:top_right_xy[1],bottom_left_xy[0]:top_right_xy[0]]
    

def Image_mod(image,threshold):
    """Image_mod

    This function low pass filters image depending on a specified threshold
    
    Args:
        image (2D array): image to be sliced 
        threshold (float): maximum value of the modified image
    Returns:
        image_mod (2D array): modified image
    """
    image_mod = (image<threshold)
    image_mod = (0<image)*image_mod*50
    return image_mod

def create_image(data,threshold,bottom_xy,top_xy):
    """Gaussian_derivative

    This function creates and saves a png image created by slicing and modifying the data
    
    Args:
        data (2D array): data from which image is created
        threshold (float): maximum value of the modified image
        bottom_left_xy ([x,y] x,y integers): x,y coordinates of bottom left point of sliced 2D array
        top_right_xy ([x,y] x,y integers): x,y coordinates of top right point of sliced 2D array

    Returns:
        None: None
    """
    data_sliced = slicing(bottom_xy,top_xy,data)
    data_mod = Image_mod(threshold,data_sliced)
    matplotlib.image.imsave('name.png', data_mod)

def Gaussian_derivative(path_to_image,size,sigma):
    """Gaussian_derivative

    This function acts a gaussian derivative kernel on a image
    
    Args:
        path_to_image (string): path to the image on which Gaussian derivative kernel is acted 
        size (postive integer): size of the Gaussian derivative kernel
        sigma (float): sigma parameter of the Gausian derivative kernerl
    Returns:
        image_gradient_x (2D array): Gradient in the x direction 
        image_gradient_y (2D array) Gradient in the y direction
        image_gradient_total (2D array): Total gradient
        gradient_angle (2D array): Angle between gradients

    """
    image = cv2.imread(path_to_image) # load image
    image = image[:,:,0]
    def convolution(derivative,gaussian,image):
        """convolution

        This function calculates the gaussian derivative kernel then convolves it with the image
        
        Args:
            derivative (2D array): X/Y direction derivative kernel 
            gaussian (2D array): gaussian smoothing kernel
            image (2D array): image on which the gaussian derivative kernel acts
        Returns:
            image_convolution (2D array): image convolved with the gaussian derivative kernel
        """
        gaussian_derivative = cv2.filter2D(derivative, -1,gaussian)
        image_convolution = cv2.filter2D(image,-1,gaussian_derivative)
        return image_convolution

    def mod_arctan(x,y):
        """mod_arctan

        This function calculates the modified inverse tan taking care of singularities
        
        Args:
            x (float): x direction magnitude
            y (float): y direction magnitude
        Returns:
            inverse tan (float): angle between x and y
        """
        if x==0 and y==0:
            return 0
        else:
            return (180/np.pi)*np.arctan(x/y)

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
    image_gradient_total = np.sqrt(image_gradient_x**2 + image_gradient_y**2)
    
    x = image_gradient_x
    y = image_gradient_y

    gradient_angle = np.ones_like(image_gradient_x)
    for i in range(gradient_angle.shape[0]):
        for j in range(gradient_angle.shape[1]):
            gradient_angle[i,j] = mod_arctan(x[i,j],y[i,j])
    
    
    return image_gradient_x,image_gradient_y,image_gradient_total,gradient_angle
    
def HRO_analysis(RA_grid,DEC_grid,polarization_file,path_to_image):
    # this is function is in the test phase
    Polarization_data = pd.read_csv(polarization_file,delimiter=',')
    polarization_RA = Polarization_data['ra']
    polarization_DEC = Polarization_data['dec']
    Polarization_angle = Polarization_data['PA']
    
    x_gradient,y_gradient,total_gradient,gradient_angle = Gaussian_derivative(path_to_image,25,1)
    mask_gradient = (total_gradient>0)
    not_mask_gradient = (total_gradient==0)

    optical_pa_mod = Polarization_angle+360
    tol = abs(DEC_grid[1,0]-DEC_grid[0,0])*10
    optical_polarization_data_map = mapping_func(RA_grid,DEC_grid,polarization_RA,polarization_DEC,optical_pa_mod,tol)
    mask_polarization_map = (optical_polarization_data_map>0)
    not_mask_polarization_map = (optical_polarization_data_map==0)

    relative_orientation = (gradient_angle*mask_polarization_map - optical_polarization_data_map*mask_gradient + (not_mask_polarization_map+not_mask_gradient)*(-720))
    # wieghted_relative_orientation = relative_orientation*(total_gradient*mask_polarization_map)
    relative_orientation_mod = (relative_orientation[-360<relative_orientation] +360)

    return relative_orientation_mod