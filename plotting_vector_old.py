import numpy as np
import matplotlib.pyplot as plt
import aplpy
import pandas as pd
from astropy.io import fits


FITs_file = 'FITS_files/Optical/Halp_rsubt_L328_323_331_new_header.fits'
datafile= "data_files_generated/whole_region_plank_polarizationv5.csv"

scuba = aplpy.FITSFigure(FITs_file)
scuba.ticks.set_xspacing(0.1) # number in degrees
# scuba.show_contour(FITs_file,colors='yellow',kernel='gauss',smooth=1,levels=10,overlap=False,linewidths=0.3,alpha=1)
scuba.show_grayscale()
# scuba.show_colorscale(cmap='default')
# scuba.set_theme('publication')

scuba.ticks.show()
scuba.ticks.show_x()
scuba.ticks.show_y()
# scuba.tick_labels.set_xformat('hh:mm:ss.s')
# scuba.tick_labels.set_yformat('dd:mm:ss.s')


# scuba.add_scalebar(0.02, "0.02 degrees", color='w')
# scuba.scalebar.set_label("0.05 pc")


# scuba.add_scalebar(0.0057, color='k', corner='top left',lw=2)
# scuba.scalebar.set_label("10 %")
# scuba.scalebar.set_font(size='xx-large')


# scuba.add_colorbar()
# scuba.colorbar.set_width(0.25)
# scuba.colorbar.set_location('right')
# scuba.colorbar.set_axis_label_font(size=14, weight='medium')
# # scuba.colorbar.set_axis_label_text('$\\rm Intensity_{850\mu m}$(mJy/beam)')
# # scuba.colorbar.set_ticks([0.0,20,40,60,80,100,120,140,160,180,200])




df = pd.read_csv(datafile, delimiter=',')
rat = np.array(df['ra'])
det = np.array(df['dec'])
prt = np.array(df['P'])
eprt = np.array(df['eP'])
pa = np.array(df['PA'])
epart = np.array(df['ePA'])
prt1 = (prt/prt)*100
prt = prt1
part=pa
#prt=10
meanp = 100 * np.mean(prt)
# meanp = 1000

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


# # plt.text(150,110,'10%', fontsize = 10,color='w')
plt.show()


# # scuba.save('POL2SN1.5_vector.png',dpi=150)
# # scuba.save('POL2SN1.5_vector_Norm.png',dpi=150)



# #-----------------save plot--------------------
# #plt.savefig('g34_POL2.eps',format='eps',dpi=300)
# #plt.savefig('vector_plot_POL2_normalised.eps',format='eps',dpi=300)
