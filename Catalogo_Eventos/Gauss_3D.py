from functions_py import *
# import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
import pandas as pd 
# import cv2
import skimage as sk
from sympy import Ellipse, Point
import random

from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from scipy.stats import multivariate_normal
from matplotlib import cbook, cm
from matplotlib.colors import LightSource

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

path = '/home/labdet/Documents/MauSan/imagenesMicrochip/17OCT23/am241_gammas/proc_skp_m-009_microchip_T_170__Vv82_NSAMP_324_NROW_50_NCOL_700_EXPOSURE_0_NBINROW_1_NBINCOL_1_img_690.fits'
hdu_list = fits.open(path)
numero_bins = 5000

extension = 4
# plt.imshow(hdu_list[0].data[:, 550:], vmin = 0, vmax = 80000, origin = 'lower')
# Overscan = hdu_list[extension - 1].data[:, 550:]
# oScan_mask=sk.measure.label(Overscan>=20000, connectivity=2)
# oScan=ma.masked_array(Overscan,mask=(oScan_mask>0))
# plt.imshow(oScan, origin='lower')
# plt.colorbar()

# plt.figure(figsize=[20,15])
# active_area = hdu_list[extension - 1].data[:, :550]
# active_area_mask=sk.measure.label(active_area>=500000, connectivity=2)
# active_area_true=ma.masked_array(active_area,mask=(active_area_mask>0))
# plt.imshow(active_area_true, origin='lower')
# plt.colorbar(location = 'bottom')



Colors = ['b', 'r', 'b', 'r', 'b', 'r', 'b'] * 11


Overscan = hdu_list[extension -1].data[:, 550:]
oScan_mask=sk.measure.label(Overscan>= np.max(Overscan) - np.mean(Overscan), connectivity=2)
oScan = ma.masked_array(Overscan, mask=(oScan_mask>0))

# active_area = hdu_list[extension].data[:, :550]
# Overscan = hdu_list[extension - 1].data[:, 550:]

active_area = hdu_list[extension -1].data[:, :550]
active_area_mask = sk.measure.label(active_area>= np.max(active_area), connectivity=2)
data = ma.masked_array(active_area, mask=(active_area_mask>0))


header = hdu_list[extension - 1].header
nsamp = float(header['NSAMP'])

del header

hist , bins_edges = np.histogram(oScan.flatten(), bins = numero_bins) #'auto'

del oScan

offset = bins_edges[np.argmax(hist)]
dataP = data-offset
# dataCal = dataP * expgain[extension] ## En keV  
dataCal = dataP ## En ADUs

del hist
del data
del dataP

bin_heights, bin_borders = np.histogram(dataCal.flatten(), bins= numero_bins) #'auto'
bin_centers=np.zeros(len(bin_heights), dtype=float)
offset_fit = bin_borders[np.argmax(bin_heights)]
for p in range(len(bin_heights)):
    bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

# xmin_fit, xmax_fit = offset_fit-(10*expgain[extension])/math.sqrt(nsamp), offset_fit+(10*expgain[extension])/math.sqrt(nsamp)			# Define fit range
xmin_fit, xmax_fit = -abs(offset), abs(offset)

bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

del nsamp

try:
    popt,_ = curve_fit(gaussian, bin_centers, bin_heights, maxfev=100000, p0 = [10, 10, 500])		# Fit histogram with gaussian
    fondo_value = 3 * abs(popt[2])

except:
    print('Fit error in image ')

del bin_heights
del bin_centers
del offset
del xmin_fit
del xmax_fit

label_img, nlabels_img = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
prop = sk.measure.regionprops(label_img, dataCal)
# list_totalEvents.append(n_events)
# print(nlabels_img)
# list_labels.append(label_img)
# list_EventsNumber.append(n_events)

## Obteniendo el valor promedio del fondo
fondo_mask = np.invert(label_img == 0)
fondo = ma.masked_array(dataCal,fondo_mask)
valor_promedio_fondo = fondo.mean()

Gammas_Events = []

for event in range(1, nlabels_img):
    mask = np.invert(label_img == event)
    loc = ndimage.find_objects(label_img == event)[0]
    
    data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                         mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

    coordX_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[1])
    coordY_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[0])

    # MaxValue_Event = data_maskEvent.max()
    MinValue_Event = data_maskEvent.min()
    MeanValue_Event = data_maskEvent.mean()
    # MeanValue_Event = (MaxValue_Event - MinValue_Event)/2
    Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]
    try:
        differval = abs(Barycentercharge - MinValue_Event) 
    except:
        differval = 0 

    rM = prop[event-1].axis_major_length
    rm = prop[event-1].axis_minor_length
    Solidity = prop[event-1].solidity
    charge = data_maskEvent.sum()
    miny, minx, maxy, maxx = prop[event-1].bbox
    Long_y = maxy - miny
    Long_x = maxx - minx 

    if Long_x < 4 or Long_y < 4 :
        continue

    if differval < MeanValue_Event:
        continue

    if Solidity < 0.8:
        continue 

    if  rM <= 1.5 * rm:
        # print(Barycentercharge, differval)
        Gammas_Events.append(event)

print('Gammas Detected: ', len(Gammas_Events))
# print(Gammas_Events)


# ficticio = np.array([[2,2,2,2,2,2,2,2,2,2,2],
#                     [2,2,2,2,2,2,2,2,2,2,2],
#                     [2,2,2,2,10,10,10,2,2,2,2],
#                     [2,2,2,2,10,40,10,2,2,2,2],
#                     [2,2,2,2,10,10,10,2,2,2,2], 
#                     [2,2,2,2,2,2,2,2,2,2,2], 
#                     [2,2,2,2,2,2,2,2,2,2,2]])

# ficticio = np.array([[0,0,0,0,0,0], [0,1,1,1,1,0],[0,1,0,0,1,0],[0,1,0,0,1,0], [0,1,0,0,1,0], [0,1,1,1,1,0]])
ficticio = active_area

# x, y = np.mgrid[-1.0:1.0:30j, -1.0:1.0:30j]
# xy = np.column_stack([x.flat, y.flat])
# # print(xy)
# mu = np.array([0.0, 0.0])

# sigma = np.array([.5, .5])
# covariance = np.diag(sigma**2)

# z = multivariate_normal.pdf(xy, mean=mu, cov=covariance)

# # Reshape back to a (30, 30) grid.
# z = z.reshape(x.shape)





# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(x,y,z)
# #ax.plot_wireframe(x,y,z)

# plt.show()

for index_event in range(0, len(Gammas_Events)):
    mask = np.invert(label_img == Gammas_Events[index_event])
    loc = ndimage.find_objects(label_img == Gammas_Events[index_event])[0]
    
    data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                         mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
    
    data_array = np.array(data_maskEvent)
    #
    # Create a figure for plotting the data as a 3D histogram.
    #
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #
    # Create an X-Y mesh of the same dimension as the 2D data. You can
    # think of this as the floor of the plot.
    #
    x_data, y_data = np.meshgrid( np.arange(data_array.shape[1]),
                                np.arange(data_array.shape[0]) )

    # print(data_array.shape[1])
    # print(data_array.shape[0])

    # Flatten out the arrays so that they may be passed to "ax.bar3d".
    # Basically, ax.bar3d expects three one-dimensional arrays:
    # x_data, y_data, z_data. The following call boils down to picking
    # one entry from each array and plotting a bar to from
    # (x_data[i], y_data[i], 0) to (x_data[i], y_data[i], z_data[i]).

    x_data = x_data.flatten()
    y_data = y_data.flatten()
    z_data = data_array.flatten()

    ls = LightSource(270, 45)
    # To use a custom hillshading mode, override the built-in shading and pass
    # in the rgb colors of the shaded surface calculated from "shade".
    # rgb = ls.shade(z_data, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')

    hist3D = ax.bar3d( x_data,
            y_data,
            np.zeros(len(z_data)), 1, 1, z_data, 'b')

    hist3D.set_facecolors('b')
    hist3D.set_edgecolors('k')
    # hist3D.set_cmap('hot')
    # hist3D.set_clim(vmin = 0, vmax = 100)

    #
    # Finally, display the plot.
    #
    plt.show()


### Revisar (paquetería de gráficas): PyqtGraph
### Revisar (multiprocesos): JobLib
