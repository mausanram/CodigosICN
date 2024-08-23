from functions_py import *
# import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
import pandas as pd 
import skimage as sk


path = '/home/labdet/Documents/MauSan/imagenesMicrochip/01AGO23/proc_skp_m-009_microchip_T_170__Vv82_NSAMP_1_NROW_650_NCOL_700_EXPOSURE_1200_NBINROW_1_NBINCOL_1_img_325.fits'
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




"""
This example demonstrates the use of GLBarGraphItem.

"""

import numpy as np

import pyqtgraph as pg
import pyqtgraph.opengl as gl


# ficticio = np.array([[2,2,2,2,2,2,2,2,2,2,2],
#                     [2,2,2,2,2,2,2,2,2,2,2],
#                     [2,2,2,2,10,10,10,2,2,2,2],
#                     [2,2,2,2,10,20,10,2,2,2,2],
#                     [2,2,2,2,10,10,10,2,2,2,2], 
#                     [2,2,2,2,2,2,2,2,2,2,2], 
#                     [2,2,2,2,2,2,2,2,2,2,2]])

shape_y, shape_x = oScan.shape
shape_y, shape_x = data.shape

app = pg.mkQApp("GLBarGraphItem Example")
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLBarGraphItem')
w.setCameraPosition( distance=300, elevation=50, azimuth=90)

cm = pg.colormap.get('CET-L17')
cm.reverse() 
pen = cm.getPen( span=(0.0,1.0), width=5 )

w.pan(shape_y/2, shape_x/2,0)

gx = gl.GLGridItem()
gx.rotate(90, 0, 1, 0)
gx.translate(-10, 0, 10)
w.addItem(gx)
gy = gl.GLGridItem()
gy.rotate(90, 1, 0, 0)
gy.translate(0, -10, 10)
w.addItem(gy)
gz = gl.GLGridItem()
gz.translate(0, 0, 0)
w.addItem(gz)

# regular grid of starting positions
pos = np.mgrid[0:shape_y, 0:shape_x, 0:1].reshape(3,shape_y, shape_x).transpose(1,2,0)
# fixed widths, random heights
size = np.empty((shape_y, shape_x, 3))
size[...,0:2] = 1
# size[...,2] = oScan * 0.001
size[...,2] = data * 0.001

bg = gl.GLBarGraphItem(pos, size)
# bg.setColor(')
w.addItem(bg)


if __name__ == '__main__':
    pg.exec()
