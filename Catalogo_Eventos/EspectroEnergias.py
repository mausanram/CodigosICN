# from functions_py import math
import math
from astropy.io import fits
import scipy.ndimage as ndimage
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys
import skimage as sk
import datetime

numero_bins = 10000
def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7]
    list_EventCharge_AllExtensions=[]
    Inicio = datetime.datetime.now()
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        hdu_list = fits.open(img)
        
        for extension in (0,1,3):
            StraightEvents_list = []

            data = hdu_list[extension].data
            header = hdu_list[extension].header
            oScan = hdu_list[extension].data[638:,530:]
            nsamp = float(header['NSAMP'])

            del header

            hist , bins_edges = np.histogram(oScan.flatten(), bins = numero_bins) #'auto'

            del oScan

            offset = bins_edges[np.argmax(hist)]
            dataP = data-offset
            dataCal = dataP/expgain[extension] ## En electrones  
            
            del hist
            del data
            del dataP
            
            bin_heights, bin_borders = np.histogram(dataCal.flatten(), bins= numero_bins) #'auto'
            bin_centers=np.zeros(len(bin_heights), dtype=float)
            offset_fit = bin_borders[np.argmax(bin_heights)]
            for p in range(len(bin_heights)):
                bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

            xmin_fit, xmax_fit = offset_fit-(10*expgain[extension])/math.sqrt(nsamp), offset_fit+(10*expgain[extension])/math.sqrt(nsamp)			# Define fit range
            bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
            bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

            del nsamp

            try:
                popt,_ = curve_fit(gaussian, bin_centers, bin_heights)#, p0=[np.max(bin_heights), 0, 1], maxfev=100000)		# Fit histogram with gaussian

            except:
                continue
            del bin_heights
            del bin_centers
            del offset
            del xmin_fit
            del xmax_fit

            # label, n_events = ndimage.label(dataCal>6*abs(popt[2]),structure=[[1,1,1],[1,1,1],[1,1,1]])
            label_img, n_events = sk.measure.label(dataCal>6*popt[2], connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            
            # print(nlabels_img)
            # list_labels.append(label_img)
            # list_EventsNumber.append(n_events)
            
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)
            valor_promedio_fondo = fondo.data.mean()

            for event in range(1,n_events):
                mask = np.invert(label_img == event)
                loc = ndimage.find_objects(label_img == event)[0]
                data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                         mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
                
                # del dataCal

                coordX_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[1])
                coordY_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[0])

                MeanValue_Event = data_maskEvent.mean()
                MinValue_Event = data_maskEvent.min()

                Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]

                try:
                    differval = abs(Barycentercharge - MinValue_Event) 
                except:
                    differval = 0 

                rM = prop[event-1].axis_major_length
                rm = prop[event-1].axis_minor_length

                minr, minc, maxr, maxc = prop[event-1].bbox

                if rM == 0 or rm == 0:
                    continue 

                elif maxc - minc <= 3:
                    continue

                elif not Barycentercharge:
                    continue

                elif  rM < 4.5 * rm:
                    continue

                elif differval < MeanValue_Event + 100:
                    continue
                
                charge = 0
                for element in data_maskEvent.data.flatten():
                    if element >= valor_promedio_fondo:
                        charge = charge + element

                # list_charge.append(charge)
                list_EventCharge_AllExtensions.append(charge)

                del data_maskEvent
                del Barycentercharge
                del charge

            # print(list_eventos_rectos)

        del hdu_list            

    Final = datetime.datetime.now()
    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    plt.hist(list_EventCharge_AllExtensions, bins = numero_bins)    
    plt.title('Espectro de Energía de Muones')
    plt.xlabel(r'e⁻')
    plt.ylabel('Cuentas') 
    plt.xlim(0, 200000)  
    plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)