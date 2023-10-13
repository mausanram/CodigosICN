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
import pickle

ratio_keV = 0.0037
Solidit = 0.7
Elip = 2
numero_bins = 5000

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def Landau(x,a, MP,xi):
    C1 = a/np.sqrt((2 * np.pi))
    C2 = np.exp(-(x-MP)/xi)
    C3 = np.exp((-0.5 * ((x-MP)/xi + C2 )))
    return  C1 * C3


def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7]
    list_EventCharge_AllExtensions=[]
    list_EventosRectos = []
    list_EventosCirc = []
    list_totalEvents = []
    Inicio = datetime.datetime.now()
    num_images =  'Imágenes Analizadas: ' +  str(len(argObj))
    num_event = 0
    
    
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        hdu_list = fits.open(img)
        
        for extension in (0,1,3):

            data = hdu_list[extension].data
            header = hdu_list[extension].header
            oScan = hdu_list[extension].data[638:,530:]
            nsamp = float(header['NSAMP'])

            del header

            hist , bins_edges = np.histogram(oScan.flatten(), bins = numero_bins) #'auto'

            del oScan

            offset = bins_edges[np.argmax(hist)]
            dataP = data-offset
            dataCal = (ratio_keV * dataP)/expgain[extension] ## En keV  
            
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
            list_totalEvents.append(n_events)
            # print(nlabels_img)
            # list_labels.append(label_img)
            # list_EventsNumber.append(n_events)
            
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)

            for event in range(1,n_events):
                mask = np.invert(label_img == event)
                loc = ndimage.find_objects(label_img == event)[0]
                data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                         mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

                minr, minc, maxr, maxc = prop[event-1].bbox

                if maxc - minc <= 3:
                    continue

                num_event = num_event + 1
                charge = data_maskEvent.sum()
                list_EventCharge_AllExtensions.append(charge)

                del data_maskEvent

        del hdu_list            
    
    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    print(Eventos_Totales)

    dict_to_save_pkl = {'Events_Detected' : num_event, 'charge' : list_EventCharge_AllExtensions}

    total_events = sum(list_totalEvents)
    fig, axs = plt.subplots(1,1)
    fig.canvas.manager.set_window_title('histogram_Imgs_'+str(len(argObj))+'_Sol_'+str(Solidit)+'_Elip_'+str(Elip)+'.pkl')
    
    fig.suptitle('Espectro de Energía de Muones')

    bin_heights, bin_borders, _ = axs.hist(list_EventCharge_AllExtensions, bins = numero_bins, label= num_images + '\n' + Eventos_Totales)  
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2 

    axs.legend(loc="upper right") 
    axs.set_xlabel(r'keV')
    axs.set_ylabel('Cuentas') 
    axs.set_xlim([0, 400])  

    file_object_histogram = open('data_AlEvents_Imgs_'+str(len(argObj))+'_.pkl', 'wb')
    pickle.dump(list_EventCharge_AllExtensions, file_object_histogram)
    file_object_histogram.close()

    plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
