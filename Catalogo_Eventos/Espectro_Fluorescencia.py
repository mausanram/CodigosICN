from functions_py import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
import pandas as pd 
# import cv2
import skimage as sk
from sympy import Ellipse, Point
import pickle

ratio_keV = 0.0037
CCD_depth = 725 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000
DeltaEL_range = 85
Solidit = 0.7
Elip = 4.5
numero_bins = 5000
DeltaEL_range_min, DeltaEL_range_max = 0.9, 3.55

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def Landau(x,a, MP,xi):
    C1 = a/np.sqrt((2 * np.pi))
    C2 = np.exp(-(x-MP)/xi)
    C3 = np.exp((-0.5 * ((x-MP)/xi + C2 )))
    return  C1 * C3

def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7]
    Fluorescence_Events = []
    list_totalEvents = []
    Inicio = datetime.datetime.now()
    num_images =  'Imágenes Analizadas: ' +  str(len(argObj))
    num_EventosRectos = 0
    
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        hdu_list = fits.open(img)
        for extension in [0,1]:
            # print('Estoy en las extensiones')
            data = hdu_list[extension].data
            header = hdu_list[extension].header
            oScan = hdu_list[extension].data[:,530:]
            nsamp = float(header['NSAMP'])

            del header

            hist , bins_edges = np.histogram(oScan.flatten(), bins = numero_bins) #'auto'

            del oScan

            offset = bins_edges[np.argmax(hist)]
            dataP = data - offset
            dataCal = (ratio_keV * dataP)/expgain[extension] ## En keV  
            # print('Ya aplané los datos')
            
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
                popt,_ = curve_fit(gaussian, bin_centers, bin_heights, maxfev=100000)		# Fit histogram with gaussian

            except:
                print('Error in image ' + str(img))
                continue

            del bin_heights
            del bin_centers
            del offset
            del xmin_fit
            del xmax_fit

            fondo_value = 6 * abs(popt[2])
            # print('Ya saqué el valor del fondo')
            # label, n_events = ndimage.label(dataCal>6*abs(popt[2]),structure=[[1,1,1],[1,1,1],[1,1,1]])
            label_img, n_events = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            list_totalEvents.append(n_events)
            # print(nlabels_img)
            # list_labels.append(label_img)
            # list_EventsNumber.append(n_events)
            
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)
            valor_promedio_fondo = fondo.data.mean()
            # print('Ya estoy en los labels')
            for event in range(1,n_events):
                # print('Estoy en los eventos')
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
                # Solidity = prop[event-1].solidity

                miny, minx, maxy, maxx = prop[event-1].bbox
                Longitud_y, Longitud_x = maxy - miny , maxx - minx # px
                # Diagonal_lenght= np.sqrt(Longitud_x**2 + Longitud_y**2) - np.sqrt(2) # px
                # Delta_L = np.sqrt( (Diagonal_lenght * px_to_micras)**2 + (CCD_depth)**2) # micras

                # if Longitud_y > 5 or Longitud_x > 5:
                #     continue
                
                if differval < MeanValue_Event + 0.4:
                    continue
                
                if  rM < 8 and rm < 8:
                    charge =  data_maskEvent.sum()
                    # if charge < 1:
                    Fluorescence_Events.append(charge)

                del data_maskEvent
                del Barycentercharge

        del hdu_list    

    list_All = [Fluorescence_Events]  
    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()
    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)

    fig, axs = plt.subplots(1,1)
    fig.canvas.manager.set_window_title('Histograma de eventos de fluorescencia del silicio')
    fig.suptitle('Espectro de Energía Eventos de baja energía')

    bin_heights, bin_borders, _ = axs.hist(Fluorescence_Events, bins = numero_bins, label= num_images) 
    # axs.hist(list_EventCharge_AllExtensions, bins = numero_bins, label= num_images + '\n' + eventos_rectos) 
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2      

    axs.legend(loc="upper right") 
    axs.set_xlabel(r'Energy (keV)')
    axs.set_ylabel('Events') 
    # axs.set_xlim([0, 400000])  

    file_object_histogram = open('data_fluorescense_hist_img_' + str(len(argObj)) + '_.pkl', 'wb')
    pickle.dump(list_All, file_object_histogram)
    file_object_histogram.close()

    plt.show() 

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)