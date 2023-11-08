# from functions_py import math
import math
import os
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
numero_bins = 5000

current_path = os.getcwd()

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def Landau(x,a, MP,xi):
    C1 = a/np.sqrt((2 * np.pi))
    C2 = np.exp(-(x-MP)/xi)
    C3 = np.exp((-0.5 * ((x-MP)/xi + C2 )))
    return  C1 * C3


def main(argObj):
    # expgain = [227, 220.4, 94.72, 197.7]
    # expgain = [0.004251274073061228, 0.003830082649184408, 0.003, 0.0034160708295529487] ## KeV/ADUs
    list_EventCharge_extension_1 =[]
    list_EventCharge_extension_2 =[]
    list_EventCharge_extension_4 = []
    list_totalEvents = []

    Inicio = datetime.datetime.now()

    total_images = len(argObj)
    image_in_bucle = 1

    num_images =  'Imágenes Analizadas: ' +  str(total_images)
    num_gammas = 0
    
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        hdu_list = fits.open(img)
        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end = '\r')
        image_in_bucle += 1

        for extension in (0,1,3):
            Overscan = hdu_list[extension].data[:, 550:]
            oScan_mask=sk.measure.label(Overscan>= np.max(Overscan) - np.mean(Overscan), connectivity=2)
            oScan=ma.masked_array(Overscan,mask=(oScan_mask>0))

            # active_area = hdu_list[extension].data[:, :550]
            # Overscan = hdu_list[extension - 1].data[:, 550:]

            active_area = hdu_list[extension].data[:, :550]
            active_area_mask = sk.measure.label(active_area>= np.max(active_area), connectivity=2)
            data = ma.masked_array(active_area,mask=(active_area_mask>0))


            header = hdu_list[extension].header
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
                popt,_ = curve_fit(gaussian, bin_centers, bin_heights, maxfev=200, p0 = [10, 10, 500])		# Fit histogram with gaussian
                fondo_value = 4 * abs(popt[2])
            except:
                print('Fit error in image ' + str(img))
                continue
            del bin_heights
            del bin_centers
            del offset
            del xmin_fit
            del xmax_fit

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

            for event in range(1,n_events):
                mask = np.invert(label_img == event)
                loc = ndimage.find_objects(label_img == event)[0]
                data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
                
                # del dataCal

                # coordX_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[1])
                # coordY_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[0])

                # MeanValue_Event = data_maskEvent.mean()
                # MinValue_Event = data_maskEvent.min()

                # Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]

                # try:
                #     differval = abs(Barycentercharge - MinValue_Event) 
                # except:
                #     differval = 0 

                rM = prop[event-1].axis_major_length
                rm = prop[event-1].axis_minor_length

                # miny, minx, maxy, maxx = prop[event-1].bbox
                
                # if differval < MeanValue_Event: #keV
                #     continue
        
                if  rM <= 2 * rm:
                    charge = data_maskEvent.sum()
                    num_gammas = num_gammas + 1

                    if extension == 0:
                        list_EventCharge_extension_1.append(charge)

                    elif extension == 1:
                        list_EventCharge_extension_2.append(charge)

                    elif extension == 3:
                        list_EventCharge_extension_4.append(charge)

                    del data_maskEvent
                    # del Barycentercharge

        del hdu_list            

    dict_to_save_pkl = {'All_Gammas_Detected' : num_gammas, 
                        'extension_1' : {'charge' : list_EventCharge_extension_1}, 
                        'extension_2' : {'charge' : list_EventCharge_extension_2},
                        'extension_4' : {'charge' : list_EventCharge_extension_4} }


    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    gammas = 'Gammas Detectadas: ' + str(num_gammas)
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    print(Eventos_Totales)
    print(gammas)
    # print(eventos_circulares)

    # fig, axs = plt.subplots(1,1)
    # fig.canvas.manager.set_window_title('Histogram of Am-241')
    
    # fig.suptitle('Energy Spectrum of Am-241')

    # axs.hist(list_EventCharge_extension_2, bins = 1000, label= 'Extension 2') 
    # axs.hist(list_EventCharge_extension_1, bins = 1000, label= 'Extension 1') 
    # axs.hist(list_EventCharge_extension_4, bins = 1000, label= 'Extension 4') 

    # axs.hist(list_EventCharge_AllExtensions, bins = numero_bins, label= num_images + '\n' + eventos_rectos) 
    # bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2 


    # # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
    # # xmin_fit, xmax_fit = min(list_EventCharge_AllExtensions), max(list_EventCharge_AllExtensions)
    # # bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
    # # bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]


    # popt, _ = curve_fit(Landau, bin_centers, bin_heights) # p0=[50,10000,1000])
    # print('Los parámetros del ajuste son: Altura: ', popt[0], ' EMP: ', popt[1], r' $\xi$: ',popt[2])
    # strAjuste = 'Ajuste: EMP: ' + str(popt[1]) + r' $\xi$: ' + str(popt[2])
    # x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
    # axs.plot(x_interval_for_fit, Landau(x_interval_for_fit,*popt), label=strAjuste)

    # axs.legend(loc="upper right") 
    # axs.set_xlabel(r'Energy (keV)')
    # axs.set_ylabel('Events') 
    # axs.set_xlim([0, 400])  

    file_name = 'dict_Am241_Imgs_'+ str(len(argObj)) +'.pkl'
    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in ', current_path + '/' + file_name, ' as a binary file.')
    print('To open use library "pickle". ')

    plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)