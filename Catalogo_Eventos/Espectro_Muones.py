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
import os


current_path = os.getcwd()

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
    # expgain = [227, 220.4, 94.72, 197.7] ## Ganancia para 324nsamp 

    # expgain = [110.78158608889959, 144.4661840232508, 100,  63.730700893071976] ## 


    list_EventCharge_extension_2 =[]
    list_EventCharge_extension_1 = []
    list_EventCharge_extension_4 = []

    list_totalEvents = []
    list_DeltaEL_extension_2 = []
    list_DeltaEL_extension_1 = []
    list_DeltaEL_extension_4 = []

    list_DeltaL_extension_2 = []
    list_DeltaL_extension_1 = []
    list_DeltaL_extension_4 = []


    total_images = len(argObj)
    image_in_bucle = 0

    Inicio = datetime.datetime.now()
    num_images =  'Imágenes Analizadas: ' +  str(total_images)
    num_muons = 0
    
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        hdu_list = fits.open(img)
        image_in_bucle += 1
        
        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end = '\r')
        
        for extension in (0,1,3):
            # extension = 1

            data = hdu_list[extension].data[:,:550]
            header = hdu_list[extension].header
            oScan = hdu_list[extension].data[:,550:]
            nsamp = float(header['NSAMP'])

            del header

            hist , bins_edges = np.histogram(oScan.flatten(), bins = numero_bins) #'auto'

            del oScan

            offset = bins_edges[np.argmax(hist)]
            dataP = data - offset
            # dataCal = (dataP)/expgain[extension] ## En keV  
            dataCal = dataP ## En ADUs
            
            del hist
            del data
            del dataP
            
            bin_heights, bin_borders = np.histogram(dataCal.flatten(), bins= numero_bins) #'auto'
            bin_centers=np.zeros(len(bin_heights), dtype=float)
            offset_fit = bin_borders[np.argmax(bin_heights)]
            for p in range(len(bin_heights)):
                bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

            xmin_fit, xmax_fit = -abs(offset), abs(offset)			# Define fit range
            bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
            bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

            del nsamp

            try:
                popt,_ = curve_fit(gaussian, bin_centers, bin_heights, maxfev=1000, p0 = [1,1,100])		# Fit histogram with gaussian
                fondo_value = 6 * abs(popt[2])
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
                Solidity = prop[event-1].solidity

                miny, minx, maxy, maxx = prop[event-1].bbox
                Longitud_y, Longitud_x = maxy - miny , maxx - minx # px
                Diagonal_lenght= np.sqrt(Longitud_x**2 + Longitud_y**2) - np.sqrt(2) # px
                Delta_L = np.sqrt( (Diagonal_lenght * px_to_micras)**2 + (CCD_depth)**2) * micra_to_cm # cm

                if rM == 0 or rm == 0:
                    continue 

                elif maxx - minx <= 3:
                    continue

                elif not Barycentercharge:
                    continue
                
                # if differval < MeanValue_Event: #ADUs
                #     continue
                
                elif  Solidity < Solidit: ## Ver el artículo: Recognition and classification of the cosmic-ray events in images captured by CMOS/CCD cameras
                    continue 

                if  rM >= Elip * rm: ## Muones
                    charge = data_maskEvent.sum()
                    # DeltaEL = (charge / 1000) / (Delta_L * micra_to_cm) # MeV/cm
                    DeltaEL = charge / (Delta_L) # ADUs/cm
                    num_muons = num_muons + 1

                    if extension == 0:
                        list_DeltaEL_extension_1.append(DeltaEL)
                        list_EventCharge_extension_1.append(charge)
                        list_DeltaL_extension_1.append(Delta_L)

                    elif extension == 1:
                        list_DeltaEL_extension_2.append(DeltaEL)
                        list_EventCharge_extension_2.append(charge)
                        list_DeltaL_extension_2.append(Delta_L)

                    elif extension == 3:
                        list_DeltaEL_extension_4.append(DeltaEL)
                        list_EventCharge_extension_4.append(charge)
                        list_DeltaL_extension_4.append(Delta_L)

                del data_maskEvent
                del Barycentercharge

        del hdu_list            

    # list_EventCharge_AllExtensions = list_EventCharge_extension_2 + list_EventCharge_extension_1 + list_EventCharge_extension_4

    # dict_to_save_pkl = {'Muons_Detected' : num_muons, 'charge' : list_EventCharge_AllExtensions, 'DeltaEL' : list_DeltaEL}

    dict_to_save_pkl = {'All_Muons_Detected' : num_muons, 
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'deltaEL' : list_DeltaEL_extension_1, 'deltaL' : list_DeltaL_extension_1}, 
                        'extension_2' : {'charge' : list_EventCharge_extension_2, 'deltaEL' : list_DeltaEL_extension_2, 'deltaL' : list_DeltaL_extension_2},
                        'extension_4' : {'charge' : list_EventCharge_extension_4, 'deltaEL' : list_DeltaEL_extension_4, 'deltaL' : list_DeltaL_extension_4}}

    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    eventos_rectos = 'Muones Rectos Detectados: ' + str(num_muons)
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    print(Eventos_Totales)
    print(eventos_rectos)
    # print(eventos_circulares)

    # fig, axs = plt.subplots(1,1)
    # fig.canvas.manager.set_window_title('histogram_Imgs_'+str(len(argObj))+'_Sol_'+str(Solidit)+'_Elip_'+str(Elip))
    # fig.suptitle('Espectro de Energía de Muones')

    # bin_heights, bin_borders, _ = axs.hist(list_EventCharge_AllExtensions, bins = numero_bins, label= num_images + '\n' + eventos_rectos) 
    # axs.hist(list_EventCharge_extension_2, bins = numero_bins, label= 'Extension 2') 
    # axs.hist(list_EventCharge_extension_1, bins = numero_bins, label= 'Extension 1') 
    # axs.hist(list_EventCharge_extension_4, bins = numero_bins, label= 'Extension 4') 
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
    # axs.set_xlabel(r'Energy (ADUs)')
    # axs.set_ylabel('Events') 
    # axs.set_xlim([0, 400])  

    file_name = 'dict_muons_Extensions_1_to_4_Imgs_' + str(len(argObj))+'_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_ADUs__.pkl'
    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in ', current_path + '/' + file_name, ' as a binary file.')
    print('To open use library "pickle". ')

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

