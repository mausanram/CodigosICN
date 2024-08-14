import math
from functions_MuonsNSAMP1 import *
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

## CONSTANTES ## 
current_path = os.getcwd()
units = 0



n_sigmas = 5
ratio_keV = 0.0037
CCD_depth = 725 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000
DeltaEL_range = 85
Solidit = 0.7
elip = 0.9
numero_bins = 1000

## DEFINICIONES ##
def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))


## FUNCIÓN PRINCIPAL ##
def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7] ## Ganancia para 324nsamp 

    # expgain = [110.78158608889959, 144.4661840232508, 100,  63.730700893071976] ## 


    list_EventCharge_extension_2 =[]
    list_EventCharge_extension_1 = []
    list_EventCharge_extension_4 = []

    list_totalEvents = []
    list_vertical_event_extension_1 = []
    list_vertical_event_extension_2 = []
    list_vertical_event_extension_4 = []

    list_horizontal_event_extension_1 = []
    list_horizontal_event_extension_2 = []
    list_horizontal_event_extension_4 = []


    total_images = len(argObj)
    image_in_bucle = 0

    Inicio = datetime.datetime.now()
    num_images =  'Imágenes Analizadas: ' +  str(total_images)
    num_muons = 0
    
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        try:
            hdu_list = fits.open(img)
            image_in_bucle += 1
        except:
                print('Loading error in image ' + str(img) + 'in open the image.')
                continue
        
        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end = '\r')
        
        for extension in (0,1,3):
            # extension = 1
            try :
                data = hdu_list[extension].data[:,:550]
                header = hdu_list[extension].header
                oScan = hdu_list[extension].data[:,550:]
                nsamp = float(header['NSAMP'])
            except:
                print('Loading error in image ' + str(img) + 'in load the data.')
                continue

            del header

            try:
                dict_popt = oScan_fit_NSAMP1(extensión=extension, active_area=data, oScan=oScan, Bins=numero_bins, make_figure_flag=False)

            except:
                print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                continue
                
            sig_ADUs = dict_popt['sigma']
            Offset = dict_popt['Offset']

            dataCal = data_calibrated(active_area=data, extension=extension, list_gain=expgain, ratio_keV=ratio_keV, unidades= units, offset=Offset)

            if units == 0:
                fondo_value = n_sigmas * sig_ADUs
            elif units == 1:
                fondo_value = n_sigmas * sig_electrons
            elif units == 2:
                fondo_value = n_sigmas * sig_KeV
            
            del oScan


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
                
                MinValue_Event = data_maskEvent.min()
                MeanValue_Event = data_maskEvent.mean()

                try: 
                    coordX_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[1])
                    coordY_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[0])
                    Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]

                    differval = abs(Barycentercharge - MinValue_Event) 

                except:
                    Barycentercharge = np.nan()
                    differval = 0

                rM = prop[event-1].axis_major_length/2
                rm = prop[event-1].axis_minor_length/2
                
                try:
                    Elipticity = (rM - rm)/rM 
                except: 
                    Elipticity = 0


                Solidity = prop[event-1].solidity
                miny, minx, maxy, maxx = prop[event-1].bbox
                Longitud_y, Longitud_x = maxy - miny , maxx - minx # px

                if rM == 0 or rm == 0:
                    continue 

                elif maxx - minx <= 3:
                    continue

                elif not Barycentercharge:
                    continue

                elif differval < MeanValue_Event: #keV
                    continue

                elif  Solidity < Solidit:
                    continue 

                elif  elip >= Elipticity :
                    charge = data_maskEvent.sum()

                    if charge < 1000:
                        continue
                    
                    if (Longitud_x < 9 and Longitud_y > 10): 
                        num_muons = num_muons + 1

                        if extension == 0:
                            list_vertical_event_extension_1.append(data_maskEvent)
                            list_EventCharge_extension_1.append(charge)

                        elif extension == 1:
                            list_vertical_event_extension_2.append(data_maskEvent)
                            list_EventCharge_extension_2.append(charge)

                        elif extension == 3:
                            list_vertical_event_extension_4.append(data_maskEvent)
                            list_EventCharge_extension_4.append(charge)

                    if ( Longitud_y < 9 and Longitud_x > 10):
                        num_muons = num_muons + 1
                        charge = data_maskEvent.sum()

                        if extension == 0:
                            list_horizontal_event_extension_1.append(data_maskEvent)
                            list_EventCharge_extension_1.append(charge)

                        elif extension == 1:
                            list_horizontal_event_extension_2.append(data_maskEvent)
                            list_EventCharge_extension_2.append(charge)

                        elif extension == 3:
                            list_horizontal_event_extension_4.append(data_maskEvent)
                            list_EventCharge_extension_4.append(charge)


                del data_maskEvent
                # del Barycentercharge

        del hdu_list            

    # list_EventCharge_AllExtensions = list_EventCharge_extension_2 + list_EventCharge_extension_1 + list_EventCharge_extension_4

    # dict_to_save_pkl = {'Muons_Detected' : num_muons, 'charge' : list_EventCharge_AllExtensions, 'DeltaEL' : list_DeltaEL}

    dict_to_save_pkl = {'All_Muons_Detected' : num_muons, 
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'Vertical_Events' : list_vertical_event_extension_1, 'Horizontal_Events' : list_horizontal_event_extension_1}, 
                        'extension_2' : {'charge' : list_EventCharge_extension_2, 'Vertical_Events' : list_vertical_event_extension_2, 'Horizontal_Events' : list_horizontal_event_extension_2},
                        'extension_4' : {'charge' : list_EventCharge_extension_4, 'Vertical_Events' : list_vertical_event_extension_4, 'Horizontal_Events' : list_horizontal_event_extension_4}}

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


    file_name = 'dict__straight_muons_Extensions_1_to_4_Imgs_' + str(len(argObj)) + '_Elip_'+str(elip) + '_Sol_' + str(Solidit) + '_ADUs__.pkl'
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