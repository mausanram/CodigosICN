import math
from functions_CONNIE import *
from astropy.io import fits
import scipy.ndimage as nd
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
units = 1

n_sigmas = 4
ratio_keV = 0.0037
CCD_depth = 725 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000
DeltaEL_range = 85


Solidit = 0.7
Elipticity = 0.8
# min_Charge =  100 # keV
min_Charge =  100 # electrons
n_skirts = 1

numero_bins = 500


def main(argObj):
    list_sigmas_vertical_event_extension_1 = []

    list_sigmas_horizontal_event_extension_1 = []
    list_EventCharge_extension_1 = []

    list_totalEvents = []
    list_vertical_event_extension_1 = []

    list_horizontal_event_extension_1 = []

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
        
        # for extension in (0,1,3):
        extension = 0
        try :
            # print('Voy a obtener el OsCan y el active area')
            dataCal = hdu_list[extension].data[:,:]
            header = hdu_list[extension].header

        except:
            print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
            continue

        # sig_electrons = abs((sig_ADUs) / Gain)

        sigma_eletrons = header['RD_NOISE']
        fondo_value = n_sigmas * sigma_eletrons

        # label, n_events = ndimage.label(dataCal>6*abs(popt[2]),structure=[[1,1,1],[1,1,1],[1,1,1]])
        label_img, n_events = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
        prop = sk.measure.regionprops(label_img, dataCal)
        # print(type(prop))

        list_totalEvents.append(n_events)
        # print(nlabels_img)
        # list_labels.append(label_img)
        # list_EventsNumber.append(n_events)
        
        ## Obteniendo el valor promedio del fondo
        fondo_mask = np.invert(label_img == 0)
        fondo = ma.masked_array(dataCal,fondo_mask)
        valor_promedio_fondo = fondo.data.mean()

        list_vertical, list_horizontal = muon_straight_filter(dataCal= dataCal, label_img=label_img, n_events=n_events, Solidit=Solidit, Elipticity=Elipticity, 
                                                                Prop= prop, min_Charge=min_Charge, Sigma=sigma_eletrons, skirts= n_skirts)

        for index in np.arange(0, len(list_vertical[0])):
            list_EventCharge_extension_1.append(list_vertical[2][index])
            list_sigmas_vertical_event_extension_1.append(list_vertical[0][index])
            list_vertical_event_extension_1.append(list_vertical[1][index])

        for index in np.arange(0, len(list_horizontal[0])):
            list_EventCharge_extension_1.append(list_horizontal[2][index])
            list_sigmas_horizontal_event_extension_1.append(list_horizontal[0][index])
            list_horizontal_event_extension_1.append(list_horizontal[1][index])
        
        # print('Extension ' + str(extension) + ' terminada.')

        del hdu_list         
    
    num_muons = len(list_EventCharge_extension_1)

    dict_to_save_pkl = {'All_Muons_Detected' : num_muons, 
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'vertical_sigmas' : list_sigmas_vertical_event_extension_1,
                        'Vertical_Events' : list_vertical_event_extension_1, 'horizontal_sigmas' : list_sigmas_horizontal_event_extension_1, 
                        'Horizontal_Events' : list_horizontal_event_extension_1}}

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


    if units == 0:
        file_name = 'dict_straight_muons_CONNIE_NSAMP400_RUNID_116_Images_' + str(len(argObj)) + '_Elip_'+str(Elipticity) + '_Sol_' + str(Solidit) + '_Skirts_'+str(n_skirts) + '_ADUs.pkl'

    elif units == 1:
        file_name = 'dict_straight_muons_CONNIE_NSAMP400_RUNID_116_Images_' + str(len(argObj)) + '_Elip_'+str(Elipticity) + '_Sol_' + str(Solidit) + '_Skirts_'+str(n_skirts) + '_electrons.pkl'

    elif units == 2:
        file_name = 'dict_straight_muons_CONNIE_NSAMP400_RUNID_116_Images_' + str(len(argObj)) + '_Elip_'+str(Elipticity) + '_Sol_' + str(Solidit) + '_Skirts_'+str(n_skirts) + '_KeV.pkl'


    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in ', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

