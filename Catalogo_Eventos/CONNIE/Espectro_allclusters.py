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

from functions_CONNIE import *

# from ROOT import *

## CONSTANTES ## 
current_path = os.getcwd()

## Datos de la CCd
CCD_depth = 680 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000

## Datos del filtro de muones GENERAL
Solidit = 0.7
Elip = 0.65
RUNID = 116

## Datos del filtro POR EXTENSIÓN
list_Elip = [0.65, 0.65, 0, 0.65]
list_Solidit = [0.7, 0.7, 0, 0.7]

DeltaEL_range_min, DeltaEL_range_max = 0.9, 3.55

ratio_keV = 0.0036
DeltaEL_range = 85

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
#### ==== LOS DATOS DE CONNIE YA ESTÁN CALIBRADOS EN ELECTRONES Y SE CARGAN LOS DATOS ASÍ ==== ###
units = 2
n_sigmas = 4
numero_bins = 500

def main(argObj):
    
    list_totalEvents = []   
    list_EventCharge_extension_1 = []

    total_images = len(argObj)
    image_in_bucle = 0

    Inicio = datetime.datetime.now()
    num_images =  'Imágenes Analizadas: ' +  str(total_images)
    
    print('Hora de inicio del cálculo: ', Inicio)

    for img in argObj:
        try:
            hdu_list = fits.open(img)
            image_in_bucle += 1

        except:
            print('Loading error in image ' + str(img) + 'in open the image.')
            continue
        
        # for extension in (0,1,3):
            # Elip = list_Elip[extension]
            # Solidit = list_Solidit[extension]
        extension = 0

        try :
            # print('Voy a obtener el OsCan y el active area')
            dataCal = hdu_list[extension].data[:,:] # En electrones
            header = hdu_list[extension].header
            # oScan = hdu_list[extension].data[:,550:]

        except:
            print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
            continue


        sigma_eletrons = header['RD_NOISE']     # Se lee la sigma del header de cada extensión 
        fondo_value = n_sigmas * sigma_eletrons
        print(fondo_value)

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

        list_charge = all_cluster(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, prop=prop)

        for index in np.arange(0, len(list_charge)):
            list_EventCharge_extension_1.append(list_charge[index])

        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    num_clusters = len(list_EventCharge_extension_1)

    dict_to_save_pkl = {'Num_Images' : total_images , 'All_Muons_Detected' : num_clusters, 'Energy_Units' : units,
                        'extension_1' : {'charge' : list_EventCharge_extension_1}}

    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    ext = img.split('/')[-1].split('_')[-2].split('g')[-1]
    if units == 0:
        file_name = 'dict_allclustes_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + '_NSIGMAS_'+ str(n_sigmas) +'_img' + str(ext) + '_ADUs.pkl'

    elif units == 1:
        file_name = 'dict_allclustes_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + '_NSIGMAS_'+ str(n_sigmas) +'_img' + str(ext) + '_electrons.pkl'

    elif units == 2:
        file_name = 'dict_allclustes_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + '_NSIGMAS_'+ str(n_sigmas) +'_img' + str(ext) + '_KeV.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    print(Eventos_Totales)

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')

    # eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    # img_err = 'Imágenes con error al cargar: ' + str(nerr_img)
    # ext_err = 'Error en fit de extension: ' + str(nerr_ext)
    # relacion = total_events / num_muons
    
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    # print(img_err)
    # print(ext_err)
    # print(eventos_rectos)

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

