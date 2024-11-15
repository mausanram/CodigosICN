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
CCD_depth = 725 #micras
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
#### ---- LOS DATOS DE CONNIE YA ESTÁN CALIBRADOS EN ELECTRONES Y SE CARGAN LOS DATOS ASÍ --------- ###
units = 1
n_sigmas = 4
numero_bins = 500


def main(argObj):
    list_totalEvents = []

    list_charge_of_all_extension_1 = []
    list_EventCharge_extension_1 = []
    list_DeltaEL_extension_1 = []
    list_DeltaL_extension_1 = []
    list_theta_extension_1 = []
    list_phi_extension_1 = []

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
            dataCal = hdu_list[extension].data[:,:]
            header = hdu_list[extension].header
            # oScan = hdu_list[extension].data[:,550:]

        except:
            print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
            continue


        sigma_eletrons = header['RD_NOISE']     # Se lee la sigma del header de cada extensión 
        fondo_value = n_sigmas * sigma_eletrons

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

        DeltaL, DeltaEL, list_charge, _, list_theta, list_phi, list_charge_all_events = muon_filter(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, 
                                                                                        prop=prop, Solidit=Solidit, Elipticity=Elip)

        for index in np.arange(0, len(DeltaEL)):
            list_charge_of_all_extension_1.append(list_charge_all_events[index])
            list_DeltaEL_extension_1.append(DeltaEL[index])
            list_EventCharge_extension_1.append(list_charge[index])
            list_DeltaL_extension_1.append(DeltaL[index])
            list_theta_extension_1.append(list_theta[index])
            list_phi_extension_1.append(list_phi[index])

        print('Imagen ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    # num_muons = len(list_EventCharge_extension_1) + len(list_EventCharge_extension_2) + len(list_EventCharge_extension_4)
    num_muons = len(list_EventCharge_extension_1)

    dict_to_save_pkl = {'Num_Images' : total_images , 'All_Muons_Detected' : num_muons, 'Energy_Units' : units, 'Elipcidad' : Elip, 
                        'Solidity' : Solidit,
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'deltaEL' : list_DeltaEL_extension_1,
                                         'deltaL' : list_DeltaL_extension_1, 'all_events' : list_charge_of_all_extension_1,
                                         'theta': list_theta_extension_1, 'phi': list_phi_extension_1}}

    # dict_to_save_pkl = {'Num_Images' : total_images , 'All_Muons_Detected' : num_muons, 'Energy_Units' : units, 'Elipcidad' : list_Elip, 
    #                     'Solidity' : list_Solidit,
    #                     'extension_1' : {'charge' : list_EventCharge_extension_1, 'deltaEL' : list_DeltaEL_extension_1,
    #                                      'deltaL' : list_DeltaL_extension_1, 'all_events' : list_charge_of_all_extension_1,
    #                                      'theta': list_theta_extension_1}, 
    #                     'extension_2' : {'charge' : list_EventCharge_extension_2, 'deltaEL' : list_DeltaEL_extension_2, 
    #                                     'deltaL' : list_DeltaL_extension_2, 'all_events' : list_charge_of_all_extension_2,
    #                                     'theta': list_theta_extension_2},
    #                     'extension_4' : {'charge' : list_EventCharge_extension_4, 'deltaEL' : list_DeltaEL_extension_4, 
    #                                      'deltaL' : list_DeltaL_extension_4, 'all_events' : list_charge_of_all_extension_4,
    #                                      'theta': list_theta_extension_4}}

    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    # relacion = total_events / num_muons
    
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    print(Eventos_Totales)
    print(eventos_rectos)

    ext = img.split('/')[-1].split('_')[-2].split('g')[-1]
    if units == 0:
        file_name = 'dict_muons_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + '_img' + str(ext) + '_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_ADUs.pkl'

    elif units == 1:
        file_name = 'dict_muons_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + '_img' + str(ext) + '_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_electrons.pkl'

    elif units == 2:
        file_name = 'dict_muons_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + '_img' + str(ext) + '_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_KeV.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
