# from functions_py import math
from functions_MuonsNSAMP1 import *

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

## Datos de la CCd
CCD_depth = 725 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000

## Datos del filtro de muones GENERAL
Solidit = 0.7
Elip = 0.65

## Datos del filtro POR EXTENSIÓN
list_Elip = [0.65, 0.65, 0, 0.65]
list_Solidit = [0.7, 0.7, 0, 0.7]

DeltaEL_range_min, DeltaEL_range_max = 0.9, 3.55

ratio_keV = 0.0037
DeltaEL_range = 85

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 0
n_sigmas = 20
numero_bins = 500

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def Landau(x,a, MP,xi):
    C1 = a/np.sqrt((2 * np.pi))
    C2 = np.exp(-(x-MP)/xi)
    C3 = np.exp((-0.5 * ((x-MP)/xi + C2 )))
    return  C1 * C3


def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7] ## Ganancia para 324nsamp 

    # expgain = [110.78158608889959, 144.4661840232508, 100,  63.730700893071976] ## 

    list_totalEvents = []

    list_charge_of_all_extension_1 = []
    list_charge_of_all_extension_2 = []
    list_charge_of_all_extension_4 = []

    list_EventCharge_extension_2 =[]
    list_EventCharge_extension_1 = []
    list_EventCharge_extension_4 = []

    list_DeltaEL_extension_2 = []
    list_DeltaEL_extension_1 = []
    list_DeltaEL_extension_4 = []

    list_DeltaL_extension_2 = []
    list_DeltaL_extension_1 = []
    list_DeltaL_extension_4 = []

    list_theta_extension_2 = []
    list_theta_extension_1 = []
    list_theta_extension_4 = []

    list_phi_extension_2 = []
    list_phi_extension_1 = []
    list_phi_extension_4 = []


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
        
        for extension in (0,1,3):
        # for extension in (0, 1):
            # extension = 1
            Elip = list_Elip[extension]
            Solidit = list_Solidit[extension]

            try: 
                data = hdu_list[extension].data[:,:550]
                oScan = hdu_list[extension].data[:,550:]
            except:
                print('Loading error in image ' + str(img) + 'in load the data.')
                continue
            # nsamp = float(header['NSAMP'])

            # oScan_fit(extensión, active_area, oScan, Bins, make_figure_flag = False):

            try : 
                dict_popt = oScan_fit_NSAMP1(extensión=extension, active_area=data, oScan=oScan, Bins=numero_bins, make_figure_flag=False)

            except:
                print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                continue

            sig_ADUs = dict_popt['sigma']
            sig_electrons = (sig_ADUs) / expgain[extension-1]
            sig_KeV = sig_electrons * ratio_keV
            Offset = dict_popt['Offset']

            dataCal = data_calibrated(active_area=data, extension=extension, list_gain=expgain, ratio_keV=ratio_keV, unidades= units, offset=Offset)

            if units == 0:
                fondo_value = n_sigmas * sig_ADUs
            elif units == 1:
                fondo_value = n_sigmas * sig_electrons
            elif units == 2:
                fondo_value = n_sigmas * sig_KeV
            
            del oScan
            # del bin_heights
            # del bin_centers
            # del offset
            # del xmin_fit
            # del xmax_fit

            # label, n_events = ndimage.label(dataCal>6*abs(popt[2]),structure=[[1,1,1],[1,1,1],[1,1,1]])
            label_img, n_events = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)

            list_totalEvents.append(n_events)
            # print(nlabels_img)
            # list_labels.append(label_img)
            # list_EventsNumber.append(n_events)


            # list_charge_of_all = []
            for event in range(1, n_events):
                mask = np.invert(label_img == event)
                loc = ndimage.find_objects(label_img == event)[0]
                
                data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                                    mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

                charge = data_maskEvent.sum()
                # list_charge_of_all.append(charge)
                
            
            DeltaL, DeltaEL, list_charge, _, list_theta, list_phi, list_charge_all_events = muon_filter(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, 
                                                                    prop=prop, Solidit=Solidit, Elipticity=Elip)

            if extension == 0: 
                for index in np.arange(0, len(DeltaEL)):
                    list_charge_of_all_extension_1.append(list_charge_all_events[index])
                    list_DeltaEL_extension_1.append(DeltaEL[index])
                    list_EventCharge_extension_1.append(list_charge[index])
                    list_DeltaL_extension_1.append(DeltaL[index])
                    list_theta_extension_1.append(list_theta[index])
                    list_phi_extension_1.append(list_phi[index])
                    
            if extension == 1: 
                for index in np.arange(0, len(DeltaEL)):
                    list_charge_of_all_extension_2.append(list_charge_all_events[index])
                    list_DeltaEL_extension_2.append(DeltaEL[index])
                    list_EventCharge_extension_2.append(list_charge[index])
                    list_DeltaL_extension_2.append(DeltaL[index])
                    list_theta_extension_2.append(list_theta[index])
                    list_phi_extension_2.append(list_phi[index])
            
            if extension == 3: 
                for index in np.arange(0, len(DeltaEL)):
                    list_charge_of_all_extension_4.append(list_charge_all_events[index])
                    list_DeltaEL_extension_4.append(DeltaEL[index])
                    list_EventCharge_extension_4.append(list_charge[index])
                    list_DeltaL_extension_4.append(DeltaL[index])
                    list_theta_extension_4.append(list_theta[index])
                    list_phi_extension_4.append(list_phi[index])
            # del data_maskEvent
            # del Barycentercharge

        print('Imagen ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list            

    # list_EventCharge_AllExtensions = list_EventCharge_extension_2 + list_EventCharge_extension_1 + list_EventCharge_extension_4

    # dict_to_save_pkl = {'Muons_Detected' : num_muons, 'charge' : list_EventCharge_AllExtensions, 'DeltaEL' : list_DeltaEL}

    num_muons = len(list_EventCharge_extension_1) + len(list_EventCharge_extension_2) + len(list_EventCharge_extension_4)

    dict_to_save_pkl = {'Num_Images' : total_images , 'All_Muons_Detected' : num_muons, 'Energy_Units' : units, 'Elipticity' : list_Elip, 
                        'Solidity' : list_Solidit,
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'deltaEL' : list_DeltaEL_extension_1,
                                         'deltaL' : list_DeltaL_extension_1, 'all_events' : list_charge_of_all_extension_1,
                                         'theta': list_theta_extension_1, 'phi': list_phi_extension_1}, 
                        'extension_2' : {'charge' : list_EventCharge_extension_2, 'deltaEL' : list_DeltaEL_extension_2, 
                                        'deltaL' : list_DeltaL_extension_2, 'all_events' : list_charge_of_all_extension_2,
                                        'theta': list_theta_extension_2, 'phi': list_phi_extension_2},
                        'extension_4' : {'charge' : list_EventCharge_extension_4, 'deltaEL' : list_DeltaEL_extension_4, 
                                         'deltaL' : list_DeltaL_extension_4, 'all_events' : list_charge_of_all_extension_4,
                                         'theta': list_theta_extension_4, 'phi': list_phi_extension_4}}


    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    # print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    # relacion = total_events / num_muons
    
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    print(Eventos_Totales)
    print(eventos_rectos)
    # print(eventos_circulares)  

    # file_name = 'dict_muons_Extensions_1_to_4_Imgs_' + str(len(argObj))+'_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_ADUs__all.pkl'

    if units == 0:
        file_name = 'dict_muons_Extensions_1_to_4_Imgs_' + str(len(argObj))+'_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_ADUs__.pkl'
    elif units == 1:
        file_name = 'dict_muons_Extensions_1_to_4_Imgs_' + str(len(argObj))+'_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_electrons__.pkl'
    elif units == 2:
        file_name = 'dict_muons_Extensions_1_to_4_Imgs_' + str(len(argObj))+'_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_KeV__.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

