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

# from ROOT import *

## CONSTANTES ## 
current_path = os.getcwd()
units = 2

n_sigmas = 13
ratio_keV = 0.0036
CCD_depth = 725 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000
DeltaEL_range = 85


Solidit = 0.65
Elipticity = 0.9
min_Charge =  100 # keV
n_skirts = 0

numero_bins = 500

def Gaussian2(x,m,s,g,a1,a2): #data, mean, sigma, gain, height1, heigth2
    return a1*np.exp(-1/2*((x-m)/s)**2)+a2*np.exp(-1/2*((x-m-g)/s)**2)

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def main(argObj):
    list_sigmas_vertical_event_extension_1 = []
    list_sigmas_vertical_event_extension_2 = []
    list_sigmas_vertical_event_extension_4 = []

    list_sigmas_horizontal_event_extension_1 = []
    list_sigmas_horizontal_event_extension_2 = []
    list_sigmas_horizontal_event_extension_4 = []

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
                # print('Voy a obtener el OsCan y el active area')
                data = hdu_list[extension].data[:,:550]
                oScan = hdu_list[extension].data[:,550:]

                oscan_x = oScan.shape[1]
                oscan_y = oScan.shape[0]

                header = hdu_list[extension].header
                # nsamp = float(header['NSAMP'])

                # print('Voy a obtener el valor medio de los píxeles')
                mean_rows_value = []
                for element in np.arange(0, oscan_y):
                    row = oScan[element: element +1, 0: oscan_x]
                    num_row = element + 1
                    mean_value = np.median(row)
                    mean_rows_value.append([mean_value])

                true_active_area = data - mean_rows_value

            except:
                print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
                continue

            del header

            try:
                dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=numero_bins, 
                                                    Bins_fit=numero_bins,make_figure_flag=False, range_fit=[-30, 350])

                sig_ADUs = dict_popt['sigma']
                Offset = dict_popt['Offset']
                Gain = dict_popt['Gain']
                Prob = dict_popt['Prob']

                if Prob < 0.05:
                    del_Bin = 600
                    dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                        Bins_fit=del_Bin, make_figure_flag=False, range_fit=[-30, 400])
                    sig_ADUs = dict_popt['sigma']
                    Offset = dict_popt['Offset']
                    Gain = dict_popt['Gain']
                    Prob = dict_popt['Prob']

                    if  Prob < 0.05:
                        # nerr_ext = nerr_ext + 1
                        print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                        continue

                # if Gain < 100:
                #     ### Aquí se deberá poner la ganancia promedio de cada extensión una vez que se obtenga de muchas imágenes
                #     print('Fit gain error in extension ' + str(extension) + ' of image ' + str(img))
                #         # continue

            except:
                print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                continue

            # sig_electrons = abs((sig_ADUs) / Gain)

            dataCal, sigma = data_calibrated_NSAMP(active_area=true_active_area, extension=extension, gain=Gain, ratio_keV=ratio_keV, 
                                                   unidades= units, offset=Offset, sigma_ADUs = sig_ADUs)
            
            fondo_value = n_sigmas * sigma
            
            del oScan


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

            list_vertical, list_horizontal = muon_straight_filter(dataCal= dataCal, label_img=label_img, n_events=n_events, Solidit=Solidit, 
                                                                  Elipticity=Elipticity, Prop= prop, min_Charge=min_Charge, 
                                                                  Sigma=sigma, skirts= n_skirts)

            if extension == 0:
                for index in np.arange(0, len(list_vertical[0])):
                    list_EventCharge_extension_1.append(list_vertical[2][index])
                    list_sigmas_vertical_event_extension_1.append(list_vertical[0][index])
                    list_vertical_event_extension_1.append(list_vertical[1][index])

                for index in np.arange(0, len(list_horizontal[0])):
                    list_EventCharge_extension_1.append(list_horizontal[2][index])
                    list_sigmas_horizontal_event_extension_1.append(list_horizontal[0][index])
                    list_horizontal_event_extension_1.append(list_horizontal[1][index])

            if extension == 1:
                for index in np.arange(0, len(list_vertical[0])):
                    list_EventCharge_extension_2.append(list_vertical[2][index])
                    list_sigmas_vertical_event_extension_2.append(list_vertical[0][index])
                    list_vertical_event_extension_2.append(list_vertical[1][index])

                for index in np.arange(0, len(list_horizontal[0])):
                    list_EventCharge_extension_2.append(list_horizontal[2][index])
                    list_sigmas_horizontal_event_extension_2.append(list_horizontal[0][index])
                    list_horizontal_event_extension_2.append(list_horizontal[1][index])

            if extension == 3:
                for index in np.arange(0, len(list_vertical[0])):
                    list_EventCharge_extension_4.append(list_vertical[2][index])
                    list_sigmas_vertical_event_extension_4.append(list_vertical[0][index])
                    list_vertical_event_extension_4.append(list_vertical[1][index])

                for index in np.arange(0, len(list_horizontal[0])):
                    list_EventCharge_extension_4.append(list_horizontal[2][index])
                    list_sigmas_horizontal_event_extension_4.append(list_horizontal[0][index])
                    list_horizontal_event_extension_4.append(list_horizontal[1][index])
            
            # print('Extension ' + str(extension) + ' terminada.')

        del hdu_list         
    
    num_muons = len(list_EventCharge_extension_1) + len(list_EventCharge_extension_2) + len(list_EventCharge_extension_4)

    dict_to_save_pkl = {'All_Muons_Detected' : num_muons, 
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'vertical_sigmas' : list_sigmas_vertical_event_extension_1,
                        'Vertical_Events' : list_vertical_event_extension_1, 'horizontal_sigmas' : list_sigmas_horizontal_event_extension_1, 
                        'Horizontal_Events' : list_horizontal_event_extension_1}, 

                        'extension_2' : {'charge' : list_EventCharge_extension_2, 'vertical_sigmas' : list_sigmas_vertical_event_extension_2, 
                        'Vertical_Events' : list_vertical_event_extension_2, 'horizontal_sigmas' : list_sigmas_horizontal_event_extension_2, 
                        'Horizontal_Events' : list_horizontal_event_extension_2},

                        'extension_4' : {'charge' : list_EventCharge_extension_4,  'vertical_sigmas' : list_sigmas_vertical_event_extension_4, 
                        'Vertical_Events' : list_vertical_event_extension_4, 'horizontal_sigmas' : list_sigmas_horizontal_event_extension_4, 
                        'Horizontal_Events' : list_horizontal_event_extension_4}}

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
        file_name = 'dict_straight_muons_Extensions_1_to_4_Imgs_' + str(len(argObj)) + '_Elip_'+str(Elipticity) + '_Sol_' + str(Solidit) + '_Skirts_'+str(n_skirts) + '_with_sigmas_ADU__NSAMP324.pkl'
    elif units == 1:
        file_name = 'dict_straight_muons_Extensions_1_to_4_Imgs_' + str(len(argObj)) + '_Elip_'+str(Elipticity) + '_Sol_' + str(Solidit) + '_Skirts_'+str(n_skirts) + '_with_sigmas_ele__NSAMP324.pkl'
    elif units == 2:
        file_name = 'dict_straight_muons_Extensions_1_to_4_Imgs_' + str(len(argObj)) + '_Elip_'+str(Elipticity) + '_Sol_' + str(Solidit) + '_Skirts_'+str(n_skirts) + '_with_sigmas_KeV__NSAMP324.pkl'
    
    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in ', current_path + '/' + file_name, ' as a binary file.')
    print('To open use library "pickle". ')



if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
