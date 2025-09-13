# from functions_py import math

import sys
import os
## ===== Así se importan las funciones de arhcivos en otras direcciones === ###
module_path = os.path.abspath('../.') # Adjust relative path as needed
# Or use an absolute path: module_path = '/home/your_user/projects/utils'

if module_path not in sys.path:
    sys.path.append(module_path)

from functions_MuonsNSAMP1 import * 
from functions_py import *

from astropy.io import fits
import scipy.ndimage as ndimage
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

ratio_keV = 0.0036

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
n_sigmas = 40
numero_bins = 1000

def Gaussian2(x,m,s,g,a1,a2): #data, mean, sigma, gain, height1, heigth2
    return a1*np.exp(-1/2*((x-m)/s)**2)+a2*np.exp(-1/2*((x-m-g)/s)**2)

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def main(argObj):
    list_totalEvents = []


    list_EventCharge_extension_2 =[]
    list_EventCharge_extension_1 = []
    list_EventCharge_extension_4 = []

    nerr_img = 0
    nerr_ext = 0

    total_images = len(argObj)
    image_in_bucle = 0

    n_extension_1 = 0
    n_extension_2 = 0
    n_extension_4 = 0
    n_total_img = 0
    n_total_ext = 0

    Inicio = datetime.datetime.now()
    num_images =  'Imágenes Analizadas: ' +  str(total_images)
    
    print('Hora de inicio del cálculo: ', Inicio)
    for img in argObj:
        try:
            hdu_list = fits.open(img)
            image_in_bucle += 1

        except:
            nerr_img = nerr_img + 1
            print('Loading error in image ' + str(img) + 'in open the image.')
            continue
        
        for extension in (0,1,3):            
            try :
                max_y = 300
                max_x = 539

                data = hdu_list[extension].data[:max_y,10:max_x]
                oScan = hdu_list[extension].data[:max_y,max_x:]

                oscan_x = oScan.shape[1]
                oscan_y = oScan.shape[0]

                mean_rows_value = []
                for element in np.arange(0, oscan_y):
                    row = oScan[element: element +1, 0: oscan_x]
                    num_row = element + 1
                    mean_value = np.median(row)
                    mean_rows_value.append([mean_value])

                true_active_area = data - mean_rows_value

            except:
                print('Loading error in extension ' + str(extension + 1) + ' of image ' + str(img) + 'in load the data.')
                continue
            
            ### Change the range for any kind of image
            Range_fit = [-50, 350]  # FOr Fe-55
            # Range_fit = [-100, 270] # For Fe-55 & Cs-137
            try:
                dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=numero_bins, 
                                                    Bins_fit=numero_bins,make_figure_flag=False, range_fit=[Range_fit[0], Range_fit[1]])

                sig_ADUs = dict_popt['sigma']
                Offset = dict_popt['Offset']
                Gain = dict_popt['Gain']
                Prob = dict_popt['Prob']
                
                if Prob < 0.05:
                    del_Bin = 500
                    dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                        Bins_fit=del_Bin, make_figure_flag=False, range_fit=[Range_fit[0], Range_fit[1]])
                    sig_ADUs = dict_popt['sigma']
                    Offset = dict_popt['Offset']
                    Gain = dict_popt['Gain']
                    Prob = dict_popt['Prob']
                    
                    if Prob < 0.05:
                        del_Bin = 400
                        dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                            Bins_fit=del_Bin, make_figure_flag=False, range_fit=[Range_fit[0], Range_fit[1]])
                        
                        sig_ADUs = dict_popt['sigma']
                        Offset = dict_popt['Offset']
                        Gain = dict_popt['Gain']
                        Prob = dict_popt['Prob']
                    
                        
                        if Prob < 0.05:
                            del_Bin = 300
                            dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                                Bins_fit=del_Bin, make_figure_flag=False, range_fit=[Range_fit[0], Range_fit[1]])
                            
                            sig_ADUs = dict_popt['sigma']
                            Offset = dict_popt['Offset']
                            Gain = dict_popt['Gain']
                            Prob = dict_popt['Prob']
                    

                            if  Prob < 0.05:
                                nerr_ext = nerr_ext + 1
                                if extension == 0:
                                    nerr_ext1 += 1
                                elif extension == 1:
                                    nerr_ext2 += 1
                                elif extension == 3:
                                    nerr_ext4 += 1

                                print('Fit error in extension ' + str(extension + 1) + ' of image ' + str(img))
                                continue

            except:
                print('Fit error in extension ' + str(extension + 1) + ' of image ' + str(img))
                continue
            
            dataCal, sigma = data_calibrated_NSAMP(active_area=true_active_area, gain=Gain, ratio_keV=ratio_keV, unidades= units, sigma_ADUs = sig_ADUs)
            
            fondo_value = n_sigmas * sigma
            
            del oScan
            

            label_img, n_events = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            
            list_totalEvents.append(n_events)
            
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)

            list_charge = all_cluster(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, prop=prop)

            if extension == 0: 
                n_extension_1 = n_extension_1 + 1
                for index in np.arange(0, len(list_charge)):
                    list_EventCharge_extension_1.append(list_charge[index])

            if extension == 1: 
                n_extension_2 = n_extension_2 + 1
                for index in np.arange(0, len(list_charge)):
                    list_EventCharge_extension_2.append(list_charge[index])

            if extension == 3: 
                n_extension_4 = n_extension_4 + 1
                for index in np.arange(0, len(list_charge)):
                    list_EventCharge_extension_4.append(list_charge[index])

            n_total_ext = n_total_ext + 1

        n_total_img = n_total_img + 1
        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    num_clusters = len(list_EventCharge_extension_1) + len(list_EventCharge_extension_2) + len(list_EventCharge_extension_4)

    dict_to_save_pkl = {'Num_Images' : total_images , 'All_clusters_detected' : num_clusters, 'Energy_Units' : units,
                        'extension_1' : {'charge' : list_EventCharge_extension_1}, 
                        'extension_2' : {'charge' : list_EventCharge_extension_2},
                        'extension_4' : {'charge' : list_EventCharge_extension_4}}

    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    ### Change the file name when the images changes 
    # in_path = 'dict_energy_allclusters_Fe55_Cs137_NSAMP200_Extensions_1_to_4_Imgs_' # For Fe-55 & Cs-137
    in_path = 'dict_energy_allclusters_Fe55_NSAMP200_Extensions_1_to_4_Imgs_' # For Fe-55 

    if units == 0:
        file_name = in_path + str(total_images) + '_SIZE_' + str(max_y)+ 'x' + str(max_x) + '_NSIGMAS_' + str(n_sigmas)  + '_ADUs.pkl'
    elif units == 1:
        file_name = in_path + str(total_images) + '_SIZE_' + str(max_y)+ 'x' + str(max_x) + '_NSIGMAS_' + str(n_sigmas)  + '_electrons.pkl'
    elif units == 2:
        file_name = in_path + str(total_images) + '_SIZE_' + str(max_y)+ 'x' + str(max_x) + '_NSIGMAS_' + str(n_sigmas)  + '_KeV.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')
    
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    # eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    img_err = 'Imágenes con error al cargar: ' + str(nerr_img)
    ext_err = 'Error en fit de extension: ' + str(nerr_ext)
    # relacion = total_events / num_muons
    
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    print(img_err)
    print(ext_err)
    print(Eventos_Totales)
    print('Number of 1ext: ', n_extension_1)
    print('Number of 2ext: ', n_extension_2)
    print('Number of 4ext: ', n_extension_4)
    print('Number of total img: ', n_total_img)
    print('Number of total ext 1+2: ', n_extension_1 + n_extension_2)
    print('Number of total ext: ', n_total_ext)

    # print(eventos_rectos)

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

