from astropy.io import fits
import numpy as np
import numpy.ma as ma
import sys
import skimage as sk
import datetime
import pickle
import os

from functions_MuonsNSAMP1 import *

## CONSTANTES ## 
current_path = os.getcwd()

ratio_keV = 0.0036

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
n_sigmas = 5
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
                # print('Voy a obtener el OsCan y el active area')
                data = hdu_list[extension].data[:250,10:539]
                oScan = hdu_list[extension].data[:250,539:]

                oscan_x = oScan.shape[1]
                oscan_y = oScan.shape[0]


                # print('Voy a obtener el valor medio de los píxeles')
                mean_rows_value = []
                for element in np.arange(0, oscan_y):
                    row = oScan[element: element +1, 0: oscan_x]
                    mean_value = np.median(row)
                    mean_rows_value.append([mean_value])

                true_active_area = data - mean_rows_value

            except:
                print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
                continue
            
            try:
                dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=numero_bins, 
                                                    Bins_fit=numero_bins,make_figure_flag=False, range_fit=[-50, 350])

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
                        nerr_ext = nerr_ext + 1
                        print('Fit error in extension ' + str(extension) + ' of image ' + str(img))


            except:
                print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                continue
            
            
            dataCal, sigma = data_calibrated_NSAMP(active_area=true_active_area, extension=extension, gain=Gain, 
                                                   ratio_keV=ratio_keV, unidades= units, offset=Offset, sigma_ADUs = sig_ADUs)
            
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

    if units == 0:
        file_name = 'dict_energy_allclusters_NSAMP324_Extensions_1_to_4_Imgs_' + str(total_images) + '_SIZE_250x529_' + '_NSIGMAS_' + str(n_sigmas)  + '_ADUs.pkl'
    elif units == 1:
        file_name = 'dict_energy_allclusters_NSAMP324_Extensions_1_to_4_Imgs_' + str(total_images) + '_SIZE_250x529_' + '_NSIGMAS_' + str(n_sigmas)  + '_electrons.pkl'
    elif units == 2:
        file_name = 'dict_energy_allclusters_NSAMP324_Extensions_1_to_4_Imgs_' + str(total_images) + '_SIZE_250x529_' + '_NSIGMAS_' + str(n_sigmas)  + '_KeV.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')
    
    print('Eventos Detectados en Total: ' +  str(total_events))
    print('Imágenes con error al cargar: ' + str(nerr_img))
    print('Error en fit de extension: ' + str(nerr_ext))
    

    print('Number of ext1: ', n_extension_1)
    print('Number of ext2: ', n_extension_2)
    print('Number of ext4: ', n_extension_4)
    print('Number of total img: ', n_total_img)
    print('Number of total ext 1+2: ', n_extension_1 + n_extension_2)
    print('Number of total ext: ', n_total_ext)



if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

