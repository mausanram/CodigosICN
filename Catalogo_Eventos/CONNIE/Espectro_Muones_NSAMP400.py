from astropy.io import fits
import numpy as np
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
imgsize_y = 1022
imgsize_x = 420

## Datos del filtro de muones GENERAL
Solidit = 0.65
Elip = 0.65
RUNID = 116

dedl_value_min = 1300 #KeV/cm

ratio_keV = 0.00368

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
#### ==== LOS DATOS DE CONNIE YA ESTÁN CALIBRADOS EN ELECTRONES Y SE CARGAN LOS DATOS ASÍ ==== ###
units = 2
n_sigmas = 4
numero_bins = 500


def main(argObj):
    list_totalEvents = []

    list_charge_of_all_extension_1 = []
    list_elip_of_all_extension_1 = []
    list_sol_of_all_extension_1 = []


    list_EventCharge_extension_1 = []
    list_DeltaEL_extension_1 = []
    list_DeltaL_extension_1 = []
    list_theta_extension_1 = []
    list_phi_extension_1 = []
    list_elip_extension_1 = []
    list_sol_extension_1 = []
    list_datamasked_extension_1 = []

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
        
        extension = 0

        try :
            data = hdu_list[extension].data[:imgsize_y,:imgsize_x]
            header = hdu_list[extension].header

            sigma_eletrons = header['RD_NOISE']     # Se lee la sigma del header de cada extensión 

            if units == 1:
                dataCal = data## En electrones
                sigma = abs(sigma_eletrons)

            elif units == 2:
                dataCal = ratio_keV * data## En keV
                sigma= sigma_eletrons * ratio_keV


        except:
            print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
            continue

        fondo_value = n_sigmas * sigma

        label_img, n_events = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
        prop = sk.measure.regionprops(label_img, dataCal)
        
        list_totalEvents.append(n_events)
        # print(nlabels_img)
        # list_labels.append(label_img)
        # list_EventsNumber.append(n_events)
        
        ## Obteniendo el valor promedio del fondo
        fondo_mask = np.invert(label_img == 0)
        fondo = ma.masked_array(dataCal,fondo_mask)


        dict_lists = muon_filter(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, prop=prop, Solidit=Solidit, Elipticity=Elip, dedl_min= dedl_value_min)

        DeltaL = dict_lists["muons"]["l"]
        DeltaEL = dict_lists["muons"]["dedl"]
        list_charge = dict_lists["muons"]["charge_muons"]
        list_theta = dict_lists["muons"]["theta"]
        list_phi = dict_lists["muons"]["phi"]
        list_elip = dict_lists["muons"]["elip"]
        list_sol = dict_lists["muons"]["sol"]
        list_datamasked = dict_lists["muons"]["image"]

        list_charge_all_events = dict_lists["non_muons"]["charge"]
        list_elip_all = dict_lists["non_muons"]["elip"]
        list_sol_all = dict_lists["non_muons"]["sol"]

        for index in np.arange(0, len(DeltaEL)):
            list_DeltaEL_extension_1.append(DeltaEL[index])
            list_EventCharge_extension_1.append(list_charge[index])
            list_DeltaL_extension_1.append(DeltaL[index])
            list_theta_extension_1.append(list_theta[index])
            list_phi_extension_1.append(list_phi[index])
            list_elip_extension_1.append(list_elip[index])
            list_sol_extension_1.append(list_sol[index])
            list_datamasked_extension_1.append(list_datamasked[index])

        for index in np.arange(0, len(list_charge_all_events)):
            list_charge_of_all_extension_1.append(list_charge_all_events[index])
            list_elip_of_all_extension_1.append(list_elip_all[index])
            list_sol_of_all_extension_1.append(list_sol_all[index])


        print('Imagen ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    # num_muons = len(list_EventCharge_extension_1) + len(list_EventCharge_extension_2) + len(list_EventCharge_extension_4)
    num_muons = len(list_EventCharge_extension_1)

    dict_to_save_pkl = {'Num_Images' : total_images , 'All_Muons_Detected' : num_muons, 'Energy_Units' : units, 'Elipcidad' : Elip, 
                        'Solidity' : Solidit,
                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'deltaEL' : list_DeltaEL_extension_1,
                                         'deltaL' : list_DeltaL_extension_1, 'all_events' : list_charge_of_all_extension_1,
                                         'theta': list_theta_extension_1, 'phi': list_phi_extension_1,
                                         'elip' : list_elip_extension_1, 'sol' : list_sol_extension_1,
                                         'all_events_elip' : list_elip_of_all_extension_1, 'all_events_sol' : list_sol_of_all_extension_1,
                                         'datamasked' : list_datamasked_extension_1}
                        }


    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    # relacion = total_events / num_muons
    

    print(Eventos_Totales)
    print(eventos_rectos)

    ext = img.split('/')[-1].split('_')[-2].split('g')[-1]
    if units == 0:
        file_name = 'dict_muons_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_Images_' + str(len(argObj)) + \
            '_img' + str(ext) + '_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_ADUs.pkl'

    elif units == 1:
        file_name = 'dict_muons_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_NIMG_' + str(len(argObj)) + \
            '_SOL_' + str(Solidit) + '_ELIP_'+str(Elip) + '_DEDL_' + str(dedl_value_min) + '_SIZE_' + str(imgsize_y) + 'x' +str(imgsize_x) + '_electrons_new.pkl'

    elif units == 2:
        file_name = 'dict_muons_NSAMP400_CONNIE_RUNID_' + str(RUNID) + '_NIMG_' + str(len(argObj)) + \
            '_SOL_' + str(Solidit) + '_ELIP_'+str(Elip) + '_DEDL_' + str(dedl_value_min) + '_SIZE_' + str(imgsize_y) + 'x' +str(imgsize_x) + '_KeV_new.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
