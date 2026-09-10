from astropy.io import fits
import numpy as np
import numpy.ma as ma
import sys
import skimage as sk
import datetime
import pickle as pkl
import os

from functions_MuonsNSAMP1 import *

## CONSTANTES ## 
current_path = os.getcwd()

## Datos del filtro de muones GENERAL
Solidit = 0.65
Elip = 0.65
dedl_value_min = 1400

ratio_keVtoe = 0.00367

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
n_sigmas = 13
numero_bins = 600

## === Active Area range
x_min, x_max  = 10, 529
y_min, y_max  = 10, 250

def main(argObj):
    total_images = len(argObj)
    image_in_bucle = 0

    start = datetime.datetime.now()
    
    print(f'=== START TIME: {start} === ')

    # path= './dict_mean_gains_Muons_NSAMP324.pkl'
    path_gains= './dict_mean_gains_NSAMP324.pkl'

    try:
        dict_gain = open(path_gains, 'rb')
        data_dict_gain = pkl.load(dict_gain)
        dict_gain.close()

        print(f"Gains file LOADED: {path_gains}. Analizing images...")

        dict_to_save_pkl = dict()
        list_CCD_array = []
        for element in list(data_dict_gain.keys()):
            list_CCD_array.append(int(element.split("_")[1]))
            dict_to_save_pkl[f"extension_{int(element.split('_')[1])}"] = {"tot_images": total_images, 
                                                                           "elip_used": Elip, 
                                                                           "sol_used" : Solidit,
                                                                           "units" : units,
                                                                           "nsigmas" : n_sigmas,
                                                                           "muons" : {"tot_events": 0, "charge": [], "deltaL": [], 
                                                                                      "deltaEL": [], "theta" : [], "phi": [], 
                                                                                      "elip": [], "sol": [], "gain": [], 
                                                                                      "datamasked": [], "run": []},
                                                                           "all_events" : {"tot_events": 0,"charge": [], "elip": [],
                                                                                           "sol": []}
                                                                        }
    except:
        print(f"Gain file wasn't found in {path}. Aborting ... ")
        exit()

    for img in argObj:
        try:
            hdu_list = fits.open(img)
            image_in_bucle += 1

            path = img.split('/')
            run = path[2]
            # print(run)

        except:
            nerr_img = nerr_img + 1
            print('Loading error in image ' + str(img) + ' in open the image.')
            continue
        
        for extension in list_CCD_array:
            extension -= 1
            # Elip = list_Elip[extension]
            # Solidit = list_Solidit[extension]

            try :
                data = hdu_list[extension].data[y_min:y_max, x_min:x_max]
                oScan = hdu_list[extension].data[y_min:y_max, x_max:]
                oscan_shape = oScan.shape

                true_active_area = cleaning_actArea(activeArea=data, OvScan=oScan, x_range=[0, oscan_shape[1]], y_range=[0, oscan_shape[0]])
            except:
                print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
                continue

            Gain = data_dict_gain[f"extension_{extension+1}"]["Gain"] # ADU/e-
            sig_ADUs = data_dict_gain[f"extension_{extension+1}"]["Sigma"] # ADUs  ### CHANGE THE KEY FOR "Sig"
            
            dataCal, sigma = data_calibrated(active_area=true_active_area, gain=Gain, 
                                             ratio_keVtoe=ratio_keVtoe, units= units, sigma_ADU=sig_ADUs)

            threshold = n_sigmas * sigma
            del oScan

            label_img, n_events = sk.measure.label(dataCal > threshold, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
           
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)

            dict_lists = muon_filter(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, 
                                     prop=prop, Solidit=Solidit, Elipticity=Elip, dedl_min= dedl_value_min)
            
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
                ### ===== Muons ===== ###
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["deltaEL"].append(DeltaEL[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["charge"].append(list_charge[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["deltaL"].append(DeltaL[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["theta"].append(list_theta[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["phi"].append(list_phi[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["elip"].append(list_elip[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["sol"].append(list_sol[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["gain"].append(Gain)
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["datamasked"].append(list_datamasked[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["muons"]["run"].append(run[index])

            for index in np.arange(0, len(list_charge_all_events)):
                ### ==== All events ==== ###
                dict_to_save_pkl[f"extension_{extension+1}"]["all_events"]["charge"].append(list_charge_all_events[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["all_events"]["elip"].append(list_elip_all[index])
                dict_to_save_pkl[f"extension_{extension+1}"]["all_events"]["sol"].append(list_sol_all[index])
                    
        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    total_events_allext = 0
    muons_detected = 0
    for extension in list_CCD_array:
        dict_to_save_pkl[f"extension_{extension}"]["muons"]["tot_events"] = len(dict_to_save_pkl[f"extension_{extension}"]["muons"]["deltaEL"])
        dict_to_save_pkl[f"extension_{extension}"]["all_events"]["tot_events"] = len(dict_to_save_pkl[f"extension_{extension}"]["all_events"]["charge"]) 

        total_events_allext += dict_to_save_pkl[f"extension_{extension}"]["muons"]["tot_events"]
        total_events_allext += dict_to_save_pkl[f"extension_{extension}"]["all_events"]["tot_events"]

        muons_detected += dict_to_save_pkl[f"extension_{extension}"]["muons"]["tot_events"]

    End = datetime.datetime.now()

    print(f'=== End time: {End}')
    print(f'=== Ellapsed time: {End - start} === \n' )
    print(f'Analized Images: {total_images}')

    print(f"Total events detected: {total_events_allext}")
    print(f"Muons detected: {muons_detected}")

    init_path = 'dict_muons_NSAMP324_Extensions_1_2_4_NIMGS_' + str(len(argObj)) + \
                '_SOL_' + str(Solidit) + '_ELIP_'+str(Elip) + '_NSIGMAS_' + str(n_sigmas) + \
                '_DEDL_' + str(dedl_value_min)+'_SIZE_' + str(x_max) + 'x' + str(y_max)
    
    if units == 0: end_path = '_ADU.pkl'
    elif units == 1: end_path = '_electron.pkl'
    elif units == 2: end_path = '_keV.pkl'

    file_name = init_path + end_path

    file_object_histogram = open(file_name, 'wb')
    pkl.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

