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
ratio_keVtoe = 0.00367

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
n_sigmas = 35
nbins = 200

## === Active Area range
x_min, x_max  = 10, 529
y_min, y_max  = 10, 250

## ===== WRITE WORKED EXTENSIONS ARRAY ====== ##
list_CCD_array = [1,2,4] # Use the real extension number
Nsamp = 324

## === CHECK THE FILE NAME AT THE END OF main() function
path_gains= './dict_mean_gains_Muons_NSAMP300.pkl'

def main(argObj):
    start = datetime.datetime.now()
    dict_to_save_pkl = {'Energy_Units': units}

    Bins = nbins
    Bins_fit = Bins

    total_images = len(argObj)
    image_in_bucle = 0
    
    print(f'=== START TIME: {start} === ')

    try:
        dict_gain = open(path_gains, 'rb')
        data_dict_gain = pkl.load(dict_gain)
        dict_gain.close()

        print(f"Gains file LOADED: {path_gains}. Analizing images...")

        dict_to_save_pkl = dict()
        list_CCD_array = []
        for element in list(data_dict_gain.keys()):
            list_CCD_array.append(int(element.split("_")[1]))
            dict_to_save_pkl[f"extension_{int(element.split('_')[1])}"] = {"tot_images": total_images, "units" : units,
                                                                           "all_events" : {"tot_events": 0,"charge": []}}
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

            ### ====  CLUSTERING PROCESS ===== ###
            label_img, n_events = sk.measure.label(dataCal > threshold, connectivity=2, return_num=True)
            # print('Events detected:', n_events)
            list_charge = all_cluster(dataCal=dataCal, label_img=label_img, nlabels_img=n_events)
            
            for index in np.arange(0, len(list_charge)):
                dict_to_save_pkl[f"extension_{extension+1}"]["all_events"]["charge"].append(list_charge[index])
                    
        print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    total_events_allext = 0
    for extension in list_CCD_array:
        dict_to_save_pkl[f"extension_{extension}"]["all_events"]["tot_events"] = len(dict_to_save_pkl[f"extension_{extension}"]["all_events"]["charge"])
        total_events_allext += dict_to_save_pkl[f"extension_{extension}"]["all_events"]["tot_events"]

    End = datetime.datetime.now()
    print(f'=== End time: {End}')
    print(f'=== Ellapsed time: {End - start} === \n' )
    print(f'Analized Images: {total_images}')
    print(f"Total events detected: {total_events_allext}")

    init_path = 'dict_energy_allclusters_NSAMP' + str(Nsamp) + '_Extensions_1_2_4_NIMGS_' + str(len(argObj)) + \
                '_NSIGMAS_' + str(n_sigmas) + '_SIZE_' + str(x_max) + 'x' + str(y_max)
    
    if units == 0: end_path = '_ADU.pkl'
    elif units == 1: end_path = '_electron.pkl'
    elif units == 2: end_path = '_keV.pkl'

    file_name = init_path + end_path # FILE NAME

    file_object_histogram = open(file_name, 'wb')
    pkl.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

