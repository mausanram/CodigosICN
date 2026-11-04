from astropy.io import fits
import numpy as np
import numpy.ma as ma
import sys
import skimage as sk
import datetime
import pickle as pkl
import os
import pandas as pd

module_path = os.path.abspath('../functions') # Adjust relative path as needed
# # Or use an absolute path: module_path = '/home/your_user/projects/utils'

if module_path not in sys.path:
    sys.path.append(module_path)

from functionsNEW import *


## Write paths ## 
current_path = os.getcwd()
folder_to_save_path = "./data/" ## ADD THE  "/" AT THE END

## CCD arrays ##
CCD_array = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
# CCD_arrayOfficial = (1,2,5,6,7,9,14,15,16)
CCD_bad = (4, 7, 9, 11, 13, 14, 15)

CCD_depth = 675 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000

## MUON FILTER PARAMETERS
Solidit = 0.65
Elip = 0.65
dedl_value_min = 0

ratio_keV = 0.00367

## UNITS AND #SIGMAS (for units: 0 = ADUs, 1 = e-, 2 = KeV)
units = 1 #this images are in e-
n_sigmas = 10


def main(argObj):
    total_images = len(argObj)
    Inicio = datetime.datetime.now()
    dict_to_save_pkl = {'Energy_Units' : units, 'Num_Images' : total_images,}

    for extension in CCD_array:
        if extension in CCD_bad:
            continue
        dict_to_save_pkl[f"extension_{extension}"] = {"muons" : {}}
        dict_to_save_pkl[f"extension_{extension}"]["muons"] = {"l": [], "dedl" : [], "charge_muons" : [], "mean_charge" : [], 
                                                               "median_charge": [], "elip": [], "sol": [], "image" : [], "label" : [],
                                                               "theta": [], "phi": []}

    list_totalEvents = []
    image_in_bucle = 0
    total_clusters = 0

    num_images =  'Analyzed images: ' +  str(total_images)
    
    print('Start time: ', Inicio)
    for img in argObj:
        try:
            list_path_img = img.split("/")
            # print(list_path_img)
            hdu_list = fits.open(img)
            image_in_bucle += 1
            runID = list_path_img[-2]
            imgType = list_path_img[-1].split("_")[0]
            # print(imgType)
            folder_path = f"{folder_to_save_path}"
            os.makedirs(folder_to_save_path, exist_ok=True)

            img_path = ""
            for element in list_path_img[:-1]:
                img_path  = img_path + element + "/"

            # print(f"Current path: {current_path}")
            # print(f"RunID: {runID}")
            # print(f"ImgType: {imgType}")

        except:
            nerr_img = nerr_img + 1
            print('Loading error in image ' + str(img) + ' in open the image.')
            continue
        
        for extension in CCD_array:
            if extension in CCD_bad:
                continue

            try :
                x_max = 430
                y_max = 500
                border = 10
                data = hdu_list[extension].data[border:y_max-border, border:x_max-border] # Without the border
                # oScan = hdu_list[extension].data[border:y_max, x_max:]
                Gain = hdu_list[extension].header["GAIN_RUN"]   # ADU/e
                NSAMP = hdu_list[extension].header["NSAMP"]
            except:
                print('Loading error in extension ' + str(extension) + ' of image ' + str(img) + 'in load the data.')
                continue

            try:
                path_csvfile = f"{img_path}PeaksParams_{runID}_CCD_{extension}.csv"
                df = pd.read_csv(path_csvfile)
                Sigma = np.mean(df['sigma']) # ADU
            except:
                Sigma = Gain * 0.2 # ADU
                print("PeaksParams NOT founded")
            
            # dataCal, sigma = data_calibrated_NSAMP(active_area=data, gain=Gain, ratio_keV=ratio_keV, unidades= units, sigma_ADUs = Sigma)
            dataCal = data # If the images are calibrated
            sigma = Sigma/Gain
            threshold = n_sigmas * sigma

            label_img, n_events = sk.measure.label(dataCal > threshold, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            list_totalEvents.append(n_events)
            total_clusters += n_events
            
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)

            dict_lists = muon_filter(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, 
                                     prop=prop, Solidit=Solidit, Elipticity=Elip, dedl_min= dedl_value_min)
            
            for index in range(0, len(dict_lists["muons"]["l"])):
                l_val = np.array(dict_lists["muons"]["l"][index], dtype=np.float32)
                dedl_val = np.array(dict_lists["muons"]["dedl"][index], dtype=np.float32)
                charge_muon_val = np.array(dict_lists["muons"]["charge_muons"][index], dtype=np.float32)
                mean_charge_val = np.array(dict_lists["muons"]["mean_charge"][index], dtype=np.float32)
                median_charge_val = np.array(dict_lists["muons"]["median_charge"][index], dtype=np.float32)
                elip_val = np.array(dict_lists["muons"]["elip"][index], dtype=np.float32)
                sol_val = np.array(dict_lists["muons"]["sol"][index], dtype=np.float32)
                theta_val = np.array(dict_lists["muons"]["theta"][index], dtype=np.float32)
                phi_val = np.array(dict_lists["muons"]["phi"][index], dtype=np.float32)

                # img_patch = np.array(dict_lists["muons"]["image"][index], dtype=np.float32)
                img_patch = dict_lists["muons"]["image"][index].astype(np.float32)

                dict_to_save_pkl[f"extension_{extension}"]["muons"]["l"].append(l_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["dedl"].append(dedl_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["charge_muons"].append(charge_muon_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["mean_charge"].append(mean_charge_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["median_charge"].append(median_charge_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["theta"].append(theta_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["phi"].append(phi_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["elip"].append(elip_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["sol"].append(sol_val)
                dict_to_save_pkl[f"extension_{extension}"]["muons"]["image"].append(img_patch)

        print('Imagen ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    Final = datetime.datetime.now()
    print('End time: ', Final)
    print('Elapsed time: ', Final-Inicio)

    dict_to_save_pkl['Detected_Muons'] = total_clusters
    dict_to_save_pkl[f"RUNID"] = {"muons" : {}}

    print(num_images)
    print(f'Total Detected Muons: {total_clusters}')

    if units == 0:
        file_name = f'dict_muons_RUNID_{runID}_IMGTYPE_{imgType}_NSAMP{NSAMP}_NIMGS_{total_images}_SIZE_{x_max-border}x{y_max-border}_NSIGMAS_{n_sigmas}_ADUs.pkl'
    elif units == 1:
        file_name = f'dict_muons_RUNID_{runID}_IMGTYPE_{imgType}_NSAMP{NSAMP}_NIMGS_{total_images}_SIZE_{x_max-border}x{y_max-border}_NSIGMAS_{n_sigmas}_electrons.pkl'
    elif units == 2:
        file_name = f'dict_muons_RUNID_{runID}_IMGTYPE_{imgType}_NSAMP{NSAMP}_NIMGS_{total_images}_SIZE_{x_max-border}x{y_max-border}_NSIGMAS_{n_sigmas}_keV.pkl'

    file_object_histogram = open(folder_path + file_name, 'wb')
    pkl.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', folder_path + file_name, ' as a binary file.\nTo open use library "pickle". ')

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

