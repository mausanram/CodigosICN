from astropy.io import fits
import numpy as np
import sys
import skimage as sk
import datetime
import pickle as pkl
import os

module_path = os.path.abspath('../functions') # Adjust relative path as needed
# # Or use an absolute path: module_path = '/home/your_user/projects/utils'

if module_path not in sys.path:
    sys.path.append(module_path)

from functionsNEW import *

## CONSTANTES ## 
current_path = os.getcwd()
ratio_keV = 0.00368 # Search the CONNIE's official value

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = keV)
units = 1
n_sigmas = 5
numero_bins = 600

def main(argObj):
    CCD_array = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
    # CCD_array = (1,2,5,6,7,9,14,15,16)
    CCD_bad = (4, 7, 9, 11, 13, 14, 15)
    

    # Gain_list = [676.84, 683.27, 0, 0, 700.82, 693.]

    nerr_img = 0
    nerr_ext = 0

    total_images = len(argObj)
    total_clusters = 0

    image_in_bucle = 0

    Inicio = datetime.datetime.now()
    dict_to_save_pkl = {'Energy_Units' : units, 'Num_Images' : total_images,}

    for extension in CCD_array:
        if extension in CCD_bad:
            continue
        dict_to_save_pkl[f"extension_{extension}"] = {"charge" : []}

    # print(dict_to_save_pkl["extension_1"]["charge"])

    print('Start time: ', Inicio)

    for img in argObj:
        try:
            hdu_list = fits.open(img)
            image_in_bucle += 1
            runId = img[0].split("/")[4].split("_")[1]

        except:
            nerr_img = nerr_img + 1
            print('Loading error in image ' + str(img) + 'in open the image.')
            continue
        
        for extension in CCD_array:
            print(f'Image {image_in_bucle}/{total_images} Extension {extension}/16  ', end='\r')
            
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
                path_csvfile = f"../../images/CONNIE_MCM/RUN_214/PeaksParams_214_CCD_{extension}.csv"
                df = pd.read_csv(path_csvfile)
                Sigma = np.mean(df['sigma']) # ADU

            except:
                Sigma = Gain * 0.2 # ADU

            # dataCal, sigma = data_calibrated_NSAMP(active_area=data, gain=Gain, ratio_keV=ratio_keV, unidades= units, sigma_ADUs = Sigma)
            dataCal = data # If the images are calibrated
            sigma = Sigma/Gain
            threshold = n_sigmas * sigma

            # print(f"sigma:{sigma}")
            # del oScan

            label_img, n_events = sk.measure.label(dataCal > threshold, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            # print(n_events)
            total_clusters += n_events

            list_charge = all_cluster(dataCal=dataCal, label_img=label_img, nlabels_img=n_events, prop=prop)

            print(f"Clusters detected: {len(list_charge)} \n")
            for index in np.arange(0, len(list_charge)):
                dict_to_save_pkl[f"extension_{extension}"]["charge"].append(list_charge[index])    


    for extension in CCD_array:
        if extension in CCD_bad:
            continue
        raw_list = dict_to_save_pkl[f"extension_{extension}"]["charge"]
        dict_to_save_pkl[f"extension_{extension}"]["charge"] = np.array(raw_list, dtype=np.float32)        

    # print(total_clusters)
    dict_to_save_pkl['All_clusters_detected'] = total_clusters
    Final = datetime.datetime.now()

    print('End time: ', Final)
    print('Elapsed time: ', Final-Inicio)

    if units == 0:
        file_name = f'dict_energy_ALLclusters_RUNID{runId}_NSAMP{NSAMP}_NImgs_{total_images}_SIZE_{x_max-border}x{y_max-border}_NSIGMAS_{n_sigmas}_ADUs.pkl'
    elif units == 1:
        file_name = f'dict_energy_ALLclusters_RUNID{runId}_NSAMP{NSAMP}_NImgs_{total_images}_SIZE_{x_max-border}x{y_max-border}_NSIGMAS_{n_sigmas}_electrons.pkl'
    elif units == 2:
        file_name = f'dict_energy_ALLclusters_RUNID{runId}_NSAMP{NSAMP}_NImgs_{total_images}_SIZE_{x_max-border}x{y_max-border}_NSIGMAS_{n_sigmas}_keV.pkl'

    folder_path = "./RUN214/"
    file_object_histogram = open(folder_path+file_name, 'wb')
    pkl.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print(f'Dictionary saved in {folder_path+file_name} as a binary file. To open use the library "pickle". ')
    
    Eventos_Totales = f'Total Clusters Detected: {total_clusters}'
    # eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    img_err = f'Imágenes con error al cargar: ' + str(nerr_img)
    ext_err = f'Error en fit de extension: ' + str(nerr_ext)
    # relacion = total_events / num_muons
    

    # print(img_err)
    # print(ext_err)
    print(Eventos_Totales)
    # print('Number of total img: ', n_total_img)
    # print('Number of total ext: ', n_total_ext)


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

