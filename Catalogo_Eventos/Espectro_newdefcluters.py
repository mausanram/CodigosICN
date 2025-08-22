from functions_py import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
import pandas as pd 
import skimage as sk
import scipy.ndimage as nd
from array import array
from functions_MuonsNSAMP1 import *
import pickle as pkl
import time

## CONSTANTES ## 
current_path = os.getcwd()

ratio_keV = 0.00367

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
nsigmas_for_seed = 13
nsigmas_for_skirts = 10
 
numero_bins = 1000

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
                max_y = 150
                max_x = 539

                data = hdu_list[extension].data[:max_y,10:max_x]
                oScan = hdu_list[extension].data[:max_y,max_x:]

                oscan_x = oScan.shape[1]
                oscan_y = oScan.shape[0]

                mean_rows_value = []
                for element in range(0, oscan_y):
                    row = oScan[element: element +1, 0: oscan_x]
                    mean_value = np.median(row)
                    mean_rows_value.append([mean_value])

                true_active_area = data - mean_rows_value

            except:
                print('Loading error in extension ' + str(extension + 1) + ' of image ' + str(img) + 'in load the data.')
                continue

            if extension == 0:
                Gain = 187.898 # ADU/e-
                sig_ADUs = 81.0868 # ADUs
            if extension == 1:
                Gain = 192.728
                sig_ADUs = 61.626
            if extension == 3:
                Gain = 189.728
                sig_ADUs = 3081
            
            dataCal, sigma = data_calibrated_NSAMP(active_area=true_active_area, gain=Gain, ratio_keV=ratio_keV, unidades= units, sigma_ADUs = sig_ADUs)
            
            seed_charge = nsigmas_for_seed * sigma
            skirts_charge = nsigmas_for_skirts * sigma
            
            del oScan

            # dataCal = example
            xshape = dataCal.shape[1]
            yshape = dataCal.shape[0]
            n_seeds = 0
            list_clusters = []

            copy_dataCal = dataCal.copy()

            n_clusters = 0
            for pix_Y in range(0, yshape):
                for pix_X in range(0, xshape):
                    px_charge = copy_dataCal[pix_Y][pix_X]

                    if px_charge > seed_charge:
                        data_zeros = np.zeros((yshape, xshape), dtype=bool)
                        n_seeds += 1
                        set_neighbor = set()

                        neighbor_matrix = [(pix_X-1, pix_Y+1), (pix_X, pix_Y+1), (pix_X+1, pix_Y+1),
                                            (pix_X-1, pix_Y),   (pix_X, pix_Y),   (pix_X+1, pix_Y), 
                                        (pix_X-1, pix_Y-1), (pix_X, pix_Y-1), (pix_X+1, pix_Y-1)]
                            

                        for px in neighbor_matrix:
                            try:
                                if copy_dataCal[px[1]][px[0]] > skirts_charge:
                                    set_neighbor.add(px)
                                else:
                                    continue
                            except:
                                continue

                        copy_set = set_neighbor.copy()

                        set_size = len(copy_set)
                        flag_loop = True
                        ncomp = 0
                        aux_set = set()

                        while flag_loop:
                            ncomp += 1
                            other_set = set()
                            # print(len(copy_set))

                            for coord_px in copy_set:
                                # print(coord_px)
                                other_set.add(coord_px)
                                neighbor_matrix = [(coord_px[0]-1, coord_px[1]+1), (coord_px[0], coord_px[1]+1), (coord_px[0]+1, coord_px[1]+1),
                                                    (coord_px[0]-1, coord_px[1]),     (coord_px[0], coord_px[1]),    (coord_px[0]+1, coord_px[1]), 
                                                    (coord_px[0]-1, coord_px[1]-1),  (coord_px[0], coord_px[1]-1),  (coord_px[0]+1, coord_px[1]-1)]
                                
                                for px in neighbor_matrix:
                                    try:
                                        if copy_dataCal[px[1]][px[0]] > skirts_charge:
                                            copy_dataCal[px[1]][px[0]] = 0
                                            set_neighbor.add(px)
                                            aux_set.add(px)

                                    except:
                                        continue

                            for element in copy_set:
                                if element in other_set:
                                    aux_set.remove(element)


                            copy_set = aux_set.copy()

                            set_size = len(set_neighbor)

                            if ncomp == set_size:
                                flag_loop = False
                    
                        for coord_px in set_neighbor:
                            data_zeros[coord_px[1]][coord_px[0]] = bool(1)

                        data_zeros = nd.binary_dilation(data_zeros, structure=[[1,1,1],[1,1,1],[1,1,1]])
                        list_clusters.append(data_zeros)

                        copy_dataCal = ma.masked_array(copy_dataCal, mask=data_zeros)

                        del set_neighbor
                    
                    else:
                        continue
            
            if extension == 0:
                for index in range(0, len(list_clusters)):
                    sample_mask = list_clusters[index]
                    inv_mask = np.invert(sample_mask)
                    datamask = ma.masked_array(dataCal, inv_mask)

                    charge = datamask.sum()
                    list_EventCharge_extension_1.append(charge)
                
                del list_clusters
                del copy_dataCal

            if extension == 1:
                for index in range(0, len(list_clusters)):
                    sample_mask = list_clusters[index]
                    inv_mask = np.invert(sample_mask)
                    datamask = ma.masked_array(dataCal, inv_mask)

                    charge = datamask.sum()
                    list_EventCharge_extension_2.append(charge)

                del list_clusters
                del copy_dataCal

            if extension == 3:
                for index in range(0, len(list_clusters)):
                    sample_mask = list_clusters[index]
                    inv_mask = np.invert(sample_mask)
                    datamask = ma.masked_array(dataCal, inv_mask)

                    charge = datamask.sum()
                    list_EventCharge_extension_4.append(charge)

                del list_clusters
                del copy_dataCal

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
    in_path = 'dict_energy_allclusters_Fe55_NSAMP200_Extensions_1_to_4_Imgs_NEWCLUSTERIZATION_' # For Fe-55 

    if units == 0:
        file_name = in_path + str(total_images) + '_SIZE_' + str(max_y)+ 'x' + str(max_x) + '_NSIGMAS_SEED&SKIRTS_' + str(nsigmas_for_seed)  + '_' + str(nsigmas_for_skirts) + '_ADUs.pkl'
    elif units == 1:
        file_name = in_path + str(total_images) + '_SIZE_' + str(max_y)+ 'x' + str(max_x) + '_NSIGMAS_SEED&SKIRTS_' + str(nsigmas_for_seed)  + '_' + str(nsigmas_for_skirts) + '_electrons.pkl'
    elif units == 2:
        file_name = in_path + str(total_images) + '_SIZE_' + str(max_y)+ 'x' + str(max_x) + '_NSIGMAS_SEED&SKIRTS_' + str(nsigmas_for_seed)  + '_' + str(nsigmas_for_skirts) + '_KeV.pkl'

    file_object_histogram = open(file_name, 'wb')
    pkl.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(num_clusters)

    print(Eventos_Totales)


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)