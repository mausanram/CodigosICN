from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys
import skimage as sk
import datetime
import pickle
import os

from functions_MuonsNSAMP1 import *

# from ROOT import *

## CONSTANTES ## 
current_path = os.getcwd()

## Datos de la CCd
CCD_depth = 725 #micras
px_to_cm = 0.0015
px_to_micras = 15
micra_to_cm = 1 / 10000

## Datos del filtro de muones GENERAL
Solidit = 0.65
Elip = 0.65
dedl_value_min = 1400

## Datos del filtro POR EXTENSIÓN
list_Elip = [0.65, 0.65, 0, 0.65]
list_Solidit = [0.65, 0.65, 0, 0.65]

ratio_keV = 0.0036

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
n_sigmas = 13
numero_bins = 600

def Gaussian2(x,m,s,g,a1,a2): #data, mean, sigma, gain, height1, heigth2
    return a1*np.exp(-1/2*((x-m)/s)**2)+a2*np.exp(-1/2*((x-m-g)/s)**2)

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7] ## Ganancia para 324nsamp 

    # expgain = [110.78158608889959, 144.4661840232508, 100,  63.730700893071976] ## 


    list_totalEvents = []

    ### ===== List og all events === ##
    list_charge_of_all_extension_1 = []
    list_charge_of_all_extension_2 = []
    list_charge_of_all_extension_4 = []

    list_elip_of_all_extension_1 = []
    list_elip_of_all_extension_2 = []
    list_elip_of_all_extension_4 = []

    list_sol_of_all_extension_1 = []
    list_sol_of_all_extension_2 = []
    list_sol_of_all_extension_4 = []

    ### ============================ ###


    ### ===== List for muons ===== ###
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

    list_elip_extension_2 =[]
    list_elip_extension_1 = []
    list_elip_extension_4 = []

    list_sol_extension_2 =[]
    list_sol_extension_1 = []
    list_sol_extension_4 = []

    list_fit_gain_2 = []
    list_fit_gain_1 = []
    list_fit_gain_4 = []

    list_datamasked_extension_2 =[]
    list_datamasked_extension_1 = []
    list_datamasked_extension_4 = []

    ### ========================= ###

    nerr_img = 0
    nerr_ext = 0

    nerr_ext1 = 0
    nerr_ext2 = 0
    nerr_ext4 = 0


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
            nerr_img = nerr_img + 1
            print('Loading error in image ' + str(img) + 'in open the image.')
            continue
        
        for extension in (0,1,3):
            # extension = 3
            # extension = 1
            Elip = list_Elip[extension]
            Solidit = list_Solidit[extension]

            try :
                # print('Voy a obtener el OsCan y el active area')
                data = hdu_list[extension].data[:250,10:539]
                oScan = hdu_list[extension].data[:250,539:]

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
            
            try:
                dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=numero_bins, 
                                                    Bins_fit=numero_bins,make_figure_flag=False, range_fit=[-50, 360])

                sig_ADUs = dict_popt['sigma']
                Offset = dict_popt['Offset']
                Gain = dict_popt['Gain']
                Prob = dict_popt['Prob']
                
                if Prob < 0.05:
                    del_Bin = 500
                    dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                        Bins_fit=del_Bin, make_figure_flag=False, range_fit=[-30, 390])
                    sig_ADUs = dict_popt['sigma']
                    Offset = dict_popt['Offset']
                    Gain = dict_popt['Gain']
                    Prob = dict_popt['Prob']
                    
                    if Prob < 0.05:
                        del_Bin = 400
                        dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                            Bins_fit=del_Bin, make_figure_flag=False, range_fit=[-30, 390])
                        
                        sig_ADUs = dict_popt['sigma']
                        Offset = dict_popt['Offset']
                        Gain = dict_popt['Gain']
                        Prob = dict_popt['Prob']
                    
                        
                        if Prob < 0.05:
                            del_Bin = 300
                            dict_popt = oScan_fit_NSAMP324_ROOT(extensión=extension, active_area=true_active_area, oScan=oScan, Bins=del_Bin, 
                                                                Bins_fit=del_Bin, make_figure_flag=False, range_fit=[-50, 400])
                            
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

                                print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                                continue

                # if Gain < 100 or Gain > 240:
                #     ### Aquí se deberá poner la ganancia promedio de cada extensión una vez que se obtenga de muchas imágenes
                #     print('Fit gain error in extension ' + str(extension) + ' of image ' + str(img))

            except:
                print('Fit error in extension ' + str(extension) + ' of image ' + str(img))
                continue
            
            
            dataCal, sigma = data_calibrated_NSAMP(active_area=true_active_area, extension=extension, gain=Gain, ratio_keV=ratio_keV, 
                                                   unidades= units, offset=Offset, sigma_ADUs = sig_ADUs)
            
            fondo_value = n_sigmas * sigma
            
            del oScan
            

            label_img, n_events = sk.measure.label(dataCal > fondo_value, connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            
            list_totalEvents.append(n_events)
            # print(nlabels_img)
            # list_labels.append(label_img)
            # list_EventsNumber.append(n_events)
            
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)
            valor_promedio_fondo = fondo.data.mean()

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

            if extension == 0: 
                for index in np.arange(0, len(DeltaEL)):
                    ### ===== Muons ===== ###
                    list_DeltaEL_extension_1.append(DeltaEL[index])
                    list_EventCharge_extension_1.append(list_charge[index])
                    list_DeltaL_extension_1.append(DeltaL[index])
                    list_theta_extension_1.append(list_theta[index])
                    list_phi_extension_1.append(list_phi[index])
                    list_elip_extension_1.append(list_elip[index])
                    list_sol_extension_1.append(list_sol[index])
                    list_fit_gain_1.append(Gain)
                    list_datamasked_extension_1.append(list_datamasked[index])

                for index in np.arange(0, len(list_charge_all_events)):
                    ### ==== All events ==== ###
                    list_charge_of_all_extension_1.append(list_charge_all_events[index])
                    list_elip_of_all_extension_1.append(list_elip_all[index])
                    list_sol_of_all_extension_1.append(list_sol_all[index])
                    ### ==================== ###
                    
            if extension == 1: 
                for index in np.arange(0, len(DeltaEL)):
                    list_DeltaEL_extension_2.append(DeltaEL[index])
                    list_EventCharge_extension_2.append(list_charge[index])
                    list_DeltaL_extension_2.append(DeltaL[index])
                    list_theta_extension_2.append(list_theta[index])
                    list_phi_extension_2.append(list_phi[index])
                    list_elip_extension_2.append(list_elip[index])
                    list_sol_extension_2.append(list_sol[index])
                    list_fit_gain_2.append(Gain)
                    list_datamasked_extension_2.append(list_datamasked[index])

                for index in np.arange(0, len(list_charge_all_events)):
                    list_charge_of_all_extension_2.append(list_charge_all_events[index])
                    list_elip_of_all_extension_2.append(list_elip_all[index])
                    list_sol_of_all_extension_2.append(list_sol_all[index])
            
            if extension == 3: 
                for index in np.arange(0, len(DeltaEL)):
                    list_DeltaEL_extension_4.append(DeltaEL[index])
                    list_EventCharge_extension_4.append(list_charge[index])
                    list_DeltaL_extension_4.append(DeltaL[index])
                    list_theta_extension_4.append(list_theta[index])
                    list_phi_extension_4.append(list_phi[index])
                    list_elip_extension_4.append(list_elip[index])
                    list_sol_extension_4.append(list_sol[index])
                    list_fit_gain_4.append(Gain)
                    list_datamasked_extension_4.append(list_datamasked[index])

                for index in np.arange(0, len(list_charge_all_events)):
                    list_charge_of_all_extension_4.append(list_charge_all_events[index])
                    list_elip_of_all_extension_4.append(list_elip_all[index])
                    list_sol_of_all_extension_4.append(list_sol_all[index])

        print('Imagen ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
        del hdu_list              

    num_muons = len(list_EventCharge_extension_1) + len(list_EventCharge_extension_2) + len(list_EventCharge_extension_4)

    dict_to_save_pkl = {'Num_Images' : total_images , 'All_Muons_Detected' : num_muons, 'Energy_Units' : units, 
                        'Elipticity' : list_Elip, 'Solidity' : list_Solidit, 'Fit_errors' : (nerr_ext1, nerr_ext2, nerr_ext4),

                        'extension_1' : {'charge' : list_EventCharge_extension_1, 'deltaEL' : list_DeltaEL_extension_1,
                                         'deltaL' : list_DeltaL_extension_1, 'all_events' : list_charge_of_all_extension_1,
                                         'theta': list_theta_extension_1, 'phi': list_phi_extension_1, 'gain' : list_fit_gain_1,
                                         'elip' : list_elip_extension_1, 'sol' : list_sol_extension_1,
                                         'all_events_elip' : list_elip_of_all_extension_1, 'all_events_sol' : list_sol_of_all_extension_1,
                                         'datamasked' : list_datamasked_extension_1},

                        'extension_2' : {'charge' : list_EventCharge_extension_2, 'deltaEL' : list_DeltaEL_extension_2, 
                                        'deltaL' : list_DeltaL_extension_2, 'all_events' : list_charge_of_all_extension_2,
                                        'theta': list_theta_extension_2, 'phi': list_phi_extension_2,'gain' : list_fit_gain_2, 
                                        'elip' : list_elip_extension_2, 'sol' : list_sol_extension_2,
                                        'all_events_elip' : list_elip_of_all_extension_2, 'all_events_sol' : list_sol_of_all_extension_2,
                                        'datamasked' : list_datamasked_extension_2},

                        'extension_4' : {'charge' : list_EventCharge_extension_4, 'deltaEL' : list_DeltaEL_extension_4, 
                                         'deltaL' : list_DeltaL_extension_4, 'all_events' : list_charge_of_all_extension_4,
                                         'theta': list_theta_extension_4, 'phi': list_phi_extension_4, 'gain' : list_fit_gain_4, 
                                         'elip' : list_elip_extension_4, 'sol' : list_sol_extension_4, 
                                         'all_events_elip' : list_elip_of_all_extension_4, 'all_events_sol' : list_sol_of_all_extension_4,
                                         'datamasked' : list_datamasked_extension_4}}

    total_events = sum(list_totalEvents)
    Final = datetime.datetime.now()

    print('Hora del final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print(num_images)
    Eventos_Totales = 'Eventos Detectados en Total: ' +  str(total_events)
    eventos_rectos = 'Muones Detectados: ' + str(num_muons)
    img_err = 'Imágenes con error al cargar: ' + str(nerr_img)
    ext_err = 'Error en fit de extensiones: ' + str(nerr_ext)
    # relacion = total_events / num_muons
    
    # eventos_circulares = 'Muones Circulares Detectados: ' + str(len(list_EventosCirc))
    # print('Número de elementos de la lista "list_EventCharge_AllExtensions": ', len(list_EventCharge_AllExtensions))
    # print('elementos de la lista "list_EventCharge_AllExtensions":', list_EventCharge_AllExtensions)
    print(img_err)
    print(ext_err)
    print(Eventos_Totales)
    print(eventos_rectos)
    print('Error in extension 1, 2, 4 fits: ', nerr_ext1, nerr_ext2, nerr_ext4)
    

    if units == 0:
        file_name = 'dict_muons_NSAMP324_Extensions_1_to_4_Imgs_' + str(len(argObj)) + \
            '_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_NSIGMAS_' + str(n_sigmas) + '_ADUs.pkl'
    
    elif units == 1:
        file_name = 'dict_muons_NSAMP324_Extensions_1_to_4_Imgs_' + str(len(argObj)) + \
            '_Sol_' + str(Solidit) + '_Elip_'+str(Elip) + '_NSIGMAS_' + str(n_sigmas) + '_electrons.pkl'
    
    elif units == 2:
        file_name = 'dict_muons_NSAMP324_Extensions_1_2_4_NIMGS_' + str(len(argObj)) + \
            '_SOL_' + str(Solidit) + '_ELIP_'+str(Elip) + '_NSIGMAS_' + str(n_sigmas) + \
            '_DEDL_' + str(dedl_value_min) + '_SIZE_250x539_KeV_n.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')

    # plt.show() 


if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)

