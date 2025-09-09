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

ratio_keV = 0.00368

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
nsigmas_for_seed = 10
nsigmas_for_skirts = 8
 
numero_bins = 1000
# numero_bins = 600

def main(argObj):
    list_gain_extension_2 = []
    list_gain_extension_1 = []
    list_gain_extension_4 = []

    list_gainerr_extension_2 = []
    list_gainerr_extension_1 = []
    list_gainerr_extension_4 = []

    list_sig_extension_2 = []
    list_sig_extension_1 = []
    list_sig_extension_4 = []

    nused_img_ext1 = 0
    nused_img_ext2 = 0
    nused_img_ext4 = 0


    nerr_img = 0
    nerr_ext = 0

    total_images = len(argObj)
    image_in_bucle = 0

    set_blacklist = set()


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
                max_y = 250
                max_x = 539

                true_active_area = hdu_list[extension].data[:max_y,10:max_x]
                oScan = hdu_list[extension].data[:max_y,max_x:]

            except:
                print('Loading error in extension ' + str(extension + 1) + ' of image ' + str(img) + 'in load the data.')
                continue

            ### Change the range for any kind of image
            # Range_fit = [-50, 350]  # FOr Fe-55
            # Range_fit = [-100, 270] # For Fe-55 & Cs-137
            Range_fit = [-50, 350] # For muons

            # file_name = 'dict_mean_gains_NSAMP324.pkl'
            file_name = 'dict_mean_gains_NSAMP200.pkl'

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
                print('Fit double gaussian error in extension ' + str(extension + 1) + ' of image ' + str(img))
                set_blacklist.add(str(img))
                continue

            Bins = 300
            Bins_fit = Bins
            # Range_fit = [-100, 400]

            Range_fit_1 = [-70, 115]
            Range_fit_2 = [195, 350]

            hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins,  range=(oScan.min(), 18000))
            offset = bins_edges[np.argmax(hist)]
            Overscan_plane = oScan - offset

            fgaus_fir = TF1("gaus1","gaus", Range_fit_1[0], Range_fit_1[1],3) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)
            fgaus_sec = TF1("gaus2","gaus", Range_fit_2[0], Range_fit_2[1],3)


            h3=TH1F("histogram", r"Distribucion del Overscan",Bins_fit, -200, 400)
            for pixel_value in Overscan_plane.flatten():
                h3.Fill(pixel_value)

            fgaus_fir.SetParameters(200, 20, 60) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
            fgaus_sec.SetParameters(100, 200, 60)

            h3.Fit(fgaus_fir, "RNQ")
            h3.Fit(fgaus_sec, "RNQ")

            true_gain = fgaus_sec.GetParameters()[1] - fgaus_fir.GetParameters()[1]
            err_true_gain = fgaus_sec.GetParError(1) + fgaus_fir.GetParError(1)
            sigma = fgaus_fir.GetParameters()[2]
            # print('Ext ' + str(extension) + ':', true_gain)

            if 170 < true_gain < 210:
                # print('Fit done')
                if extension == 0:
                    nused_img_ext1+=1
                    list_gain_extension_1.append(true_gain)
                    list_gainerr_extension_1.append(err_true_gain)
                    list_sig_extension_1.append(sigma)
                if extension == 1:
                    nused_img_ext2+=1
                    list_gain_extension_2.append(true_gain)
                    list_gainerr_extension_2.append(err_true_gain)
                    list_sig_extension_2.append(sigma)
                if extension == 3:
                    nused_img_ext4+=1
                    list_gain_extension_4.append(true_gain)
                    list_gainerr_extension_4.append(err_true_gain)
                    list_sig_extension_4.append(sigma)
                
                print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')

            else:
                print('Error individual gaussians fit')
                print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
                continue


    gain_mean_ext1 = 0
    err_gain_mean_ext1 = 0
    sig_mean_ext1 = 0
    for index in range(0, len(list_gain_extension_1)):
        gain = list_gain_extension_1[index]
        err_gain = list_gainerr_extension_1[index]
        gain_mean_ext1 += (gain * (1 / err_gain**2))/((1 / err_gain**2))

        err_gain_mean_ext1 += np.sqrt(1 / (1 / err_gain**2))
        sig_mean_ext1 += list_sig_extension_1[index]



    gain_mean_ext2 = 0
    err_gain_mean_ext2 = 0
    sig_mean_ext2 = 0
    for index in range(0, len(list_gain_extension_2)):
        gain = list_gain_extension_2[index]
        err_gain = list_gainerr_extension_2[index]
        gain_mean_ext2 += (gain * (1 / err_gain**2))/((1 / err_gain**2))

        err_gain_mean_ext2 += np.sqrt(1 / (1 / err_gain**2))
        sig_mean_ext2 += list_sig_extension_2[index]

    gain_mean_ext4 = 0
    err_gain_mean_ext4 = 0
    sig_mean_ext4 = 0
    for index in range(0, len(list_gain_extension_4)):
        gain = list_gain_extension_4[index]
        err_gain = list_gainerr_extension_4[index]
        gain_mean_ext4 += (gain * (1 / err_gain**2))/((1 / err_gain**2))

        err_gain_mean_ext4 += np.sqrt(1 / (1 / err_gain**2))
        sig_mean_ext4 += list_sig_extension_4[index]


    print('Number of elements per ext: ', nused_img_ext1, nused_img_ext2, nused_img_ext4)

    dict_gains = {'extension_1' : {'Gain' : gain_mean_ext1/len(list_gain_extension_1), 'Err_gain' : err_gain_mean_ext1/len(list_gain_extension_1), 
                                   'Sigma' : sig_mean_ext1/len(list_gain_extension_1)}, 
                  'extension_2' : {'Gain' : gain_mean_ext2/len(list_gain_extension_2), 'Err_gain' : err_gain_mean_ext2/len(list_gain_extension_2),
                                   'Sigma' : sig_mean_ext1/len(list_gain_extension_2)},
                  'extension_4' : {'Gain' : gain_mean_ext4, 'Err_gain' : err_gain_mean_ext4, 'Sigma' : sig_mean_ext1} }
    
    print('The main gain of extension 1 is: ', dict_gains['extension_1']['Gain'], ' +- ', dict_gains['extension_1']['Err_gain'], ' & Sigma: ', dict_gains['extension_1']['Sigma'], ' ADU/e-')
    print('The main gain of extension 2 is: ', dict_gains['extension_2']['Gain'], ' +- ', dict_gains['extension_2']['Err_gain'], ' & Sigma: ', dict_gains['extension_2']['Sigma'], ' ADU/e-')
    print('The main gain of extension 4 is: ', dict_gains['extension_4']['Gain'], ' +- ', dict_gains['extension_4']['Err_gain'], ' & Sigma: ', dict_gains['extension_4']['Sigma'], ' ADU/e-')



    file_object_histogram = open(file_name, 'wb')
    pkl.dump(dict_gains, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    file_object = open('set_blacklist.pkl', 'wb')
    pkl.dump(dict_gains, file_object) ## Save the dictionary with all info 
    file_object.close()
    print('The blacklist is in set_blacklist.pkl, use the pickle library to open')
        

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)