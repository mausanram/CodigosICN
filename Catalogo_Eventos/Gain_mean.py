import datetime
import numpy as np
import os
import pickle as pkl
from ROOT import TF1, TH1F
import sys
from astropy.io import fits


## CONSTANTES ## 
current_path = os.getcwd()
ratio_keV = 0.00368

## Unidades, número de sigmas y número de bins (en las unidades 0 = ADUs, 1 = e-, 2 = KeV)
units = 2
nsigmas_for_seed = 10
nsigmas_for_skirts = 8
 
# numero_bins = 900
numero_bins = 400

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

    list_sigerr_extension_2 = []
    list_sigerr_extension_1 = []
    list_sigerr_extension_4 = []

    nused_img_ext1 = 0
    nused_img_ext2 = 0
    nused_img_ext4 = 0

    nerr_img = 0

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
                oScan = hdu_list[extension].data[:max_y,max_x:]
                active_area = hdu_list[extension].data[:max_y,:max_x]
            except:
                print('Loading error in extension ' + str(extension + 1) + ' of image ' + str(img) + 'in load the data.')
                continue

            ## Media extraction per row
            medi_rows_value = []
            for element in range(0, max_y):
                row = oScan[element: element +1, 0: max_x]
                num_row = element + 1
                medi_value = np.median(row)
                medi_rows_value.append([medi_value])
            active_area = active_area - medi_rows_value    

            ### Change the range for any kind of image
            # Range_fit = [-65, 300]  # FOr Fe-55
            # Range_fit = [-100, 270] # For Fe-55 & Cs-137
            # Range_fit = [-50, 350] # For muons

            # file_name = 'dict_mean_gain_NSAMP324.pkl'
            # file_name = 'dict_mean_gains_NSAMP200.pkl'
            file_name = 'dict_mean_gains_Fe55_NSAMP300.pkl'

            Bins = numero_bins
            Bins_fit = Bins

            min_oScan = np.min(oScan)
            min_active_area = np.min(active_area)

            hist, bins_edges = np.histogram(oScan.flatten(), bins = Bins,  range=(min_active_area, 18000))
            offset = bins_edges[np.argmax(hist)]
            # print('Offset Value: ', offset, ' ADUs')

            Overscan_plane = oScan - offset
            acArea_plane = active_area

            # Bins = numero_bins
            # Range_fit = [-100, 400]

            ### FOr Fe-55
            # if extension == 0:
            #     Range_fit_1 = [-70, 115]
            #     Range_fit_2 = [165, 300]

            # elif extension == 1:
            #     Range_fit_1 = [-70, 115]
            #     Range_fit_2 = [140, 300]

            ### FOr Muons
            if extension == 0:
                Range_fit_1 = [-100, 110]
                Range_fit_2 = [120, 320]

            elif extension == 1:
                Range_fit_1 = [-100, 110]
                Range_fit_2 = [100, 320]

            fgaus_fir = TF1("gaus1","gaus", Range_fit_1[0], Range_fit_1[1],3) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)
            fgaus_sec = TF1("gaus2","gaus", Range_fit_2[0], Range_fit_2[1],3)


            h3=TH1F("h3", r"Distribucion del Overscan",Bins_fit, -200, 400)
            for pixel_value in Overscan_plane.flatten():
                h3.Fill(pixel_value)

            # for pixel_value in acArea_plane.flatten():
            #     h3.Fill(pixel_value)

            fgaus_fir.SetParameters(200, 20, 60) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
            fgaus_sec.SetParameters(100, 200, 60)

            h3.Fit(fgaus_fir, "RNQ")
            h3.Fit(fgaus_sec, "RNQ")

            true_gain = fgaus_sec.GetParameters()[1] - fgaus_fir.GetParameters()[1]
            err_true_gain = fgaus_sec.GetParError(1) + fgaus_fir.GetParError(1)
            sigma = fgaus_fir.GetParameters()[2]
            sigma_err = fgaus_fir.GetParError(2)
            del h3

            if 180 < true_gain < 215:
                # print('Fit done')
                if extension == 0:
                    nused_img_ext1+=1
                    list_gain_extension_1.append(true_gain)
                    list_gainerr_extension_1.append(err_true_gain)
                    list_sig_extension_1.append(sigma)
                    list_sigerr_extension_1.append(sigma_err)
                if extension == 1:
                    nused_img_ext2+=1
                    list_gain_extension_2.append(true_gain)
                    list_gainerr_extension_2.append(err_true_gain)
                    list_sig_extension_2.append(sigma)
                    list_sigerr_extension_2.append(sigma_err)
                if extension == 3:
                    nused_img_ext4+=1
                    list_gain_extension_4.append(true_gain)
                    list_gainerr_extension_4.append(err_true_gain)
                    list_sig_extension_4.append(sigma)
                    list_sigerr_extension_4.append(sigma_err)
                print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')

            else:
                if extension == 0:
                    set_blacklist.add(img)
                    print('Error individual gaussians fit in ext ' + str(extension+1) + ' of image ' + str(img))
                    # print('Gain:', true_gain)
                    print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
                    continue
                continue

    sum_gain_weight_ext1 = 0
    sum_wheith_gain_ext1 = 0
    sum_sig_weight_ext1 = 0
    sum_wheith_sig_ext1 = 0
    for index in range(0, len(list_gain_extension_1)):
        gain = list_gain_extension_1[index]
        err_gain = list_gainerr_extension_1[index]
        weight = (1 / (err_gain**2))

        sum_gain_weight_ext1 += gain*weight 
        sum_wheith_gain_ext1 += weight

        sig = list_sig_extension_1[index]
        err_sig = list_sigerr_extension_1[index]
        sig_weight = (1/(err_sig**2))

        sum_sig_weight_ext1 += sig*sig_weight
        sum_wheith_sig_ext1 += sig_weight

    if sum_wheith_sig_ext1>0 and sum_wheith_gain_ext1>0:
        mean_gain_ext1 = sum_gain_weight_ext1/sum_wheith_gain_ext1
        true_gain_error_ext1 = np.sqrt(sum_wheith_gain_ext1)

        mean_sig_ext1 = sum_sig_weight_ext1/sum_wheith_sig_ext1
        true_sig_error_ext1 = np.sqrt(sum_wheith_sig_ext1)

    sum_gain_weight_ext2 = 0
    sum_wheith_gain_ext2 = 0
    sum_sig_weight_ext2 = 0
    sum_wheith_sig_ext2 = 0
    for index in range(0, len(list_gain_extension_2)):
        gain = list_gain_extension_2[index]
        err_gain = list_gainerr_extension_2[index]
        weight = (1 / (err_gain**2))

        sum_gain_weight_ext2 += gain*weight 
        sum_wheith_gain_ext2 += weight

        sig = list_sig_extension_2[index]
        err_sig = list_sigerr_extension_2[index]
        sig_weight = (1/(err_sig**2))

        sum_sig_weight_ext2 += sig*sig_weight
        sum_wheith_sig_ext2 += sig_weight

    if sum_wheith_sig_ext2>0 and sum_wheith_gain_ext2>0:
        mean_gain_ext2 = sum_gain_weight_ext2/sum_wheith_gain_ext2
        true_gain_error_ext2 = np.sqrt(sum_wheith_gain_ext2)

        mean_sig_ext2 = sum_sig_weight_ext2/sum_wheith_sig_ext2
        true_sig_error_ext2 = np.sqrt(sum_wheith_sig_ext2)

    sum_gain_weight_ext4 = 0
    sum_wheith_gain_ext4 = 0
    sum_sig_weight_ext4 = 0
    sum_wheith_sig_ext4 = 0
    for index in range(0, len(list_gain_extension_4)):
        gain = list_gain_extension_4[index]
        err_gain = list_gainerr_extension_4[index]
        weight = (1 / (err_gain**2))

        sum_gain_weight_ext4 += gain*weight 
        sum_wheith_gain_ext4 += weight

        sig = list_sig_extension_4[index]
        err_sig = list_sigerr_extension_4[index]
        sig_weight = (1/(err_sig**2))

        sum_sig_weight_ext4 += sig*sig_weight
        sum_wheith_sig_ext4 += sig_weight
    if sum_wheith_sig_ext4>0 and sum_wheith_gain_ext4>0:
        mean_gain_ext4 = sum_gain_weight_ext4/sum_wheith_gain_ext4
        true_gain_error_ext4 = np.sqrt(sum_wheith_gain_ext4)

        mean_sig_ext4 = sum_sig_weight_ext4/sum_wheith_sig_ext4
        true_sig_error_ext4 = np.sqrt(sum_wheith_sig_ext4)

    
    with open("black_list.txt", "w") as f:
        for item in set_blacklist:
            f.write(item + "\n")
    print(f"Black List saved in black_list.txt file")

    print('Number of elements per ext: ', nused_img_ext1, nused_img_ext2, nused_img_ext4)

    dict_gains = {'extension_1' : {'Gain' : mean_gain_ext1, 'Err_gain' : true_gain_error_ext1, 'Sigma' : mean_sig_ext1, 'Err_sig' : true_sig_error_ext1}, 
                  'extension_2' : {'Gain' : mean_gain_ext2, 'Err_gain' : true_gain_error_ext2, 'Sigma' : mean_sig_ext2, 'Err_sig' : true_sig_error_ext2},
                  'extension_4' : {'Gain' : mean_gain_ext4, 'Err_gain' : true_gain_error_ext4, 'Sigma' : mean_sig_ext4, 'Err_sig' : true_sig_error_ext4} }
    
    print('The main gain of extension 1 is: ', dict_gains['extension_1']['Gain'], ' +- ', dict_gains['extension_1']['Err_gain'], ' & Sigma: ', dict_gains['extension_1']['Sigma'], ' +- ', dict_gains['extension_1']['Err_sig'], ' ADU/e-')
    print('The main gain of extension 2 is: ', dict_gains['extension_2']['Gain'], ' +- ', dict_gains['extension_2']['Err_gain'], ' & Sigma: ', dict_gains['extension_2']['Sigma'], ' +- ', dict_gains['extension_2']['Err_sig'], ' ADU/e-')
    print('The main gain of extension 4 is: ', dict_gains['extension_4']['Gain'], ' +- ', dict_gains['extension_4']['Err_gain'], ' & Sigma: ', dict_gains['extension_4']['Sigma'], ' +- ', dict_gains['extension_4']['Err_sig'], ' ADU/e-')



    file_object_histogram = open(file_name, 'wb')
    pkl.dump(dict_gains, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()
        

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
