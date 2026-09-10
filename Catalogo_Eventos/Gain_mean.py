import datetime
import numpy as np
import os
import pickle as pkl
from ROOT import TF1, TH1F
import sys
from astropy.io import fits
# import matplotlib.pyplot as plt

from functions_MuonsNSAMP1 import cleaning_actArea


## FIXES VARIABLES ## 
current_path = os.getcwd()
dict_type_experiment = {0: "Muons", 1: "Fe55", 2:"Fe55pCs137"}

## ===== WRITE WORKED EXTENSIONS ARRAY ====== ##
list_CCD_array = [1,2,4] # Use the real extension number

##  ============= SELECT IMAGES' TYPE ============ ##
##  ==(Muons == 0, Fe55 == 1, Fe55+Cs137 == 2,) == ##
type_experiment = 0

## ======= WRITE NSAMP & NBINS ======= ##
Nsamp = 300
nbins = 200

## ====== WRITE THE ACTIVE AREA's PIXELS ====== ##
min_y, max_y = 0, 250
min_x, max_x = 10, 539

try:
    int_type_exp = int(type_experiment)
except:
    print("You must enter a valid number. Finishing program...")
    exit()

if int_type_exp not in dict_type_experiment.keys():
    print(f"No images' type specified")
    exit()
else:
    string_typeE = dict_type_experiment[int_type_exp]

file_name =f"dict_mean_gains_{string_typeE}_NSAMP{Nsamp}.pkl"

def main(argObj):

    start_time = datetime.datetime.now()
    print('Start time: ', start_time)

    dict_gains = dict()
    dict_auxiliar = dict()

    Bins = nbins
    Bins_fit = Bins

    for extension in list_CCD_array:
        dict_auxiliar[f"extension_{extension}"] = {"gain": [], "gain_err": [], "sig":[], "sig_err":[], "images_used":0} 

    nerr_img = 0

    total_images = len(argObj)
    image_in_bucle = 0

    set_blacklist = set()

    for img in argObj:
        try:
            hdu_list = fits.open(img)
            image_in_bucle += 1
        except:
            nerr_img = nerr_img + 1
            print('Loading error in image ' + str(img) + 'in open the image.')
            continue
        
        for extension in list_CCD_array:                 
            extension -= 1
            try :
                actArea = hdu_list[extension].data[min_y:max_y, min_x:max_x]
                oScan = hdu_list[extension].data[min_y:max_y, max_x:]
            except:
                print('Loading error in extension ' + str(extension + 1) + ' of image ' + str(img) + 'in load the data.')
                continue

            oScan = oScan.flatten()

            lower_bound = np.percentile(oScan, 1)   # Cuts off bottom 1%
            upper_bound = np.percentile(oScan, 99)  # Cuts off top 1%

            filtered_oScan = oScan[(oScan >= lower_bound) & (oScan <= upper_bound)]
            hist, bins_edges = np.histogram(filtered_oScan, bins='auto')
            max_bin_idx = np.argmax(hist)

            lower_boundary = bins_edges[max_bin_idx]
            upper_boundary = bins_edges[max_bin_idx + 1]


            tallest_bin_data = filtered_oScan[(filtered_oScan >= lower_boundary) & (filtered_oScan <= upper_boundary)]
            # print(tallest_bin_data)

            offset = tallest_bin_data.max()
            # print(f"Offset: {offset}")
            Overscan_plane = oScan - offset
            # plt.hist(Overscan_plane, range = (-500, 400), bins=100)
            # plt.show()

            ## FOr Fe-55
            if extension == 0:
                Range_fit_1 = [-100, 60]
                Range_fit_2 = [120, 300]

            elif extension == 1:
                Range_fit_1 = [-120, 45]
                Range_fit_2 = [90, 250]

            # ### For Muons
            # if extension == 0:
            #     Range_fit_1 = [-100, 110]
            #     Range_fit_2 = [120, 320]

            # elif extension == 1:
            #     Range_fit_1 = [-100, 110]
            #     Range_fit_2 = [100, 320]

            # print("ending program...")
            # exit()

            fgaus_fir = TF1("gaus1","gaus", Range_fit_1[0], Range_fit_1[1],3) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)
            fgaus_sec = TF1("gaus2","gaus", Range_fit_2[0], Range_fit_2[1],3)


            h3=TH1F("h3", r"Distribucion del Overscan", Bins_fit, -200, 400)
            for pixel_value in Overscan_plane.flatten():
                h3.Fill(pixel_value)

            fgaus_fir.SetParameters(200, 20, 60) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
            fgaus_sec.SetParameters(100, 200, 60)

            h3.Fit(fgaus_fir, "RNQ")
            Prob_1 = fgaus_fir.GetProb()

            h3.Fit(fgaus_sec, "RNQ")
            Prob_2 = fgaus_sec.GetProb()

            true_gain = fgaus_sec.GetParameters()[1] - fgaus_fir.GetParameters()[1]
            err_true_gain = fgaus_sec.GetParError(1) + fgaus_fir.GetParError(1)
            sigma = fgaus_fir.GetParameters()[2]
            err_sigma = fgaus_fir.GetParError(2)
            del h3

            # print(f"Probabilities: {Prob_1}, {Prob_2}")
            # print(f"Gain: {true_gain}, +- {err_true_gain}")
            # exit()

            if 180 < true_gain < 215:
                dict_auxiliar[f"extension_{extension+1}"]["images_used"] = dict_auxiliar[f"extension_{extension+1}"]["images_used"] + 1
                dict_auxiliar[f"extension_{extension+1}"]["gain"].append(true_gain)
                dict_auxiliar[f"extension_{extension+1}"]["gain_err"].append(err_true_gain)
                dict_auxiliar[f"extension_{extension+1}"]["sig"].append(sigma)
                dict_auxiliar[f"extension_{extension+1}"]["sig_err"].append(err_sigma)
                print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
            else:
                if extension == 0:
                    set_blacklist.add(img)
                    print('Error individual gaussians fit in ext ' + str(extension+1) + ' of image ' + str(img))
                    # print('Gain:', true_gain)
                    print('Image ' + str(image_in_bucle) + '/' + str(total_images), end='\r')
                    continue
                continue

    for extension in list_CCD_array:
        sum_gain_weight = 0
        sum_weight_gain = 0

        sum_sig_weight = 0
        sum_weight_sig = 0

        nimages = dict_auxiliar[f"extension_{extension}"]["images_used"]
        list_gain = dict_auxiliar[f"extension_{extension}"]["gain"]
        list_gainerr =  dict_auxiliar[f"extension_{extension}"]["gain_err"]
        list_sig = dict_auxiliar[f"extension_{extension}"]["sig"]
        list_sigerr = dict_auxiliar[f"extension_{extension}"]["sig_err"]

        for index in range(0, len(list_gain)):
            gain = list_gain[index]
            err_gain = list_gainerr[index]
            weight = (1 / (err_gain**2))

            sum_gain_weight += gain*weight 
            sum_weight_gain += weight

            sig = list_sig[index]
            err_sig = list_sigerr[index]
            sig_weight = (1/(err_sig**2))

            sum_sig_weight += sig*sig_weight
            sum_weight_sig += sig_weight

        if sum_weight_sig>0 and sum_weight_gain>0:
            mean_gain = sum_gain_weight/sum_weight_gain
            true_gain_error = np.sqrt(sum_weight_gain)

            mean_sig = sum_sig_weight/sum_weight_sig
            true_sig_error = np.sqrt(sum_weight_sig)
            # print("Sig_err: ", true_sig_error)

            dict_gains[f"extension_{extension}"] = {'NImages': nimages,
                                                    'Gain' : mean_gain, 'Err_gain' : true_gain_error, 
                                                    'Sigma' : mean_sig, 'Err_sig' : true_sig_error}

    with open("black_list.txt", "w") as f:
        for item in set_blacklist:
            f.write(item + "\n")
    print(f"Black List saved in black_list.txt file")

    del dict_auxiliar

    end_time = datetime.datetime.now()
    print('End time: ', end_time)

    ellapsed_time = end_time - start_time
    print(f"Ellapsed time: {ellapsed_time} \n" )
    
    for extension in list_CCD_array:
        try:
            print(f"=========== EXTENSION { extension} ===========")
            gain = dict_gains[f'extension_{extension}']['Gain']
            err_gain = dict_gains[f'extension_{extension}']['Err_gain']
            sigma = dict_gains[f'extension_{extension}']['Sigma']
            err_sigma = dict_gains[f'extension_{extension}']['Err_sig']
            nused_imgs = dict_gains[f"extension_{extension}"]["NImages"]

            print(f"Gain: {gain} +- {err_gain} ADU/e-")
            print(f"Sigma: {sigma} +- {err_sigma} ADU")
            print(f'Extension used: {nused_imgs}')
        except:
            print(f"---- ERROR: GAIN OR SIGMA NOT DETECTED")

    file_object = open(file_name, 'wb')
    pkl.dump(dict_gains, file_object) ## Save the dictionary with all info 
    file_object.close()
        

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
