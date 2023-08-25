from functions_py import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
import pandas as pd 
import sys
import skimage as sk
from scipy.stats import linregress
import datetime

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def recta(x, m, b):
    return m * x + b

def main(argObj):
    expgain = [227, 220.4, 94.72, 197.7]
    list_labels = []
    list_EventsNumber = []
    list_EventosRectos_AllExtensions = []
    list_runid = []

    fig_all, axs_all = plt.subplots(1, 1, figsize=(10, 10))
    Inicio = datetime.datetime.now()

    for img in argObj:
        hdu_list = fits.open(img)
        hdu_path = img.split('_')
        list_eventos_rectos = []

        for extension in range(0,4):
            data = hdu_list[extension].data
            header=hdu_list[extension].header
            oScan=hdu_list[extension].data[638:,530:]
            nsamp = float(header['NSAMP'])
            Runid = int(header['RUNID'])
            Exposure = float(hdu_path[15])/60 ## En minutos

            del header

            hist , bins_edges = np.histogram(oScan.flatten(), bins = 'auto')

            del oScan

            offset = bins_edges[np.argmax(hist)]
            dataP = data-offset
            dataCal = dataP/expgain[extension] ## En electrones  
            
            del hist
            del data
            del dataP

            bin_heights, bin_borders = np.histogram(dataCal.flatten(), bins= 'auto') 
            bin_centers=np.zeros(len(bin_heights), dtype=float)
            offset_fit = bin_borders[np.argmax(bin_heights)]
            for p in range(len(bin_heights)):
                bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

            xmin_fit, xmax_fit = offset_fit-(10*expgain[extension])/math.sqrt(nsamp), offset_fit+(10*expgain[extension])/math.sqrt(nsamp)			# Define fit range
            bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
            bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

            del nsamp

            try:
                popt,_ = curve_fit(gaussian, bin_centers, bin_heights)#, p0=[np.max(bin_heights), 0, 1], maxfev=100000)		# Fit histogram with gaussian
            except:
                continue
            
            del bin_heights
            del bin_centers
            del offset
            del xmin_fit
            del xmax_fit
            
            # label, n_events =ndimage.label(dataCal>6*abs(popt[2]),structure=[[1,1,1],[1,1,1],[1,1,1]])
            label_img, n_events = sk.measure.label(dataCal>6*popt[2], connectivity=2, return_num=True)
            prop = sk.measure.regionprops(label_img, dataCal)
            
            # print(nlabels_img)
            list_labels.append(label_img)
            list_EventsNumber.append(n_events)
            
            ## Obteniendo el valor promedio del fondo
            fondo_mask = np.invert(label_img == 0)
            fondo = ma.masked_array(dataCal,fondo_mask)
            valor_promedio_fondo = fondo.data.mean()
            cuentas = 0

            for event in range(1,n_events):
                mask = np.invert(label_img == event)
                loc = ndimage.find_objects(label_img == event)[0]
                data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                         mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
                
                coordX_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[1])
                coordY_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[0])
                # print(Runid, extension+1, event+1, coordX_centerCharge, coordY_centerCharge)

                MeanValue_Event = data_maskEvent.mean()
                MinValue_Event = data_maskEvent.min()

                Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]

                try:
                    differval = abs(Barycentercharge - MinValue_Event) 
                except:
                    differval = 0 

                rM = prop[event-1].axis_major_length
                rm = prop[event-1].axis_minor_length
                _, minc, _, maxc = prop[event-1].bbox

                if rM == 0 or rm == 0:
                    continue

                elif maxc - minc <= 3:
                    continue

                elif not Barycentercharge:
                    continue

                elif rM > 4.5 * rm:
                    cuentas = cuentas +1
                
                elif differval < MeanValue_Event + 100:
                    continue


                del MeanValue_Event
                del MinValue_Event
                del data_maskEvent
                del Barycentercharge
                
            list_eventos_rectos.append(cuentas)
    
        list_runid.append(Exposure)
        list_EventosRectos_AllExtensions.append(sum(list_eventos_rectos))

        del list_eventos_rectos
                

            # print(list_eventos_rectos)

    # poptR, pcovR = curve_fit(recta, list_runid, list_EventosRectos_AllExtensions)
    # x_line = np.arange(int(list_runid[0]), int(list_runid[-1])+1, 1)
    # y_line = recta(x_line, poptR[0], poptR[1])

    # residuals = list_EventosRectos_AllExtensions - recta(x_line, *poptR)
    # ss_res = np.sum(residuals**2)
    # ss_tot = np.sum((list_EventosRectos_AllExtensions-np.mean(list_EventosRectos_AllExtensions))**2)

    # RSquare = 1-(ss_res/ss_tot)
    # print('R²: ', RSquare)

    # slope, intercept, r_value, p_value, std_err = linregress(list_runid,list_EventosRectos_AllExtensions)
    # print("R²: ", r_value**2)
    Fin = datetime.datetime.now()
    process_time = Fin - Inicio

    print('Tiempo que se tardó el proceso: {0}'.format(process_time))

    axs_all.plot(list_runid, list_EventosRectos_AllExtensions, 'ok') 
    # axs_all.plot(list_runid, y_line, label = 'R^2: '+str(r_value**2))
    axs_all.set_title('Eventos Rectos Detectados por imagen')
    axs_all.set_xlabel('Exposure (min)')
    axs_all.set_ylabel('Cuentas') 
    # plt.xlim(0, 200000)  
    plt.show() 


if __name__ == "__main__":
    argObj=sys.argv[1:]
    exitcode=main(argObj)
    exit(code=exitcode)