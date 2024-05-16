from functions_py import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
import pandas as pd 
import skimage as sk
import random
import time

### Distribución Gaussiana ###
def gaussian(x, a, mean, sigma): 
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

### Distribución de Doble Gaussiana ###
def Gaussian2(x,m,s,g,a1,a2): #data, mean, sigma, gain, height1, heigth2
    return a1*np.exp(-1/2*((x-m)/s)**2)+a2*np.exp(-1/2*((x-m-g)/s)**2)


### Distribución aproximada de Landau ###
def Landau(x,a, MP,xi):
    C1 = a/np.sqrt((2 * np.pi))
    C2 = np.exp(-((x-MP)-1)/xi)
    C3 = np.exp((-0.5 * (((x-MP)-1)/xi + C2 )))

    # C1 = 1/np.sqrt((2 * np.pi))
    # C3 = np.exp((-0.5 * (x + np.exp(-x))))
    return  C1 * C3

def oScan_fit_NSAMP1(extensión, active_area, oScan, Bins, make_figure_flag = False) -> dict:
    Maxfev = 1000000
    P0=[10, 2000, 900]

    if make_figure_flag:
        fig_all, axs_all = plt.subplots(1, 1, figsize=(10, 10))
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins)
        offset = bins_edges[np.argmax(hist)]
        print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset
        median_oScan = np.median(Overscan_plane.flatten())
        max_oScan = np.max(Overscan_plane.flatten())
        min_oScan = np.min(Overscan_plane.flatten())

        diff = (max_oScan - abs(median_oScan)) / 2
        # print(median_oScan, diff, min_oScan)

        # Range = (min_oScan, max_oScan)
        Range = (-3000, 3000)
        # if 30 * abs(median_oScan) <  diff:
        #     Range = (min_oScan, diff)
        #     # print(Range)

        # else:
        #     Range = (min_oScan, diff * 2)
        #     # print(Range)

        bin_heights, bin_borders, _ = axs_all.hist(Overscan_plane.flatten(), bins = Bins,range = Range , label="Pixeles del Overscan")
        bin_centers = np.zeros(len(bin_heights), dtype=float)
        offset_fit = bin_borders[np.argmax(bin_heights)]

        for p in range(len(bin_heights)):
            bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

        # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
        xmin_fit, xmax_fit = bin_centers[0], bin_centers[-1]
        # print(xmin_fit, xmax_fit)

        bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
        bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
 
        popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, maxfev=Maxfev, p0 = [1,100,100])		# Fit histogram with gaussian
        axs_all.plot(bin_centers, gaussian(bin_centers, *popt), 'k', label = 'Ajuste Gaussiano')	

        dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]), 'Offset' : offset, 'Pcov' : pcov, 'Bins' : [bin_centers, bin_heights]}
        print('Centroide: ',popt[1], ' Amplitud: ', popt[0], 'sigma: ', abs(popt[2])) #gaussian(x, a, mean, sigma)

        axs_all.set_title("Distribución de pixeles del Overscan")
        axs_all.legend()
        plt.show()
        
    else:
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins)
        offset = bins_edges[np.argmax(hist)]

        Overscan_plane = oScan - offset
        median_oScan = np.median(Overscan_plane.flatten())
        max_oScan = np.max(Overscan_plane.flatten())
        min_oScan = np.min(Overscan_plane.flatten())

        diff = (max_oScan - abs(median_oScan)) / 2
        # print(median_oScan, diff, min_oScan)


        # if 40 * abs(median_oScan) <  diff:
        #     Range = (min_oScan, diff)
        #     # print(Range)

        # else:
        #     Range = (min_oScan, diff * 2)
        #     # print(Range)

        Range = (min_oScan, max_oScan)

        bin_heights, bin_borders = np.histogram(Overscan_plane.flatten(), bins= Bins, range = Range) #'auto'
        bin_centers = np.zeros(len(bin_heights), dtype=float)
        offset_fit = bin_borders[np.argmax(bin_heights)]

        for p in range(len(bin_heights)):
            bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

        # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
        xmin_fit, xmax_fit = bin_centers[0], bin_centers[-1]
        bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
        bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

        popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, maxfev = Maxfev, p0 = [1,100,100])		# Fit histogram with gaussiano')
        dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]),'Offset' : offset, 'Pcov' : pcov}
    
    return dict_popt

def oScan_fit_NSAMP324(extensión, active_area, oScan, Bins, make_figure_flag = False) -> dict:
    Maxfev = 100000
    P0=[10, 2000, 900]

    if make_figure_flag:
        fig_all, axs_all = plt.subplots(1, 1, figsize=(10, 10))
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins)
        offset = bins_edges[np.argmax(hist)]
        print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset
        median_oScan = np.median(Overscan_plane.flatten())
        max_oScan = np.max(Overscan_plane.flatten())
        min_oScan = np.min(Overscan_plane.flatten())

        diff = (max_oScan - abs(median_oScan)) / 2
        # print(median_oScan, diff, min_oScan)

        # Range = (min_oScan, max_oScan)
        if 30 * abs(median_oScan) <  diff:
            Range = (min_oScan, diff)
            # print(Range)

        else:
            Range = (min_oScan, diff * 2)
            # print(Range)

        bin_heights, bin_borders, _ = axs_all.hist(Overscan_plane.flatten(), bins = Bins,range = Range , label="Pixeles del Overscan")
        bin_centers = np.zeros(len(bin_heights), dtype=float)
        offset_fit = bin_borders[np.argmax(bin_heights)]

        for p in range(len(bin_heights)):
            bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

        # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
        xmin_fit, xmax_fit = bin_centers[0], bin_centers[-1]
        # print(xmin_fit, xmax_fit)

        bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
        bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
 
        popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, maxfev=Maxfev, p0 = [1,100,100])		# Fit histogram with gaussian
        axs_all.plot(bin_centers, gaussian(bin_centers, *popt), 'k', label = 'Ajuste Gaussiano')	

        dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]), 'Offset' : offset}
        print('Centroide: ',popt[1], ' Amplitud: ', popt[0], 'sigma: ', abs(popt[2])) #gaussian(x, a, mean, sigma)

        axs_all.set_title("Distribución de pixeles del Overscan")
        axs_all.legend()
        plt.show()
        
    else:
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins)
        offset = bins_edges[np.argmax(hist)]

        Overscan_plane = oScan - offset
        median_oScan = np.median(Overscan_plane.flatten())
        max_oScan = np.max(Overscan_plane.flatten())
        min_oScan = np.min(Overscan_plane.flatten())

        diff = (max_oScan - abs(median_oScan)) / 2
        # print(median_oScan, diff, min_oScan)


        if 40 * abs(median_oScan) <  diff:
            Range = (min_oScan, diff)
            # print(Range)

        else:
            Range = (min_oScan, diff * 2)
            # print(Range)


        bin_heights, bin_borders = np.histogram(Overscan_plane.flatten(), bins= Bins, range = Range) #'auto'
        bin_centers = np.zeros(len(bin_heights), dtype=float)
        offset_fit = bin_borders[np.argmax(bin_heights)]

        for p in range(len(bin_heights)):
            bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

        # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
        xmin_fit, xmax_fit = bin_centers[0], bin_centers[-1]
        bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
        bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

        popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, maxfev = Maxfev, p0 = [1,100,100])		# Fit histogram with gaussiano')
        dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]),'Offset' : offset}
    
    return dict_popt

def data_calibrated(active_area, extension, offset, list_gain, ratio_keV, unidades):
    dataP = active_area - offset

    if unidades == 0:
        data = dataP

    elif unidades == 1:
        data = dataP / list_gain[extension - 1]

    elif unidades == 2:
        data = (ratio_keV * dataP) / list_gain[extension - 1]

    return data

def event_DataFrame(dataCal, label_img, nlabels_img, prop, header, extension, unidades) -> pd.DataFrame:
    list_Runid = []
    list_ext = []
    list_Matrix_Slice_Event = []
    list_Size_Matrix_Event = []
    list_n_events = []

    #Listas a mano
    list_event_size = []
    list_charge = []
    list_mean_charge = []
    list_Barycenter= []
    list_Barycenter_charge = []
    list_n_events = []
    
    list_event_size_sk = []
    list_mean_charge_sk = []
    list_Barycenter_sk= []
    list_Barycenter_charge_sk = []
    list_n_events_sk = []

    extra = 0
    # data = hdu_list[extension-1].data
    Runid = str(int(header['RUNID']))

    fondo_mask = np.invert(label_img==0)
    fondo = ma.masked_array(dataCal,fondo_mask)
    # print(fondo)
    valor_promedio_fondo = fondo.data.mean()

    for i in range(0,nlabels_img):
        list_n_events.append(i+1)   

    for event in range(1, nlabels_img + 1):
        mask = np.invert(label_img == event)
        loc = ndimage.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

        if len(data_maskEvent)<1:
            list_Size_Matrix_Event.append('NaN')

        else:

            ## Número de imagen
            list_Runid.append(Runid)

            ## Número de Extensión
            list_ext.append(extension)

            ## Lista de Coordenadas de cada evento
            list_Size_Matrix_Event.append(str(data_maskEvent.shape[1])+'x'+str(data_maskEvent.shape[0])) ## La dimensión de la matriz del evento en pixeles
            
            event_size = 0

            # Obtiene los pixeles que componen al evento
            num_pixels = prop[event-1].num_pixels
            list_event_size_sk.append(num_pixels)

            # Obtiene la carga total del evento en electrones
            charge = data_maskEvent.sum()
            list_charge.append(charge)

            # Carga promedio en electrones
            mean_ch = prop[event-1].intensity_mean
            list_mean_charge_sk.append(round(mean_ch,3))
            

            ## Baricentro (con skmeasure)
            coordY_centerCharge, coordX_centerCharge = round(prop[event-1].centroid_local[0],4), round(prop[event-1].centroid_local[1],4)
            list_coordCenterCharge = [coordX_centerCharge, coordY_centerCharge]
            # print('Barycenter: ', prop[n_label-1].centroid_local)
            list_Barycenter_sk.append(list_coordCenterCharge)
            # print(centerMass)


            ## Carga del Baricentro (con skmeasure)
            BarycenterChage = prop[event-1].centroid_weighted_local
            # if BarycenterChage:
            list_Barycenter_charge_sk.append(BarycenterChage)
            # else: 
            #     list_Barycenter_charge_sk.append('NaN')
            # list_centerCharge.append(centerCharge)
            # print(list_centerCharge)

    ## DataFrame de Cada evento
    print('Events: '+ str(list_n_events[-1]))

    RunidFrame = pd.DataFrame(list_Runid, columns = ['Image ID'])
    ExtensionFrame = pd.DataFrame(list_ext, columns = ['Extension'])
    Event_IDFrame= pd.DataFrame(list_n_events, columns = ['Event ID'])
    Matrix_Size_EventFrame = pd.DataFrame(list_Size_Matrix_Event, columns = ['Matrix Size (px)'])
    EventSK_SizeFrame = pd.DataFrame(list_event_size_sk, columns = ['Event Size (px)'])


    if unidades == 0:
        ChargeFrame = pd.DataFrame(list_charge, columns = ['Total Charge (ADUs)'])
        MeanChargeSKFrame = pd.DataFrame(list_mean_charge_sk, columns = ['Mean Charge (ADUs)'])
        totalFrame = pd.concat([Event_IDFrame, RunidFrame, ExtensionFrame, Matrix_Size_EventFrame, EventSK_SizeFrame, ChargeFrame, MeanChargeSKFrame], axis = 1 )

        totalFrame['Barycenter (px)'] = pd.Series(list_Barycenter_sk)

        # totalFrame["Barycenter Charge (keV)"]=pd.Series(list_Barycenter_charge)
        # totalFrame["Barycenter Charge SK (e-)"]=pd.Series(list_Barycenter_charge_sk)
        TF = totalFrame.set_index('Event ID')

    elif unidades == 1: 
        ChargeFrame = pd.DataFrame(list_charge, columns = ['Total Charge (e-)'])
        MeanChargeSKFrame = pd.DataFrame(list_mean_charge_sk, columns = ['Mean Charge (e-)'])
        totalFrame = pd.concat([Event_IDFrame, RunidFrame, ExtensionFrame, Matrix_Size_EventFrame, EventSK_SizeFrame, ChargeFrame, MeanChargeSKFrame], axis = 1 )

        totalFrame['Barycenter (px)'] = pd.Series(list_Barycenter_sk)

        # totalFrame["Barycenter Charge (keV)"]=pd.Series(list_Barycenter_charge)
        # totalFrame["Barycenter Charge SK (e-)"]=pd.Series(list_Barycenter_charge_sk)
        TF = totalFrame.set_index('Event ID')

    elif unidades == 2: 
        ChargeFrame = pd.DataFrame(list_charge, columns = ['Total Charge (KeV)'])
        MeanChargeSKFrame = pd.DataFrame(list_mean_charge_sk, columns = ['Mean Charge (KeV)'])
        totalFrame = pd.concat([Event_IDFrame, RunidFrame, ExtensionFrame, Matrix_Size_EventFrame, EventSK_SizeFrame, ChargeFrame, MeanChargeSKFrame], axis = 1 )

        totalFrame['Barycenter (px)'] = pd.Series(list_Barycenter_sk)

        # totalFrame["Barycenter Charge (keV)"]=pd.Series(list_Barycenter_charge)
        # totalFrame["Barycenter Charge SK (e-)"]=pd.Series(list_Barycenter_charge_sk)
        TF = totalFrame.set_index('Event ID')

    return TF 

def muon_filter(dataCal, label_img, nlabels_img, prop, Solidit, Elipticity):
    CCD_depth = 725 ## micras
    px_to_micras = 15 ## micras
    px_to_cm = 0.0015 ## cm/px
    micra_to_cm = 1 / 10000 ## micras/cm

    DeltaEL_range_min, DeltaEL_range_max = 1, 4

    list_Muon_labels = []
    list_DeltaEL = []
    list_DeltaL = []
    list_charge = []
    list_theta = []

    for event in range(1, nlabels_img):
        mask = np.invert(label_img == event)
        loc = ndimage.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

        coordX_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[1])
        coordY_centerCharge = round(ndimage.center_of_mass(data_maskEvent)[0])

        coordY_centerCharge, coordX_centerCharge = int(prop[event-1].centroid_local[0]), int(prop[event-1].centroid_local[1])
        # print(type(coordY_centerCharge))
        # MaxValue_Event = data_maskEvent.max()
        MinValue_Event = data_maskEvent.min()
        MeanValue_Event = data_maskEvent.mean()
        # MeanValue_Event = (MaxValue_Event - MinValue_Event)/2
        Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]
        
        try:
            differval = abs(Barycentercharge - MinValue_Event) 
        except:
            differval = 0 

        rM = prop[event-1].axis_major_length
        rm = prop[event-1].axis_minor_length
        Solidity = prop[event-1].solidity
        miny, minx, maxy, maxx = prop[event-1].bbox
        Longitud_y, Longitud_x = maxy - miny , maxx - minx # px

        Diagonal_lenght= np.sqrt(Longitud_x**2 + Longitud_y**2) - np.sqrt(2) # px
        Delta_L = np.sqrt( (Diagonal_lenght * px_to_micras)**2 + (CCD_depth)**2) * micra_to_cm # cm

        if rM == 0 or rm == 0:
            continue 

        elif maxx - minx <= 3:
            continue

        elif not Barycentercharge:
            continue

        # elif differval < MeanValue_Event: #keV
        #     continue

        elif  Solidity < Solidit:
            continue 

        elif  rM >= Elipticity * rm:
            charge = data_maskEvent.sum()

            if charge > 100:
                Delta_EL = (charge)/ (Delta_L) 
                theta = np.arctan((Diagonal_lenght * px_to_cm)/(CCD_depth * micra_to_cm)) *(180 /np.pi)

                list_DeltaL.append(Delta_L)
                list_DeltaEL.append(Delta_EL)
                list_charge.append(charge)
                list_theta.append(theta)
                # print(charge, DeltaEL)

                # if DeltaEL_range_min <= DeltaEL <= DeltaEL_range_max:
                list_Muon_labels.append(event)

    return list_DeltaL, list_DeltaEL, list_charge, list_Muon_labels, list_theta 

    



        
