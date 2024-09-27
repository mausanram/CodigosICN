from functions_py import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd 
import skimage as sk
import scipy.ndimage as nd
import random
import time

from ROOT import TF1, TH1F, Fit, TCanvas, gStyle

### Distribución Gaussiana ###
def gaussian(x, a, mean, sigma): 
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

### Distribución de Doble Gaussiana ###
def Gaussian2(x,m,s,g,a1,a2): #data, mean, sigma, gain, height1, heigth2
    return a1*np.exp(-1/2*((x-m)/s)**2)+a2*np.exp(-1/2*((x-m-g)/s)**2)

# ### Distribución aproximada de Landau ###
# def Landau(x,a, MP,xi):
#     C1 = a/np.sqrt((2 * np.pi))
#     C2 = np.exp(-((x-MP)-1)/xi)
#     C3 = np.exp((-0.5 * (((x-MP)-1)/xi + C2 )))

#     # C1 = 1/np.sqrt((2 * np.pi))
#     # C3 = np.exp((-0.5 * (x + np.exp(-x))))
#     return  C1 * C3


##### ------------------- Catálogo de Muones ----------------- #####
def oScan_fit_NSAMP1(extensión, active_area, oScan, Bins, make_figure_flag = False) -> dict:
    Maxfev = 10000000
    P0=[10, 2000, 900]

    if make_figure_flag:
        fig_all, axs_all = plt.subplots(1, 1, figsize=(10, 10))
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins, )
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
    Range_x_label = [-200, 400]

    if make_figure_flag:
        fig_all, axs_all = plt.subplots(1, 1, figsize=(10, 10))
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins,)
        offset = bins_edges[np.argmax(hist)]
        print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset
        median_oScan = np.median(Overscan_plane.flatten())
        max_oScan = np.max(Overscan_plane.flatten())
        min_oScan = np.min(Overscan_plane.flatten())

        diff = (max_oScan - abs(median_oScan)) / 2
        # print(median_oScan, diff, min_oScan)

        # Range = (min_oScan, max_oScan)
        if 40 * abs(median_oScan) <  diff:
            Range = (min_oScan, diff)
            print(Range)

        else:
            Range = (min_oScan, diff * 2)
            print(Range)

        bin_heights, bin_borders, _ = axs_all.hist(Overscan_plane.flatten(), bins = Bins,range = Range_x_label, label="Pixeles del Overscan")

        bin_centers = np.zeros(len(bin_heights), dtype=float)
        offset_fit = bin_borders[np.argmax(bin_heights)]

        for p in range(len(bin_heights)):
            bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

        # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
        xmin_fit, xmax_fit = bin_centers[0], bin_centers[-1]
        # print(xmin_fit, xmax_fit)

        bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
        bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
 
        popt, pcov = curve_fit(Gaussian2, bin_centers, bin_heights, maxfev=Maxfev, p0 = [3, 200, 200, 900, 100], bounds = ((-10, 10), (150, 200), (150, 210)))		# Fit histogram with gaussian
        axs_all.plot(bin_centers, Gaussian2(bin_centers, *popt), 'k', label = 'Ajuste Gaussiano')	

        dict_popt = {'Mean' : popt[0], 'sigma' : abs(popt[1]), 'Gain' : popt[2], 'Offset' : offset}

        print('Centroide: ',dict_popt['Mean'], ' Sigma: ', dict_popt['sigma'], 'Offset: ', dict_popt['Offset'], 'Gain: ', dict_popt['Gain']) 

        axs_all.set_title("Distribución de pixeles del Overscan")
        # axs_all.set_xlim(xmin_fit, xmax_fit)
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


        if 30 * abs(median_oScan) <  diff:
            Range = (min_oScan, diff)
            # print(Range)

        else:
            Range = (min_oScan, diff * 2)
            print(Range)


        bin_heights, bin_borders = np.histogram(Overscan_plane.flatten(), bins= Bins, range = Range) #'auto'
        bin_centers = np.zeros(len(bin_heights), dtype=float)
        offset_fit = bin_borders[np.argmax(bin_heights)]

        for p in range(len(bin_heights)):
            bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

        # xmin_fit, xmax_fit = offset_fit-(10*expgain[extension-1])/math.sqrt(nsamp), offset_fit+(10*expgain[extension-1])/math.sqrt(nsamp)			# Define fit range
        xmin_fit, xmax_fit = bin_centers[0], bin_centers[-1]
        bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
        bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

        popt, pcov = curve_fit(Gaussian2, bin_centers, bin_heights, maxfev = Maxfev, p0 = [0,100,200, 900, 100])		# Fit histogram with gaussiano')
        dict_popt = {'Mean' : popt[0], 'sigma' : abs(popt[1]), 'Gain' : popt[2], 'Offset' : offset}
    
    return dict_popt


def oScan_fit_NSAMP324_ROOT(extensión, active_area, oScan, Bins, make_figure_flag = False) -> dict:
    Maxfev = 100000
    Range_x_label = [-200, 400]

    min_oScan = np.min(oScan)

    if make_figure_flag == True:
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins,  range=(min_oScan, 18000))
        offset = bins_edges[np.argmax(hist)]
        # print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset 

        fgaus2 = TF1("fgauss2","[3]*exp(-0.5*((x-[0])/[1])^2)+[4]*exp(-0.5*((x-[0]-[2])/[1])^2)",-300,600,5) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)

        h3=TH1F("histogram", "Distribution of OsCan",1500,-300,500)
        for pixel_value in Overscan_plane.flatten():
            # if not np.ma.is_masked(pixel_value):
            h3.Fill(pixel_value)
            #print(pixel_value)

        fgaus2.SetParameters(0,10,100, 100, 100) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
        h3.Fit(fgaus2)

        c1=TCanvas()
        h3.Draw()
        c1.Draw()


        fgaus2.Draw("same")
        gStyle.SetOptFit(1100)
        gStyle.SetPadGridX (True)
        # fgaus2.Draw('Quiet')


        print('Parameters of the Doble-Gaussian Fit')
        print('Mean: ', fgaus2.GetParameters()[0],  ' +- ', fgaus2.GetParError(0))
        print('Sigma: ', fgaus2.GetParameters()[1],  ' +- ', fgaus2.GetParError(1))
        print('Gain: ', fgaus2.GetParameters()[2],  ' +- ', fgaus2.GetParError(2), '\n')

        print("chiSquare: " + str(fgaus2.GetChisquare()))
        print("NDegrees of Freedom: " + str(fgaus2.GetNDF()))
        print("chiSquare / NDF :", fgaus2.GetChisquare() / fgaus2.GetNDF(), '\n')

        dict_popt = {'Mean' :fgaus2.GetParameters()[0], 'sigma' : abs(fgaus2.GetParameters()[1]), 'Gain' : abs(fgaus2.GetParameters()[2]), 'Offset' : offset}
        
    if make_figure_flag == False:
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins,  range=(min_oScan, 18000))
        offset = bins_edges[np.argmax(hist)]
        # print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset

        fgaus2 = TF1("fgauss2","[3]*exp(-0.5*((x-[0])/[1])^2)+[4]*exp(-0.5*((x-[0]-[2])/[1])^2)",-300,600,5) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)
        
        h3=TH1F("histogram", "Distribution of OsCan",1500,-300,500)
        for pixel_value in Overscan_plane.flatten():
            # if not np.ma.is_masked(pixel_value):
            h3.Fill(pixel_value)
            #print(pixel_value)
        fgaus2.SetParameters(0,10,100, 100, 100) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
        h3.Fit(fgaus2)

        # c3=TCanvas()

        # h3.Draw()
        # c3.Draw()
        # fgaus2.Draw("same")
        # gStyle.SetOptFit(1100)
        # gStyle.SetPadGridX (True)
        # # fgaus2.Draw('Quiet')

        dict_popt = {'Mean' :fgaus2.GetParameters()[0], 'sigma' : abs(fgaus2.GetParameters()[1]), 'Gain' : abs(fgaus2.GetParameters()[2]), 'Offset' : offset}
    
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

def data_calibrated_NSAMP(active_area, extension, offset, gain, ratio_keV, unidades, sigma_ADUs):
    ## Se aplica el offset ##
    dataP = active_area - offset ## En ADUs

    if unidades == 0:
        data = dataP ## En ADUs
        sigma = sigma_ADUs

    elif unidades == 1:
        data = dataP / gain ## En electrones
        sigma = abs(sigma_ADUs / gain)

    elif unidades == 2:
        data = ratio_keV * (dataP / gain) ## En keV
        sigma = abs( ratio_keV *  (sigma_ADUs/ gain))

    return data, sigma


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

    list_charge_all_events = []

    for event in np.arange(1, nlabels_img):
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

        ## Aquí se obtiene los radios de la elipse del cluster
        rM = prop[event-1].axis_major_length / 2
        rm = prop[event-1].axis_minor_length / 2


        ## Aquí se calcula la elipcidad del cluster ##
        try:
            elip = (rM - rm) / rM
        except:
            elip = 0


        ## Aquí se obtiene el solidity del cluster
        Solidity = prop[event-1].solidity


        ## Aquí se obtiene las medidas de la matriz que contiene al cluster
        miny, minx, maxy, maxx = prop[event-1].bbox
        Longitud_y, Longitud_x = maxy - miny , maxx - minx # px


        ## Aquí se calcula la diagonal en el plano XY
        Diagonal_lenght= np.sqrt(Longitud_x**2 + Longitud_y**2) - np.sqrt(2) # px


        ## Aquí se calcula el Delta L del muon
        Delta_L = np.sqrt( (Diagonal_lenght * px_to_micras)**2 + (CCD_depth)**2) * micra_to_cm # cm


        ## Aquí se calcula la carga total del cluster
        charge = data_maskEvent.sum()

        ## Aquí comienza el filtro de muones 
        if rM == 0 or rm == 0:
            list_charge_all_events.append(charge)
            continue 

        elif maxx - minx <= 3:
            list_charge_all_events.append(charge)
            continue

        elif maxy - miny <= 3:
            list_charge_all_events.append(charge)
            continue

        elif not Barycentercharge:
            list_charge_all_events.append(charge)
            continue

        elif differval < MeanValue_Event: 
            continue

        elif  Solidity < Solidit:
            list_charge_all_events.append(charge)
            continue 

        # elif rM <= Elipticity * rm:
        #     list_charge_all_events.append(charge)
        #     continue

        elif elip < Elipticity:
            list_charge_all_events.append(charge)
            continue

        elif  elip >= Elipticity :
            # charge = data_maskEvent.sum()

            # if charge > 100:
            Delta_EL = (charge)/ (Delta_L) 
            theta = np.arctan((Diagonal_lenght * px_to_cm)/(CCD_depth * micra_to_cm)) *(180 /np.pi)

            list_DeltaL.append(Delta_L)
            list_DeltaEL.append(Delta_EL)
            list_charge.append(charge)
            list_theta.append(theta)
            # print(charge, DeltaEL)

            # if DeltaEL_range_min <= DeltaEL <= DeltaEL_range_max:
            list_Muon_labels.append(event)

    return list_DeltaL, list_DeltaEL, list_charge, list_Muon_labels, list_theta, list_charge_all_events

##### ------------------- Catálogo de Muones Rectos (verticales/horizontales) ---------------------- ###

def muon_straight_filter(dataCal, label_img, n_events, Solidit, Elipticity, Prop, min_Charge, Sigma, skirts):
    list_sigmas_vertical_event = []
    list_vertical_event = []
    list_charge_vertical_event = []

    list_sigmas_horizontal_event = []
    list_horizontal_event = []
    list_charge_horizontal_event = []

    list_vertical_events = []
    list_horizontal_events = []

    # num_muons = 0

    for event in np.arange(1,n_events):

        if skirts > 0:
            mask = np.invert(nd.binary_dilation(label_img == event, iterations= skirts))
            loc = nd.find_objects(label_img == event)[0]
        
        elif skirts == 0:
            mask = np.invert(label_img == event)
            loc = nd.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
        
        MinValue_Event = data_maskEvent.min()
        MeanValue_Event = data_maskEvent.mean()

        try: 
            coordX_centerCharge = round(nd.center_of_mass(data_maskEvent)[1])
            coordY_centerCharge = round(nd.center_of_mass(data_maskEvent)[0])
            Barycentercharge = data_maskEvent[coordY_centerCharge, coordX_centerCharge]

            differval = abs(Barycentercharge - MinValue_Event) 

        except:
            Barycentercharge = np.nan()
            differval = 0

        rM = Prop[event-1].axis_major_length/2
        rm = Prop[event-1].axis_minor_length/2
        
        try:
            elip = (rM - rm)/rM 
        except: 
            elip = 0


        Solidity = Prop[event-1].solidity
        miny, minx, maxy, maxx = Prop[event-1].bbox
        Longitud_y, Longitud_x = maxy - miny , maxx - minx # px

        if rM == 0 or rm == 0:
            continue 

        elif maxx - minx <= 3:
            continue

        elif not Barycentercharge:
            continue

        elif differval < MeanValue_Event: #keV
            continue

        elif  Solidity < Solidit:
            continue 

        elif  elip >= Elipticity :
            charge = data_maskEvent.sum()

            if charge < min_Charge:
                continue
            
            if (Longitud_x < 7 and Longitud_y > 10): 
                # num_muons = num_muons + 1

                list_sigmas_vertical_event.append(Sigma)
                list_vertical_event.append(data_maskEvent)
                list_charge_vertical_event.append(charge)

            if ( Longitud_y < 7 and Longitud_x > 10):
                # num_muons = num_muons + 1

                list_sigmas_horizontal_event.append(Sigma)
                list_horizontal_event.append(data_maskEvent)
                list_charge_horizontal_event.append(charge)



        del data_maskEvent
        # del Barycentercharge


    list_vertical_events.append(list_sigmas_vertical_event)
    list_vertical_events.append(list_vertical_event)
    list_vertical_events.append(list_charge_vertical_event)

    list_horizontal_events.append(list_sigmas_horizontal_event)
    list_horizontal_events.append(list_horizontal_event)
    list_horizontal_events.append(list_charge_horizontal_event)


    return list_vertical_events, list_horizontal_events


################### ---------------------- Funciones para el Modelo de Difusión ---------------------------- ###############

def check_flip_vertical_muon(dict, label_muon, Delta_in, Delta_fin, extension):

    Delta_inicial = Delta_in    # px
    Delta_final = Delta_fin     # px

    event =  dict['extension_' + str(extension)]['Vertical_Events'][label_muon]

    label_verticalMuon, nlabels_verticalMuon = ndimage.label(event,structure=[[0,0,0],[1,1,1],[0,0,0]])

    ### Parte de abajo de la imagen ##
    line = label_verticalMuon == Delta_inicial
    # print(Delta_inicial)
    loc = ndimage.find_objects(label_verticalMuon == Delta_inicial)[0]
    mask_35 = np.invert(label_verticalMuon == Delta_inicial)
    data_mask = ma.masked_array(event[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop], mask_35[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

    Longitud_linea = len(data_mask[0])
    Carga_renglon = data_mask[0].sum()
    Mean_carga = np.mean(data_mask[0])

    Mean_in_1 = 0
    var_1 = 0
    carga_cuadrada = 0
    NaN_pixel_1 = 0

    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask[0][pixel]:
            element_pixel = (pixel * data_mask[0][pixel]) / Carga_renglon
            Mean_in_1 = Mean_in_1 + element_pixel
        else:
            element_pixel = 0
        
        Mean_in_1 = Mean_in_1 + element_pixel

    ## Calcula la suma de las cargas al cuadrado ##
    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask[0][pixel]:
            element_pixel = data_mask[0][pixel]**2
        else:
            NaN_pixel_1 = NaN_pixel_1 + 1
            element_pixel = 0
        carga_cuadrada = carga_cuadrada + element_pixel 

    Mean_carga_cuadrada_1 = carga_cuadrada/Longitud_linea


    for pixel in np.arange(0, Longitud_linea, 1):
        element_pixel =(1 / (Longitud_linea - 1)) * (pixel - Mean_in_1)**2
        var_1 = var_1 + element_pixel 

    var_1_true = var_1 * (Mean_carga_cuadrada_1 / (Mean_carga**2))
    # var_1_true = var_1

    sigma_in = np.sqrt(var_1_true)


    ### Parte de arriba de la imagen ###
    line = label_verticalMuon ==  nlabels_verticalMuon - Delta_final
    # print( nlabels_verticalMuon - Delta_final)
    loc = ndimage.find_objects(label_verticalMuon == nlabels_verticalMuon - Delta_final)[0]
    mask_35 = np.invert(label_verticalMuon == nlabels_verticalMuon - Delta_final)
    data_mask = ma.masked_array(event[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop], mask_35[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
    # print(data_mask[0])

    Longitud_linea = len(data_mask[0])
    Carga_renglon = data_mask[0].sum()
    Mean_carga = np.mean(data_mask[0])

    Mean_in_2 = 0
    var_2 = 0
    carga_cuadrada = 0
    NaN_pixel_2 = 0

    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask[0][pixel]:
            element_pixel = (pixel * data_mask[0][pixel]) / Carga_renglon
        else:
            NaN_pixel_2 = NaN_pixel_2 + 1
            element_pixel = 0
        Mean_in_2 = Mean_in_2 + element_pixel
        # print('Valor mean: ', Mean_in_2)

    ## Calcula la suma de las cargas al cuadrado ##
    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask[0][pixel]:
            element_pixel = data_mask[0][pixel]**2
        else:
            element_pixel = 0
        carga_cuadrada = carga_cuadrada + element_pixel 

    Mean_carga_cuadrada_2 = carga_cuadrada / Longitud_linea

    for pixel in np.arange(0, Longitud_linea, 1):
        element_pixel = (1 / (Longitud_linea - 1)) * (pixel - Mean_in_2)**2
        var_2 = var_2 + element_pixel

    var_2_true = var_2 * (Mean_carga_cuadrada_2 / (Mean_carga**2))
    # var_2_true = var_2

    sigma_fn = np.sqrt(var_2_true)

    # if sigma_in > sigma_fn:
    #     turn_event = np.flip(event, 0)
    #     event = turn_event

    # if sigma_in < sigma_fn:
    #     n = 0


    if NaN_pixel_1 >= NaN_pixel_2:
        Event = event
        flag_turn = False

    elif NaN_pixel_1 < NaN_pixel_2:
        turn_event = np.flip(event, 0)
        Event = turn_event
        flag_turn = True
    
    return Event, flag_turn

def check_flip_horizontal_muon(dict, label_muon, Delta_in, Delta_fin, extension):

    Delta_inicial = Delta_in    # px
    Delta_final = Delta_fin     # px

    event =  dict['extension_' + str(extension)]['Horizontal_Events'][label_muon]

    label_verticalMuon, nlabels_verticalMuon = ndimage.label(event,structure=[[0,1,0],[0,1,0],[0,1,0]])

    ### Parte de abajo de la imagen ##
    line = label_verticalMuon == Delta_inicial
    # print(Delta_inicial)
    loc = ndimage.find_objects(label_verticalMuon == Delta_inicial)[0]
    mask_35 = np.invert(label_verticalMuon == Delta_inicial)
    data_mask = ma.masked_array(event[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop], mask_35[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

    Longitud_linea = len(data_mask.flatten())
    Carga_renglon = data_mask.flatten().sum()
    Mean_carga = np.mean(data_mask.flatten())

    Mean_in_1 = 0
    var_1 = 0
    carga_cuadrada = 0
    NaN_pixel_1 = 0

    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask.flatten()[pixel]:
            element_pixel = (pixel *data_mask.flatten()[pixel]) / Carga_renglon
            Mean_in_1 = Mean_in_1 + element_pixel
        else:
            element_pixel = 0
        
        Mean_in_1 = Mean_in_1 + element_pixel

    ## Calcula la suma de las cargas al cuadrado ##
    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask.flatten()[pixel]:
            element_pixel = data_mask.flatten()[pixel]**2
        else:
            NaN_pixel_1 = NaN_pixel_1 + 1
            element_pixel = 0
        carga_cuadrada = carga_cuadrada + element_pixel 

    Mean_carga_cuadrada_1 = carga_cuadrada/Longitud_linea


    for pixel in np.arange(0, Longitud_linea, 1):
        element_pixel =(1 / (Longitud_linea - 1)) * (pixel - Mean_in_1)**2
        var_1 = var_1 + element_pixel 

    var_1_true = var_1 * (Mean_carga_cuadrada_1 / (Mean_carga**2))
    # var_1_true = var_1

    sigma_in = np.sqrt(var_1_true)


    ### Parte de arriba de la imagen ###
    line = label_verticalMuon ==  nlabels_verticalMuon - Delta_final
    # print( nlabels_verticalMuon - Delta_final)
    loc = ndimage.find_objects(label_verticalMuon == nlabels_verticalMuon - Delta_final)[0]
    mask_35 = np.invert(label_verticalMuon == nlabels_verticalMuon - Delta_final)
    data_mask = ma.masked_array(event[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop], mask_35[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
    # print(data_mask[0])

    Longitud_linea = len(data_mask.flatten())
    Carga_renglon = data_mask.flatten().sum()
    Mean_carga = np.mean(data_mask.flatten())

    Mean_in_2 = 0
    var_2 = 0
    carga_cuadrada = 0
    NaN_pixel_2 = 0

    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask.flatten()[pixel]:
            element_pixel = (pixel * data_mask.flatten()[pixel]) / Carga_renglon
        else:
            NaN_pixel_2 = NaN_pixel_2 + 1
            element_pixel = 0
        Mean_in_2 = Mean_in_2 + element_pixel
        # print('Valor mean: ', Mean_in_2)

    ## Calcula la suma de las cargas al cuadrado ##
    for pixel in np.arange(0, Longitud_linea, 1):
        if data_mask.flatten()[pixel]:
            element_pixel = data_mask.flatten()[pixel]**2
        else:
            element_pixel = 0
        carga_cuadrada = carga_cuadrada + element_pixel 

    Mean_carga_cuadrada_2 = carga_cuadrada / Longitud_linea

    for pixel in np.arange(0, Longitud_linea, 1):
        element_pixel = (1 / (Longitud_linea - 1)) * (pixel - Mean_in_2)**2
        var_2 = var_2 + element_pixel

    var_2_true = var_2 * (Mean_carga_cuadrada_2 / (Mean_carga**2))
    # var_2_true = var_2

    sigma_fn = np.sqrt(var_2_true)

    # if sigma_in > sigma_fn:
    #     turn_event = np.flip(event, 0)
    #     event = turn_event

    # if sigma_in < sigma_fn:
    #     n = 0


    if NaN_pixel_1 >= NaN_pixel_2:
        Event = event
        flag_turn = False

    elif NaN_pixel_1 < NaN_pixel_2:
        turn_event = np.flip(event, 1)
        Event = turn_event
        flag_turn = True
    
    return Event, flag_turn

def diffution_vertical_muon(dict, list_vertical_labels, Delta_in, Delta_fin, extension):

    list_all_sigmas = []
    list_deep = []

    CCD_depth = 725 # micras
    Delta_inicial = Delta_in    # px
    Delta_final = Delta_fin     # px

    for label_muon in list_vertical_labels:
        list_sigmas = []
        # event = data_histogram['extension_' + str(extension)]['Vertical_Events'][label_muon]
        event, _ = check_flip_vertical_muon(dict  = dict, label_muon = label_muon, Delta_in= Delta_inicial, Delta_fin=Delta_final, extension=extension)
        

        label_verticalMuon, nlabels_verticalMuon = ndimage.label(event,structure=[[0,0,0],[1,1,1],[0,0,0]])

        size_x = event.shape[1]
        size_y = event.shape[0]

        lines = 0
        Longitud_XY = size_y 

        Z_inicial = (Delta_inicial * CCD_depth) / (Longitud_XY - Delta_final)

        for lable_line in np.arange(Delta_inicial, nlabels_verticalMuon - Delta_final):
            ## Enmascara la linea en turno
            line = label_verticalMuon == lable_line
            loc = ndimage.find_objects(label_verticalMuon == lable_line)[0]
            mask_35 = np.invert(label_verticalMuon == lable_line)
            data_mask = ma.masked_array(event[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop], mask_35[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

            Longitud_linea = len(data_mask[0])
            Carga_renglon = data_mask[0].sum()
            Mean_carga = np.mean(data_mask[0])

            Mean_in  = 0
            var = 0
            carga_cuadrada = 0

            ### Se calcula el X promedio ##
            for pixel in np.arange(0, Longitud_linea, 1):
                if data_mask[0][pixel]:
                    element_pixel = (pixel * data_mask[0][pixel]) / Carga_renglon
                    Mean_in = Mean_in + element_pixel
                else:
                    element_pixel = 0
                
                Mean_in = Mean_in + element_pixel

            ## Calcula la suma de las cargas al cuadrado ##
            for pixel in np.arange(0, Longitud_linea, 1):
                if data_mask[0][pixel]:
                    element_pixel = data_mask[0][pixel]**2
                else:
                    element_pixel = 0

                carga_cuadrada = carga_cuadrada + element_pixel 

            Mean_carga_cuadrada = carga_cuadrada/Longitud_linea

            ### Se calcula la varianza ##
            for pixel in np.arange(0, Longitud_linea, 1):
                element_pixel =(1 / (Longitud_linea - 1)) * (pixel - Mean_in)**2
                var = var + element_pixel  ### COreggir la varianza con otro estimados

            ### Se corrige la varianza con la carga ###
            var_true = var * (Mean_carga_cuadrada / (Mean_carga**2))

            ### Se calcula la sigma ###
            sigma_in = np.sqrt(var_true)
            

            ##Se crea un arreglo para usarlo en el plot de los datos, y se realiza el juste ##
            list_xlabel = np.arange(0.5, len(data_mask[0]), 1)

            list_xlabel_long = np.linspace(-Longitud_linea + int(Longitud_linea/2) , Longitud_linea + int(Longitud_linea/2), Longitud_linea)
            
            popt, pcov = curve_fit(gaussian, list_xlabel, data_mask[0], maxfev=100000, p0 = [1000, Mean_in, sigma_in])		# Fit histogram with gaussian

            ## Se guardan lo parámetros del ajuste en un diccionario ##
            dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]), 'Pcov' : pcov}
            Centroide = popt[1]
            Sigma = abs(popt[2])

            # print(Sigma)

            # if Sigma > 6:
            #     continue

            # ## Se grafican los puntos experimentales ##
            # axs_all.scatter(list_xlabel_long, data_mask[0], lable_line, '.', color = 'k')

            # ## Se crea otro arreglo para el plot del ajusto y se dibuja ##
            # list_xlabel_long = np.linspace( Centroide - 4 , Centroide + 4)
            
            # axs_all.plot(list_xlabel_long, gaussian(list_xlabel_long, *popt), lable_line, 'k')	
            # axs_all.legend()

            # Se guarda la sigma de la distribución en una lista ##
            list_all_sigmas.append(Sigma)
            list_sigmas.append(Sigma)
            # list_all_sigmas.append(Sigma)
            # print('Centroide: ',popt[1], ' Amplitud: ', popt[0], 'sigma: ', abs(popt[2]))  #gaussian(x, a, mean, sigma)
            lines = lines + 1

        list_xlabel_sigmas = np.linspace(Z_inicial, CCD_depth, len(list_sigmas))

        for deep in list_xlabel_sigmas:
            list_deep.append(deep)

    
    return list_all_sigmas, list_deep

def diffution_horizontal_muon(dict, list_horizontal_labels, Delta_in, Delta_fin, extension):

    list_all_sigmas = []
    list_deep = []

    CCD_depth = 725 # micras
    Delta_inicial = Delta_in    # px
    Delta_final = Delta_fin     # px

    for label_muon in list_horizontal_labels:
        list_sigmas = []
        # event = data_histogram['extension_' + str(extension)]['Vertical_Events'][label_muon]
        event, _ = check_flip_horizontal_muon(dict  = dict, label_muon = label_muon, Delta_in= Delta_inicial, Delta_fin=Delta_final, extension=extension)
        

        label_verticalMuon, nlabels_verticalMuon = ndimage.label(event,structure=[[0,1,0],[0,1,0],[0,1,0]])

        size_x = event.shape[1]
        size_y = event.shape[0]

        lines = 0
        Longitud_XY = size_x 

        Z_inicial = (Delta_inicial * CCD_depth) / (Longitud_XY - Delta_final)

        for lable_line in np.arange(Delta_inicial, nlabels_verticalMuon - Delta_final):
            ## Enmascara la linea en turno
            line = label_verticalMuon == lable_line
            loc = ndimage.find_objects(label_verticalMuon == lable_line)[0]
            mask_35 = np.invert(label_verticalMuon == lable_line)
            data_mask = ma.masked_array(event[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop], mask_35[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

            Longitud_linea = len(data_mask.flatten())
            Carga_renglon = data_mask.flatten().sum()
            Mean_carga = np.mean(data_mask.flatten())

            Mean_in  = 0
            var = 0
            carga_cuadrada = 0

            ### Se calcula el X promedio ##
            for pixel in np.arange(0, Longitud_linea, 1):
                if data_mask.flatten()[pixel]:
                    element_pixel = (pixel * data_mask.flatten()[pixel]) / Carga_renglon
                    Mean_in = Mean_in + element_pixel
                else:
                    element_pixel = 0
                
                Mean_in = Mean_in + element_pixel

            ## Calcula la suma de las cargas al cuadrado ##
            for pixel in np.arange(0, Longitud_linea, 1):
                if data_mask.flatten()[pixel]:
                    element_pixel = data_mask.flatten()[pixel]**2
                else:
                    element_pixel = 0

                carga_cuadrada = carga_cuadrada + element_pixel 

            Mean_carga_cuadrada = carga_cuadrada/Longitud_linea

            ### Se calcula la varianza ##
            for pixel in np.arange(0, Longitud_linea, 1):
                element_pixel =(1 / (Longitud_linea - 1)) * (pixel - Mean_in)**2
                var = var + element_pixel  ### COreggir la varianza con otro estimados

            ### Se corrige la varianza con la carga ###
            var_true = var * (Mean_carga_cuadrada / (Mean_carga**2))

            ### Se calcula la sigma ###
            sigma_in = np.sqrt(var_true)
            

            ##Se crea un arreglo para usarlo en el plot de los datos, y se realiza el juste ##
            list_xlabel = np.arange(0.5, len(data_mask.flatten()), 1)

            list_xlabel_long = np.linspace(-Longitud_linea + int(Longitud_linea/2) , Longitud_linea + int(Longitud_linea/2), Longitud_linea)
            
            try:
                popt, pcov = curve_fit(gaussian, list_xlabel, data_mask.flatten(), maxfev=100000, p0 = [1000, Mean_in, sigma_in])		# Fit histogram with gaussian

            except: 
                print('Fit error in label muon: ' + str(label_muon))
                continue

            ## Se guardan lo parámetros del ajuste en un diccionario ##
            dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]), 'Pcov' : pcov}
            Centroide = popt[1]
            Sigma = abs(popt[2])

            # print(Sigma)

            # if Sigma > 6:
            #     continue

            # ## Se grafican los puntos experimentales ##
            # axs_all.scatter(list_xlabel_long, data_mask[0], lable_line, '.', color = 'k')

            # ## Se crea otro arreglo para el plot del ajusto y se dibuja ##
            # list_xlabel_long = np.linspace( Centroide - 4 , Centroide + 4)
            
            # axs_all.plot(list_xlabel_long, gaussian(list_xlabel_long, *popt), lable_line, 'k')	
            # axs_all.legend()

            # Se guarda la sigma de la distribución en una lista ##
            list_all_sigmas.append(Sigma)
            list_sigmas.append(Sigma)
            # list_all_sigmas.append(Sigma)
            # print('Centroide: ',popt[1], ' Amplitud: ', popt[0], 'sigma: ', abs(popt[2]))  #gaussian(x, a, mean, sigma)
            lines = lines + 1

        list_xlabel_sigmas = np.linspace(Z_inicial, CCD_depth, len(list_sigmas))

        for deep in list_xlabel_sigmas:
            list_deep.append(deep)

    
    return list_all_sigmas, list_deep


        
