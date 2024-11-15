import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import numpy.ma as ma
import pandas as pd 
import skimage as sk
import scipy.ndimage as nd
import random
import time

# from ROOT import TF1, TH1F, Fit, TCanvas, gStyle


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
        loc = nd.find_objects(label_img == event)[0]
        
        data_maskEvent = np.ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
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
    list_phi = []

    list_charge_all_events = []

    for event in np.arange(1, nlabels_img):
        mask = np.invert(label_img == event)
        loc = nd.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

        coordX_centerCharge = round(nd.center_of_mass(data_maskEvent)[1])
        coordY_centerCharge = round(nd.center_of_mass(data_maskEvent)[0])

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

        elif elip < Elipticity:
            list_charge_all_events.append(charge)
            continue

        elif  elip >= Elipticity :
            # charge = data_maskEvent.sum()

            # if charge > 100:
            Delta_EL = (charge)/ (Delta_L) 
            theta = np.arctan((Diagonal_lenght * px_to_cm)/(CCD_depth * micra_to_cm)) *(180 /np.pi)

            
            #### ------------------------ CÁLCULO DEL ÁNGULO PHI ---------------------------- ### 
            len_y, len_x = data_maskEvent.shape
            flag_ld, flag_rd, flag_ru, flag_lu = False, False, False, False

            try:
                # print('Estoy intentando callcular phi')
                left_down = data_maskEvent[ 0:int(len_y/2), 0:int(len_x/2)]
                right_down = data_maskEvent[0:int(len_y/2), int(len_x/2):len_x]
                right_up = data_maskEvent[int(len_y/2):len_y, int(len_x/2):len_x]
                left_up = data_maskEvent[int(len_y/2):len_y, 0:int(len_x/2)]

                charge_ld = left_down.sum()
                charge_rd = right_down.sum()
                charge_ru = right_up.sum()
                charge_lu = left_up.sum()

                diff_deltas = 20
                if charge_ld > charge_rd and charge_ld > charge_ru and charge_ld > charge_lu: 
                    if charge_ld - charge_ru < diff_deltas:
                        flag_ru = True
                        # print('La cola está arriba derecha (correccion)')
                    else: 
                        flag_ld = True
                        # print('La cola está abajo izquierda')
                    
                elif charge_rd > charge_ld and charge_rd > charge_ru and charge_rd > charge_lu:
                    if charge_rd - charge_lu < diff_deltas:
                        flag_lu = True
                        # print('La cola está arriba izquierda (correccion)')
                    else:
                        flag_rd = True
                        # print('La cola está abajo derecha')

                elif charge_ru > charge_ld and charge_ru >charge_rd and charge_ru > charge_lu:
                    if charge_ru - charge_ld < diff_deltas:
                        flag_ld = True
                        # print('La cola está abajo izquierda (correccion)')
                    else:
                        flag_ru = True
                        # print('La cola está arriba derecha')

                elif charge_lu > charge_ld and charge_lu > charge_rd and charge_lu > charge_ru:
                    if charge_lu - charge_rd < diff_deltas:
                        flag_rd = True
                        # print('La cola está abajo derecha (correccion)') 
                    else:
                        flag_lu = True
                        # print('La cola está arriba izquierda')
                else:
                    # print('No sirvió para calcular phi')
                    continue

                if len_x < 7:
                    len_x = len_x  / 2
                    if flag_ld: # Cuadrante I
                        phi = np.arctan(len_y/len_x) # En radianes

                    elif flag_rd: # Cuadrante II
                        phi_comp = np.arctan(len_x/len_y)
                        # print(phi_comp)
                        phi = phi_comp + np.pi/2

                    elif flag_ru: # Cuadrante III
                        phi_comp = np.arctan(len_y/len_x) 
                        phi = phi_comp + np.pi

                    elif flag_lu: # Cuadrante IV
                        phi_comp = np.arctan(len_x/len_y)
                        phi = phi_comp + 3 * np.pi/2
            
                elif len_y < 7:
                    len_y = len_y  / 2
                    if flag_ld: # Cuadrante I
                        phi = np.arctan(len_y/len_x) # En radianes

                    elif flag_rd: # Cuadrante II
                        phi_comp = np.arctan(len_x/len_y)
                        # print(phi_comp)
                        phi = phi_comp + np.pi/2

                    elif flag_ru: # Cuadrante III
                        phi_comp = np.arctan(len_y/len_x) 
                        phi = phi_comp + np.pi

                    elif flag_lu: # Cuadrante IV
                        phi_comp = np.arctan(len_x/len_y)
                        phi = phi_comp + 3 * np.pi/2

                else:
                    if flag_ld: # Cuadrante I
                        phi = np.arctan(len_y/len_x) # En radianes

                    elif flag_rd: # Cuadrante II
                        phi_comp = np.arctan(len_x/len_y)
                        # print(phi_comp)
                        phi = phi_comp + np.pi/2

                    elif flag_ru: # Cuadrante III
                        phi_comp = np.arctan(len_y/len_x) 
                        phi = phi_comp + np.pi

                    elif flag_lu: # Cuadrante IV
                        phi_comp = np.arctan(len_x/len_y)
                        phi = phi_comp + 3 * np.pi/2
            except:
                # print('No pude callcular phi')
                phi = -1
            # print(phi)



            list_DeltaL.append(Delta_L)
            list_DeltaEL.append(Delta_EL)
            list_charge.append(charge)
            list_theta.append(theta)
            list_phi.append(phi)
            # print(charge, DeltaEL)

            # if DeltaEL_range_min <= DeltaEL <= DeltaEL_range_max:
            list_Muon_labels.append(event)

    return list_DeltaL, list_DeltaEL, list_charge, list_Muon_labels, list_theta, list_phi, list_charge_all_events


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
            Barycentercharge = np.nan
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

