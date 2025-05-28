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
from ROOT import TMath, TF1, TH1F, TH2F, TCanvas, gStyle, TProfile, TGraphErrors


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

def phi_angle_ROOT(data_mask):
    ## ======== Dimensiones del evento a analizar ============ ###
    NBX = data_mask.shape[1]
    xmin = 0
    xmax = NBX
    # print(NBX, len(data_mask[:, 1]))

    NBY = data_mask.shape[0]
    ymin = 0
    ymax = NBY

    ### ========== Se crean y llenan los TProfile de X y Y ===== ###
    YProf = TProfile("YProf", "", NBX, 0, xmax)
    XProf = TProfile("XProf", "", NBY, 0, ymax)

    for i in np.arange(0, NBX):
    # print(i)
        for j in np.arange(0, NBY):
            # print(j)
            if data_mask[j][i]:
                XProf.Fill(j, i, data_mask[j][i])
                YProf.Fill(i, j, data_mask[j][i])

    ### ======= Se obtienen las medias y las sigmas de cada TProfile ====== ###
    list_yprofile_mean = []
    list_yprofile_sigma = []
    for index in np.arange(0, NBX):
        mean_y = YProf.GetBinContent(int(index+1))
        sigma_y = YProf.GetBinError(int(index+1))

        list_yprofile_mean.append(mean_y)
        list_yprofile_sigma.append(sigma_y)
    # print(list_yprofile_sigma)

    list_xprofile_mean = []
    list_xprofile_sigma = []
    for index in np.arange(0, NBY):
        mean_x = XProf.GetBinContent(int(index+1))
        sigma_x = XProf.GetBinError(int(index+1))

        list_xprofile_mean.append(mean_x)
        list_xprofile_sigma.append(sigma_x)

    ### ====== Se llenan los TGrapsErrors para X y Y ====== ###
    GRprofX = TGraphErrors()
    GRprofY = TGraphErrors()

    for i in np.arange(0, NBX):
        GRprofY.SetPoint(int(i), i, list_yprofile_mean[i])
        GRprofY.SetPointError(int(i), 1/2, list_yprofile_sigma[i])  ### Se colocan los errores en X y Y
    # print('Se llenó el GRprofY')
        
    for ip in np.arange(0, NBY):
        GRprofX.SetPoint(int(ip),  list_xprofile_mean[ip], ip)
        GRprofX.SetPointError(int(ip), list_xprofile_sigma[ip], 1/2)
    # print('Se llenó el GRprofY')
    
    if NBX > 3 * NBY: ### El muon es muy horizontal y se tomará la componente X para hacer los ajustes
        lox = 0
        hix = NBX

        pend = (NBY)/(NBX)
        fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
        fitline.SetParameters(0, pend)
        GRprofY.Fit("fitline", "W")

        pendiente = fitline.GetParameters()[1]
        Prob_fitline = fitline.GetProb()

        if Prob_fitline > 0.01:
            First = list_yprofile_sigma[3]
            Last = list_yprofile_sigma[-4]

            if pendiente > 0:
                if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                    phi = np.arctan((NBY/2)/NBX)
                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                    phi = np.arctan((NBY/2)/NBX) + np.radians(180)
            else:
                if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                    phi = np.arctan(NBX/(NBY/2)) + np.radians(90)
                else: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 2
                    phi = np.arctan(NBX/(NBY/2)) + np.radians(270)

                
        else:
            fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
            fitline.SetParameters(NBY-1, -pend)
            GRprofY.Fit("fitline", "W")
            Prob_fitline = fitline.GetProb()
            pendiente = fitline.GetParameters()[1]

            if Prob_fitline > 0.01:
                First = list_yprofile_sigma[3]
                Last = list_yprofile_sigma[-4]

                if pendiente > 0:
                    if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                        phi = np.arctan((NBY/2)/NBX)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                        phi = np.arctan((NBY/2)/NBX) + np.radians(180)
                else:
                    if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/(NBY/2)) + np.radians(90)
                    else: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 2
                        phi = np.arctan(NBX/(NBY/2)) + np.radians(270)
            else:
                phi = -1

    elif NBY > 3 * NBX: ### El muon es muy vertical y se tomará la componente Y para hacer los ajustes
        lox = 0
        hix = NBX

        pend = (NBY)/(NBX)

        fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
        fitline.SetParameters(0, pend)
        # GRprofX.Fit("fitline", "W")
        GRprofX.Fit("fitline")
        Prob_fitline = fitline.GetProb()
        pendientefit = fitline.GetParameters()[1]

        if Prob_fitline > 0.01:
            First = list_xprofile_sigma[4]
            Last = list_xprofile_sigma[-3]
            # print(First, Last)
            if pendientefit > 0:
                if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                    phi = np.arctan(NBY/(NBX/2))
                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                    phi = np.arctan(NBY/(NBX/2)) + np.radians(180)
            else: 
                if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                    phi = np.arctan((NBX/2)/NBY) + np.radians(270)
                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                    phi = np.arctan((NBX/2)/NBY) + np.radians(90)
        
        else:
            fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
            fitline.SetParameters(NBY-1, -pend)
            # GRprofX.Fit("fitline", "W")
            GRprofX.Fit("fitline")
            Prob_fitline = fitline.GetProb()
            pendientefit = fitline.GetParameters()[1]

            if Prob_fitline > 0.01:
                First = list_xprofile_sigma[4]
                Last = list_xprofile_sigma[-3]
                # print(First, Last)
                if pendientefit > 0:
                    if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                        phi = np.arctan(NBY/(NBX/2))
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                        phi = np.arctan(NBY/(NBX/2)) + np.radians(180)
                else:
                    if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan((NBX/2)/NBY) + np.radians(270)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan((NBX/2)/NBY) + np.radians(90)
            else:
                fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                fitline.SetParameters(0, pend)
                GRprofX.Fit("fitline")
                Prob_fitline = fitline.GetProb()
                pendientefit = fitline.GetParameters()[1]

                if Prob_fitline > 0.01:
                    First = list_xprofile_sigma[4]
                    Last = list_xprofile_sigma[-3]
                    # print(First, Last)
                    if pendientefit > 0:
                        if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                            phi = np.arctan(NBY/(NBX/2))
                        else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                            phi = np.arctan(NBY/(NBX/2)) + np.radians(180)
                    else: 
                        if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                            phi = np.arctan((NBX/2)/NBY) + np.radians(270)
                        else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                            phi = np.arctan((NBX/2)/NBY) + np.radians(90)
                
                else:
                    fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                    fitline.SetParameters(NBY-1, -pend)
                    GRprofX.Fit("fitline")
                    Prob_fitline = fitline.GetProb()
                    pendientefit = fitline.GetParameters()[1]

                    if Prob_fitline > 0.01:
                        First = list_xprofile_sigma[4]
                        Last = list_xprofile_sigma[-3]
                        # print(First, Last)
                        if pendientefit > 0:
                            if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                                phi = np.arctan(NBY/(NBX/2))
                            else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                                phi = np.arctan(NBY/(NBX/2)) + np.radians(180)
                        else:
                            if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                phi = np.arctan((NBX/2)/NBY) + np.radians(270)
                            else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                phi = np.arctan((NBX/2)/NBY) + np.radians(90)
                    else:
                        phi = -2
        
    else:   ### El muon es oblicuo y se tomarán las componentes X y Y para hacer los ajustes
        #### ================ MEAN XY =============== ###
        list_x_XY = list(np.arange(0, NBX))
        for element in list_xprofile_mean:
            list_x_XY.append(element)
        # print(list_x_XY)

        list_y_XY = list(list_yprofile_mean)
        for element in np.arange(0, NBY):
            list_y_XY.append(element)
        # print(list_y_XY)

        #### =============== SIGMA XY =============== ###
        zerosXY = list(np.zeros(NBX))
        list_xsigma_XY = []
        for index in np.arange(0, len(zerosXY)):
            list_xsigma_XY.append(1/2)
        # print(list_xsigma_XY)

        for element in list_xprofile_sigma:
            list_xsigma_XY.append(element)
        # print(list_xsigma_XY)

        list_ysigma_XY = list(list_yprofile_sigma)
        for element in list(np.zeros(NBY)):
            list_ysigma_XY.append(1/2)
        # print(list_ysigma_XY)

        ### ===== Se llena el TGraphErrors de las componentes X y Y JUNTAS ====== ###
        GRprofXY = TGraphErrors()
        for index in np.arange(0, len(list_y_XY)):
            GRprofXY.SetPoint(int(index), list_x_XY[index], list_y_XY[index])
            GRprofXY.SetPointError(int(index), list_xsigma_XY[index], list_ysigma_XY[index])

        lox = 0
        hix = NBX
        pend = (NBY)/(NBX)

        if NBY < NBX:
            fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
            fitline.SetParameters(0, pend)
            GRprofX.Fit("fitline", "W")
            Prob_fitline = fitline.GetProb()
            pendientefit = fitline.GetParameters()[1]

            if Prob_fitline > 0.01: 
                First = list_yprofile_sigma[3]
                Last = list_yprofile_sigma[-4]
                # print(len(list_xprofile_sigma), len(list_yprofile_sigma), First, Last)
                # print(list_yprofile_mean[3], list_yprofile_mean[-4], First, Last)

                if pendientefit > 0:
                    if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                        phi = np.arctan(NBY/NBX)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                        phi = np.arctan(NBY/NBX) + np.radians(180)
                else: 
                    if First < Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(270)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(90)
            else:
                fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                fitline.SetParameters(NBY-1, -pend)
                GRprofX.Fit("fitline", "W")
                Prob_fitline = fitline.GetProb()
                pendientefit = fitline.GetParameters()[1]

                if Prob_fitline > 0.01: 
                    First = list_yprofile_sigma[3]
                    Last = list_yprofile_sigma[-4]

                    if pendientefit > 0:
                        if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                            phi = np.arctan(NBY/NBX)
                        else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                            phi = np.arctan(NBY/NBX) + np.radians(180)
                    else: 
                        # print('pendiente negativa')
                        if First < Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                            phi = np.arctan(NBX/NBY) + np.radians(270)
                        else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                            phi = np.arctan(NBX/NBY) + np.radians(90)
                else: 
                    # print('No fit- fail in fit -pend (NBY<NBX)')
                    fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                    fitline.SetParameters(0, pend)
                    GRprofX.Fit("fitline")
                    Prob_fitline = fitline.GetProb()
                    pendientefit = fitline.GetParameters()[1]

                    if Prob_fitline > 0.01: 
                        First = list_yprofile_sigma[3]
                        Last = list_yprofile_sigma[-4]
                        # print(len(list_xprofile_sigma), len(list_yprofile_sigma), First, Last)
                        # print(list_yprofile_mean[3], list_yprofile_mean[-4], First, Last)

                        if pendientefit > 0:
                            if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                                phi = np.arctan(NBY/NBX)
                            else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                                phi = np.arctan(NBY/NBX) + np.radians(180)
                        else: 
                            if First < Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                phi = np.arctan(NBX/NBY) + np.radians(270)
                            else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                phi = np.arctan(NBX/NBY) + np.radians(90)
                    else:
                        fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                        fitline.SetParameters(NBY-1, -pend)
                        GRprofX.Fit("fitline")
                        Prob_fitline = fitline.GetProb()
                        pendientefit = fitline.GetParameters()[1]

                        if Prob_fitline > 0.01: 
                            First = list_yprofile_sigma[3]
                            Last = list_yprofile_sigma[-4]

                            if pendientefit > 0:
                                if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                                    phi = np.arctan(NBY/NBX)
                                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                                    phi = np.arctan(NBY/NBX) + np.radians(180)
                            else: 
                                # print('pendiente negativa')
                                if First < Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                    phi = np.arctan(NBX/NBY) + np.radians(270)
                                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                    phi = np.arctan(NBX/NBY) + np.radians(90)
                        else: 
                            # print('No fit- fail in fit -pend (NBY<NBX)')
                            phi = -3
                        
        else:
            fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
            fitline.SetParameters(0, pend)
            GRprofY.Fit("fitline", "W")
            Prob_fitline = fitline.GetProb()
            pendientefit = fitline.GetParameters()[1]

            if Prob_fitline > 0.01: 
                First = list_xprofile_sigma[3]
                Last = list_xprofile_sigma[-4]

                if pendientefit > 0:
                    if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                        phi = np.arctan(NBY/NBX)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                        phi = np.arctan(NBY/NBX) + np.radians(180)
                else: 
                    if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(270)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(90)
            else:
                fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                fitline.SetParameters(NBY-1, -pend)
                GRprofY.Fit("fitline", "W")
                Prob_fitline = fitline.GetProb()
                pendientefit = fitline.GetParameters()[1]

                if Prob_fitline > 0.01: 
                    First = list_xprofile_sigma[3]
                    Last = list_xprofile_sigma[-4]
                    # print(First, Last)
                    if pendientefit > 0:
                        if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                            phi = np.arctan(NBY/NBX)
                        else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                            phi = np.arctan(NBY/NBX) + np.radians(180)
                    else: 
                        if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                            phi = np.arctan(NBX/NBY) + np.radians(270)
                        else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                            phi = np.arctan(NBX/NBY) + np.radians(90)
                else: 
                    fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                    fitline.SetParameters(0, pend)
                    GRprofY.Fit("fitline")
                    Prob_fitline = fitline.GetProb()
                    pendientefit = fitline.GetParameters()[1]

                    if Prob_fitline > 0.01: 
                        First = list_xprofile_sigma[3]
                        Last = list_xprofile_sigma[-4]

                        if pendientefit > 0:
                            if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                                phi = np.arctan(NBY/NBX)
                            else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                                phi = np.arctan(NBY/NBX) + np.radians(180)
                        else: 
                            if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                phi = np.arctan(NBX/NBY) + np.radians(270)
                            else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                phi = np.arctan(NBX/NBY) + np.radians(90)
                    else:
                        fitline = TF1("fitline", "[0] + [1]*x",lox,hix) 
                        fitline.SetParameters(NBY-1, -pend)
                        GRprofY.Fit("fitline")
                        Prob_fitline = fitline.GetProb()
                        pendientefit = fitline.GetParameters()[1]

                        if Prob_fitline > 0.01: 
                            First = list_xprofile_sigma[3]
                            Last = list_xprofile_sigma[-4]
                            # print(First, Last)
                            if pendientefit > 0:
                                if First < Last: ## La "cola" está en la parte de abajo, y el muon está en el cuadrante 1
                                    phi = np.arctan(NBY/NBX)
                                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 3
                                    phi = np.arctan(NBY/NBX) + np.radians(180)
                            else: 
                                if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                    phi = np.arctan(NBX/NBY) + np.radians(270)
                                else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                                    phi = np.arctan(NBX/NBY) + np.radians(90)
                        else:
                            phi = -3
                            
    return phi


def phi_angle_pixels(data_mask):
    NBX = data_mask.shape[1]
    NBY = data_mask.shape[0]
    XY_proportion = 2


    nbx = NBX
    lox = 0
    hix = data_mask.shape[1]

    nby = NBY
    loy = 0
    hiy = data_mask.shape[0]

    hist2d = TH2F("hist2d", "", nbx,lox,hix, nby,loy,hiy)
    for i in np.arange(0, nbx):
        for j in np.arange(0, nby):
            # cont = data_mask[i][j]
            if data_mask[j][i]:
                hist2d.SetBinContent(int(i+1),int(j+1), data_mask[j][i])
            else:
                # n = n+1
                # print(n)
                cont = 0
                hist2d.SetBinContent(int(i+1),int(j+1), cont)

    fitline = TF1("fitline", "[0] + [1]*x",lox,hix)  
    hist2d.Fit("fitline", "R")

    ### ===== Parámetros del Ajuste ===== ###
    ordenada = fitline.GetParameters()[0]
    pendiente = fitline.GetParameters()[1]
    Prob = fitline.GetProb()
    Chi2 = fitline.GetChisquare()
    ### ================================== ###

    flag_hor = False
    flag_ver = False

    ### ====== Esta sección determina si es un muon horizontal, vertical u oblicuo ===== ###
    if NBX > NBY * XY_proportion:
        flag_hor = True
        c_l = data_mask[:, 0:int(NBX/2)]
        c_r = data_mask[:, int(NBX/2):NBX] 

    elif NBY > NBX * XY_proportion:
        flag_ver = True
        c_d = data_mask[:int(NBY/2), :]
        c_u = data_mask[int(NBY/2):NBY, :]
        
    else: 
        c_ld = data_mask[0:int(NBY/2), 0:int(NBX/2)]
        c_lu = data_mask[int(NBY/2):NBY, 0:int(NBX/2)]
        c_rd = data_mask[0:int(NBY/2),int(NBX/2):NBX]
        c_ru = data_mask[int(NBY/2):NBY, int(NBX/2):NBX]


    ### ====== Esta sección determina el ángulo phi con la pendiente del ajuste === ###
    if flag_hor: ## Es un muon horizontal
        if pendiente > 0:
            ## El muon puede estar en el cuadrante 1 o 3
            # print('pendiente positiva')

            n_cl = 0
            n_cr = 0
            for index_y in np.arange(0, len(c_l)):
                # print(index)
                for index_x in np.arange(0, len(c_l[0])):
                    # print(index_x)
                    if c_l[index_y][index_x] != 0:
                        n_cl = n_cl +1

            for index_y in np.arange(0, len(c_r)):
                # print(index)
                for index_x in np.arange(0, len(c_r[0])):
                    # print(index_x)
                    if c_r[index_y][index_x] != 0:
                        n_cr = n_cr +1

            if n_cl < n_cr:
                # print('El muon está en el sector 1')
                phi = np.arctan(pendiente)

            elif n_cl > n_cr:
                # print('El muon está en el sector 3')
                phi = np.arctan(pendiente) + np.pi
        else: 
            # print('pendiente negativa')
            ### El muon puede estar en el cuadrante 2 o 4

            n_cl = 0
            n_cr = 0
            for index_y in np.arange(0, len(c_l)):
                # print(index)
                for index_x in np.arange(0, len(c_l[0])):
                    # print(index_x)
                    if c_l[index_y][index_x] != 0:
                        n_cl = n_cl +1

            for index_y in np.arange(0, len(c_r)):
                # print(index)
                for index_x in np.arange(0, len(c_r[0])):
                    # print(index_x)
                    if c_r[index_y][index_x] != 0:
                        n_cr = n_cr +1

            if n_cl > n_cr:
                # print('El muon está en el sector 2')
                phi = np.arctan(pendiente)  + np.pi/2

            elif n_cl < n_cr:
                # print('El muon está en el sector 4')
                phi = np.arctan(pendiente) + 3 * np.pi/ 2

    if flag_ver: ## Es un muon horizontal
        if pendiente > 0:
            ## El muon puede estar en el cuadrante 1 o 3
            # print('pendiente positiva')

            n_cu = 0
            n_cd = 0
            for index_y in np.arange(0, len(c_d)):
                # print(index)
                for index_x in np.arange(0, len(c_d[0])):
                    # print(index_x)
                    if c_d[index_y][index_x] != 0:
                        n_cd = n_cd + 1

            for index_y in np.arange(0, len(c_u)):
                # print(index)
                for index_x in np.arange(0, len(c_u[0])):
                    # print(index_x)
                    if c_u[index_y][index_x] != 0:
                        n_cu = n_cu + 1

            if n_cu > n_cd:
                # print('El muon está en el sector 1')
                phi = np.arctan(pendiente)

            elif n_cu < n_cd:
                # print('El muon está en el sector 3')
                phi = np.arctan(pendiente) + np.pi
        else: 
            # print('pendiente negativa')
            ### El muon puede estar en el cuadrante 2 o 4

            n_cu = 0
            n_cd = 0
            for index_y in np.arange(0, len(c_u)):
                # print(index)
                for index_x in np.arange(0, len(c_u[0])):
                    # print(index_x)
                    if c_u[index_y][index_x] != 0:
                        n_cu = n_cu + 1

            for index_y in np.arange(0, len(c_d)):
                # print(index)
                for index_x in np.arange(0, len(c_d[0])):
                    # print(index_x)
                    if c_d[index_y][index_x] != 0:
                        n_cd = n_cd + 1

            if n_cu > n_cd:
                # print('El muon está en el sector 2')
                phi = np.arctan(pendiente)  + np.pi/2

            elif n_cu < n_cd:
                # print('El muon está en el sector 4')
                phi = np.arctan(pendiente) + 3 * np.pi/ 2


    if not flag_ver and not flag_hor: ## Otros casos
        if pendiente > 0:
            ## El muon puede estar en el cuadrante 1 o 3
            # print('pendiente positiva')

            n_cu = 0
            n_cd = 0
            for index_y in np.arange(0, len(c_ld)):
                # print(index)
                for index_x in np.arange(0, len(c_ld[0])):
                    # print(index_x)
                    if c_ld[index_y][index_x] != 0:
                        n_cd = n_cd + 1

            for index_y in np.arange(0, len(c_ru)):
                # print(index)
                for index_x in np.arange(0, len(c_ru[0])):
                    # print(index_x)
                    if c_ru[index_y][index_x] != 0:
                        n_cu = n_cu + 1

            if n_cu > n_cd:
                # print('El muon está en el sector 1')
                phi = np.arctan(pendiente)

            elif n_cu < n_cd:
                # print('El muon está en el sector 3')
                phi = np.arctan(pendiente) + np.pi
        else: 
            # print('pendiente negativa')
            ### El muon puede estar en el cuadrante 2 o 4

            n_cu = 0
            n_cd = 0
            for index_y in np.arange(0, len(c_lu)):
                # print(index)
                for index_x in np.arange(0, len(c_lu[0])):
                    # print(index_x)
                    if c_lu[index_y][index_x] != 0:
                        n_cu = n_cu + 1

            for index_y in np.arange(0, len(c_rd)):
                # print(index)
                for index_x in np.arange(0, len(c_rd[0])):
                    # print(index_x)
                    if c_rd[index_y][index_x] != 0:
                        n_cd = n_cd + 1

            if n_cu > n_cd:
                # print('El muon está en el sector 2')
                phi = np.arctan(pendiente)  + np.pi/2

            elif n_cu < n_cd:
                # print('El muon está en el sector 4')
                phi = np.arctan(pendiente) + 3 * np.pi/ 2
        
    return phi


def muon_filter(dataCal, label_img, nlabels_img, prop, Solidit, Elipticity):
    CCD_depth = 680 ## micras
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
    list_elip = []
    list_sold = []

    list_charge_all_events = []
    list_elip_all_events = []
    list_sol_all_events = []

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
            list_sol_all_events.append(Solidity)
            list_elip_all_events.append(elip)
            continue 

        elif maxx - minx <= 3:
            list_charge_all_events.append(charge)
            list_sol_all_events.append(Solidity)
            list_elip_all_events.append(elip)
            continue

        elif maxy - miny <= 3:
            list_charge_all_events.append(charge)
            list_sol_all_events.append(Solidity)
            list_elip_all_events.append(elip)
            continue

        elif not Barycentercharge:
            list_charge_all_events.append(charge)
            list_sol_all_events.append(Solidity)
            list_elip_all_events.append(elip)
            continue

        elif differval < MeanValue_Event: 
            continue

        elif  Solidity < Solidit:
            list_charge_all_events.append(charge)
            list_sol_all_events.append(Solidity)
            list_elip_all_events.append(elip)
            continue 

        elif elip < Elipticity:
            list_charge_all_events.append(charge)
            list_sol_all_events.append(Solidity)
            list_elip_all_events.append(elip)
            continue

        elif  elip >= Elipticity :
            # charge = data_maskEvent.sum()

            # if charge > 100:
            Delta_EL = (charge)/ (Delta_L) 
            theta = np.arctan((Diagonal_lenght * px_to_cm)/(CCD_depth * micra_to_cm))

            try:
                phi = phi_angle_pixels(data_maskEvent)
            except:
                phi = -4



            list_DeltaL.append(Delta_L)
            list_DeltaEL.append(Delta_EL)
            list_charge.append(charge)
            list_theta.append(theta)
            list_phi.append(phi)
            list_elip.append(elip)
            list_sold.append(Solidity)
            # print(charge, DeltaEL)

            # if DeltaEL_range_min <= DeltaEL <= DeltaEL_range_max:
            list_Muon_labels.append(event)

    return list_DeltaL, list_DeltaEL, list_charge, list_Muon_labels, list_theta, list_phi, list_charge_all_events, list_elip, list_sold, list_elip_all_events, list_sol_all_events


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


#### =================== FUNCIONES DE CLUSTERIZACIÓN Y CREACCIÓN DE PDFs ================== ###

def all_cluster(dataCal, label_img, nlabels_img, prop):
    list_charge = []

    for event in np.arange(1, nlabels_img):
        mask = np.invert(label_img == event)
        loc = nd.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

        ## Aquí se calcula la carga total del cluster
        charge = data_maskEvent.sum()
        list_charge.append(charge)

    return list_charge
