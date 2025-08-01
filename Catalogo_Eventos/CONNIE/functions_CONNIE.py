import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import numpy.ma as ma
import pandas as pd 
import skimage as sk
import scipy.ndimage as nd
import scipy.ndimage as ndimage

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

def phi_angle_pixels(data_mask, pendiente, flag_rot):
    NBX = data_mask.shape[1]
    NBY = data_mask.shape[0]
    XY_proportion = 2

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
    flag_control = False
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
                flag_control = True

            elif n_cl > n_cr:
                # print('El muon está en el sector 3')
                phi = np.pi + np.arctan(pendiente) 
                flag_control = True
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
                phi = np.pi + np.arctan(pendiente) 
                flag_control = True 

            elif n_cl < n_cr:
                # print('El muon está en el sector 4')
                phi = 2 * np.pi + np.arctan(pendiente) 
                flag_control = True

    elif flag_ver: ## Es un muon vertical
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
                flag_control = True

            elif n_cu < n_cd:
                # print('El muon está en el sector 3')
                phi = np.arctan(pendiente) + np.pi
                flag_control = True
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
                phi = np.pi + np.arctan(pendiente)
                flag_control = True

            elif n_cu < n_cd:
                # print('El muon está en el sector 4')
                phi =  2* np.pi + np.arctan(pendiente)
                flag_control = True

    elif not flag_ver and not flag_hor: ## Otros casos
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
                flag_control = True

            elif n_cu < n_cd:
                # print('El muon está en el sector 3')
                phi =  np.pi + np.arctan(pendiente)
                flag_control = True
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
                phi = np.pi + np.arctan(pendiente)
                flag_control = True

            elif n_cu < n_cd:
                # print('El muon está en el sector 4')
                phi = 2 * np.pi + np.arctan(pendiente)
                flag_control = True

    if flag_rot and flag_control:
        if np.pi/2 < phi < np.pi:
            # print('Estoy en el sector 1 (rotado)')
            phi = phi - np.pi/2
        elif np.pi < phi < 3*np.pi/2:
            # print('Estoy en el sector 2 (rotado)')
            phi = phi - np.pi/2
        elif 3*np.pi/2 < phi < 2 * np.pi:
            # print('Estoy en el sector 3 (rotado)')
            phi = phi - np.pi/2
        elif 0 < phi < np.pi/2:
            # print('Estoy en el sector 4 (rotado)')
            phi = phi + 3*np.pi/2 
    
    return phi

def pixel_rot(x_bin, x0, y_bin, y0, theta):
    diff_x = x_bin - x0
    diff_y = y_bin - y0

    new_x = diff_x * np.cos(theta) - diff_y * np.sin(theta) + x0
    new_y = diff_x * np.sin(theta) + diff_y * np.cos(theta) + y0

    # return int(np.around(new_x, 0)), int(np.around(new_y, 0))
    return int(new_x), int(new_y)

def linear_fit(data_mask):
    NBX = data_mask.shape[1]
    NBY = data_mask.shape[0]
    flag_rot = False

    nbx = NBX
    lox = 0
    hix = data_mask.shape[1]

    nby = NBY
    loy = 0
    hiy = data_mask.shape[0]

    Qs = 0

    hist2d = TH2F("hist2d", "", nbx,lox,hix, nby,loy,hiy)
    for i in np.arange(0, nbx):
        for j in np.arange(0, nby):
            q = data_mask[j][i]
            if q != 0:
                Qs += q
                hist2d.SetBinContent(int(i+1),int(j+1), q)
            else:
                hist2d.SetBinContent(int(i+1),int(j+1), 0)

    fitline = TF1("fitline", "[0] + [1]*x",lox,hix)  
    hist2d.Fit("fitline", "RNQ")

    ### ===== Parámetros del Ajuste ===== ###
    ordenada = fitline.GetParameters()[0]
    pendiente = fitline.GetParameters()[1]
    Chi_2 = fitline.GetChisquare()
    ### ================================== ###

    if Chi_2/Qs > 2:
        # print('Intenté el ajuste pero fallé.')
        flag_rot = True
        long_y = data_mask.shape[0] - 1

        data_mask_zeros = np.empty((NBX, NBY))
        # data_mask_rot[:] = np.nan
        data_mask_zeros[:] = 0

        angle_rot = TMath.Pi()/2

        for y_bin in range(0, NBY):
            for x_bin in range(0, NBX):
                if data_mask[y_bin][x_bin] != 0:
                    # nx, ny = pixel_rot(x_bin=x_bin, x0=0, y_bin=y_bin, y0=0, theta= angle_rot)
                    data_mask_zeros[x_bin][long_y - y_bin] = data_mask[y_bin][x_bin]

        label_img, nlabels_img = sk.measure.label(data_mask_zeros > 0, connectivity=2, return_num=True)
        # loc_rot = ndimage.find_objects(label_img == 1)[0]
        mask_rot = np.invert(label_img==1)
        data_mask_rot = ma.masked_array(data_mask_zeros, mask_rot)

        
        NBX = data_mask_rot.shape[1]
        NBY = data_mask_rot.shape[0]

        nbx = NBX
        lox = 0
        hix = data_mask_rot.shape[1]

        nby = NBY
        loy = 0
        hiy = data_mask_rot.shape[0]

        Qs = 0

        hist2d_rot = TH2F("hist2d_rot", "", nbx,lox,hix, nby,loy,hiy)
        for j in range(0, nby):
            for i in range(0, nbx):
                q = data_mask_rot[j][i]
                if q != 0:
                    Qs += q
                    hist2d_rot.SetBinContent(int(i+1),int(j+1), q)
                else:
                    hist2d_rot.SetBinContent(int(i+1),int(j+1), 0)

        fitline = TF1("fitline", "[0] + [1]*x",lox,hix)  
        hist2d_rot.Fit("fitline", "RNQ")

        ### ===== Parámetros del Ajuste ===== ###
        ordenada = fitline.GetParameters()[0]
        pendiente = fitline.GetParameters()[1]
        Chi_2 = fitline.GetChisquare()
        data_mask = data_mask_rot
        # print('Ord: ', ordenada, 'Pend: ', pendiente, 'Chi2: ', Chi_2/Qs, 'flag_rot: ', flag_rot)
        ### ================================== ###

    return ordenada, pendiente, Chi_2/Qs, flag_rot, data_mask

def muon_filter(dataCal, label_img, nlabels_img, prop, Solidit, Elipticity, dedl_min):
    CCD_depth = 725 ## micras
    px_to_micras = 15 ## micras
    px_to_cm = 0.0015 ## cm/px
    micra_to_cm = 1 / 10000 ## micras/cm

    list_Muon_labels = []
    list_DeltaEL = []
    list_DeltaL = []
    list_charge = []
    list_theta = []
    list_phi = []
    list_elip = []
    list_sold = []
    list_datamasked = []

    list_charge_all_events = []
    list_elip_all_events = []
    list_sol_all_events = []

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

        ## Aquí se obtiene los radios de la elipse del cluster
        rM = prop[event-1].axis_major_length / 2
        rm = prop[event-1].axis_minor_length / 2

        ## Tensor de inercia ##
        # inercia_tensor = prop[event-1].inertia_tensor
        # non_diag_inercia_tensor = inercia_tensor[0][1]

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

            if Delta_EL < dedl_min:
                list_charge_all_events.append(charge)
                list_sol_all_events.append(Solidity)
                list_elip_all_events.append(elip)
                continue
                

            # print('Estoy por entrar a linear_fit')
            ord, pend, chi2_Qs, flag_rot, data_mask = linear_fit(data_mask=data_maskEvent)
            # print('Ord: ', ord, 'Pend: ', pend, 'Chi2: ', chi2_Qs, 'flag_rot: ', flag_rot)

            if chi2_Qs > 2:
                list_charge_all_events.append(charge)
                list_sol_all_events.append(Solidity)
                list_elip_all_events.append(elip)
                continue
            
            #### ======================== CÁLCULO DEL ÁNGULO THETA ============================ ###
            #### ---------  Se toma que TODOS los muones atravezaron por completo la CCD ------ ###
            theta = np.arctan((Diagonal_lenght * px_to_cm)/(CCD_depth * micra_to_cm)) 

            ### ============ CÁLCULO DEL ÁNGULO PHI (pixels) ===================== ###
            try:
                # print('Phi angle function')
                phi = phi_angle_pixels(data_mask, pend, flag_rot)

            except:
                phi = -4
                
            if phi < 0:
                continue

            list_phi.append(phi)    
            list_DeltaL.append(Delta_L)
            list_DeltaEL.append(Delta_EL)
            list_charge.append(charge)
            list_theta.append(theta)
            list_elip.append(elip)
            list_sold.append(Solidity)
            list_datamasked.append(data_maskEvent)

            list_Muon_labels.append(event)

    dict_lists = {"muons" : {"l": list_DeltaL, "dedl" : list_DeltaEL, "charge_muons" : list_charge, "theta": list_theta, 
                             "phi": list_phi, "elip": list_elip, "sol": list_sold, "image" : list_datamasked,
                             "label" : list_Muon_labels}, 

                 "non_muons" : {"charge": list_charge_all_events, "elip": list_elip_all_events, "sol": list_sol_all_events}}
    
    return dict_lists

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

        elif maxy - miny <= 3:
            continue

        elif not Barycentercharge:
            continue

        elif differval < MeanValue_Event: #keV
            continue

        elif  Solidity < Solidit:
            continue 

        elif  elip >= Elipticity :
            charge = data_maskEvent.sum()

            # if charge < min_Charge:
            #     continue
            
            if (Longitud_x < 12 and Longitud_y > 12): 
                # num_muons = num_muons + 1

                list_sigmas_vertical_event.append(Sigma)
                list_vertical_event.append(data_maskEvent)
                list_charge_vertical_event.append(charge)

            if ( Longitud_y < 12 and Longitud_x > 12):
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

def DataFrame_muons(dict_muons, extension):
    KeV_elec_ratio = 0.00368    # KeV/e-
    if extension == 1:
        dict_extension = dict_muons['extension_1']
    elif extension == 2:
        dict_extension = dict_muons['extension_2']
    elif extension == 4:
        dict_extension = dict_muons['extension_4']

    list_datamask = dict_extension['datamasked']
    DF_charge = pd.DataFrame(dict_extension['charge'] * np.array(KeV_elec_ratio), columns=['Charge (KeV)'])
    DF_sol = pd.DataFrame(dict_extension['sol'], columns=['Solidity'])
    DF_eli = pd.DataFrame(dict_extension['elip'], columns=['Elipticity'])
    DF_thet = pd.DataFrame(np.degrees(dict_extension['theta']), columns=['Theta (Deg)'])
    DF_phi = pd.DataFrame(np.degrees(dict_extension['phi']), columns=['Phi (Deg)'])
    DF_dedl = pd.DataFrame(dict_extension['deltaEL'] * np.array(KeV_elec_ratio), columns=['dEdL (KeV/cm)'])
    DF_l = pd.DataFrame(dict_extension['deltaL'], columns=['l (cm)'])

    list_muonid =[]
    for index in range(0, len(list_datamask)):
        list_muonid.append(index)
    DF_muonid = pd.DataFrame(list_muonid, columns=['Muon ID'])

    TotalFrame = pd.concat([DF_muonid, DF_sol, DF_eli, DF_thet, DF_phi, DF_charge, DF_l, DF_dedl], axis=1)
    Frame = TotalFrame.set_index('Muon ID')

    return Frame, list_datamask

def all_cluster(dataCal, label_img, nlabels_img, prop):
    list_charge = []
    list_xsize = []
    list_ysize = []

    for event in np.arange(1, nlabels_img):
        mask = np.invert(label_img == event)
        loc = nd.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])
        
        ysize, xsize = data_maskEvent.shape[0], data_maskEvent.shape[1]

        ## Aquí se calcula la carga total del cluster
        charge = data_maskEvent.sum()
        list_charge.append(charge)
        list_xsize.append(xsize)
        list_ysize.append(ysize)

    return list_charge, list_xsize, list_ysize

def DataFrame_allclusters(dict_all, units):
    KeV_elec_ratio = 0.00368    # KeV/e-
    dict_extension = dict_all['extension_1']
    Len = len(dict_extension['charge'])

    if units == 1:
        DF_charge = pd.DataFrame(dict_extension['charge'], columns=['Charge (e-)'])
    elif units == 2:
        DF_charge = pd.DataFrame(dict_extension['charge'] * np.array(KeV_elec_ratio), columns=['Charge (KeV)'])
        
    DF_xsize = pd.DataFrame(dict_extension['xsize (px)'], columns=['Xsize (px)'])
    DF_ysize = pd.DataFrame(dict_extension['ysize (px)'], columns=['Ysize (px)'])

    list_clusterid =[]
    for index in range(0, Len):
        list_clusterid.append(index)
    DF_clusterid = pd.DataFrame(list_clusterid, columns=['Cluster ID'])

    TotalFrame = pd.concat([DF_clusterid, DF_xsize, DF_ysize, DF_charge], axis=1)
    Frame = TotalFrame.set_index('Cluster ID')
    return Frame

### =================== 