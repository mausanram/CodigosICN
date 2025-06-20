from functions_py import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd 
import skimage as sk
import scipy.ndimage as nd
import random as rand
import time

from ROOT import TMath, TF1, TH1F, TH2F, TCanvas, gStyle, TProfile, TGraphErrors

### Distribución Gaussiana ###
def gaussian(x, a, mean, sigma): 
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

### Distribución de Doble Gaussiana ###
def Gaussian2(x,m,s,g,a1,a2): #data, mean, sigma, gain, height1, heigth2
    return a1*np.exp(-1/2*((x-m)/s)**2)+a2*np.exp(-1/2*((x-m-g)/s)**2)

### ============================== Funciones de calibración de imágenes ============================== ###
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

def oScan_fit_NSAMP324_ROOT(extensión, active_area, oScan, Bins, Bins_fit, make_figure_flag, range_fit):
    Range_x_label = range_fit

    min_oScan = np.min(oScan)

    if make_figure_flag == True:
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins,  range=(min_oScan, 18000))
        offset = bins_edges[np.argmax(hist)]
        # print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset 

        fgaus2 = TF1("fgauss2","[3]*exp(-0.5*((x-[0])/[1])^2)+[4]*exp(-0.5*((x-[0]-[2])/[1])^2)",-300,600,5) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)

        h3=TH1F("histogram", "Distribution of OsCan", Bins_fit, Range_x_label[0], Range_x_label[1])
        for pixel_value in Overscan_plane.flatten():
            # if not np.ma.is_masked(pixel_value):
            h3.Fill(pixel_value)
            #print(pixel_value)

        fgaus2.SetParameters(0,10,100, 100, 100) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
        h3.Fit(fgaus2, "R")

        canv=TCanvas()
        canv.SetLogy()
        h3.SetStats(0)
        h3.Draw()
        fgaus2.Draw("same")
        canv.Draw()

        gStyle.SetOptFit(1100)
        gStyle.SetPadGridX (True)
        # fgaus2.Draw('Quiet')


        print('Parameters of the Doble-Gaussian Fit')
        print('Mean: ', fgaus2.GetParameters()[0],  ' +- ', fgaus2.GetParError(0))
        print('Sigma: ', fgaus2.GetParameters()[1],  ' +- ', fgaus2.GetParError(1))
        print('Gain: ', fgaus2.GetParameters()[2],  ' +- ', fgaus2.GetParError(2), '\n')

        print("chiSquare: " + str(fgaus2.GetChisquare()))
        print("NDegrees of Freedom: " + str(fgaus2.GetNDF()))
        print("chiSquare / NDF :", fgaus2.GetChisquare() / fgaus2.GetNDF())
        print("Prob:", fgaus2.GetProb(), '\n')

        dict_popt = {'Mean' :fgaus2.GetParameters()[0], 'sigma' : abs(fgaus2.GetParameters()[1]), 'Gain' : abs(fgaus2.GetParameters()[2]), 
                     'Offset' : offset, 'Prob' :  fgaus2.GetProb()}
        
    elif make_figure_flag == False:
        hist , bins_edges = np.histogram(oScan.flatten(), bins = Bins,  range=(min_oScan, 18000))
        offset = bins_edges[np.argmax(hist)]
        # print('Offset Value: ', offset, ' ADUs')

        Overscan_plane = oScan - offset

        fgaus2 = TF1("fgauss2","[3]*exp(-0.5*((x-[0])/[1])^2)+[4]*exp(-0.5*((x-[0]-[2])/[1])^2)", range_fit[0], range_fit[1],5) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)
        
        h3=TH1F("histogram", "Distribution of OsCan",Bins_fit, range_fit[0], range_fit[1])
        for pixel_value in Overscan_plane.flatten():
            # if not np.ma.is_masked(pixel_value):
            h3.Fill(pixel_value)
            #print(pixel_value)
        fgaus2.SetParameters(0,10,100, 100, 100) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
        h3.Fit(fgaus2, "RQN")

        dict_popt = {'Mean' :fgaus2.GetParameters()[0], 'sigma' : abs(fgaus2.GetParameters()[1]), 'Gain' : abs(fgaus2.GetParameters()[2]), 
                     'Offset' : offset, 'Prob' :  fgaus2.GetProb()}
    
    return dict_popt

def data_calibrated(active_area, extension, offset, list_gain, ratio_keV, unidades):
    dataP = active_area

    if unidades == 0:
        data = dataP

    elif unidades == 1:
        data = dataP / list_gain[extension - 1]

    elif unidades == 2:
        data = (ratio_keV * dataP) / list_gain[extension - 1]

    return data

def data_calibrated_NSAMP(active_area, extension, offset, gain, ratio_keV, unidades, sigma_ADUs):
    ## NO se aplica el offset porque ya se le aplicó la mediana del OsCan##
    dataP = active_area ## En ADUs

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
### =================================================================================================== ###

### ================ Gaus-Poisson convolution fit (only for NSAMP324) ================= ###
def gauss_comppoisson_fit(x, par):
    k =5
    #  m = 4
    #  ydata = 0;
    xval = x[0]
    a     = par[0]
    mu    = par[1]
    sigma = par[2]
    lambda_poisson = par[3]
    pgeom = par[4]
    gain  = par[5]

    fitval = 0.0
    for p in np.arange(0, k):
        fitval = fitval + a * TMath.Gaus(xval*gain,p-mu,sigma,1) * TMath.PoissonI(p,lambda_poisson)

    return fitval

def fit_gausCONVcomppois(oScan_data,):
    data = oScan_data.flatten() ## get oScan data and turn 1D array
    long_data = len(data)
    NBins = long_data - 1

    histo = TH1F("histo","",NBins, -100, 350)
    histo.GetXaxis().SetTitle("charge (e^{-})")

    # for ibin in np.arange(0, NBins):
    #     histo.SetBinContent(ibin,BinCont[ibin])

    # nevents = histo.Integral()
    # print('N_events: ', nevents)

    npar = 2
    norm = 3091.7   
    offs = 0.020    # offset
    sigm = 0.211    # sigma
    lamb = 0.163    # lambda de Poisson
    pgeo = 0.083    # probabilidad (se ignora en el cálculo)
    gain = 1.014

    lofit = -0.5
    hifit =  4.5

    fitf = TF1("fitf",gauss_comppoisson_fit,lofit,hifit,npar)
    fitf.SetParameter(0,norm)
    fitf.SetParameter(1,offs)
    fitf.SetParameter(2,sigm)
    fitf.SetParameter(3,lamb)
    fitf.SetParameter(4,pgeo)
    fitf.SetParameter(5,gain)

    # fitf.SetNpx(400)
    # fitf.SetMinimum(1e-3)
    fitf.SetLineWidth(1)

    # histo.Fit("fitf","R")
    norm = fitf.GetParameter(0)
    offs = fitf.GetParameter(1)
    sigm = fitf.GetParameter(2)
    lamb = fitf.GetParameter(3)
    pgeo = fitf.GetParameter(4)
    gain = fitf.GetParameter(5)

    norme = fitf.GetParError(0)
    offse = fitf.GetParError(1)
    sigme = fitf.GetParError(2)
    lambe = fitf.GetParError(3)
    pgeoe = fitf.GetParError(4)
    gaine = fitf.GetParError(5)

    chisq = fitf.GetChisquare()
    ndegf = fitf.GetNDF()
    proba = fitf.GetProb()

    canv = TCanvas("canv","",600,400)
    histo.Draw("h same")
    fitf.Draw("same")

    canv.Draw()

    dict_info = {'par' : {'Offset':offs, 'Sigma':sigm, 'Gain':gain}, 'err' :  {'Offset':offse, 'Sigma':sigme, 'Gain':gaine},
                 'fit_quality': {'Chiq' : chisq, 'Ndegf' : ndegf, 'Prob': proba}}
                 
    return dict_info
### ==================================================================================== ###

### =========== Funciones para determinar ángulo Phi ================== ###
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

def phi_angle_ROOT_pendpos(data_mask):
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
        GRprofY.Fit("fitline")

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


    elif NBY > 3 * NBX: ### El muon es muy vertical y se tomará la componente Y para hacer los ajustes
        lox = 0
        hix = NBX

        pend = (NBY)/(NBX)
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
            fitline.SetParameters(0, pend)
            GRprofY.Fit("fitline")
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
                    if First > Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(270)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(90)

    return phi

def phi_angle_ROOT_pendneg(data_mask):
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
        fitline.SetParameters(NBY-1, -pend)
        GRprofY.Fit("fitline")
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


    elif NBY > 3 * NBX: ### El muon es muy vertical y se tomará la componente Y para hacer los ajustes
        lox = 0
        hix = NBX

        pend = (NBY)/(NBX)
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
            fitline.SetParameters(NBY-1, -pend)
            GRprofX.Fit("fitline")
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
                    # print('pendiente negativa')
                    if First < Last: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(90)
                    else: ## La "cola" está en la parte de arriba, y el muon está en el cuadrante 4
                        phi = np.arctan(NBX/NBY) + np.radians(270)

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
                phi = np.arctan(pendiente) + np.pi
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

    elif flag_ver: ## Es un muon horizontal
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
                phi = phi = 2 * np.pi + np.arctan(pendiente)
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
                phi = np.arctan(pendiente) + np.pi
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

### =================================================================== ###

def pixel_rot(x_bin, x0, y_bin, y0, theta):
    diff_x = x_bin - x0
    diff_y = y_bin - y0

    new_x = diff_x * np.cos(theta) - diff_y * np.sin(theta) + x0
    new_y = diff_x * np.sin(theta) + diff_y * np.cos(theta) + y0

    # return int(np.around(new_x, 0)), int(np.around(new_y, 0))
    return int(new_x), int(new_y)

### ================================ Filtro de Muones General ============================================ ###
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
        flag_rot = True
        long_y = data_mask.shape[0] - 1

        data_mask_rot = np.empty((NBX, NBY))
        data_mask_rot[:] = np.nan

        angle_rot = TMath.Pi()/2

        for y_bin in range(0, NBY):
            for x_bin in range(0, NBX):
                if data_mask[y_bin][x_bin] != 0:
                    # nx, ny = pixel_rot(x_bin=x_bin, x0=0, y_bin=y_bin, y0=0, theta= angle_rot)
                    data_mask_rot[x_bin][long_y - y_bin] = data_mask[y_bin][x_bin]

        
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
        print(ordenada, pendiente, Chi_2/Qs, flag_rot)
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

            _, pend, chi2_Qs, flag_rot, data_mask = linear_fit(data_mask=data_maskEvent)

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
                phi = phi_angle_pixels(data_mask, pend, flag_rot)

            except:
                phi = -4
                
            # if phi < 0:
            #     continue

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
### ====================================================================================================== ###

### ========================= Catálogo de Muones Rectos (verticales/horizontales) ========================= ###
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
### ======================================================================================================= ###

### ================ Funciones para el Modelo de Difusión ========================= ###
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

            try:
                popt, pcov = curve_fit(gaussian, list_xlabel, data_mask[0], maxfev=100000, p0 = [1000, Mean_in, sigma_in], method={'lm'})		# Fit histogram with gaussian

                ## Se guardan lo parámetros del ajuste en un diccionario ##
                dict_popt = {'Mean' : popt[1], 'Hight' : popt[0], 'sigma' : abs(popt[2]), 'Pcov' : pcov}
                Centroide = popt[1]
                Sigma = abs(popt[2])

                # print(Sigma)

                # if Sigma > 6:
                #     continue

            except: 
                continue

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

def diffution_vertical_muon_ROOT(dict, list_vertical_labels, Delta_in, Delta_fin, extension):

    list_all_sigmas = []
    list_deep = []

    CCD_depth = 725 # micras
    Delta_inicial = Delta_in    # px
    Delta_final = Delta_fin     # px

    for label_muon in list_vertical_labels:
        list_sigmas = []
        # event = data_histogram['extension_' + str(extension)]['Vertical_Events'][label_muon]
        event, flag = check_flip_vertical_muon(dict  = dict, label_muon = label_muon, Delta_in= Delta_inicial, Delta_fin=Delta_final, extension=extension)
        # print('Turn flag: ', flag)

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
            # list_xlabel = np.arange(0.5, len(data_mask[0]), 1)

            # list_xlabel_long = np.linspace(-Longitud_linea + int(Longitud_linea/2) , Longitud_linea + int(Longitud_linea/2), Longitud_linea)

            try:
                fgaus2 = TF1("fgauss","gaus", -7,  7, 3) # TF1("nombre", "funcion escrita como en root", min, max, #parametros)
                h3 = TH1F("histogram", "Distribution of Line", 50, -7,  7)
                
                for pixel_value in data_mask[0].flatten():
                    # if not np.ma.is_masked(pixel_value):
                    # print(pixel_value)
                    h3.Fill(pixel_value)
                        #print(pixel_value)

                fgaus2.SetParameters(1, Mean_in, sigma_in) # Establecer parametros iniciales del fit, de manera visual es posible determinarlos como una primera aproximacion
                h3.Fit(fgaus2)

                # dict_popt = {'Mean' :fgaus2.GetParameters()[1], 'sigma' : abs(fgaus2.GetParameters()[2])}
                Sigma = abs(fgaus2.GetParameters()[2])

            except: 
                continue

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
### =============================================================================== ###

#### =================== FUNCIONES DE CLUSTERIZACIÓN Y CREACCIÓN DE PDFs ================== ###

def all_cluster(dataCal, label_img, nlabels_img, prop):
    list_charge = []

    for event in np.arange(1, nlabels_img):
        mask = np.invert(label_img == event)
        loc = ndimage.find_objects(label_img == event)[0]
        
        data_maskEvent = ma.masked_array(dataCal[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop],
                                            mask[loc[0].start:loc[0].stop, loc[1].start:loc[1].stop])

        ## Aquí se calcula la carga total del cluster
        charge = data_maskEvent.sum()
        list_charge.append(charge)

        # if DeltaEL_range_min <= DeltaEL <= DeltaEL_range_max:
        # list_Muon_labels.append(event)

    return list_charge


