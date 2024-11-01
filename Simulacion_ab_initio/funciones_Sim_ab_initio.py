import numpy as np
import mpmath as mp
import random as rand
import datetime
import os
import subprocess
import matplotlib.pyplot as plt
from array import array

from ROOT import TMath, TF1

def dis_probability(theta, I_0):
    return I_0 * np.cos(theta)

def dis_angular(theta): ## Distribucion angular
    return 1 * np.cos(theta)**2 * np.sin(theta)

def dis_energy(E_mu, theta, units): ### Modelo de Smith-Duller (units: 0 == MeV, units: 1 == GeV)
    ## Constantes físicas
    k = 8 / 3
    b = 0.771
    lambda_pi = 120     ## g/cm^2
    y_0 = 1000      ## g/cm^2
    r = 0.76
    b_mu = 0.8      

    if units == 0:
        a = 2.5         ## MeV cm^2/g
        m_mu = 105.7    ## MeV/c^2 
        m_pi = 139.6    ## MeV/c^2

    elif units == 1:
        a = 0.0025         ## GeV cm^2/g
        m_mu = 0.1057    ## GeV/c^2 
        m_pi = 0.1396    ## GeV/c^2

    tau_mu_0 = 2.2 * 10**(-6)   ## s
    tau_0 = 2.6 * 10 **(-8)     ## s
    rho_0 = 0.00129 ## g/cm^3
    c = 3 * 10 ** 10 ## cm/s

    ### Parámetros
    E_pi = (1 / r) * (E_mu + a * y_0 * ((1/np.cos(theta)) - 0.1))
    B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c)
    P_mu = ((0.1 * np.cos(theta)) * (1 - (a * (y_0 *(1/np.cos(theta)) - 100))/( r * E_pi)) ) ** ((B_mu)/((r * E_pi + 100 * a) * np.cos(theta)))
    j_pi = (m_pi * y_0)/(c * tau_0 * rho_0)

    # E_pi = (1 / r) * (E_mu + a * y_0 * (mp.sec(theta) - 0.1))
    # B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c)
    # P_mu = ((0.1 * np.cos(theta)) * (1 - (a * (y_0 *mp.sec(theta) - 100))/( r * E_pi)) ) ** ((B_mu)/((r * E_pi + 100 * a) * np.cos(theta)))
    # j_pi = (m_pi * y_0)/(c * tau_0 * rho_0)

    ## Intensidad diferencial
    # C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    C_2 = E_pi * np.cos(theta)
    C_3 = b * j_pi

    # return (C_1 * np.sin(theta)) / (C_2 + C_3)
    return (C_1 ) / (C_2 + C_3)

def coord_cartesian(Thet, Phi): ##### Coordenadas Cartesianas
    coord_X = np.sin(Thet) * np.cos(Phi)
    coord_Y = np.sin(Thet) * np.sin(Phi)
    coord_Z = np.cos(Thet)

    Vec_sph = [float(coord_X), float(coord_Y), float(coord_Z)]
    return Vec_sph 

def norma_vec(Vector):
    norma = np.sqrt(Vector[0] ** 2 + Vector[1] ** 2  + Vector[2] ** 2 )
    return norma

def Energy_list_log(lim_inf, lim_sup, N):
    list_Energy = []
    Delta_log = (np.log10(lim_sup) - np.log10(lim_inf))/ (N - 1)

    for i in np.arange(0, N):
        En = 10 ** (np.log10(lim_inf) + i * Delta_log)
        list_Energy.append(En)

    return list_Energy

def dimension_x(long_x):
    step = 0.0001

    list_long_x = [-long_x]

    while long_x:
        x = np.round(list_long_x[-1] + step, 4)
        list_long_x.append(x)

        if list_long_x[-1] == long_x:
            break

    return list_long_x

def dimension_y(long_y):
    step = 0.0001

    list_long_y = [-long_y]

    while long_y:
        y = np.round(list_long_y[-1] + step, 4)
        list_long_y.append(y)

        if list_long_y[-1] == long_y:
            break

    return list_long_y

def dimension_z(long_z):
    step = 0.0001

    list_long_z = [-long_z]

    while long_z:
        z = np.round(list_long_z[-1] + step, 4)
        list_long_z.append(z)

        if list_long_z[-1] == long_z:
            break

    return list_long_z

# def intersection_CCD(flags_CCD, list_z, medida_z, Random_th ):
#     # list_delta_L = []
#     delta_L = None

#     ## Caras 1 y 2
#     if flags_CCD[0] and flags_CCD[1]:
#         delta_L = np.abs(medida_z / np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         n_muons_in_CCD = 1
#         # list_delta_L.append(delta_L)
#         # continue

#     ## Caras 1 y 3
#     if flags_CCD[0] and flags_CCD[2]:
#         h = medida_z - list_z[0]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 1 y 4
#     if flags_CCD[0] and flags_CCD[3]:
#         h = medida_z - list_z[1]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 1 y 5
#     if flags_CCD[0] and flags_CCD[4]:
#         h = medida_z - list_z[2]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 1 y 6
#     if flags_CCD[0] and flags_CCD[5]:
#         h = medida_z - list_z[3]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 3 y 2
#     if flags_CCD[2] and flags_CCD[1]:
#         h = list_z[0]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 3 y 4
#     if flags_CCD[2] and flags_CCD[3]:
#         h = list_z[0]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 3 y 5
#     if flags_CCD[2] and flags_CCD[4]:
#         h = list_z[0]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 3 y 6
#     if flags_CCD[2] and flags_CCD[5]:
#         h = list_z[0]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 4 y 2
#     if flags_CCD[3] and flags_CCD[1]:
#         h = list_z[1]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 4 y 3
#     if flags_CCD[3] and flags_CCD[2]:
#         h = list_z[1]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 4 y 5
#     if flags_CCD[3] and flags_CCD[4]:
#         h = list_z[1]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 4 y 6
#     if flags_CCD[3] and flags_CCD[5]:
#         h = list_z[1]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 5 y 2
#     if flags_CCD[4] and flags_CCD[1]:
#         h = list_z[2]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 5 y 3
#     if flags_CCD[4] and flags_CCD[2]:
#         h = list_z[2]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 5 y 4
#     if flags_CCD[4] and flags_CCD[3]:
#         h = list_z[2]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 5 y 6
#     if flags_CCD[4] and flags_CCD[5]:
#         h = list_z[2]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 6 y 2
#     if flags_CCD[5] and flags_CCD[1]:
#         h = list_z[3]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 6 y 3
#     if flags_CCD[5] and flags_CCD[2]:
#         h = list_z[3]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue
#     ## Caras 6 y 4
#     if flags_CCD[5] and flags_CCD[3]:
#         h = medida_z - list_z[3]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     ## Caras 6 y 5
#     if flags_CCD[5] and flags_CCD[4]:
#         h = list_z[3]
#         delta_L = np.abs(h /  np.abs(np.cos(Random_th)))  ## cm
#         # n_muons_in_CCD = n_muons_in_CCD + 1
#         # list_delta_L.append(delta_L)
#         n_muons_in_CCD = 1
#         # continue

#     # if (flags_CCD[0] == False) and  (flags_CCD[1] == False) and (flags_CCD[2] == False) and (flags_CCD[3] == False) and (flags_CCD[4] == False) and (flags_CCD[5] == False):
#     if delta_L is None:    
#         delta_L, n_muons_in_CCD = 0, 0

#     # if delta_L < 0.0725 and delta_L != 0:
#     #     print(delta_L)

#     return delta_L, n_muons_in_CCD


def intersection_CCD(flags_CCD, list_z, medida_z, Random_th ):
    # list_delta_L = []
    delta_L = None

    ## Caras 1 y 2
    if flags_CCD[0] and flags_CCD[1]:
        delta_L = medida_z / np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        n_muons_in_CCD = 1
        # list_delta_L.append(delta_L)
        # continue

    ## Caras 1 y 3
    if flags_CCD[0] and flags_CCD[2]:
        h = medida_z - list_z[0]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 4
    if flags_CCD[0] and flags_CCD[3]:
        h = medida_z - list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 5
    if flags_CCD[0] and flags_CCD[4]:
        h = medida_z - list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 6
    if flags_CCD[0] and flags_CCD[5]:
        h = medida_z - list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 2
    if flags_CCD[2] and flags_CCD[1]:
        h = list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 4
    if flags_CCD[2] and flags_CCD[3]:
        h = list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 5
    if flags_CCD[2] and flags_CCD[4]:
        h = list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 6
    if flags_CCD[2] and flags_CCD[5]:
        h = list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 2
    if flags_CCD[3] and flags_CCD[1]:
        h = list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 3
    if flags_CCD[3] and flags_CCD[2]:
        h = list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 5
    if flags_CCD[3] and flags_CCD[4]:
        h = list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 6
    if flags_CCD[3] and flags_CCD[5]:
        h = list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 2
    if flags_CCD[4] and flags_CCD[1]:
        h = list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 3
    if flags_CCD[4] and flags_CCD[2]:
        h = list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 4
    if flags_CCD[4] and flags_CCD[3]:
        h = list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 6
    if flags_CCD[4] and flags_CCD[5]:
        h = list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 2
    if flags_CCD[5] and flags_CCD[1]:
        h = list_z[3]
        delta_L =  (h /np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 3
    if flags_CCD[5] and flags_CCD[2]:
        h = list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue
    ## Caras 6 y 4
    if flags_CCD[5] and flags_CCD[3]:
        h = medida_z - list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 5
    if flags_CCD[5] and flags_CCD[4]:
        h = list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    # if (flags_CCD[0] == False) and  (flags_CCD[1] == False) and (flags_CCD[2] == False) and (flags_CCD[3] == False) and (flags_CCD[4] == False) and (flags_CCD[5] == False):
    if delta_L is None:    
        delta_L, n_muons_in_CCD = 0, 0

    # if delta_L < 0.0725 and delta_L != 0:
    #     print(delta_L)

    return delta_L, n_muons_in_CCD


## ----------------------- Funciones de Landau (ROOT) ------------------------------ ##
def LandV(lx, lpar):
    Delta = lx[0]	# Energy loss in absorber
    L = lpar[0]		# Thickness of absorber (Distance crossed by the particle)
    p = lpar[1]		# Momentum (in MeV/c)

    K = 0.1535	# K coefficient = 2*pi*N*r^2*m*c^2 (in MeV g^-1 cm^2)

    z = -1;		# Charge number of incident particle
    ZA = 0.498487	# Atomic number over Atomic mass of absorber (for Si)
    c = TMath.C()	# Speed of light
    me = 0.510998928	#  Electron mass in MeV/c^2
    M = 105.65839	# Muon mass in MeV/c^2
    I = 0.000173		# Mean excitation energy (for Si)
    bg = p/M
    beta = bg/np.sqrt(1+bg**2)	# Beta factor
    gamma = 1/np.sqrt(1-beta**2)	# Gamma factor
    pi = TMath.Pi()
    rho = 2.33	# Density of material (for Si)

    a = 0.1492	# Parameters (taken from W.R. Leo for SI)
    k = 3.25		
    X0 = 0.2014		
    X1 = 2.87		
    C = -4.44		
    # // double d0 = 0.0;		//
    X = np.log10(bg)
    # d;	// Variable for the Density effect
    
    if X>=X1:
        d = 2 * np.log(10.0) * X - C
    elif X0<=X and X<X1:
        d = 2 * np.log(10.0) * X - C + a * ((X1-X)**k)
    elif X<X0:
        # d = d0*(pow(10,(2*(X-X0))));
        d = 0

    # WM = 2 * me * ((beta*gamma)**2)/(1+(2*me*gamma/M)+(me/M)**2)   #Maximum energy tranfer
    WM = (2 * me * (beta*gamma)**2)/(1+(2*me/M)*np.sqrt(1 + (beta*gamma)**2) + (me/M)**2)   #Maximum energy tranfer

    loge = np.log((1-beta**2)*(I**2)/(2*M*(beta**2))) + (beta**2) # log epsilon variable

    EC = 0.577	# Euler's constant

    DeltaAv = K*rho*L*(z**2)*ZA*(1.0/(beta**2))*(np.log(2*me*(gamma**2)*(p**2)*WM/(I**2))-(2*beta**2)-(d)) # Mean energy loss (Bethe-Bloch)

    xi = (K)*rho*ZA*L*(z/beta)**2		# Xi variable 

    Lambda = (Delta-xi*(np.log(xi)-loge+1-EC))/xi # Lambda parameter

    Deltamp = xi*(np.log(xi/np.exp(loge))+0.198-d)		# Most probable energy loss
    lambdamp = (Deltamp-xi*(np.log(xi)-loge+1-EC))/xi 

    kappa = xi/WM		# Kappa ratio
    beta2 = beta**2
    sigma2 = (xi**2)*(1-beta2/2)/kappa		# Standard deviation for relativistic particles

    if kappa<=0.01:
        phi = TMath.Landau(Lambda, lambdamp, 1.0)
        return phi/xi

    elif 0.01<kappa and kappa<10:
        vav = TMath.Vavilov(Delta-Deltamp, kappa, beta2)
        return vav

    else:
        # gauss = exp(((Delta-DeltaAv)**2)/(2*sigma2));
        gauss = TMath.Gaus(Delta, DeltaAv, np.sqrt(sigma2))
        return gauss

def random_LV(s, p):
    f = TF1("f", LandV, 0.0, 10.0,2)
    # f = TF1("f", LandV, 0, 10, 2)
    f.SetParameter(0, s)
    f.SetParameter(1, p)

    Edep = f.GetRandom()
    return Edep * 1000


def muon_generator_1(Radio, long_a, long_b, number_thet, number_points_per_angle, Theta, 
                                  Theta_true, Phi, Energy):
    list_rand_thet = []
    list_rand_phi = []

    list_rand_thet_deg = []
    list_rand_phi_deg = []

    # list_delta_L = []
    list_random_energy = []

    for i in np.arange(0,number_thet):
        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        Random_th_deg = np.degrees(Random_th[0])
        Random_phi_deg = np.degrees(Random_phi)
        # print(Random_th[0])

        list_dis_Energy = []
        for energy in Energy:
            dis_Energy = dis_energy(energy, Random_th[0])
            list_dis_Energy.append(dis_Energy)

        # Vec = coord_cartesian(Random_th, Random_phi)
        # Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        # Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # # print('Vector sobre la esfera: ', Point)
        # # print('Norma del Vector:', norma)
        # # Points.append(Point)

        # # normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
        #                                                                                     ##  hacia el centro de coordenadas
        # # normal_Norma_Vec = norma_vec(normal_Vec)
        # # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # # print(len(normal_Vec))
        # # Vectors.append(normal_Vec)

        # vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        # vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        for i in np.arange(0,number_points_per_angle):
            # random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            # random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            # # list_random_a.append(random_a)
            # # list_random_b.append(random_b)

            # P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
            #             random_a * vec_thet[1] + random_b * vec_phi[1], 
            #             random_a * vec_thet[2] + random_b * vec_phi[2]]
            
            # list_P_vector.append(P_vector)

            # random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

            # print(random_plane_point)
            # list_random_point.append(random_plane_point)

            list_rand_thet.append(Random_th[0]) 
            list_rand_phi.append(Random_phi) ## En radianes

            list_rand_thet_deg.append(Random_th_deg)
            list_rand_phi_deg.append(Random_phi_deg) ## En grados

            Random_energy = rand.choices(Energy, list_dis_Energy)
            list_random_energy.append(Random_energy[0])


    dict_muons =  {'Theta(Rad)': list_rand_thet, 'Theta(Deg)': list_rand_thet_deg, 
                   'Phi(Rad)' : list_rand_phi, 'Phi(Deg)' : list_rand_phi_deg, 
                   'Energy(MeV)' : list_random_energy} 

    return dict_muons

def muon_generator(Energy, number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, 
                  long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z):
    
    list_rand_thet = []
    list_rand_phi = []
    
    list_rand_thet_deg = []
    list_rand_phi_deg = []

    list_P_vector = []
    list_random_point = []
    list_delta_L = []
    list_random_energy = []
    list_energy_Landau = []

    m_mu = 105.7

    n_muons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 1

    for i in np.arange(0,number_thet):
        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi en radianes
        # print(Random_th[0])
        Random_th_deg = np.degrees(Random_th[0]) ## El ángulo theta en grados
        Random_phi_deg = np.degrees(Random_phi) ## El ángulo phi en grados

        list_rand_thet.append(Random_th)
        list_rand_phi.append(Random_phi)
        list_rand_thet_deg.append(Random_th_deg)
        list_rand_phi_deg.append(Random_phi_deg)

        list_dis_Energy = []
        for energy in Energy:   ## Aquí se crea la distribución de Smith-Duller en MeV
            dis_Energy = dis_energy(energy, Random_th[0])
            list_dis_Energy.append(dis_Energy)
            
        Random_energy = rand.choices(Energy, list_dis_Energy) ## Escoje una energía segun la distribución de Smith-Duller en 
        list_random_energy.append(Random_energy[0])

        ### Momento del muon ###
        # momentum = np.sqrt(Random_energy[0]**2 - m_mu**2)
        momentum = Random_energy[0]

        os.environ["EN_SMITH"] = str(momentum)
        
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
                                                                                            ##  hacia el centro de coordenadas
        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
        # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
        # list_random_phi.append(Random_phi)


        random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
        random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

        # list_random_a.append(random_a)
        # list_random_b.append(random_b)

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]
        
        # list_P_vector.append(P_vector)

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
        # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

        # print(random_plane_point)
        # list_random_point.append(random_plane_point)

        #### Intersecciones con cada cara   ####
        
        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

        #### Cara Inferior ###
        t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3

        ### Cara 4 ###
        t_4 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4

        #### Caras en Y ###
        ### Cara 3 ###
        t_5 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5

        ### Cara 4 ###
        t_6 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_4 
        x_6 = random_plane_point[0] + normal_Vec[0] * t_4

        list_z = [z_3, z_4, z_5, z_6]
        
        if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
            flag_cara_1 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
            flag_cara_2 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
            flag_cara_3 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
            flag_cara_4 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
            flag_cara_5 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
            flag_cara_6 = True
            # print('Bandera 2: ', flag_cara_2)

        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
        
        Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        if Delta_L != 0:

            if Delta_L > 0 and Delta_L < 2.1:
                list_delta_L.append(Delta_L)
                # print('Estoy agregando el deltaL')

                n_muons_in_CCD = n_muons_in_CCD + muon

                ## Para la laptop en el ICN  ##
                #new_env = subprocess.run(["root", "-l", "-b", "/home/labdet/Documents/MauSan/Programas/Repositorio_Git/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"],
                #                     capture_output=True)

                ## Para la computadora de casa ##
                # new_env = subprocess.run(["root", "-l", "-b", "/home/bruce/Documents/Programas/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"], 
                #                             capture_output=True)

                ## Para el CLUSTER ##
                new_env = subprocess.run(["root", "-l", "-b", "/home/icn/mausanram/Software/CodigosICN/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"], 
                                             capture_output=True)

                # print('Energía de SMith-Duller: ', os.getenv("EN_SMITH"))
                # print(new_env.stdout)
                # print(new_env.stderr)

                # print(os.getenv('PATH'))
                # subprocess.run()
                # print(new_env.stdout)
                Random_energy_Landau = float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1])
                print(Random_energy_Landau)

                # print(float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1]))
                list_energy_Landau.append(Random_energy_Landau)
                
                muon_in_bucle += 1
                # print("El valor de EDEP", str(os.getenv('USERNAME')))

                print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet * number_points_per_angle), end = '\r')

            else:
                n_negative_long = n_negative_long + 1
                continue

        else:
                continue
        
        # print(os.environ)

        

    # dict_muons =  {'Random_Thet': list_rand_thet, 'Random_Phi' : list_rand_phi, 'Random_Energy' : list_random_energy, 'DeltaL' : list_delta_L} 

    dict_muons =  {'Theta(Rad)': list_rand_thet, 'Theta(Deg)': list_rand_thet_deg, 
                   'Phi(Rad)' : list_rand_phi, 'Phi(Deg)' : list_rand_phi_deg, 
                   'Energy-SD(MeV)' : list_random_energy, 'Energy_Landau' : list_energy_Landau} 

    return dict_muons, n_muons_in_CCD, n_negative_long

def muon_generator_3(Energy, number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, 
                  long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z):
    ### Primary lists ###
    list_rand_thet = []
    list_rand_phi = []
    list_random_energy = []
    
    list_rand_thet_deg = []
    list_rand_phi_deg = []

    ## Another lists ###
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []
    list_nmuons = []

    list_P_vector = []
    list_random_point = []
    list_delta_L = []
    list_energy_Landau = []

    m_mu = 105.7

    n_muons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    for i in np.arange(0,number_thet):
        n_muon = i
        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi en radianes
        # print(Random_th[0])
        Random_th_deg = np.degrees(Random_th[0]) ## El ángulo theta en grados
        Random_phi_deg = np.degrees(Random_phi) ## El ángulo phi en grados

        # list_rand_thet.append(Random_th[0])
        # list_rand_phi.append(Random_phi)
        # list_rand_thet_deg.append(Random_th_deg)
        # list_rand_phi_deg.append(Random_phi_deg)

        # list_dis_Energy = []
        # for energy in Energy:   ## Aquí se crea la distribución de Smith-Duller en MeV
        #     dis_Energy = dis_energy(energy, Random_th[0])
        #     list_dis_Energy.append(dis_Energy)

        # In = datetime.datetime.now()
        list_dis_Energy = dis_energy(Energy, Random_th[0], units=1)
        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para distribucion_En: ', Fin-In)
            
        Random_energy = rand.choices(Energy, list_dis_Energy)[0] * 1000 ## Escoje una energía segun la distribución de Smith-Duller en MeV
        # list_random_energy.append(Random_energy[0])

        list_rand_thet_deg.append(Random_th_deg)
        list_rand_phi_deg.append(Random_phi_deg)
        list_random_energy.append(Random_energy)

        ### Momento del muon ###
        # momentum = np.sqrt(Random_energy[0]**2 - m_mu**2)
        momentum = Random_energy

        # os.environ["EN_SMITH"] = str(momentum)
        

        # In = datetime.datetime.now()

        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        # normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)
        normal_Vec =  [-1 * np.sin(Random_th) * np.cos(Random_phi), -1 * np.sin(Random_th) * np.sin(Random_phi), -1 * np.cos(Random_th)]     ## Es un vector normal unitario apuntando 
                                                                                            ##  hacia el centro de coordenadas
        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), -np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
        # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
        # list_random_phi.append(Random_phi)


        random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
        random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

        # list_random_a.append(random_a)
        # list_random_b.append(random_b)

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]
        
        # list_P_vector.append(P_vector)

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
        # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

        # print(random_plane_point)
        # list_random_point.append(random_plane_point)

        #### Intersecciones con cada cara   ####
        
        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1[0] 
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1[0] 

        #### Cara Inferior ###
        t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2[0] 
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2[0]

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3[0] 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3[0]

        ### Cara 4 ###
        t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4[0] 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4[0]

        #### Caras en Y ###
        ### Cara 3 ###
        t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5[0] 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5[0]

        ### Cara 4 ###
        t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_4[0]
        x_6 = random_plane_point[0] + normal_Vec[0] * t_4[0]

        list_z = [z_3, z_4, z_5, z_6]
        
        if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
            flag_cara_1 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
            flag_cara_2 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
            flag_cara_3 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
            flag_cara_4 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
            flag_cara_5 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
            flag_cara_6 = True
            # print('Bandera 2: ', flag_cara_2)

        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
        
        Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para Delta L: ', Fin-In)

        if Delta_L > 0:

            # if Delta_L > 0 and Delta_L < 2.1:
                # list_delta_L.append(Delta_L)
                # print('Estoy agregando el deltaL')
                # print(Delta_L[0])
                n_muons_in_CCD = n_muons_in_CCD + 1

                Random_energy_Landau = random_LV(s = Delta_L, p = momentum) # En KeV

                list_nmuons.append(n_muon)
                list_energy_Landau.append(Random_energy_Landau)
                list_thet_in_CCD.append(Random_th[0])
                list_phi_in_CCD.append(Random_phi)
                list_energy_pri_in_CCD.append(Random_energy)
                list_delta_L.append(Delta_L[0])

                # list_rand_thet.append(Random_th[0])
                # list_rand_phi.append(Random_phi)
                
                muon_in_bucle += 1

                print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet * number_points_per_angle), end = '\r')

            # else:
            #     n_negative_long = n_negative_long + 1
            #     continue

        else:
            Random_energy_Landau = 0 # En KeV
            Delta_L = 0

            list_nmuons.append(n_muon)
            list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th[0])
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)

            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet * number_points_per_angle), end = '\r')
        
        # print(os.environ)

        
    Len_Landau = len(list_energy_Landau)
    Len_Theta = len(list_rand_thet)

    # if Len

    # dict_muons =  {'Random_Thet': list_rand_thet, 'Random_Phi' : list_rand_phi, 'Random_Energy' : list_random_energy, 'DeltaL' : list_delta_L} 

    # dict_all_muons =  {'Theta(Deg)': list_rand_thet_deg, 'Phi(Deg)' : list_rand_phi_deg, 'Energy-SD(MeV)' : list_random_energy} 

    dict_muons =  {'NMuon': list_nmuons, 'Theta(Rad)': list_thet_in_CCD, 'Phi(Rad)' : list_phi_in_CCD, 'Energy-SD(MeV)' : list_energy_pri_in_CCD, 
                        'Delta_L(cm)' : list_delta_L, 'Energy_Landau(KeV)' : list_energy_Landau} 

    return dict_muons, n_muons_in_CCD 

def muon_generator_CLUSTER(Energy, number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, 
                  long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z):
    
    list_rand_thet = []
    list_rand_phi = []
    list_rand_thet_deg = []
    list_rand_phi_deg = []
    list_random_energy = []

    ## Another lists ###
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []

    list_P_vector = []
    list_random_point = []
    list_delta_L = []
    list_energy_Landau = []

    m_mu = 105.7

    n_muons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    for i in np.arange(0,number_thet):
        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi en radianes
        # print(Random_th[0])
        Random_th_deg = np.degrees(Random_th[0]) ## El ángulo theta en grados
        Random_phi_deg = np.degrees(Random_phi) ## El ángulo phi en grados

        # list_rand_thet.append(Random_th[0])
        # list_rand_phi.append(Random_phi)
        # list_rand_thet_deg.append(Random_th_deg)
        # list_rand_phi_deg.append(Random_phi_deg)

        # list_dis_Energy = []
        # for energy in Energy:   ## Aquí se crea la distribución de Smith-Duller en MeV
        #     dis_Energy = dis_energy(energy, Random_th[0])
        #     list_dis_Energy.append(dis_Energy)
        
        list_dis_Energy = dis_energy(Energy, Random_th[0], units=1)

        Random_energy = rand.choices(Energy, list_dis_Energy)[0] * 1000 ## Escoje una energía segun la distribución de Smith-Duller en MeV
        # list_random_energy.append(Random_energy[0])

        list_rand_thet_deg.append(Random_th_deg)
        list_rand_phi_deg.append(Random_phi_deg)
        list_random_energy.append(Random_energy)

        ### Momento del muon ###
        # momentum = np.sqrt(Random_energy[0]**2 - m_mu**2)
        momentum = Random_energy

        # os.environ["EN_SMITH"] = str(momentum)
        
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
                                                                                            ##  hacia el centro de coordenadas
        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
        # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
        # list_random_phi.append(Random_phi)


        random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
        random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

        # list_random_a.append(random_a)
        # list_random_b.append(random_b)

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]
        
        # list_P_vector.append(P_vector)

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
        # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

        # print(random_plane_point)
        # list_random_point.append(random_plane_point)

        #### Intersecciones con cada cara   ####
        
        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

        #### Cara Inferior ###
        t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3

        ### Cara 4 ###
        t_4 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4

        #### Caras en Y ###
        ### Cara 3 ###
        t_5 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5

        ### Cara 4 ###
        t_6 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_4 
        x_6 = random_plane_point[0] + normal_Vec[0] * t_4

        list_z = [z_3, z_4, z_5, z_6]
        
        if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
            flag_cara_1 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
            flag_cara_2 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
            flag_cara_3 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
            flag_cara_4 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
            flag_cara_5 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
            flag_cara_6 = True
            # print('Bandera 2: ', flag_cara_2)

        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
        
        Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        if Delta_L != 0:

            if Delta_L > 0 and Delta_L < 2.1:
                # list_delta_L.append(Delta_L)
                # print('Estoy agregando el deltaL')

                n_muons_in_CCD = n_muons_in_CCD + muon


                # os.environ["DELTA_L"] = str(Delta_L)

                # print(os.environ["DELTA_L"])
                # print('Entrando al proceso de ROOT')
                
                # new_env = subprocess.run(["root", "-l", "-b", "/home/icn/mausanram/Software/CodigosICN/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"], 
                #                              capture_output=True)

                # print('Energía de SMith-Duller: ', os.getenv("EN_SMITH"))
                # print(new_env.stdout)
                # print(new_env.stderr)

                # print(os.getenv('PATH'))
                # subprocess.run()
                # print(new_env.stdout)
                # Random_energy_Landau = float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1]) # En KeV
                # print(Random_energy_Landau)

                # print(float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1]))

                Random_energy_Landau = random_LV(s = Delta_L, p = momentum) # En KeV


                list_rand_thet.append(Random_th[0])
                list_rand_phi.append(Random_phi)
                list_rand_thet_deg.append(Random_th_deg)
                list_rand_phi_deg.append(Random_phi_deg)
                list_random_energy.append(Random_energy)
                list_delta_L.append(Delta_L)
                list_energy_Landau.append(Random_energy_Landau)
                
                muon_in_bucle += 1
                # print("El valor de EDEP", str(os.getenv('USERNAME')))

                # print(muon_in_bucle, Random_th, Random_phi, Random_energy[0], Random_energy_Landau)

                print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet * number_points_per_angle), end = '\r')

            else:
                n_negative_long = n_negative_long + 1
                continue

        else:
                continue
        
        # print(os.environ)

        
    Len_Landau = len(list_energy_Landau)
    Len_Theta = len(list_rand_thet)

    # if Len

    # dict_muons =  {'Random_Thet': list_rand_thet, 'Random_Phi' : list_rand_phi, 'Random_Energy' : list_random_energy, 'DeltaL' : list_delta_L} 

    dict_muons =  {'Theta(Rad)': list_rand_thet, 'Theta(Deg)': list_rand_thet_deg, 
                   'Phi(Rad)' : list_rand_phi, 'Phi(Deg)' : list_rand_phi_deg, 
                   'Energy-SD(MeV)' : list_random_energy, 'Delta_L' : list_delta_L,
                   'Energy_Landau' : list_energy_Landau} 

    return dict_muons, n_muons_in_CCD, n_negative_long

def func_longitud(number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, 
                  long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z):
    
    list_P_vector = []
    list_random_point = []
    list_delta_L = []

    list_random_th = []
    list_random_phi = []

    n_muons_in_CCD = 0
    n_negative_long = 0

    for i in np.arange(0,number_thet):
        list_random_energy_per_plane = []
        list_points_per_plane = []

        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        # print('Random Thet: ', np.degrees(Random_th[0]))
        # print('Random Phi: ', np.degrees(Random_phi))
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))

        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        # normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma) 
        normal_Vec =  [-1 * np.sin(Random_th) * np.cos(Random_phi), -1 * np.sin(Random_th) * np.sin(Random_phi), -1 * np.cos(Random_th)]  ## Es un vector normal unitario apuntando hacia el centro de coordenadas
        # print('Direction Vec: ', normal_Vec)

        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), -np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        for i in np.arange(0,number_points_per_angle):
            flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
            # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
            # list_random_phi.append(Random_phi)


            random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            # list_random_a.append(random_a)
            # list_random_b.append(random_b)

            P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                        random_a * vec_thet[1] + random_b * vec_phi[1], 
                        random_a * vec_thet[2] + random_b * vec_phi[2]]
            
            list_P_vector.append(P_vector)

            random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

            # print(random_plane_point)
            list_random_point.append(random_plane_point)

            #### Intersecciones con cada cara   ####
            # print('Coord. del punto sobre el plano: (', random_plane_point[0], random_plane_point[1], random_plane_point[2],')')

            #### Cara Superior ###

            t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
            x_1 = random_plane_point[0] + normal_Vec[0] * t_1[0] 
            y_1 = random_plane_point[1] + normal_Vec[1] * t_1[0] 

            #### Cara Inferior ###
            t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
            x_2 = random_plane_point[0] + normal_Vec[0] * t_2[0] 
            y_2 = random_plane_point[1] + normal_Vec[1] * t_2[0]

            ### Caras en X ###
            ### Cara 3 ###
            t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0]
            z_3 = random_plane_point[2] + normal_Vec[2] * t_3[0] 
            y_3 = random_plane_point[1] + normal_Vec[1] * t_3[0]

            ### Cara 4 ###
            t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
            z_4 = random_plane_point[2] + normal_Vec[2] * t_4[0] 
            y_4 = random_plane_point[1] + normal_Vec[1] * t_4[0]

            #### Caras en Y ###
            ### Cara 5 ###
            t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1]
            z_5 = random_plane_point[2] + normal_Vec[2] * t_5[0] 
            x_5 = random_plane_point[0] + normal_Vec[0] * t_5[0]

            ### Cara 4 ###
            t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
            z_6 = random_plane_point[2] + normal_Vec[2] * t_6[0] 
            x_6 = random_plane_point[0] + normal_Vec[0] * t_6[0]

            # print(random_plane_point[0][0])
            # x, y, z = np.indices((10, 10, 10))

            # draw cuboids in the top left and bottom right corners, and a link between
            # them
            # cube1 = (x > 0) & (y > 0) & (z > 0)
            # cube2 = (x >= 5) & (y >= 5) & (z >= 5)
            # link = abs(x - y) + abs(y - z) + abs(z - x) <= 2

            # combine the objects into a single boolean array
            # voxelarray = cube1 | cube2 | link
            # voxelarray = cube1 | cube2 
            # voxelarray = cube1

            # set the colors of each object
            # colors = np.empty(voxelarray.shape, dtype=object)
            # colors[link] = 'red'
            # colors[cube1] = 'blue'
            # colors[cube2] = 'green'

            # and plot everything
            # ax = plt.figure().add_subplot(projection='3d')
            # # ax.voxels(filled = voxelarray, facecolors=colors)
            # VecStart_x = [-50, -50, -50, -50]
            # VecStart_y = [-5, -5, 5, 5]
            # VecStart_z = [0, 10, 0, 10]
            # VecEnd_x = [50, 50, 50, 50]
            # VecEnd_y = [-5, -5, 5, 5]
            # VecEnd_z  =[0, 10, 0, 10]

            # for i in range(0, len(VecStart_x)):
            #     ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'b')

            # VecStart_x = [-50, -50, -50, -50]
            # VecStart_y = [-5, -5, 5, 5]
            # VecStart_z = [0, 0, 0, 10]
            # VecEnd_x = [-50, -50, -50, -50]
            # VecEnd_y = [-5, 5, 5, -5]
            # VecEnd_z  =[10, 0, 10, 10]
            # for i in range(0, len(VecStart_x)):
            #     ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'b')

            # # t = 100
            # # t = (medida_z - random_plane_point[2][0]) / normal_Vec[2][0]
            # t = (-10 - random_plane_point[2][0]) / normal_Vec[2][0]
            # VecStart_x = [random_plane_point[0][0]]
            # VecStart_y = [random_plane_point[1][0]]
            # VecStart_z = [random_plane_point[2][0]]
            # VecEnd_x = [random_plane_point[0][0] + normal_Vec[0][0] * t]
            # VecEnd_y = [random_plane_point[1][0] + normal_Vec[1][0] * t]
            # VecEnd_z  =[random_plane_point[2][0] + normal_Vec[2][0] * t]

            # ax.scatter(VecStart_x[0], VecStart_y[0], VecStart_z[0], color = 'k')

            # for i in range(1):
            #     ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'k')
            
            # ax.scatter(0,0,0, color = 'k')
            # ax.scatter(x_1[0], y_1[0], medida_z, color = 'r')
            # ax.scatter(x_2[0], y_2[0], 0, color = 'r')
            # ax.scatter(medida_x, y_3[0], z_3[0], color = 'g')
            # ax.scatter(-medida_x, y_4[0], z_4[0], color = 'g')
            # ax.scatter(x_5[0], medida_y, z_5[0], color = 'purple')
            # # ax.scatter(random_plane_point[0], random_plane_point[1], random_plane_point[2], 'k')
            # # ax.plot([random_plane_point[0], 0], [random_plane_point[1],0], [random_plane_point[2],0])
            # # ax.scatter(x_2, y_2,0, 'k')

            # plt.show()


            list_z = [z_3, z_4, z_5, z_6]
            # print('Valores de t: ', t_1[0], t_2[0], t_3[0], t_4[0], t_5[0], t_6[0])
            # print('Valores de x: ', x_1[0], x_2[0], x_5[0], x_6[0])
            # print('Valores de y: ', y_1[0], y_2[0], y_3[0], y_4[0])
            # print('Valores de z: ', z_3[0], z_4[0], z_5[0], z_6[0])
            
            if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
                flag_cara_1 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
                flag_cara_2 = True
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
                flag_cara_3 = True
                # print('Bandera 3: ', flag_cara_3) 

            if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
                flag_cara_4 = True
                # print('Bandera 4: ', flag_cara_4)

            if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
                flag_cara_5 = True
                # print('Bandera 5: ', flag_cara_5) 

            if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
                flag_cara_6 = True
                # print('Bandera 6: ', flag_cara_6)

            list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
            
            Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

            if Delta_L > 0:
                # if Delta_L > 0 and Delta_L < 2.1:
                    list_delta_L.append(Delta_L)
                    list_random_th.append(np.degrees(Random_th[0]))
                    list_random_phi.append(np.degrees(Random_phi))

                    n_muons_in_CCD = n_muons_in_CCD + muon

                # else:
                #     n_negative_long = n_negative_long + 1
                #     continue

            else:
                    continue
            
    return list_delta_L, n_muons_in_CCD, n_negative_long, list_random_th, list_random_phi

def deffuse_model(number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, long_a, long_b, 
                        medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z, dict_diffuse_model):
    
    list_P_vector = []
    list_random_point = []
    list_delta_L = []
    n_muons_in_CCD = 0
    n_negative_long = 0

    for i in np.arange(0,number_thet):
        list_random_energy_per_plane = []
        list_points_per_plane = []

        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
                                                                                            ##  hacia el centro de coordenadas
        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        for i in np.arange(0,number_points_per_angle):
            flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
            # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
            # list_random_phi.append(Random_phi)


            random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            # list_random_a.append(random_a)
            # list_random_b.append(random_b)

            P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                        random_a * vec_thet[1] + random_b * vec_phi[1], 
                        random_a * vec_thet[2] + random_b * vec_phi[2]]
            
            list_P_vector.append(P_vector)

            random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

            # print(random_plane_point)
            list_random_point.append(random_plane_point)

            #### Intersecciones con cada cara   ####


            #### Cara Superior ###
            t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
            x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
            y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

            #### Cara Inferior ###
            t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
            x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
            y_2 = random_plane_point[1] + normal_Vec[1] * t_2

            ### Caras en X ###
            ### Cara 3 ###
            t_3 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
            z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
            y_3 = random_plane_point[1] + normal_Vec[1] * t_3

            ### Cara 4 ###
            t_4 = (medida_x - random_plane_point[0]) / normal_Vec[0]
            z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
            y_4 = random_plane_point[1] + normal_Vec[1] * t_4

            #### Caras en Y ###
            ### Cara 3 ###
            t_5 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
            z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
            x_5 = random_plane_point[0] + normal_Vec[0] * t_5

            ### Cara 4 ###
            t_6 = (medida_y - random_plane_point[1]) / normal_Vec[1]
            z_6 = random_plane_point[2] + normal_Vec[2] * t_4 
            x_6 = random_plane_point[0] + normal_Vec[0] * t_4

            list_z = [z_3, z_4, z_5, z_6]
            
            if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
                flag_cara_1 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
                flag_cara_2 = True
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
                flag_cara_3 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
                flag_cara_4 = True
                # print('Bandera 2: ', flag_cara_2)

            if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
                flag_cara_5 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
                flag_cara_6 = True
                # print('Bandera 2: ', flag_cara_2)

            list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
            
            Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

            if Delta_L != 0:
                if Delta_L > 0 and Delta_L < 2.1:
                    list_delta_L.append(Delta_L) 

                    ### -------------------- Modelo de Difusión ----------------- ##

                    # CCD_deep = 725 # micras
                    # Alpha = dict_diffuse_model['Alpha']
                    # Beta = dict_diffuse_model['Beta']

                    # z_deeps = np.arange()





                    n_muons_in_CCD = n_muons_in_CCD + muon

                else:
                    n_negative_long = n_negative_long + 1
                    continue

            else:
                    continue
            
    return list_delta_L, n_muons_in_CCD, n_negative_long


