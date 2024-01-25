import numpy as np
import mpmath as mp

def dis_probability(theta, I_0):
    return I_0 * np.cos(theta)

def dis_angular(theta): ## Distribucion angular
    return 1 * np.cos(theta)**2 * np.sin(theta)

def dis_energy(E_mu, theta): ### Poner la distribución de energías en función del angulo theta
    ## Constantes físicas
    k = 8 / 3
    b = 0.771
    lambda_pi = 120     ## g/cm^2
    y_0 = 1000      ## g/cm^2
    r = 0.76
    a = 2.5         ## MeV cm^2/g
    b_mu = 0.8      
    m_mu = 105.7    ## MeV/c^2 
    m_pi = 139.6    ## MeV/c^2
    tau_mu_0 = 2.2 * 10**(-6)   ## s
    tau_0 = 2.6 * 10 **(-8)     ## s
    rho_0 = 0.00129 ## g/cm^3
    c = 3 * 10 ** 10 ## cm/s

    ### Parámetros
    E_pi = (1 / r) * (E_mu + a * y_0 * (mp.sec(theta) - 0.1))
    B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c)
    P_mu = ((0.1 * np.cos(theta)) * (1 - (a * (y_0 * mp.sec(theta) - 100))/( r * E_pi)) ) ** ((B_mu)/((r * E_pi + 100 * a) * np.cos(theta)))
    j_pi = (m_pi * y_0)/(c * tau_0 * rho_0)


    ## Intensidad diferencial
    C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    C_2 = E_pi * np.cos(theta)
    C_3 = b * j_pi

    return (C_1 * np.sin(theta)) / (C_2 + C_3)

##### Coordenadas Cartesianas Unitarias######
def coord_cartesian(Thet, Phi):
    coord_X = np.sin(Thet) * np.cos(Phi)
    coord_Y = np.sin(Thet) * np.sin(Phi)
    coord_Z = np.cos(Thet)

    Vec_sph = (float(coord_X), float(coord_Y), float(coord_Z))
    return Vec_sph 
