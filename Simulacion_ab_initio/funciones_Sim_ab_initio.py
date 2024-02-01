import numpy as np
import mpmath as mp
import random as rand
import datetime
import os

def dis_probability(theta, I_0):
    return I_0 * np.cos(theta)

def dis_angular(theta): ## Distribucion angular
    return 1 * np.cos(theta)**2 * np.sin(theta)

def dis_energy(E_mu, theta): ### Pone la distribución de energías (en MeV) en función del angulo theta
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

    # a = 0.0025         ## GeV cm^2/g
    # m_mu = 0.1057    ## GeV/c^2 
    # m_pi = 0.1396    ## GeV/c^2

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

def muon_generator(E, Theta, Theta_true, Phi, Radio, long_a, long_b, n_thet, n_points):
    list_random_th = []
    list_random_phi = []
    list_random_a = []
    list_random_b = []
    list_random_point = []
    list_random_energy = []
    list_delta_L = []

    Vectors = []
    Points = []
    Inicio = datetime.datetime.now()
    print('Hora de inicio del cálculo: ', Inicio)

    for i in np.arange(0,n_thet):
        list_points_per_plane = []
        list_random_energy_per_plane = []

        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true

        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        Energy_dist = dis_energy(E, Random_th[0])

        Vec = coord_cartesian(Random_th, Random_phi)
        # print(type(Vec[0]))
        Point = (Radio * Vec[0], Radio * Vec[1], Radio * Vec[2])  ## Genera un punto sobre la esfera.
        Points.append(Point)

        normal_Vec = (-1 * Vec[0], -1 * Vec[1], -1 * Vec[2])     ## Es un vector apuntando hacia el centro de coordenadas
        # print(len(normal_Vec))
        Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]

        for i in np.arange(0,n_points):
            list_random_th.append(Random_th[0])    ## Lo anexa en una lista
            list_random_phi.append(Random_phi)


            random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            list_random_a.append(random_a)
            list_random_b.append(random_b)

            delta_L = 0.0725 / np.cos(Random_th)  ## cm
            list_delta_L.append(delta_L)

            P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                        random_a * vec_thet[1] + random_b * vec_phi[1], 
                        random_a * vec_thet[2] + random_b * vec_phi[2]]

            random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            list_random_point.append(random_plane_point)

            Random_energy = rand.choices(E, Energy_dist)
            list_random_energy.append(Random_energy[0])

        # list_random_energy_per_plane.append(list_random_energy)
        # list_points_per_plane.append(list_random_plane_point)

    # dict_simulation = {'Theta_Radianes': list_random_th, 'Phi_radianes': list_random_phi, 'Points_per_plane' : list_planes,
    #                    'Energy_per_muon' : list_random_energy }

    dict_simulation = {'Number_of_angles' : n_thet, 'Points_per_angle' : n_points,  'Theta_Radianes': list_random_th, 
                       'Phi_Radianes': list_random_phi, 'list_random_a' : list_random_a, 'list_random_b' : list_random_b, 
                       'Points' : list_random_point, 'Energy_per_muon' : list_random_energy, 'list_Delta_L' : list_delta_L }
    
    # del list_random_th, list_random_phi, list_points_per_plane, list_random_energy, list_random_energy_per_plane, Vectors, Points

    Final = datetime.datetime.now()

    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    return dict_simulation
