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

def norma_vec(Vector):
    norma = np.sqrt(Vector[0] ** 2 + Vector[1] ** 2  + Vector[2] ** 2 )
    return norma

def dimension_x():
    long_x = 1.197 / 2
    step = 0.0001

    list_long_x = [-long_x]

    while long_x:
        x = np.round(list_long_x[-1] + step, 4)
        list_long_x.append(x)

        if list_long_x[-1] == long_x:
            break

    return list_long_x

def dimension_y():
    long_y = 1.587 / 2
    step = 0.0001

    list_long_y = [-long_y]

    while long_y:
        y = np.round(list_long_y[-1] + step, 4)
        list_long_y.append(y)

        if list_long_y[-1] == long_y:
            break

    return list_long_y

def dimension_z():
    long_z = 0.0725
    step = 0.0001

    list_long_z = [-long_z]

    while long_z:
        z = np.round(list_long_z[-1] + step, 4)
        list_long_z.append(z)

        if list_long_z[-1] == long_z:
            break

    return list_long_z

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
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 5
    if flags_CCD[0] and flags_CCD[4]:
        h = medida_z - list_z[2]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 6
    if flags_CCD[0] and flags_CCD[5]:
        h = medida_z - list_z[3]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 2
    if flags_CCD[2] and flags_CCD[1]:
        h = list_z[0]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 4
    if flags_CCD[2] and flags_CCD[3]:
        h = list_z[0]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 5
    if flags_CCD[2] and flags_CCD[4]:
        h = list_z[0]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 6
    if flags_CCD[2] and flags_CCD[5]:
        h = list_z[0]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 2
    if flags_CCD[3] and flags_CCD[1]:
        h = list_z[1]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 3
    if flags_CCD[3] and flags_CCD[2]:
        h = list_z[1]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 5
    if flags_CCD[3] and flags_CCD[4]:
        h = list_z[1]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 6
    if flags_CCD[3] and flags_CCD[5]:
        h = list_z[1]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 2
    if flags_CCD[4] and flags_CCD[1]:
        h = list_z[2]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 3
    if flags_CCD[4] and flags_CCD[2]:
        h = list_z[2]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 4
    if flags_CCD[4] and flags_CCD[3]:
        h = list_z[2]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 6
    if flags_CCD[4] and flags_CCD[5]:
        h = list_z[2]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 2
    if flags_CCD[5] and flags_CCD[1]:
        h = list_z[3]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 3
    if flags_CCD[5] and flags_CCD[2]:
        h = list_z[3]
        delta_L = h /  np.cos(Random_th)  ## cm
        n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue
    ## Caras 6 y 4
    if flags_CCD[5] and flags_CCD[3]:
        h = medida_z - list_z[3]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 5
    if flags_CCD[5] and flags_CCD[4]:
        h = list_z[3]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    # if (flags_CCD[0] == False) and  (flags_CCD[1] == False) and (flags_CCD[2] == False) and (flags_CCD[3] == False) and (flags_CCD[4] == False) and (flags_CCD[5] == False):
    if delta_L is None:    
        delta_L, n_muons_in_CCD = 0, 0

    return delta_L, n_muons_in_CCD


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
