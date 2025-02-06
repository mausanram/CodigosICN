print('Se están importando todas las paqueterías')
import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from funciones_Sim_ab_initio import *
import datetime
import os
# import pickle as pkl
# import ROOT 

from ROOT import TFile, TTree
from array import array

def main():
    current_path = os.getcwd()
    Inicio = datetime.datetime.now()
    print('Hora de inicio del cálculo: ', Inicio)

    ### Número de muones a simular ### 
    number_thet = 10000    ## Valores de un ángulo Theta.

    print('Se simularán ' + str(number_thet) + ' muones.')

    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons = muon_generator_1(number_thet)


    Thet_Rad = array('f', [-9999])
    Phi_Rad = array('f', [-9999])
    Energy_array = array('f', [-9999])

    file_root_name = 'MuonGen_NMUONS_' + str(number_thet) + '_.root'
    
    file = TFile.Open(file_root_name, "RECREATE")
    tree = TTree('tree', 'tree')

    tree.Branch('thet', Thet_Rad, 'thet/F')
    tree.Branch('phi', Phi_Rad, 'phi/F')
    tree.Branch('epri', Energy_array, 'epri/F')

    for i in np.arange(0, len(dict_muons['Theta(Rad)'])):
        Thet_Rad[0] = dict_muons['Theta(Rad)'][i]
        Phi_Rad[0] = dict_muons['Phi(Rad)'][i]
        Energy_array[0] =  dict_muons['Energy-SD(MeV)'][i] 
        tree.Fill()

    tree.Write()
    file.Close()

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    print('Se guardó la información de los muones simulados en el archivo ' + file_root_name)

if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)


