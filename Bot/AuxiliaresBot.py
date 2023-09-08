#######################################################################################
###                                                                                 ###
###                                                                                 ###
###          ARCHIVO DE FUNCIONES AUXILIARES UTILIZADAS POR EL BOT PARA             ###
###             MONITOREO DEL  LABORATORIO DE DETECTORES DEL ICN                    ###
###                                                                                 ###
###          * ULTIMA ACTUALIZACIÓN: 04 DE SEPTIEMBRE DEL 2023                      ###  
###                                                                                 ###                            
####################################################################################### 

import pandas as pd
global path_userinfo_csv

def ReadTemp(historial_Temp):
    with open(historial_Temp, 'r') as fs:
        lines = fs.readlines()
        last_line = lines[-1]
        fs.close()

    last_line_list = last_line.split(' ')
    clean_line_list = []

    for line in last_line_list:
        if line.startswith('\n'):
            continue    
        if len(line) == 0:
            continue
        clean_line_list.append(line)

    return clean_line_list

def dictConfigFile_335(configFileName):
    Config335_dict = {}
    with open(configFileName, 'r') as fl:
        lines = fl.readlines()
        fl.close()

    for line in lines:
        if line.startswith("#"):
            continue
        if line.startswith('\n'):
            continue
        list = line.split(':')
        Config335_dict.setdefault(list[0], list[1].strip())
        
    return Config335_dict

def Users_DataFrame(path_userinfo_csv):
    User_dataframe = pd.read_csv(path_userinfo_csv)
    return User_dataframe

def AddUser_to_csv(path_userinfo_csv, DataFrame, user_id, number_of_jobs):
    user_ID_Flag = True
    for user in list(DataFrame['User_ID']):
        if int(user) == user_id:
            #print('El usuario ya está registrado')
            user_ID_Flag = False
        
    if user_ID_Flag:
        if number_of_jobs == 1:
            DataFrame.loc[len(DataFrame.index)] = [str(user_id), False]
            DataFrame.to_csv(path_or_buf = path_userinfo_csv, index_label=False)
            DataFrame = pd.read_csv(path_userinfo_csv)
        
        elif number_of_jobs == 2:
            DataFrame.loc[len(DataFrame.index)] = [str(user_id), False, True]
            DataFrame.to_csv(path_or_buf = path_userinfo_csv, index_label=False)
            DataFrame = pd.read_csv(path_userinfo_csv)
    
    else: 
        DataFrame = DataFrame

    return DataFrame

def UpdateValue_to_csv(path_userinfo_csv, user_id, DataFrame, option = 0):
    for index_user in range(0,len(DataFrame['User_ID'])):
        print('Indice de usuario: ', index_user)
        print(DataFrame['User_ID'][index_user])
        print(user_id)
        if str(DataFrame['User_ID'][index_user]) == str(user_id):
            Index_user = index_user

            if option == 0:
                DataFrame.loc[Index_user,'Temp_Alarm'] = not DataFrame['Temp_Alarm'][Index_user]
                DataFrame.to_csv(path_or_buf = path_userinfo_csv, index_label=False)
                DataFrame = pd.read_csv(path_userinfo_csv)

        else:
            continue

    return DataFrame
    
    ## Agregar mas valores para mas alertas o tareas.

def AddJob_to_csv(path_userinfo_csv, DataFrame, job_name, job_value_default):
    DataFrame[str(job_name)] = job_value_default
    DataFrame.to_csv(path_or_buf = path_userinfo_csv, index_label=False)
    DataFrame = pd.read_csv(path_userinfo_csv)
    return DataFrame

def Search_User(DataFrame, user_id):
    for index_user in range(0,len(DataFrame['User_ID'])):
        if DataFrame['User_ID'][index_user] == user_id: return index_user