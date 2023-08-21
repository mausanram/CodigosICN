#######################################################################################
###                                                                                 ###
###                                                                                 ###
###          ARCHIVO DE FUNCIONES AUXILIARES UTILIZADAS POR EL BOT PARA             ###
###             MONITOREO DEL  LABORATORIO DE DETECTORES DEL ICN                    ###
###                                                                                 ###
###          * ULTIMA ACTUALIZACIÓN: 18 DE AGOSTO DEL 2023                          ###  
###                                                                                 ###                            
####################################################################################### 

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
