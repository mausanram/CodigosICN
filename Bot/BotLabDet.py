#######################################################################################
###                                                                                 
###                                                                                 
###      BOT PARA MONITOREO DE PARÁMETROS DEL LABORATORIO DE DETECTORES DEL ICN     
###                                                                                 
###      * ULTIMA ACTUALIZACIÓN: 08 DE SEPTIEMBRE DEL 2023                          
###                                                                                 
###                                                                                 
#######################################################################################     

import logging, os, signal
from telegram import __version__ as TG_VER
from telegram import Update, ReplyKeyboardMarkup
from telegram.ext import Application, CommandHandler, ContextTypes, MessageHandler, filters
from AuxiliaresBot import ReadTemp, dictConfigFile_335, Users_DataFrame, AddJob_to_csv, AddUser_to_csv, UpdateValue_to_csv, Search_User
import asyncio
    
try:
    from telegram import __version_info__
except ImportError:
    __version_info__ = (0, 0, 0, 0, 0)  # type: ignore[assignment]

if __version_info__ < (20, 0, 0, "alpha", 1):
    raise RuntimeError(
        f"This example is not compatible with your current PTB version {TG_VER}. To view the "
        f"{TG_VER} version of this example, "
        f"visit https://docs.python-telegram-bot.org/en/v{TG_VER}/examples.html")

logging.basicConfig(format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s", level = logging.INFO)
logger = logging.getLogger(__name__)
Token = " " #Token del Bot. Nombre de Bot en Telegram: PruebainBot
intervalo_de_rutinaTemp = 45 # Tiempo para repetir la alarma de temperatura
toloreancia_Temp = 1.5 # K 
restored_conection = True # Variable auxiliar para reestablecer las tareas que se hayan quedado activas.

# # Computadora de Laboratorio # 
# historyTemp_path = '/home/detectores/Software/LakeshoreControl/history_2023Aug11.txt' # Dirección del historial de temperatuas
# ConfigFile335_path = '/home/detectores/Software/LakeshoreControl/ConfigFile_M335' # Dirección del archivo de configuración del 335
# UserInfo_path = '/home/detectores/Software/SSocial/MauSan/CodigosICN/Bot/UsersBot.csv' # Dirección del archivo de información de tareas de todos los usuarios.

# # Computadora de Casa # 
# historyTemp_path = '/home/bruce/Documents/Programas_Pruebas/Bot/history_2023Sep04.txt' # Dirección del historial de temperatuas
# ConfigFile335_path = '/home/bruce/Documents/Programas/Bot/ConfigFile335.txt' # Dirección del archivo de configuración del 335
# UserInfo_path = '/home/bruce/Documents/Programas/Bot/UsersBot.csv'

# # Computadora de Oficina
historyTemp_path = '/home/labdet/Documents/MauSan/Programas/Bot/history_2023Sep04.txt' # Dirección del historial de temperatuas
ConfigFile335_path = '/home/labdet/Documents/MauSan/Programas/Repositorio_Git/Bot/ConfigFile335.txt' # Dirección del archivo de configuración del 335
UserInfo_path = '/home/labdet/Documents/MauSan/Programas/Repositorio_Git/Bot/UsersBot.csv'


try:
    df_Users = Users_DataFrame(UserInfo_path)
    df_Users_Flag = True

except: 
    df_Users_Flag = False

async def start(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None: ## Manda un mensaje de inicio y despliega en pantalla los comandos que se pueden utilizar
    """Send a message when the command /start is issued."""
    # logger.info('Estoy en la función start')
    global df_Users
    global df_Users_Flag

    user = update.effective_user
    chat_id = update.message.chat_id
    
    if df_Users_Flag:
        reply_keyboard = [["/Start"],["/WatchTemp", "/StartAlarm"], ["/Help", "/OffAlarm"] ]

        text = "¡Hola, "+user.first_name+"! 👋👋👋\n\nSoy el bot 🤖 del Lab. de Detectores del ICN.\
        \n\nTe puedo dar las últimas mediciones de los sensores de temperatura"

        await Restored_AlarmsStates(context)
        df_Users = AddUser_to_csv(UserInfo_path, df_Users, chat_id, 1) ## Añade al nuevo usuario al archivo de usuario
        await update.message.reply_text(text, reply_markup=ReplyKeyboardMarkup(
            reply_keyboard, one_time_keyboard = True, input_field_placeholder = "What do you want?"))

    else:
        text = "¡Hola, "+user.first_name+"! 👋👋👋\n\nSoy el bot 🤖 del Lab. de Detectores del ICN.\
        \n\nTe puedo dar las últimas mediciones de los sensores de temperatura"
        await update.message.reply_text(text)

        text= 'El archivo de las tareas de usuarios NO se puede leer. \nNo se puede monitorear las alarmas de los usuarios así que no podrás activar ninguna.'

        reply_keyboard = [["/Start"],["/WatchTemp"], ["/Help"] ]
        await update.message.reply_text(text, reply_markup=ReplyKeyboardMarkup(
                    reply_keyboard, one_time_keyboard = True, input_field_placeholder = "What do you want?"))

    # logger.info('This is the dictionary:')

async def help_command(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None: ## Despliega todos los comandos que se pueden utilizar y su descripción.
    """Send a message when the command /help is issued."""
    logger.info('Estoy en la función help_command')
    text = "Los comandos válidos son los siguientes: \
    \n\n/Start - Inicia el bot. \
    \n\n/Temperatura, /Temp, /WatchTemp - Regresa el último valor de temperatura de cada sensor.  \
    \n\n/Ayuda ,/Help - Regresa la lista de los comandos y su descripción. \
    \n\n/Detener, /Kill - Detiene el bot completamente, si se utiliza deberá volver a correr el bot en la terminal del laboatorio. \
    \n\n/Activar_Alarma, /StartAlarm - Activa una rutina que te avisa cuando la temperatura de la CCD cambia demasiado. \
    \n\n/Desactivar_Alarma, /OffAlarm - Desactiva la alarma de temperatura."
    await update.message.reply_text(text)

async def echo(update: Update, context: ContextTypes.DEFAULT_TYPE) -> None: ## Identifica las palabras que NO son un comando. 
    """Echo the user message."""
    await update.message.reply_text('Creo que "' + update.message.text + '" no es un comando, quieres intentar utilizando el menú?')

async def stop(update: Update, context: ContextTypes.DEFAULT_TYPE): ## Detiene el bot completamente. 
    logger.info('He recibido un comando stop')
    name = update.effective_chat.first_name
    text = "¡Hasta pronto, " + name + "! 👋👋👋"
    
    await update.message.reply_text(text)
    os.kill(os.getpid(), signal.SIGINT)

async def Restored_AlarmsStates(context: ContextTypes.DEFAULT_TYPE): ## Reestablece todas las tareas despues de una interrupción en el bot. 
    global df_Users
    global df_Users_Flag
    global restored_conection

    begin_time = 1

    for index_user in range(0,len(df_Users['User_ID'])):
        user_id = df_Users['User_ID'][index_user]

        if df_Users['Temp_Alarm'][index_user]  and restored_conection:
            context.job_queue.run_repeating(Temp_alarm, interval = intervalo_de_rutinaTemp, first = begin_time, chat_id = user_id, name = str(user_id) + '_ReadTemp_job')
            logger.info('Restableciendo las alarmas')
            # print(df_Users)
        elif df_Users['Temp_Alarm'][index_user] == False and restored_conection == True: continue

        begin_time = begin_time + 0.5
    restored_conection = False
    
async def ReadTemperature(update: Update, context: ContextTypes.DEFAULT_TYPE): ## Obtiene la última temperatura del historial y la manda al chat.  
    try: 
        tempList = ReadTemp(historyTemp_path)
        # logger.info('Si pude leer el archivo. ')
        # print(tempList)

        measureTime = tempList[0]
        Date = tempList[1][0:9]
        Temp = tempList[3][0:6]
        heater = tempList[-2]
        power = tempList[-1]
        
        text =  "Temperatura: "+ Temp +" K"+ \
                "\nHeater Power: "+ heater +" "+ power + \
                "\n\nFecha y hora de última actualización: \n" + Date +", "+ measureTime 

        await update.message.reply_text(text)

    except:
        Text = 'Archivo de historial de temperatura NO detectado. Por favor verifica la dirección del archivo y reinicia el bot. \nDirección actual: '+ historyTemp_path
        await update.message.reply_text(Text)

        Text = 'El bot se apagará automáticamente ya que no se puede obtener información del laboratorio.'
        await update.message.reply_text(Text)
        os.kill(os.getpid(), signal.SIGINT)

def Refresh_File(ReferenceValue, NewValue): ## Función que se encarga de checar si el historial de temperatura se está actualizando. 
    if ReferenceValue[0] == NewValue[0] and ReferenceValue[1] == NewValue[1]:
        return False
    else:
        return True

async def Temp_alarm(context: ContextTypes.DEFAULT_TYPE): ## Alarma para cambio de temperatura. 
    # logger.info('Estoy en alarma')
    global df_Users
    global df_Users_Flag

    job = context.job

    # print(job)
    # user_id = job.chat_id 

    #print('Class:', type(context.job.name.capitalize()))
    # print('Este es el ID del usuario: ', user_id)
    # print('User_ID:', user_id)
    list_job_name = context.job.name.title().split('_')

    User_id = list_job_name[0]

    try:
        # logger.info('Estoy en el try del historial')
        medsRef = ReadTemp(historyTemp_path)
        historyTemp_path_Flag = True
        
    except: 
        # logger.info('Estoy en el except del historial')
        historyTemp_path_Flag = False
        context.job_queue.get_jobs_by_name(str(User_id)+ '_ReadTemp_job')[0].schedule_removal()
        df_Users = UpdateValue_to_csv(UserInfo_path, User_id, df_Users, 0)

        Text = 'Archivo de historial de temperatura NO detectado. Por favor verifica la dirección del archivo.\nDirección actual: '+ historyTemp_path
        await context.bot.send_message(chat_id = User_id,text = Text)
    
    if historyTemp_path_Flag:
        try:
            # logger.info('Estoy en el try del ConfigFile335')
            Config335_dict = dictConfigFile_335(ConfigFile335_path)
            Config335_dict_Flag = True

        except: 
            # logger.info('Estoy en el except del ConfigFile335') 
            Config335_dict_Flag = False
            df_Users = UpdateValue_to_csv(UserInfo_path, User_id, df_Users, 0)
            Text = 'Archivo de configuración del 335 NO detectado. Por favor verifica la dirección del archivo.\nDirección actual: ' + ConfigFile335_path
            await context.bot.send_message(chat_id = User_id,text = Text)
        
        if Config335_dict_Flag:
            # logger.info('Estoy esperando')
            await asyncio.sleep(31)
            # time.sleep(31)
            medNew = ReadTemp(historyTemp_path)
            # logger.info(medsRef)
            # logger.info(medNew)
            # logger.info('Estoy en Refresh')
            Refresh_File_Flag = Refresh_File(medsRef, medNew)

        

            if Refresh_File_Flag:
                # logger.info('Se levantó la bandera Refresh')

                job = context.job
                SensorTemperature_Control = float(medNew[3][1:6])
                SetPointValue = float(Config335_dict['Set_Point Channel 1'])
                date = medNew[1][0:9]
                measureTime = medNew[0]

                difVal = abs(SetPointValue - SensorTemperature_Control)
                # logger.info('Ya obtuve las nuevas mediciones')
                # logger.info('medNew: '+str(SensorTemperature_Control))
                # logger.info('SetPoint: '+str(SetPointValue))
                # logger.info('difVal: '+str(difVal))

                if difVal >= 1.5:
                    text = "‼️ ‼️ ‼️ PRECAUCIÓN ‼️ ‼️ ‼️" 
                    await context.bot.send_message(chat_id = User_id,text = text)

                    text = "La temperatura se está alejando del valor establecido. \nTemperatura actual: " + str(SensorTemperature_Control) + ' K'\
                    + '\nTemperatura del SetPoint:' + str(SetPointValue) + ' K' + '\nHora y fecha de última lectura: ' + measureTime + ', ' + date
                    await context.bot.send_message(chat_id = User_id,text = text)

            else: 
                # logger.info('No se levantó la bandera Refresh')
                df_Users = UpdateValue_to_csv(UserInfo_path, User_id, df_Users, 0)
                Text = 'PRECUACIÓN: El historial de temperatura NO se está actualizando. Hora y fecha de la ultima actualización: '
                await context.bot.send_message(chat_id = User_id,text = Text + str(medNew[0])+', '+str(medNew[2][0:9]))

                await context.bot.send_message(chat_id = User_id,text='No se pueden realizar el monitoreo de la temperatura. Las alarmas se apagarán automáticamente.')
                context.job_queue.get_jobs_by_name(str(User_id)+ '_ReadTemp_job')[0].schedule_removal()

        else: 
            context.job_queue.get_jobs_by_name(str(User_id) + '_ReadTemp_job')[0].schedule_removal()
            Text = 'Las alarmas se desactivarán automáticante ya que no se puede obtener el valor de Set Point actual.'
            await context.bot.send_message(chat_id = User_id,text = Text)
        
    else:
        Text = 'Las alarmas y el bot se apagarán automáticamente ya que no se puede obtener información del laboratorio.'
        await context.bot.send_message(chat_id = User_id,text = Text)
        os.kill(os.getpid(), signal.SIGINT)

async def startAlarm(update: Update, context: ContextTypes.DEFAULT_TYPE): ## Enciende la alarma para detectar cambios en la temperatura.
    global df_Users
    global df_Users_Flag
    # logger.info('He recibido un comando startalerts')  
    user_id = update.effective_user.id

    if df_Users_Flag: 
        index_user_in_DF = Search_User(df_Users, user_id)

        if df_Users['Temp_Alarm'][index_user_in_DF]:
            Text = 'Las alarmas ya están activas. 🚨'
            await context.bot.send_message(chat_id = user_id ,text = Text)

        else: 
            context.job_queue.run_repeating(Temp_alarm, interval = intervalo_de_rutinaTemp, first = 1, chat_id= user_id, name = str(user_id) + '_ReadTemp_job')
            df_Users = UpdateValue_to_csv(UserInfo_path, user_id, df_Users, 0)

            text = "🚨 Alarma de temperatura activada 🚨 \n Si la temperatura cambia más de {0} K se le notificará. 🚨".format(toloreancia_Temp) 
            # logger.info(ReadTemp_job.from_aps_job())
            # logger.info(psutil.pids())
            await context.bot.send_message(chat_id = user_id ,text = text)
            print(df_Users)
    
    else:
        Text='El archivo de las tareas de usuarios NO se detectó. \n No se puede monitorear las alarmas de los usuarios así que no podrás activar ninguna.'
        await context.bot.send_message(chat_id = user_id ,text = Text)

async def stopAlarm(update: Update, context: ContextTypes.DEFAULT_TYPE): ## Apaga la alarma de la temperatura.
    # logger.info('He recibido un comando stopalerts')
    user_id = update.effective_user.id
    global df_Users
    global df_Users_Flag

    index_user_in_DF = Search_User(df_Users, user_id)

    if df_Users_Flag: 
        if df_Users['Temp_Alarm'][index_user_in_DF] == True: 
            df_Users = UpdateValue_to_csv(UserInfo_path, user_id, df_Users, 0)
            JobReadAutomatic = context.job_queue.get_jobs_by_name( str(user_id) + '_ReadTemp_job')
            # logger.info('Ya detecté el trabajo.')
            
            JobReadAutomatic[0].schedule_removal()
            # logger.info('Ya eliminé el trabajo.')
            
            # print(dict_Users)

            text = "Las alarmas se han desactivado. 🔕" 
            await update.effective_message.reply_text(text)
            print(df_Users)

        else:
            text = "NO hay alarmas activadas. 🔕"
            await update.effective_message.reply_text(text)
    
    else: 
        Text='El archivo de las tareas de usuarios NO se detectó. \n No se puede monitorear las alarmas de los usuarios.'
        await context.bot.send_message(chat_id = user_id ,text = Text)

def main():
    """Start the bot."""
    # Create the Application and pass it your bot's token.
    application = Application.builder().token(Token).build()

    # on different commands - answer in Telegram
    # application.add_handler(CommandHandler())
    application.add_handler(CommandHandler("Start", start))
    application.add_handler(CommandHandler("Help", help_command))
    application.add_handler(CommandHandler("Ayuda", help_command))
    application.add_handler(CommandHandler('Kill', stop))
    application.add_handler(CommandHandler('Detener', stop))

    ## Comandos para checar la temperatura ##
    application.add_handler(CommandHandler("WatchTemp", ReadTemperature))
    application.add_handler(CommandHandler("Temp", ReadTemperature))
    application.add_handler(CommandHandler('Temperatura', ReadTemperature))
    application.add_handler(CommandHandler('T', ReadTemperature))

    ## Comandos de la alarma ##
    application.add_handler(CommandHandler("Activar_Alarma",startAlarm))
    application.add_handler(CommandHandler("StartAlarm",startAlarm))

    application.add_handler(CommandHandler("Desactivar_Alarma",stopAlarm))
    application.add_handler(CommandHandler("OffAlarm",stopAlarm))

    ## Futuros comandos ##
    # application.add_handler(CommandHandler("getLog",getLog))
    # application.add_handler(CommandHandler("getConfig",getConfig))
    # application.add_handler(CommandHandler("setHeater",setHeater))
    # application.add_handler(CommandHandler("setPower", setPower))
    # application.add_handler(CommandHandler("setRamp",setRamp))
                            
    
    
        # on non command i.e message - echo the message on Telegram
    application.add_handler(MessageHandler(filters.TEXT & ~filters.COMMAND, echo))

    # Run the bot until the user presses Ctrl-C
    application.run_polling()

if __name__ == "__main__":
    main()
    



