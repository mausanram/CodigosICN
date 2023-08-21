import paramiko
import time

active = True
print('El proceso de copiado está activo.')
while active:
    try:
        client = paramiko.SSHClient() 
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect("132.248.29.249", 1007,username='detectores', password='damic', timeout=30)
        # print('Successful Connection')

        sftp = client.open_sftp()
        # print("Copying file")
        sftp.get(remotepath='/home/detectores/Software/LakeshoreControl/history_2023Aug11.txt',
                localpath='/home/bruce/Documents/Programas/Bot/history_2023Aug11.txt')
        # print("File Copied")
        time.sleep(1)

        sftp.close()
        client.close()
        # print('Connection Finished')

        # print('Time Sleep \n')
        with open('/home/bruce/Documents/Programas/Bot/history_2023Aug11.txt', 'r') as fl:
            lastline = fl.readlines()[-1]
            fl.close()
        print(lastline)
        time.sleep(20)

    except:
        print('Reintentado comunicación.')
        time.sleep(3)
