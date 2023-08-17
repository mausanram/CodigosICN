import paramiko
import time

active = True

# while active:
#     client = paramiko.SSHClient() 
#     client.load_system_host_keys()
#     client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     client.connect("132.248.29.249", 1007,username='detectores', password='damic', timeout=30)
#     # print('Successful Connection')

#     sftp = client.open_sftp()
#     # print("Copying file")
#     sftp.get(remotepath='/home/detectores/Software/LakeshoreControl/history_2023Jul27.txt',
#             localpath='/home/bruce/Documents/Programas/Bot/history_2023Jul27.txt')
#     # print("File Copied")

#     sftp.close()
#     client.close()
    # print('Connection Finished')

    # print('Time Sleep \n')
    # time.sleep(5)

client = paramiko.SSHClient() 
client.load_system_host_keys()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect("132.248.29.249", 1007,username='detectores', password='damic', timeout=30)
# print('Successful Connection')

sftp = client.open_sftp()
# print("Copying file")
sftp.get(remotepath='/home/detectores/Software/LakeshoreControl/history_2023Jul27.txt',
        localpath='/home/bruce/Documents/Programas/Bot/history_2023Jul27.txt')
# print("File Copied")

sftp.close()
client.close()