import pickle as pkl
import numpy as np

path= './dict_mean_gains.pkl'

dict_gain = open(path, 'rb')
data_dict_gain = pkl.load(dict_gain)
dict_gain.close()

# print(data_dict_gain.keys())

ext1 = data_dict_gain['extension_1']
ext2 = data_dict_gain['extension_2']
ext4 = data_dict_gain['extension_4']

# print(ext1.keys())
print('Gain in ext1: ', np.around(ext1['Gain'], 3), ' +- ', np.around(ext1['Err_gain'], 3), ' ADU/e-')
print('Gain in ext2: ', np.around(ext2['Gain'], 3), ' +- ', np.around(ext2['Err_gain'], 3), ' ADU/e-')
print('Gain in ext4: ', np.around(ext4['Gain'], 3), ' +- ', np.around(ext4['Err_gain'], 3), ' ADU/e-', end='\n\n')

print('Sigma ext1: ', np.around(ext1['Sigma'], 3), ' ADU')
print('Sigma ext2: ', np.around(ext2['Sigma'], 3), ' ADU')
print('Sigma ext4: ', np.around(ext4['Sigma'], 3), ' ADU')
