import math
from astropy.io import fits
import scipy.ndimage as ndimage
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys
import skimage as sk
import datetime
import pickle

histogram = open('histogram.pkl', 'rb')
exampleObj = pickle.load(histogram)
histogram.close()

print(type(exampleObj))
plt.show()
