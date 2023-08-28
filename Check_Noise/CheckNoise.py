import os
import glob
import sys
import numpy as np
import math
import argparse
from scipy.stats import norm
from astropy.stats import median_absolute_deviation as mad
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def sqrt(x, sigma0):
    return sigma0 / np.sqrt(x)

def varsort(item):
	return int(item.split("_SSAMP")[1].split("_")[0])

def doublegaus(x, norm, offset, noise, gain, mu):
	return (1.0-mu)*norm*np.exp(-((x-offset)**2/(2*(gain*noise)**2))) + mu*norm*np.exp(-((x-offset-gain)**2/(2*(gain*noise)**2)))

def sumgaus(x, norm, offset, noise, gain, mu, npeaks):
	return sum(norm*((mu**i)/np.math.factorial(i))*np.exp(-((x-offset-(i*gain))**2/(2*(gain*noise)**2))) for i in range(npeaks))

def hist_and_DataFr(files,optionFlag = 0, graphicFlag = 0):
    list_var=[]				# List to store the computed variable

    imgDict={}
    list_StandarDev = []
    list_Noises = []
    list_Names = []
    list_Temp = []
    list_Exposure = []
    list_VoltageFile = []
    list_nSamp = []

    # Define active and overscan areas
    active_mask = np.s_[:, 9:538]		#   9 <= x < 538
    overscan_mask = np.s_[:, 538:]		# 538 <= x
    mask=np.s_[:, :]			# Area where variable will be computed
    expgain = [191, 191, 193, 158]

    fig_all, axs_all = plt.subplots(1, 4, figsize=(20, 5))	# Define figure to stack histograms of all images

    for image in files:
        hdul=fits.open(image)
#		hdul.info()
        img=os.path.basename(image)				# Get basename of image

        j=files.index(image)
        figctr=0

        var=[]
        var_fit=[]

        for i in range(0, len(hdul)):

			# Load data and header
            data=hdul[i].data
            header=hdul[i].header

            if data is not None:				# Check if data is not empty

#				data = data - np.median(data, axis=0, keepdims=True)			# Subtract median per column
#				data = data - np.median(data[overscan_mask], axis=1, keepdims=True)	# Subtract OS median per row

				# Extract info from header
                string='RUNID'
                stringval=header[string]
                nsamp=float(header['NSAMP'])
                exposure = float(header['EXPOSURE'])
    

#				hlabel = string+" "+stringval			# Define label of histogram
                hlabel = img

				#Plot histogram of data to obtain offset
                hist, bin_edges = np.histogram(data[mask].flatten(), bins=1000000)
                offset = bin_edges[np.argmax(hist)]

#				print(offset)
                data = data-offset		# Subtract offset from data
                # offset = 0
				
                bin_heights, bin_borders, _ = axs_all[figctr].hist(data[mask].flatten(), range=[offset-13000, offset+13000], bins=100, histtype='step', label=hlabel)
                offset_fit = bin_borders[np.argmax(bin_heights)]
                axs_all[figctr].set_title('ext '+str(figctr+1))
                handles_all, labels_all = axs_all[figctr].get_legend_handles_labels()
                bin_centers=np.zeros(len(bin_heights), dtype=float)			# Compute centers of each bin

                for p in range(len(bin_heights)):
                    bin_centers[p]=(bin_borders[p+1]+bin_borders[p])/2

                xmin_fit, xmax_fit = offset_fit-(10*expgain[figctr])/math.sqrt(nsamp), offset_fit+(10*expgain[figctr])/math.sqrt(nsamp)			# Define fit range
                bin_heights = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]
                bin_centers = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

#				popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, p0=[np.max(bin_heights), 0, var[figctr]], maxfev=100000)	# Fit histogram with gaussian
                popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, p0=[np.max(bin_heights), 0, expgain[figctr]], maxfev=100000)		# Fit histogram with gaussian
                axs_all[figctr].plot(bin_centers, gaussian(bin_centers, *popt))								# Plot gaussian fit
#				print(popt)

#				var_fit.append(abs(popt[2])/float(stringval))
                var_fit.append(round(abs(popt[2])/float(expgain[figctr]),5)) 

#				var.append(np.std(data[mask])/float(stringval))		# Standard deviation
                var.append(round(np.std(data[mask])/float(expgain[figctr]),5))	# Standard deviation
#				var.append(np.mean(data[mask]))		# Mean
#				var.append(np.median(data[mask]))	# Median

                figctr=figctr+1
                
                # plt.legend()
        # plt.show()			# Show histogram per image
        ImgName=str(image)
        Tempt = ImgName.split('_')
        list_Temp.append(Tempt[6][1:])
        list_VoltageFile.append(Tempt[5][-1])

        list_StandarDev.append(var)
        list_Noises.append(var_fit)
        list_Names.append(int(stringval))
        list_nSamp.append(int(header['NSAMP']))
        list_Exposure.append(int(header['EXPOSURE']))
        
        # Añade elementos al diccionario
        imgDict.setdefault('img'+stringval,{'stdv':[round(var[0],5), round(var[1],5), round(var[2],5), round(var[3],5)],
                            'noise':[round(var_fit[0],5), round(var_fit[1],5), round(var_fit[2],5), round(var_fit[3],5)]})


        datos_imgDict= pd.DataFrame([key for key in imgDict.keys()])
        datos_imgDict['Standard Deviation']=[value['stdv'] for value in imgDict.values()]
        datos_imgDict['Noise']=[value['noise'] for value in imgDict.values()]
		


#		print(float(stringval), var_fit[0], var_fit[1], var_fit[2], var_fit[3])
#		list_var.append([float(stringval), var[0], var[1], var[2], var[3]])
        list_var.append([float(stringval), var_fit[0], var_fit[1], var_fit[2], var_fit[3]])		# To use the var obtained from the gaussian fit

		# print(imgDict)

    arr_var=np.array(list_var)
    arr_var=arr_var[np.argsort(arr_var[:, 0])]	# Sort array by values on first column
#	print(arr_var)

    fig_all.legend(handles_all, labels_all, loc='upper right') # histogram figure  
    #plt.show()

#	fig_var.tight_layout()

    if graphicFlag == 1:
        fig_var, axs_var = plt.subplots(1, 4, figsize=(20, 5), sharey=True) # Noise figure
        for k in range(0, 4):
            axs_var[k].plot([row[0] for row in arr_var], [row[k+1] for row in arr_var], ".k")
            axs_var[k].plot([row[0] for row in arr_var], sqrt([row[0] for row in arr_var], arr_var[0, k+1]), "-r")		# Sqrt fit when doing noise vs nsamp
            axs_var[k].set_title('ext '+str(k+1))
            axs_var[k].set_xlabel(string)
            axs_var[k].set_ylabel('Noise (e-)')
            axs_var[k].grid(True)
    #		axs_var[k].set_yscale('log')
    #		axs_var[k].set_ylim([0, 200])
        #plt.show()				# Show plot

    DFrame = Build_DFrame(list_Temp,list_Names,list_nSamp,list_StandarDev,list_Noises,list_Exposure,list_VoltageFile, optionFlag, graphicFlag)

    return DFrame

def Build_DFrame(list_Temp,list_Names,list_nSamp,list_StandarDev,list_Noises,list_Exposure,list_VoltageFile, optionFlag, graphicFlag):
    Tempt_frame = pd.DataFrame(list_Temp, columns=["Temp"])
    names_frame = pd.DataFrame(list_Names, columns=["Image"])
    nSamp_frame = pd.DataFrame(list_nSamp, columns=["nSamp"])
    stdev_frame = pd.DataFrame(list_StandarDev, columns=['Std ext1','Std ext2','Std ext3','Std ext4'])
    noises_frame = pd.DataFrame(list_Noises, columns=['Noise ext1','Noise ext2','Noise ext3','Noise ext4'])
    exposure_frame = pd.DataFrame(list_Exposure, columns=['Exposure'] )
    voltage_frame = pd.DataFrame(list_VoltageFile, columns=['Voltage File'])

    if optionFlag == 0: #ID //
        total_frame=pd.concat([names_frame, stdev_frame,noises_frame],axis=1)
        IDFrame = total_frame.sort_values('Image')
        return IDFrame

    elif optionFlag == 1:#nSamp //
        total_frame=pd.concat([names_frame, nSamp_frame, stdev_frame,noises_frame],axis=1)
        nSampFrame = total_frame.sort_values('nSamp')

        if graphicFlag == 2:
            _, axes= plt.subplots(nrows=1,ncols=4, figsize=(20,5))
            nSampFrame.plot(kind='scatter',x='nSamp',y='Noise ext1', ax=axes[0])
            nSampFrame.plot(kind='scatter',x='nSamp',y='Noise ext2', ax=axes[1])
            nSampFrame.plot(kind='scatter',x='nSamp',y='Noise ext3', ax=axes[2])
            nSampFrame.plot(kind='scatter',x='nSamp',y='Noise ext4', ax=axes[3])

        return nSampFrame 

    elif optionFlag == 2: #Std ext1 //
        total_frame=pd.concat([names_frame, stdev_frame,noises_frame],axis=1)
        Std_ext1Frame = total_frame.sort_values('Std ext1')

        return Std_ext1Frame
    
    elif optionFlag == 3: #Noise ext1 \\
        total_frame=pd.concat([names_frame, stdev_frame,noises_frame],axis=1)
        Noise_ext1Frame = total_frame.sort_values('Noise ext1')

        return Noise_ext1Frame
    
    elif optionFlag == 4: #Temp
        total_frame=pd.concat([names_frame,Tempt_frame, stdev_frame,noises_frame],axis=1)
        TempFrame = total_frame.sort_values('Temp')

        return TempFrame
    
    elif optionFlag == 5: #Exposure
        total_frame=pd.concat([names_frame,exposure_frame, stdev_frame,noises_frame],axis=1)
        ExpFrame = total_frame.sort_values('Exposure')

        return ExpFrame
    
    elif optionFlag == 6: #Voltage
        total_frame=pd.concat([names_frame,voltage_frame, stdev_frame,noises_frame],axis=1)
        VoltFrame = total_frame.sort_values('Voltage File')

        return VoltFrame

def Noise(files, parplotFlag = 1):
	#files.sort(key=varsort)			# Sort files by key
#	files.sort(key=os.path.getmtime)

	dirname=os.path.dirname(files[0])	# Get dirname of first element in files

#	latest_file = max(files, key=os.path.getctime)
#	print(latest_file)
#	image='image.fz'

	# Define active and overscan areas
	active_mask = np.s_[:, 9:538]		#   9 <= x < 538
	overscan_mask = np.s_[:, 538:]		# 538 <= x

	mask=np.s_[:, 538:700]			# Area where variable will be computed

	list_var=[]				# List to store the computed variable

	expgain = [201.8325949210918, 194.70825464645284, 202.97945260519987, 193.2145155088731]	# Expected gain; if only fitting 1 peak, noise is divided by this number
	numpeaks = 2			# Number of peaks to fit

	varsplot = ["Constant (ADU)", "Offset (ADU)", "Noise (e-)", "Gain (ADU/e-)", "SER (e-/pix)"]	# Variables that can be chosen to be plotted
	#parplot = 5				# From varsplot, index of parameter to plot; for example if you want to plot Noise (e-), varsplot=3
	string=''				# Variable to be the x axis, as shown in the header of the image

	#	fig_all, axs_all = plt.subplots(1, 4, figsize=(20, 5))		# Define figure to stack histograms of all images
	#	fig_all.tight_layout()

		# Open datafile and write header
	#	datafile=open(dirname+"/noisevspsamp.txt","a+")
	#	datafile.write("#PSAMP\tNoise_ext1\tNoise_ext2\tNoise_ext3\tNoise_ext4\n")

	deltaH=[]
	deltaV=[]
	deltaT=[]
	deltaSW=[]
	RunID=[]

	for image in files:
		fig_all, axs_all = plt.subplots(1, 4, figsize=(20, 5), sharey=True)	# Define figure to stack histogram of each image
#		fig_all.tight_layout()

		hdul=fits.open(image)
#		hdul.info()
		img=os.path.basename(image)				# Get basename of image
		print("\nImage: "+img)

		j=files.index(image)
		figctr=0

		# Extract info from img title
#		string=image.split("img")[1].split(".")[0]		# To obtain imgID

		var=[]; var_fit=[]

		print("# of peaks to fit: "+str(numpeaks)+"\n")

#		for i in range(0,1):
		for i in range(0, len(hdul)):
			data=hdul[i].data; header=hdul[i].header	# Load data and header

			if data is not None:								# Check if data is not empty
				# data = data - np.median(data, axis=0, keepdims=True)			# Subtract median per column
				# data = data - np.median(data[overscan_mask], axis=1, keepdims=True)	# Subtract OS median per row (use when proc*fits have no baseline substracted)

				# Extract info from header
				stringval=header["RUNID"]
				#stringval=float(header["H1AH"])-float(header["H1AL"])
				nsamp=float(header['NSAMP'])

#				hlabel = string+" "+stringval			# Define label of histogram
				hlabel = img

				#Plot histogram of data to obtain offset
				hist, bin_edges = np.histogram(data[mask].flatten(), bins=1000000)
				offset = bin_edges[np.argmax(hist)]
#				print(offset)
#				data = data-offset		# Subtract offset from data
#				offset = 0			# Offset to plot
				
				bin_heights, bin_borders, _ = axs_all[figctr].hist(data[mask].flatten(), range=[offset-3000, offset+3000], bins=200, histtype='step', label=hlabel)	# Plot histogram
				offset_fit = bin_borders[np.argmax(bin_heights)]	# Offset to fit
#				offset_fit = 0						# Offset to fit
				axs_all[figctr].set_title('ext '+str(figctr+1))
				handles_all, labels_all = axs_all[figctr].get_legend_handles_labels()

				bin_centers=np.array([(bin_borders[p+1]+bin_borders[p])/2 for p in range(len(bin_heights))])	# Compute centers of each bin

#				if nsamp>100: numpeaks=2
#				else: numpeaks=1
#				print("# of peaks to fit: "+str(numpeaks))

				# Fit gaussians
				aux_arr=np.array([])
				for npeak in range(numpeaks):
					xmin_fit = offset_fit+(npeak*expgain[figctr])-(5*expgain[figctr])/math.sqrt(nsamp) 	# Define fit range
#					xmin_fit = offset_fit+(npeak*expgain[figctr])-(0.25*expgain[figctr])
					xmax_fit = offset_fit+(npeak*expgain[figctr])+(5*expgain[figctr])/math.sqrt(nsamp)
#					xmax_fit = offset_fit+(npeak*expgain[figctr])+(0.25*expgain[figctr])
					print("Ext"+str(figctr)+" trying to fit peak of "+str(npeak)+"e- between: xmin = "+str(xmin_fit)+" & xmax = "+str(xmax_fit))

					bin_heights_peak = bin_heights[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]		# Constrain histogram to given range
					bin_centers_peak = bin_centers[(bin_centers>xmin_fit) & (bin_centers<xmax_fit)]

					try:	# Try to fit, pass if error
						popt, pcov = curve_fit(gaussian, bin_centers_peak, bin_heights_peak, p0=[np.max(bin_heights_peak), bin_centers_peak[np.argmax(bin_heights_peak)], 0.5*expgain[figctr]], maxfev=100000, bounds=([0, xmin_fit, 0.01*expgain[figctr]], [1.5*np.max(bin_heights_peak), xmax_fit, 5*expgain[figctr]]))	# Fit histogram with gaussian
#						print(popt)
						axs_all[figctr].plot(bin_centers_peak, gaussian(bin_centers_peak, *popt))			# Plot gaussian fit
						par_fit = np.append(aux_arr, popt)
						aux_arr = par_fit
						print("Successful fit :)")
					except: print("ERROR in fit :(")

				try:	# Enters here if numpeaks>=2
					norm = par_fit[0]; offset = par_fit[1]; gain = par_fit[4]-par_fit[1]; noise = par_fit[2]/gain; mu = par_fit[3]/par_fit[0];
					print("Trying complex fit with initial parameters:")
					print(norm, offset, gain, noise, mu)

					bin_heights_doublegaus = bin_heights[(bin_centers>(offset_fit-0.5*expgain[figctr])) & (bin_centers<(offset_fit+1.5*expgain[figctr]))]
					bin_centers_doublegaus = bin_centers[(bin_centers>(offset_fit-0.5*expgain[figctr])) & (bin_centers<(offset_fit+1.5*expgain[figctr]))]
					popt, pcov = curve_fit(doublegaus, bin_centers_doublegaus, bin_heights_doublegaus, p0=[norm, offset, noise, gain, mu], maxfev=100000, bounds=([0, offset_fit-0.5*expgain[figctr], 0.001, 0, 0], [np.inf, offset_fit+1.5*expgain[figctr], np.inf, 2*expgain[figctr], np.inf]))		# Fit histogram with double gaussian
					axs_all[figctr].plot(bin_centers_doublegaus, doublegaus(bin_centers_doublegaus, *popt))		# Plot gaussian fit
#					popt, pcov = curve_fit(sumgaus, bin_centers_doublegaus, bin_heights_doublegaus, p0=[norm, offset, noise, gain, mu, 2], maxfev=100000)	# Fit histogram with sum of gaussians
#					axs_all[figctr].plot(bin_centers_doublegaus, sumgaus(bin_centers_doublegaus, *popt))							# Plot gaussian fit

#					var_fit.append(abs(popt[2]*popt[3])/float(stringval))
					var_fit.append(popt[parplotFlag-1])
					print("Successful fit :)"); print("Extracting "+varsplot[parplotFlag-1])
                                        
				except:
                                        
					try:	# Enters here if numpeaks=1, try to extract parameters from gaussian fit to 0e- peak
						norm = par_fit[0]; offset = par_fit[1]; noise = par_fit[2]/expgain[figctr]
                                                
						if parplotFlag == 1: 
							var_fit.append(norm); print("Extracting "+varsplot[parplotFlag-1])

						if parplotFlag == 2: 
							var_fit.append(offset); print("Extracting "+varsplot[parplotFlag-1])

						if parplotFlag == 3: 
							var_fit.append(noise); print("Extracting "+varsplot[parplotFlag-1])

						else: 
							var_fit.append(0); print("ERROR extracting "+varsplot[parplotFlag-1])

					except: 
						var_fit.append(0); print("ERROR in fit :(")

#				var.append(np.std(data[overscan_mask])/float(stringval))	# Standard deviation
				var.append(np.std(data[overscan_mask]))				# Standard deviation
#				var.append(np.mean(data[overscan_mask]))			# Mean
#				var.append(np.median(data[overscan_mask]))			# Median

				figctr=figctr+1

#		plt.legend()
#		plt.show()			# Show histogram per image
		

		# STORE COMPUTED VARIABLE
		print("\nNoise in overscan from stddev [ADU]:")
		print(var[0], var[1], var[2], var[3])
		print(varsplot[parplotFlag-1]+" in selected area:")
		print(var_fit[0], var_fit[1], var_fit[2], var_fit[3])

#		list_var.append([float(stringval), var[0], var[1], var[2], var[3]])
		list_var.append([float(stringval), var_fit[0], var_fit[1], var_fit[2], var_fit[3]])		# To use the var obtained from the gaussian fit

		deltaH.append(float(header["H1AH"])-float(header["H1AL"]))
		deltaV.append(float(header["V1AH"])-float(header["V1AL"]))
		deltaT.append(float(header["TGAH"])-float(header["TGAL"]))
		deltaSW.append(float(header["SWAH"])-float(header["SWAL"]))
		RunID.append(float(header["RUNID"]))

		fig_all.legend(handles_all, labels_all, loc='upper right')
#		plt.legend()
		

	arr_var=np.array(list_var)
	arr_var=arr_var[np.argsort(arr_var[:, 0])]	# Sort array by values on first column
#	print(arr_var)
	plt.show()
	

	
	# PLOT
	fig_var, axs_var = plt.subplots(2, 4, figsize=(20, 5))
#	fig_var.tight_layout()

	for k in range(0, 4):
		axs_var[0][k].plot([row[0] for row in arr_var], [row[k+1] for row in arr_var], ".k")
#		axs_var[k].plot([row[0] for row in arr_var], sqrt([row[0] for row in arr_var], arr_var[0, k+1]), "-r")		# Sqrt fit when doing noise vs nsamp
		axs_var[0][k].set_title('ext '+str(k+1))
		axs_var[0][k].set_xlabel(string)
		axs_var[0][k].set_ylabel(varsplot[parplotFlag-1])
		axs_var[0][k].grid(True)
#		axs_var[k].set_xscale('log')
#		axs_var[k].set_yscale('log')
#		axs_var[k].set_ylim(ymin=1)

	axs_var[1][0].plot(RunID,deltaV,".k")
	axs_var[1][0].set_title("Delta V")
	axs_var[1][0].set_ylabel("Volts")
	axs_var[1][0].set_xlabel("RunID")
	axs_var[1][1].plot(RunID,deltaT,".k")
	axs_var[1][1].set_title("Delta T")
	axs_var[1][1].set_xlabel("RunID")
	axs_var[1][2].plot(RunID,deltaH,".k")
	axs_var[1][2].set_title("Delta H")
	axs_var[1][2].set_xlabel("RunID")
	axs_var[1][3].plot(RunID,deltaSW,".k")
	axs_var[1][3].set_title("Delta SW")
	axs_var[1][3].set_xlabel("RunID")

	#	plt.savefig(dirname+"/checknoise.png")	# Save plot
		#plt.show()				# Show plot

		# TXT FILE
	#	np.savetxt(datafile, arr_var, fmt="%s")	# Save array to datafile
	#	datafile.close()			# Close datafile'

def parser(args):

    parser = argparse.ArgumentParser(prog='Check Noise',description='Displays a dataframe of image ID, Standard Deviation and Noise for each extension and \
                                     can display extra information for each image (ex. Temperature, Voltage File, Exposure,...), and plots the histogram for each image.')
    
    parser.add_argument('path',type=str,nargs='+', help="Save the file in a list")
    parser.add_argument('--ID',type = str, nargs= '+', help = "Displays a DataFrame sorted by ID column, the histogram and the noise plot by extension" )
    parser.add_argument('--StD',type = str, nargs= '+', help = "Displays a DataFrame sorted by standard deviation column (StD ext1), the histogram and the noise plot by extension")
    parser.add_argument('--Noise', '-N','-n', type = str, nargs= '+', help = "Displays a DataFrame sorted by noise column (Noise ext1), the histogram and the noise plot by extension")
    parser.add_argument('--NSamp','--NS',type= str, nargs= '+', help = "Displays a DataFrame sorted by number of samp. column (NSamp),the histogram and graph of number of samples per image")
    parser.add_argument('--Temp','-T',type = str, nargs= '+', help = "Displays a DataFrame sorted by Temperature column (Temp) and the histogram")
    parser.add_argument('--Exposure','--Exp','-E',type = str, nargs= '+', help = "Displays a DataFrame sorted by exposure column (Exposure) and the histogram")
    parser.add_argument('--Voltage','--Volt','-V',type = str, nargs= '+', help = "Displays a DataFrame sorted by voltage file column (Voltage), and the histogram")
    
    parser.add_argument('--Const','-C','-c', action='store_true', help = " Plot the ... for each image" )
    parser.add_argument('--Off','--off', action='store_true', help = " Plot the ... for each image")
    parser.add_argument('--DeltaNoise', action='store_true', help = "Plot the noise for each image")
    parser.add_argument('--Gain','-gain','-g', action='store_true', help = "Plot the gain for each image")
    parser.add_argument('--SER',action='store_true', help = "Plot the Single Electron Rate for each image")

    argObj = parser.parse_args()
    return argObj

def main(argObj):
    value = os.fork()

    if argObj.ID:
        files = argObj.path
        graphicFlag = 1
        IDFrame = hist_and_DataFr(files, graphicFlag = 1)
        if value == 0:
            plt.show() # Show hist and noise
        else: 
            print(IDFrame.to_string(index=False))
            print('ID OK')
    
    elif argObj.NSamp:
        files = argObj.path
        optionFlag = 1
        graphicFlag = 2
        nSampFrame= hist_and_DataFr(files, optionFlag, graphicFlag)
        if value == 0:
            plt.show() # Show hist and nsamp

        else: 
            print(nSampFrame.to_string(index=False))
            print('nSamp OK')
    
    elif argObj.StD:
            files = argObj.path
            optionFlag = 2 
            graphicFlag = 1
            StD_ext1Frame= hist_and_DataFr(files,optionFlag, graphicFlag)
            if value==0:
                plt.show() # Show hist and noise 

            else:
                print(StD_ext1Frame.to_string(index=False))
                print('Std ext1 OK')
    
    elif argObj.Noise:
        files = argObj.path
        optionFlag = 3
        graphicFlag = 1
        Noise_ext1Frame = hist_and_DataFr(files,optionFlag, graphicFlag)
        if value == 0:
            plt.show()
        else:
            print(Noise_ext1Frame.to_string(index=False))
            print('Noise ext1 OK')
    
    elif argObj.Temp:
        files = argObj.path
        optionFlag = 4
        TempFrame = hist_and_DataFr(files,optionFlag)
        if value == 0:
            plt.show()
        else:    
            print(TempFrame.to_string(index=False))

    elif argObj.Exposure:
        files = argObj.path
        optionFlag = 5
        ExpFrame = hist_and_DataFr(files,optionFlag)
        if value == 0:
            plt.show()
        else:
            print(ExpFrame.to_string(index=False))
    
    elif argObj.Voltage:
        files = argObj.path
        optionFlag = 6
        VoltFrame = hist_and_DataFr(files,optionFlag)
        if value == 0:
            plt.show()
        else:
            print(VoltFrame.to_string(index=False))

#########################################################

    elif argObj.Const :
        files = argObj.path
        Noise(files)
        plt.show()
        # print('Todo OK')

    elif argObj.Off:
        files = argObj.path
        parplotFlag = 2
        Noise(files, parplotFlag)
        plt.show()
        # print('Todo OK')

    elif argObj.Noise:
        files = argObj.path
        parplotFlag = 3
        Noise(files, parplotFlag)
        plt.show()
        # print('Todo OK')

    elif argObj.Gain:
        files = argObj.path
        parplotFlag = 4
        Noise(files, parplotFlag)
        plt.show()
        # print('Todo OK')

    elif argObj.SER:
        files = argObj.path
        parplotFlag = 5
        Noise(files, parplotFlag)
        plt.show()
        # print('Todo OK')

    else:
        files = argObj.path
        Noise(files)
        plt.show()
        # print('Todo OK')

if __name__ == "__main__":
   argObj=parser(sys.argv)
   exitcode=main(argObj)
   exit(code=exitcode)