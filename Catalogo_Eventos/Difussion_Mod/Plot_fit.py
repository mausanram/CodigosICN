import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

## Plot's Configuration
plt.rcParams.update({
    "image.origin": "lower",
    "image.aspect": 0,
    #"text.usetex": True,
    "grid.alpha": .5,
    "axes.linewidth":2,
    "lines.linewidth" : 1,
    "font.size":    15.0,
    "xaxis.labellocation": 'right',  # alignment of the xaxis label: {left, right, center}
    "yaxis.labellocation": 'top',  # alignment of the yaxis label: {bottom, top, center}
    "xtick.top":           True ,  # draw ticks on the top side
    "xtick.major.size":    8    ,# major tick size in points
    "xtick.minor.size":    4      ,# minor tick size in points
    "xtick.direction":     'in',
    "xtick.minor.visible": True,
    "ytick.right":           True ,  # draw ticks on the top side
    "ytick.major.size":    8    ,# major tick size in points
    "ytick.minor.size":    4      ,# minor tick size in points
    "ytick.direction":     'in',
    "ytick.minor.visible": True,
    "ytick.major.width":   2   , # major tick width in points
    "ytick.minor.width":   1 ,
    "xtick.major.width":   2   , # major tick width in points
    "xtick.minor.width":   1 ,
    "legend.framealpha": 0 ,
    "legend.loc": 'best',

})



def diffution_curve(x, alpha, beta):
    return np.sqrt((alpha * np.log(1 - (beta * x))))/15

def data_extraction(path, extension, list_muons, flag_flip, delta_cut, muon_type):
    list_spreads_all = []
    list_depth_all = []
    list_spreads_without_deltas = []


    # if option == 2 or option == 5 or option == 7 or option == 10:
    #     flag_flip = True

    for index in range(0, len(list_muons)):
        file_path = path + f"/muon{list_muons[index]}.xye"    # For Vertical muons
        column_names = ['x', 'y', 'charge']
        df = pd.read_csv(file_path, sep='\s+', names=column_names, comment='#')

        if muon_type == 0 or muon_type == 1:
            row_charge = df.groupby('y')['charge'].sum().reset_index() # For Vertical muons
            true_list_linecharge = list(row_charge['y'])
        # print(row_charge.head())
        else:
            row_charge = df.groupby('x')['charge'].sum().reset_index() # For Horizontal muons
            true_list_linecharge = list(row_charge['x'])


        with open(path+ f"/muon{list_muons[index]}_spreads.txt", "r") as file:
            list_lines = file.readlines()
        
        list_spreads = []
        for line in list_lines:
            true_line = line.split(" ")
            line_charge = float(true_list_linecharge[list_lines.index(line)])
            # print(f"line charge: {line_charge}")

            spread = true_line[0]
            # depth = true_line[1].split("\n")[0]

            # list_spreads_without_deltas.append(spread)

            if line_charge > delta_cut:
                # print(f"line charge: {line_charge}")
                list_spreads.append(0)
                # list_depth.append(float(depth))
            else:
                list_spreads.append(float(spread))

        if flag_flip:
            list_spreads.reverse()

        thk=725.    # (um) CCD thickness.
        pixwd=15.   # (um) pixel size.
        ms = np.arange(0,len(list_spreads))+1
        muonl=ms[-4]*pixwd  # (um) Muon length
        tgang=thk/muonl # Tangent of muon angle.
        z=((ms-0.5)*pixwd)*tgang

        for jndex in range(0, len(list_spreads)-5):
            list_depth_all.append(z[jndex])
            list_spreads_all.append(list_spreads[jndex])

    return list_depth_all, list_spreads_all, list_spreads_without_deltas

def data_visualization(list_depth_all, list_spread_all, extension):
    thk = 725
    list_thk = np.arange(0, thk)

    fig, axs = plt.subplots(figsize = [10,10])
    axs.scatter(list_depth_all, list_spread_all, color="b", marker=".", s=9)

    popt_DM, _ = curve_fit(diffution_curve, list_depth_all, list_spread_all, maxfev=100000, p0= [-200, 0.001])
    dict_diffution_model = {'Alpha' : popt_DM[0], 'Beta' : popt_DM[1]}
    print('Alpha: ', dict_diffution_model['Alpha'], ' Beta: ', dict_diffution_model['Beta'])

    axs.plot(list_thk, diffution_curve(list_thk, -211.9, 0.00102), "--k", linewidth=2,
            label = r"CONNIE's paper: ($\alpha$=-211.9 $m^{2}$, $\beta$=0.00102$\mu m^{-1}$)")
    
    axs.plot(list_thk, diffution_curve(list_thk, dict_diffution_model["Alpha"], dict_diffution_model["Beta"]), "--r", linewidth=2,
            label = r"$\sqrt{\alpha \ln(1 - \beta z)}$ fit: ($\alpha$= "+ str(round(dict_diffution_model["Alpha"],1))+ 
                    r" $m^{2}$, $\beta$= " + str(round(dict_diffution_model["Beta"],5)) + r" $\mu m^{-1}$)")

    # axs.set_ylim(0,1.25)
    # axs.set_xlim(-1, 680)

    axs.set_ylabel("Spread (px)")
    axs.set_xlabel(r"Depth ($\mu$m)")
    axs.set_title(f"Size-to-Depth Relation (ICN data) Extension {extension}")

    axs.legend()
    axs.grid()
    plt.show()
    return 0

def plot_array(ext, option, muons):
    list_spreads_all = []
    list_depth_all = []

    if option == 1:
        basedir = f"./muons_ext{ext}VUP"
    if option == 2:
        basedir = f"./muons_ext{ext}VDOWN"
    if option == 3:
        basedir = f"./muons_ext{ext}HRIGHT1"
    if option == 4:
        basedir = f"./muons_ext{ext}HRIGHT2"
    if option == 5:
        basedir = f"./muons_ext{ext}HLEFT"

    if option == 6:
        basedir = f"./muons_ext{ext}VUP"
    if option == 7:
        basedir = f"./muons_ext{ext}VDOWN"
    if option == 8:
        basedir = f"./muons_ext{ext}HRIGHT1"
    if option == 9:
        basedir = f"./muons_ext{ext}HRIGHT2"
    if option == 10:
        basedir = f"./muons_ext{ext}HLEFT"

    if option == 2 or option == 5 or option == 7 or option == 10:
        flag_flip = True
    else:
        flag_flip = False

    type_muon = basedir.split(str(ext))[1]

    for index in range(0, len(muons)):
        with open(basedir + f"/muon{muons[index]}_spreads.txt", "r") as file:
            list_lines = file.readlines()
            list_spreads = []
            list_depth = []

            for line in list_lines:
                true_line = line.split(" ")
                spread = true_line[0]
                depth = true_line[1].split("\n")[0]

                list_spreads.append(float(spread))
                # list_depth.append(float(depth))

        if flag_flip:
            list_spreads.reverse()

        thk=725.    # (um) CCD thickness.
        pixwd=15.   # (um) pixel size.
        ms = np.arange(0,len(list_spreads))+1
        muonl=ms[-4]*pixwd  # (um) Muon length
        tgang=thk/muonl # Tangent of muon angle.
        z=((ms-0.5)*pixwd)*tgang

        for jndex in range(0, len(list_spreads)-5):

            list_depth_all.append(z[jndex])
            list_spreads_all.append(list_spreads[jndex])

    list_thk = np.arange(0, thk)
    # fit_values = diffution_curve(list_depth, -341.3, 0.000816)

    # print(list_depth_all[:10])

    fig, axs = plt.subplots(figsize = [10,10])
    axs.scatter(list_depth_all, list_spreads_all, color="b", marker=".", s=9)

    popt_DM, _ = curve_fit(diffution_curve, list_depth_all, list_spreads_all, maxfev=100000, p0= [-200, 0.001])
    dict_diffution_model = {'Alpha' : popt_DM[0], 'Beta' : popt_DM[1]}
    print('Alpha: ', dict_diffution_model['Alpha'], ' Beta: ', dict_diffution_model['Beta'])

    # axs.plot(list_thk, diffution_curve(list_thk, -212.15831763733, 0.001056), "--b", linewidth=2,
    #          label = r"Test")

    axs.plot(list_thk, diffution_curve(list_thk, -211.9, 0.00102), "--k", linewidth=2,
            label = r"CONNIE's paper: ($\alpha$=-211.9 $m^{2}$, $\beta$=0.00102$\mu m^{-1}$)")
    
    # axs.plot(list_thk, diffution_curve(list_thk, -305.2, 0.00094), "--c", linewidth=2,
    #         label = r"EXT1 VDOWN: ($\alpha$=-305.2 $m^{2}$, $\beta$=0.00094$\mu m^{-1}$)")
    
    # axs.plot(list_thk, diffution_curve(list_thk, -484.6, 0.00078), "--y", linewidth=2,
    #         label = r"EXT1 HRIGHT1: ($\alpha$=-484.6 $m^{2}$, $\beta$=0.00078$\mu m^{-1}$)")
    
    # axs.plot(list_thk, diffution_curve(list_thk, -443.3, 0.00081), "--g", linewidth=2,
    #         label = r"EXT1 HRIGHT2: ($\alpha$=-443.3 $m^{2}$, $\beta$=0.00081$\mu m^{-1}$)")
    
    # axs.plot(list_thk, diffution_curve(list_thk, -460.8, 0.00079), "--g", linewidth=2,
    #         label = r"EXT1 HLEFT: ($\alpha$=-460.8 $m^{2}$, $\beta$=0.00079$\mu m^{-1}$)")


    ### EXTENSION 1 ###
    # Alpha = -274.4
    # Beta = 0.00105
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--y", linewidth=2,
    #         label = r"EXT1 VUP: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -354.1
    # Beta = 0.00082
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--g", linewidth=2,
    #         label = r"EXT1 VDOWN: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -260.8
    # Beta = 0.00101
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--c", linewidth=2,
    #         label = r"EXT1 HRIGHT1: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -223.7
    # Beta = 0.00109
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--r", linewidth=2,
    #         label = r"EXT1 HRIGHT2: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -229.0
    # Beta = 0.00108
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--b", linewidth=2,
    #         label = r"EXT1 HLEFT: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")

    # axs.plot(list_thk, diffution_curve(list_thk, dict_diffution_model["Alpha"], dict_diffution_model["Beta"]), "--r", linewidth=2,
    #         label = r"$\sqrt{\alpha \ln(1 - \beta z)}$ fit: ($\alpha$= "+ str(round(dict_diffution_model["Alpha"],1))+ 
    #                 r" $m^{2}$, $\beta$= " + str(round(dict_diffution_model["Beta"],5)) + r" $\mu m^{-1}$)")
    

    ### EXTENSION 2 ###
    Alpha = -234.9
    Beta = 0.00111
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--y", linewidth=2,
    #         label = r"EXT2 VUP: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -318.2
    # Beta = 0.00088
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--g", linewidth=2,
    #         label = r"EXT2 VDOWN: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -176.9
    # Beta = 0.00109
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--c", linewidth=2,
    #         label = r"EXT2 HRIGHT1: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -229.5
    # Beta = 0.00095
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--r", linewidth=2,
    #         label = r"EXT2 HRIGHT2: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -172.9
    # Beta = 0.00112
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--b", linewidth=2,
    #         label = r"EXT2 HLEFT: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")
    
    # Alpha = -235.7
    # Beta = 0.00105
    # axs.plot(list_thk, diffution_curve(list_thk, Alpha, Beta), "--b", linewidth=2,
    #         label = r"EXT2 HLEFT: $\alpha$= " + str(Alpha) + r" $m^{2}$, $\beta$= " + str(Beta) + r" $\mu m^{-1}$)")

    axs.plot(list_thk, diffution_curve(list_thk, dict_diffution_model["Alpha"], dict_diffution_model["Beta"]), "--r", linewidth=2,
            label = r"$\sqrt{\alpha \ln(1 - \beta z)}$ fit: ($\alpha$= "+ str(round(dict_diffution_model["Alpha"],1))+ 
                    r" $m^{2}$, $\beta$= " + str(round(dict_diffution_model["Beta"],5)) + r" $\mu m^{-1}$)")

    # axs.set_ylim(0,1.25)
    # axs.set_xlim(-1, 680)

    axs.set_ylabel("Spread (px)")
    axs.set_xlabel(r"Depth ($\mu$m)")
    # axs.set_title("Size-to-Depth Relation (Horizontal Muons)")
    axs.set_title(f"Size-to-Depth Relation ({type_muon} Muons)")
    # axs.set_title(f"Size-to-Depth Relation Extension 2")

    axs.legend()
    axs.grid()
    plt.show()

def main():
    Extension = 2 # Choose the extension
    flag_wholeExt = True # Active this to analyze all the types of muons at same time

    flag_VUP = True # Active just one type of muons
    flag_VDOWN = False
    flag_HRIGHT1 = False
    flag_HRIGHT2 = False
    flag_HLEFT = False

    if flag_wholeExt:
        list_depths_all = []
        list_spreads_all = []
        list_types = ["VUP", "VDOWN", "HRIGHT1", "HRIGHT2", "HLEFT"]
        if Extension == 1:
            list_deltacut = [1000, 1000, 1000, 1000, 1000] # Delta cut threshold
            # list_deltacut = [194, 194, 194, 194, 194] # Mean Delta cut threshold
            # list_deltacut = [180, 180, 180, 180, 180]
            muons = [
                [261, 353, 760, 765, 1096, 1121],
                [119, 321, 385],
                [316, 505, 857, 948],
                [73, 157, 204, 280, 329, 429, 504, 770, 783, 1014, 1114],
                [330, 371, 773, 914, 1002, 1023]
            ]
            # muons = [
            #     [353,760,765, 1096],
            #     [119, 192, 321, 385,550], 
            #     [316,505,637,857,948],
            #     [73,110,204,280,329,448,770],
            #     [330, 371, 773, 914, 1002, 1023]
            # ]
        else: 
            list_deltacut = [1000, 1000, 1000, 1000, 1000] # Delta cut threshold
            # list_deltacut = [180, 180, 180, 180, 180] # Delta cut threshold
            muons = [
                [99, 130, 145, 389, 449, 1002, 1133, 1166, 1216, 1337],
                [20, 64, 301, 793, 1167, 1175],
                [315, 464, 709, 727, 1011, 1245],
                [395, 790, 890, 1056, 1066],
                [47, 170, 228, 246, 252, 347, 379, 388, 524, 701, 723, 839, 859, 965,1206]
            ]

        for type in range(0, 5):
            # print("Type: ", type)
            # if type == 1:
            #     print(list_types[type])
            #     continue
            basedir = f"./muons_ext{Extension}{list_types[type]}"
            threshold = list_deltacut[type]
            print(f"threshold: {threshold}")
            if type == 1 or type == 4:
                flag_flip = True
            else:
                flag_flip = False
            list_depths, list_spreads, list_spreads_deltas = data_extraction(basedir, Extension, muons[type], flag_flip, threshold, type)

            for index in range(0, len(list_depths)):
                if list_spreads[index] > 0:
                    list_depths_all.append(list_depths[index])
                    list_spreads_all.append(list_spreads[index])
                else:
                    continue
                # list_spreads_all.append(list_spreads_deltas[index])

        data_visualization(list_depths_all, list_spreads_all, Extension)
        flag_VUP, flag_VDOWN, flag_HRIGHT1, flag_HRIGHT2, flag_HLEFT = False, False, False, False, False
    
    if flag_VUP and Extension == 1:
        # list_muons = [39, 80, 81, 96, 210, 214, 233, 261, 294, 353, 406, 414, 426, 459, 532, 653, 760, 765, 915, 1096, 1121]
        list_muons = [261, 353, 760, 765, 1096, 1121]   # Whith PDF and plot
        # list_muons = [1121]

        # list_muons = [353,760,765, 1096]
        # list_muons = [915]
        flag_option = 1

    if flag_VDOWN and Extension == 1:
        # list_muons = [58, 91, 114, 119, 121, 128, 143, 168, 173, 192, 253, 277, 288, 302, 306, 309, 321, 349, 362, 383, 385, 
        #               391, 404, 434, 451, 468, 516, 550, 588, 690, 793, 831, 864, 929, 1090, 1091]

        list_muons = [119, 321, 385] # Whith PDF and plot
        # list_muons = [383]

        # list_muons = [119, 192, 321, 385, 550,1091] # With PDF and data behavior
        # list_muons = [119, 321, 385,550] 
        
        # list_muons = [119,#128,#168, #192, 302, 306, 321, $362, 383, $385, $404, #550, $793,$1090, #1091]  # latest
        # list_muons = [119,302, 306, 321,383]
        # list_muons = [362,404, 793, 1090]
        flag_option = 2

    if flag_HRIGHT1 and Extension == 1:
        # list_muons = [36, 57, 61, 152, 188, 254, 262, 265, 266, 271, 315, 316, 320, 348, 444, 454, 499, 500, 505, 573, 632, 
        #               637, 686, 703, 768, 771, 815, 818, 824, 836, 857, 909, 948, 956, 963]

        list_muons = [316, 505, 857, 948] # With PDF and plot
        # list_muons = [499]

        # list_muons = [316,505,637,857,948] # With e- (more accurate)

        # list_muons = [#36, 57, #61, #152, #188, #254, #262, #265, #266, #271, 315, #316, #320, #348, #444, #454, 499, #500, #505, #573, #632, 
        #               637, #686, #703, #768, #771, #815, #818, #824, #836, 857, #909, #948, #956, #963]
        # list_muons = [57,315,499,637,857]
        # list_muons = [315,499,637,857]
        flag_option = 3

    if flag_HRIGHT2 and Extension == 1:  
        # list_muons = [21, 73, 107, 110, 157, 199, 204, 280, 329, 373, 429, 448, 504, 587, 591, 604, 659, 738, 743, 758, 770, 
        #               783, 810, 841, 844, 855, 918, 920, 972, 974, 1014, 1061, 1084, 1114]

        # list_muons = [73, 157, 204, 280, 329, 429, 448, 504, 738, 770, 783, 918, 1014, 1114] # With PDF and plot
        list_muons = [73, 157, 204, 280, 329, 429, 504, 770, 783, 1014, 1114] # With PDF and plot
        # list_muons = [448]

        # list_muons = [73,110,204,280,329,448,770,] # With e- (more accurate)

        # list_muons = [73,110, 204,280,448,504,738,770,783,844,972,1084,1114]
        # list_muons = [1114, 1084,844,770,738,280,204,110,73] #lastest
        flag_option = 4

    if flag_HLEFT and Extension == 1:
        # list_muons = [1, 16, 34, 37, 54, 64, 162, 293, 317, 330, 333, 343, 347, 371, 415, 416, 456, 486, 493, 501, 507, 534, 
        #               545, 574, 577, 584, 594, 609, 617, 626, 701, 707, 716, 723, 726, 736, 745, 756, 759, 773, 775, 816, 819, 
        #               828, 832, 847, 867, 878, 885, 898, 914, 933, 937, 941, 946, 949, 961, 982, 993, 995, 997, 1002, 1016, 1021,
        #               1022, 1023, 1025, 1030, 1041, 1044, 1047, 1048, 1054, 1085, 1087, 1123, 1125]

        list_muons = [330, 371, 773, 914, 1002, 1023] # With PDF and plot
        # list_muons = [1054]

        # Muons with CONNIE spreads or smaller: 34,343, 416, 486, 493, 534, 594, 723, 736, 982, 997, 1030 ?????????

        # list_muons = [330, 371, 773, 1002, 1023]
        # list_muons = [330, 371, 1023, 1044] # Whit plot


        # list_muons = [16,54, 64,330,371,701,773, 775,819,898, 914,1002,1023, 1044]
        # list_muons = [64,330,371,701,775,819,898,914,1002,1023]
        # list_muons = [64,371]
        flag_option = 5

    if flag_VUP and Extension == 2:
        # list_muons = [65, 92, 99, 112, 130, 136, 145, 149, 182, 191, 358, 389, 449, 482, 545, 550, 678, 771, 795, 808, 820, 939,
        #               940, 947, 954, 956, 1002, 1133, 1164, 1166, 1216, 1337]

        # list_muons = [99, 112, 130, 136, 145, 149, 182, 389, 449, 550, 678, 820, 939, 940, 947, 1002, 1133, 1166, 1216, 1337]
        # list_muons = [99, 130, 145, 182, 389, 449, 939, 1002, 1133, 1166, 1216, 1337]
        list_muons = [99, 130, 145, 389, 449, 1002, 1133, 1166, 1216, 1337]
        # list_muons = [1337]

        # Muons with CONNIE spreads or smaller: 112, 550, 820
        flag_option = 6

    if flag_VDOWN and Extension == 2:
        # list_muons = [20, 64, 108, 128, 194, 232, 301, 349, 386, 450, 487, 593, 697, 706, 744, 793, 869, 873, 881, 903, 925, 1053,
        #               1091, 1100, 1141, 1145, 1160, 1167, 1175, 1276]
        # list_muons = [20, 64, 108, 128, 301, 386, 450, 487, 706, 744, 793, 1091, 1141, 1145, 1160, 1167, 1175, 1276]
        list_muons = [20, 64, 301, 793, 1167, 1175]
        # list_muons = [1276]
        
        # Muons with CONNIE spreads or smaller: 128,386, 487, 706, 1160
        flag_option = 7

    if flag_HRIGHT1 and Extension == 2:
        # list_muons = [37, 157, 315, 464, 513, 622, 632, 635, 636, 709, 727, 818, 819, 1011, 1221, 1245, 1252, 1334]
        # list_muons = [37, 157, 315, 464, 622, 709, 727, 818, 819, 1011, 1221, 1245, 1252, 1334]
        list_muons = [315, 464, 709, 727, 1011, 1245]
        # list_muons = [1334]

        # Muons with CONNIE spreads or smaller: 37, 819, 1011, 1221, 1334
        flag_option = 8

    if flag_HRIGHT2 and Extension == 2:
        # list_muons = [17, 29, 40, 113, 203, 258, 327, 344, 376, 395, 403, 790, 802, 833, 863, 890, 900, 975, 999, 1047, 1056, 1066, 
        #               1151, 1180, 1316]
        # list_muons = [29, 40, 203, 327, 344, 395, 403, 790, 833, 863, 890, 900, 1047, 1056, 1066, 1180]
        list_muons = [395, 790, 890, 1056, 1066]
        # list_muons = [1180]

        # Muons with CONNIE spreads or smaller: 29, 40, 327, 344, 403, 833, 863, 900, 1047
        flag_option = 9

    if flag_HLEFT and Extension == 2:
        # list_muons = [5, 43, 47, 50, 81, 96, 104, 170, 209, 219, 228, 241, 246, 252, 265, 285, 313, 347, 372, 379, 388, 399, 411, 433,
        #               442, 524, 565, 568, 701, 723, 766, 768, 785, 786, 821, 835, 838, 839, 841, 846, 859, 860, 879, 889, 899, 909, 921,
        #               945, 965, 986, 997, 1016, 1043, 1069, 1073, 1080, 1086, 1116, 1131, 1153, 1154, 1157, 1206, 1225, 1247, 1248, 1330]

        # list_muons = [5, 43, 47, 50, 81, 104, 170, 209, 228, 241, 246, 252, 347, 372, 379, 388, 399, 411, 433, 442, 524, 565, 701, 723, 766,
        #               785, 786, 835, 839, 846, 859, 899, 921, 945, 965, 997, 1043, 1069, 1131, 1154, 1206, 1225]

        # list_muons = [43, 47, 170, 228, 246, 252, 347, 379, 388, 411, 524, 701, 723, 839, 859, 945, 965, 1043, 1069, 1206]
        list_muons = [47, 170, 228, 246, 252, 347, 379, 388, 524, 701, 723, 839, 859, 965,1206]
        # list_muons = []
        
        # Muons with CONNIE spreads or smaller: 5,50, 81, 104, 209, 241, 372, 399, 433, 565, 766, 785, 786, 846, 921, 997, 1154, 1225
        flag_option = 10
    
    if flag_VUP or flag_VDOWN or flag_HRIGHT1 or flag_HRIGHT2 or flag_HLEFT:
        plot_array(Extension,flag_option,list_muons)

if __name__ == "__main__":
    main()
