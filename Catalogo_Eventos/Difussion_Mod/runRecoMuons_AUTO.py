import os

def run_automation(ext : int, option, muons : list):
    print(f"Starting processing for {len(muons)} muons...")

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

    flag_direction = basedir.split('_')[1][4]

    for i in muons:
        # Define paths
        input_xye = basedir + f"/muon{i}.xye"
        output_txt = basedir + f"/muon{i}_spreads.txt"

        # Build the command string exactly as it would appear in the terminal
        # cmdline = f"root -l -q 'RecoMuon2.C(\"{input_xye}\", \"{output_txt}\")'"
        cmdline = f"root -l -q 'RecoMuon{flag_direction}.C(\"{input_xye}\", \"{output_txt}\")'"
        
        print(f"--- Running: {cmdline}")
        
        # Execute in terminal
        os.system(cmdline)

    print("Process finished.")

def main():

    Extension = 2 # Choose the extension

    flag_VUP = False # Active just one type of muons
    flag_VDOWN = False
    flag_HRIGHT1 = False
    flag_HRIGHT2 = False
    flag_HLEFT = True

    if flag_VUP and Extension == 1:
        list_muons = [39, 80, 81, 96, 210, 214, 233, 261, 294, 353, 406, 414, 426, 459, 532, 653, 760, 765, 915, 1096, 1121]
        flag_option = 1

    if flag_VDOWN and Extension == 1:
        list_muons = [58, 91, 114, 119, 121, 128, 143, 168, 173, 192, 253, 277, 288, 302, 306, 309, 321, 349, 362, 383, 385, 
                      391, 404, 434, 451, 468, 516, 550, 588, 690, 793, 831, 864, 929, 1090, 1091]
        flag_option = 2

    if flag_HRIGHT1 and Extension == 1:
        list_muons = [36, 57, 61, 152, 188, 254, 262, 265, 266, 271, 315, 316, 320, 348, 444, 454, 499, 500, 505, 573, 632, 
                      637, 686, 703, 768, 771, 815, 818, 824, 836, 857, 909, 948, 956, 963]
        flag_option = 3

    if flag_HRIGHT2 and Extension == 1:  
        list_muons = [21, 73, 107, 110, 157, 199, 204, 280, 329, 373, 429, 448, 504, 587, 591, 604, 659, 738, 743, 758, 770, 
                      783, 810, 841, 844, 855, 918, 920, 972, 974, 1014, 1061, 1084, 1114]
        flag_option = 4

    if flag_HLEFT and Extension == 1:
        list_muons = [1, 16, 34, 37, 54, 64, 162, 293, 317, 330, 333, 343, 347, 371, 415, 416, 456, 486, 493, 501, 507, 534, 
                      545, 574, 577, 584, 594, 609, 617, 626, 701, 707, 716, 723, 726, 736, 745, 756, 759, 773, 775, 816, 819, 
                      828, 832, 847, 867, 878, 885, 898, 914, 933, 937, 941, 946, 949, 961, 982, 993, 995, 997, 1002, 1016, 1021,
                      1022, 1023, 1025, 1030, 1041, 1044, 1047, 1048, 1054, 1085, 1087, 1123, 1125]
        flag_option = 5

    if flag_VUP and Extension == 2:
        list_muons = [65, 92, 99, 112, 130, 136, 145, 149, 182, 191, 358, 389, 449, 482, 545, 550, 678, 771, 795, 808, 820, 939,
                      940, 947, 954, 956, 1002, 1133, 1164, 1166, 1216, 1337]
        flag_option = 6

    if flag_VDOWN and Extension == 2:
        list_muons = [20, 64, 108, 128, 194, 232, 301, 349, 386, 450, 487, 593, 697, 706, 744, 793, 869, 873, 881, 903, 925, 1053,
                      1091, 1100, 1141, 1145, 1160, 1167, 1175, 1276]
        flag_option = 7

    if flag_HRIGHT1 and Extension == 2:
        list_muons = [37, 157, 315, 464, 513, 622, 632, 635, 636, 709, 727, 818, 819, 1011, 1221, 1245, 1252, 1334]
        flag_option = 8

    if flag_HRIGHT2 and Extension == 2:
        list_muons = [17, 29, 40, 113, 203, 258, 327, 344, 376, 395, 403, 790, 802, 833, 863, 890, 900, 975, 999, 1047, 1056, 1066, 
                      1151, 1180, 1316]
        flag_option = 9

    if flag_HLEFT and Extension == 2:
        list_muons = [5, 43, 47, 50, 81, 96, 104, 170, 209, 219, 228, 241, 246, 252, 265, 285, 313, 347, 372, 379, 388, 399, 411, 433,
                      442, 524, 565, 568, 701, 723, 766, 768, 785, 786, 821, 835, 838, 839, 841, 846, 859, 860, 879, 889, 899, 909, 921,
                      945, 965, 986, 997, 1016, 1043, 1069, 1073, 1080, 1086, 1116, 1131, 1153, 1154, 1157, 1206, 1225, 1247, 1248, 1330]
        flag_option = 10


    run_automation(Extension, flag_option, list_muons)

if __name__ == "__main__":
    main()
