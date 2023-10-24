import matplotlib.pyplot as plt
import pickle as pck
import sys
import os

# value = os.fork()

def main(argObj):
    for path in argObj:
        # histogram = open(path, 'rb')
        # exampleObj = pickle.load(histogram)
        # histogram.close()
        data_am_241 = open(path, 'rb')
        dict_am_241 = pck.load(data_am_241)
        data_am_241.close()

        # exampleObj.canvas.manager.set_window_title(path)
        gamma_spectrum = []
        gamma_spectrum_ext2 = []
        gamma_spectrum_ext1_4 = []

        # for charge in dict_am_241['charge_All_extension']:
        #     if charge < 14000:
        #         gamma_spectrum.append(charge)

        # for charge in dict_am_241['charge_ext2']:
        #     if charge < 14000:
        #         gamma_spectrum_ext2.append(charge)

        for charge in dict_am_241['charge_ext1_4']:
            if charge < 14000:
                gamma_spectrum_ext1_4.append(charge)

        fig, axs = plt.subplots(figsize = [10,10])
        # axs.hist(gamma_spectrum, bins = 20000, color = 'k', label = 'All Extensions')

        # axs.hist(gamma_spectrum_ext1_4, bins = 20000, color = 'r', label = 'Extensions 1 and 4')
        axs.hist(gamma_spectrum_ext2, bins = 30000, color = 'b', label = 'Extension 2')

        fig.suptitle('Energy Spectrum of Am-241', size = 18, y=0.92)
        axs.set_xlim(0,14000)
        axs.set_xlabel('Energy (ADUs)', size = 18)
        axs.set_ylabel('Counts', size = 18)
        axs.legend()
        axs.grid()
        plt.show()


    # if value == 0:
    plt.show()

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
