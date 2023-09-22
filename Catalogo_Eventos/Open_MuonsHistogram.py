import matplotlib.pyplot as plt
import pickle
import sys
import os

# value = os.fork()

def main(argObj):
    for path in argObj:
        histogram = open(path, 'rb')
        exampleObj = pickle.load(histogram)
        histogram.close()

        exampleObj.canvas.manager.set_window_title(path)
    # if value == 0:
    plt.show()

if __name__ == "__main__":
    argObj = sys.argv[1:]
    exitcode = main(argObj)
    exit(code = exitcode)
