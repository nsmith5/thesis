# h5topng.py
#
# Usage:
#   python h5topng.py <filename>
#
# Convert an HDF5 Simulation record <filename> to a mountain of png images
#

import sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

def main(filename):
    try:
        fid = h5.File(filename, 'r')
    except:
        print("Error opening file, exiting.")
        return 1
    for time in fid:
        group = fid[time]
        c = group['Concentration']
        n = group['Density']
        plt.imshow(n, cmap='gray')
        plt.clim(-1, 1)
        plt.imshow(c, cmap='coolwarm', alpha=0.5)
        plt.clim(0, 1)
        plt.colorbar()

        plt.savefig(time + ".png")
        print("Write figure " + time +".png")

        plt.clf()

    return 0

if __name__ == '__main__':
    main(sys.argv[1])
