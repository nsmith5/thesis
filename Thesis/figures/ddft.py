import sys
import matplotlib.pyplot as plt 
import numpy as np

ALPHA=0.01
N = 1000
xmin = 0
xmax = 5


def make_plot(axis, x, y):
    axis.plot(x, y)


def gaussian(x, x0=0.0, sigma=1.0):
    return (np.sqrt(1/(2 * np.pi * sigma ** 2)) * 
        np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)))


def main(path):
    # Make gaussian peaks 
    x = np.linspace(xmin, xmax, N)
    peaks = (gaussian(x, x0=x0, sigma=ALPHA) for x0 in range(6))
    y = sum(peaks, np.zeros(N))
   
    # Make figure
    fig, ax = plt.subplots()
    make_plot(ax, x, y)
    
    # Save figure
    fig.savefig(path)


if __name__ == '__main__':
    if not sys.argv[1]:
        raise Exception("Requires one pathname to save figure to")
    path = sys.argv[1]
    main(path)
