import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.ndimage.filters import gaussian_filter1d
from matplotlib.style import use

def load_and_split(filename):
    data = np.loadtxt(filename)
    liq = data[data[:, 2] == 0]     # Liquid droplets
    sol = data[data[:, 2] == 1]     # Solid droplets
    return liq, sol

def means_and_std(droplets):
    times = np.unique(droplets[:, 0])
    output = np.empty((len(times), 3))  # time | mean size | std size
    output[:, 0] = times                # set time column
    
    for i, t in enumerate(times):
        droplets_t = droplets[droplets[:, 0] == t]
        output[i, 1] = np.mean(droplets_t)
        output[i, 2] = np.std(droplets_t)
    
    return output

def smooth(stats):
    gaussian_filter1d(stats[:, 1], 15, output=stats[:, 1])    

def make_errorbars(liquids, solids):
    plt.errorbar(solids[:, 0], 
                 solids[:, 1], 
                 yerr=solids[:, 2], 
                 fmt='o',
                 label="Solid")
    plt.errorbar(liquids[:, 0],
                 liquids[:, 1],
                 yerr=liquids[:, 2],
                 fmt='o',
                 label="Liquid")
    plt.show()

def make_loglog(liquids, solids):
    use('bmh')
    plt.loglog(solids[:, 0], 
               np.sqrt(solids[:, 1]),
               label="Solid")
    plt.loglog(liquids[:, 0],
               np.sqrt(liquids[:, 1]),
               label="Liquid")
    plt.xlabel("Time (dt)")
    plt.ylabel("Average droplet radius (R)")
    plt.title("Droplet growth rates")
    plt.legend()
    plt.show()

def make_plot(liquids, solids):
    plt.plot(solids[:, 0], 
             np.sqrt(solids[:, 1]),
             label="Solid")
    plt.plot(liquids[:, 0],
             np.sqrt(liquids[:, 1]),
             label="Liquid")
    plt.xlabel("Time (dt)")
    plt.ylabel(r'Average droplet radius $\langle R \rangle$')
    plt.title("Droplet growth rates")
    plt.legend()
    plt.show()

def make_slope_plot(liquids, solids): 
    use('bmh')
    sol_logtime = np.log(solids[:, 0])
    sol_slope = np.gradient(np.log(np.sqrt(solids[:, 1])), np.ediff1d(sol_logtime, to_end=100))
    plt.plot(sol_logtime, sol_slope, label="Solid")
    liq_logtime = np.log(liquids[:, 0])
    liq_slope = np.gradient(np.log(np.sqrt(liquids[:, 1])), np.ediff1d(liq_logtime, to_end=100))
    plt.plot(liq_logtime, liq_slope, label="Liquid")
    plt.legend()
    plt.xlabel(r'$\log(t)$', fontsize=18)
    plt.ylabel(r'$\frac{d \log(\langle R \rangle)}{d \log(t)}$', fontsize=18)
    plt.title("Growth exponent of droplets")
    plt.show()

def main(filename):
    liquids, solids = load_and_split(filename)
    liqstats = means_and_std(liquids)
    solstats = means_and_std(solids)
    #smooth(solstats)
    #smooth(liqstats)
    #make_loglog(liqstats, solstats)
    make_slope_plot(liqstats, solstats)
    make_plot(liqstats, solstats)
    #slope, intercept, r_value, p_value, std_err = linregress(np.log(solstats[:, 0]), np.log(np.sqrt(solstats[:, 1])))
    #print slope, r_value, std_err
    return

if __name__ == '__main__':
    main(sys.argv[1])
