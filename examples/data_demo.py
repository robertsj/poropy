# examples/data_demo.py

# This demonstrates use of the simple data model.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

from poropy.nucleardata import *

def main():

    # All data corresponds to a WH 17x17.
    
    IFBA = 0 # No IFBA
    WABA = 0 # No WABA
    GAD  = 0 # No GAD
    enrichment = 4.5 # w/o U-235
    
    # We want to get some data as a function of burnup to
    # plot.  
    number_points = 10
    burnup = np.linspace(0,50,number_points)
    D1 = np.zeros(number_points)
    D2 = np.zeros(number_points)
    A1 = np.zeros(number_points)
    A2 = np.zeros(number_points)
    F1 = np.zeros(number_points)
    F2 = np.zeros(number_points)
    S12= np.zeros(number_points)
    for b in range(0, number_points):
        print burnup[b]
        data = get_2g_parms(burnup[b], enrichment, IFBA, WABA, GAD)
        D1[b] = data.DIFF1
        D2[b] = data.DIFF2
        A1[b] = data.ABS1
        A2[b] = data.ABS2
        F1[b] = data.NUFISS1
        F2[b] = data.NUFISS2
        S12[b]= data.REMOV1 - data.ABS1


    # Plot the diffusion coefficients
    plt.plot( burnup, D1, 'b',  \
              burnup, D2, 'g--',\
              lw=2) 
    plt.title('Group Diffusion Coefficients')
    plt.xlabel('burnup [MWd/kg]')
    plt.ylabel('D [cm]')
    plt.legend(('group 1', 'group 2'),loc=5, shadow=True)
    plt.grid(True)
    plt.show()
    


if __name__ == "__main__":
    main()
