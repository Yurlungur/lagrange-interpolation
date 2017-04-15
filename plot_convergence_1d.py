#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys

NX = [11,21,41]
PLOT_NAME = "images/convergence_slice."
linewidth=3
fontsize=12
fontsize_labels=16
mpl.rcParams.update({'font.size': fontsize})

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("{} ORDER_1D FILENAME".format(sys.argv[0]))
        quit()
    order = float(sys.argv[1])
    filename = sys.argv[2]

    data = np.loadtxt(filename)
    x = data[...,0]
    errors = data[...,1:].transpose()
    dx = map(lambda r: (np.max(x)-np.min(x))/(r-1),NX)
    for i in range(errors.shape[0]):
        e = errors[i]
        plt.plot(x,e/(dx[i]**order),
                 label='nx = {}'.format(NX[i]),
                 lw=linewidth)
    plt.legend()
    plt.xlabel('x',fontsize=fontsize_labels)
    plt.ylabel(r'error$/dx^{%f}$' % order,
               fontsize=fontsize_labels)
    for extension in ["png","pdf"]:
        plt.savefig(PLOT_NAME+extension,
                    bbox_inches='tight')

