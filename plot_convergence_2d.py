#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import scipy as sp
from scipy import optimize
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys

ORDERS = [2,4,6]
PLOT_NAME = "images/convergence_norm."
linewidth=3
markersize=12
fontsize=12
fontsize_labels=16
styles_data=['bo-','go-','ro-']
mpl.rcParams.update({'font.size': fontsize})

def f(dx,alpha,o):
    return alpha + o*dx

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("{} FILENAME".format(sys.argv[0]))
        quit()
    filename = sys.argv[1]

    data = np.loadtxt(filename)
    nx = data[...,0]
    ny = data[...,1]
    errors = data[...,2:].transpose()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    for i in range(errors.shape[0]):
        dx = 2*np.pi/(nx[i]-1)
        ax1.semilogy(nx,errors[i],styles_data[i],
                     lw=linewidth,ms=markersize,
                     label="order = {}".format(ORDERS[i]))
    ax1.legend()
    ax1.set_xlabel(r"nx",fontsize=fontsize_labels)
    ax1.set_ylabel(r'$|error|_{2}$',
                   fontsize=fontsize_labels)
    ax1.set_xticks(nx)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(nx)
    ax2.set_xticklabels(map(int,ny))
    ax2.set_xlabel(r"ny",fontsize=fontsize_labels)
    for extension in ["png","pdf"]:
        plt.savefig(PLOT_NAME+extension,
                    bbox_inches='tight')
