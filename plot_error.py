from __future__ import print_function
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys

try:
    fname = sys.argv[1]
    nx = int(sys.argv[2])
    ny = int(sys.argv[3])
except:
    print("Usage: python plot_error.py filename nx ny")
    quit()

data = np.loadtxt(fname,unpack=True)
data = map(lambda x: np.reshape(x,(nx,ny)),data)
x,y,z = data[2],data[3],data[4]
z_interp = data[-1]
dz = z_interp - z

f, axes = plt.subplots(2,sharex=True)
ax1,ax2 = axes
trueplot=ax1.pcolor(x,y,z)
true_cbar = plt.colorbar(trueplot,ax=ax1)
true_cbar.set_label(r'$x^4 - y^4 + \frac{1}{2} (xy)^3$',rotation=90)
errorplot=ax2.pcolor(x,y,dz)
error_cbar = plt.colorbar(errorplot,ax=ax2)
error_cbar.set_label(r'$\Delta z$',rotation=90)
for ax in axes:
    ax.set(adjustable='box-forced',aspect='equal')
    ax.set_xlabel(r'$x$',fontsize=16)
    ax.set_ylabel(r'$y$',fontsize=16)
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)

plotname=format(reduce(lambda x,y: "{}.{}".format(x,y),
                       fname.split('.')[:-1])+'.')
plt.savefig('images/'+plotname+'pdf',
            bbox_inches='tight')
plt.savefig('images/'+plotname+'png',
            bbox_inches='tight')
#plt.show()
