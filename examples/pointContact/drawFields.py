from ROOT import TFile, TTree
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

f = TFile("ppc.root")
t = f.pcdz.GetTree()
r = t.AsMatrix(['c1'])
z = t.AsMatrix(['c2'])
v = t.AsMatrix(['v']) # fixme: v must be a 2D array

cmap = plt.get_cmap('PiYG')
levels = MaxNLocator(nbins=50).tick_values(v.min(), v.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
fig, ax = plt.subplots()
img = ax.pcolormesh(r,z,v, cmap=cmap, norm=norm)
fig.colorbar(img, ax=ax)

plt.ylabel('Axial position [cm]')
plt.xlabel('Radial position [cm]')
plt.show()
