from ROOT import TFile
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

f = TFile("ppc.root")
t = f.pcdz.GetTree()
r = t.AsMatrix(['c1'])
z = t.AsMatrix(['c2'])
v = t.AsMatrix(['v'])

# https://root-forum.cern.ch/t/pydoublebuffer-seems-to-have-incorrect-size/28118/4
# still does not work:
#r = np.array(tuple(f.pcdz.GetC1s()))
#z = np.array(tuple(f.pcdz.GetC2s()))
#v = np.array(tuple(f.pcdz.GetVs()))

n1 = f.pcdz.GetN1() # number of grid points along r
n2 = f.pcdz.GetN2() # number of grid points along z
v.shape=(n2,n1) # mind the order!

# imshow is better than pcolormesh if x, y, z are one to one mapped:
# https://stackoverflow.com/questions/24119920/how-to-plot-a-density-map-in-python
plt.imshow(v,
        extent=(np.amin(r), np.amax(r), np.amin(z), np.amax(z)),
        origin='lower', # (0,0) is at bottom left corner
        cmap=plt.get_cmap('gnuplot'), # color map
        norm=LogNorm() # log scale
        )
plt.colorbar()
plt.ylabel('Axial position [cm]')
plt.xlabel('Radial position [cm]')
plt.show()
