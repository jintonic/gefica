import sys
if len(sys.argv) < 2:
    fin = "ppc.root"
else :
    fin = sys.argv[1]

from ROOT import TFile
f = TFile(fin)
t = f.pcdz.GetTree()
r = t.AsMatrix(['c1'])
z = t.AsMatrix(['c2'])
v = t.AsMatrix(['v'])
e = t.AsMatrix(['e'])

n1 = f.pcdz.GetN1() # number of grid points along r
n2 = f.pcdz.GetN2() # number of grid points along z
v.shape=(n2,n1) # mind the order!
e.shape=(n2,n1) # mind the order!

# does not work:
# https://root-forum.cern.ch/t/pydoublebuffer-seems-to-have-incorrect-size/28118/4
# import numpy as np
#r = np.array(tuple(f.pcdz.GetC1s()))
#z = np.array(tuple(f.pcdz.GetC2s()))
#v = np.array(tuple(f.pcdz.GetVs()))
# also doesn't work:
# https://stackoverflow.com/questions/23930671/how-to-create-n-dim-numpy-array-from-a-pointer
#r = np.ctypeslib.as_array(f.pcdz.GetC1s(),shape=(n2*n1,))
#z = np.ctypeslib.as_array(f.pcdz.GetC2s(),shape=(n2*n1,))
#v = np.ctypeslib.as_array(f.pcdz.GetVs(),shape=(n2,n1))
#print(v.shape)

from matplotlib import rcParams
rcParams['font.size'] = 12
rcParams['font.family'] = 'FreeSerif'
import matplotlib.pyplot as plt
import matplotlib.colors as clr
plt.figure(1) # potential
# imshow is better than pcolormesh if x, y, z are one to one mapped:
# https://stackoverflow.com/questions/24119920/how-to-plot-a-density-map-in-python
plt.imshow(v,
        extent=(r.min(), r.max(), z.min(), z.max()),
        origin='lower', # (0,0) is at bottom left corner
        cmap=plt.get_cmap('gnuplot'), # color map
        )
plt.ylabel('Axial position [cm]')
plt.xlabel('Radial position [cm]')
plt.colorbar().ax.set_ylabel('Potential [V]')

plt.figure(2) # E field
plt.imshow(e,
        extent=(r.min(), r.max(), z.min(), z.max()),
        origin='lower', # (0,0) is at bottom left corner
        cmap=plt.get_cmap('gnuplot'), # color map
        norm=clr.LogNorm(1) # log scale
        )
plt.ylabel('Axial position [cm]')
plt.xlabel('Radial position [cm]')
plt.colorbar().ax.set_ylabel('Electric field [V/cm]')

# the vector field is different from a few nicely spaced field lines
#e1 = t.AsMatrix(['e1'])
#e2 = t.AsMatrix(['e2'])
#plt.quiver(r,z,e1,e2)

plt.show()
