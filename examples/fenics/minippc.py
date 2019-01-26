"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.
  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary
  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
from dolfin import *
from mshr import *

resolution = 100.0
step = 1.0/resolution #steps in cm

# Create mesh and define function space
domain = Rectangle(Point(0,0),Point(1,1))
mesh = generate_mesh(domain,resolution)
V = FunctionSpace(mesh, 'P', 1)

#~~~~~~~~~~Calculate Weighting Potential~~~~~~~~~~~~~#

#Define boundary condition
u_PC = Constant(1.0)

def boundary_PC(x, on_boundary):
  tol = 1E-14
  return on_boundary and near(x[1],0,tol) and x[0]>-0.00 and x[0]<0.10

bc_PC = DirichletBC(V,u_PC,boundary_PC)

u_NC = Constant(0.0)

def boundary_NC(x, on_boundary):
  tol = 1E-14
  return on_boundary and near(x[0],1,tol) or near(x[1],1,tol)

bc_NC = DirichletBC(V,u_NC,boundary_NC)

bcs = [bc_PC,bc_NC]

r = Expression('x[0]',degree=1)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
a = dot(grad(u), grad(v))*r*dx
L = f*v*r*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution and mesh
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.set_title('miniPPC Weighting Potential')
plt.ylabel('z (cm)')
plt.xlabel('r (cm)')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)
color_map = plt.cm.RdYlBu_r
p1 = plot(u, cmap=color_map)
#plot(mesh)
fig.colorbar(p1, ax=ax)
#p1.set_cmap('plasma')
#plt.show()
plt.savefig('wp.png')

#Save values in the plot to a file
wpfile = open('wp_data.dat','w')
i = 0
j = 0
while i<=1:
  while j<=1:
    #coolfile.write(f'{i:.1f}\t{j:.1f}\t{u(i,j):.1f}\n')
    wpfile.write('%.2f %.2f %.6f\n' % (i*10,j*10,u(i,j))) #r and z reported in mm for siggen
    j += step 
  i += step
  j = 0
wpfile.close()
#~~~~~~~~~~Calculate Electric~~~~~~~~~~~~~#

# Variable for bias voltage
bias = 500

# Define boundary condition
u_PC_e = Constant(0.0)

bc_PC_e = DirichletBC(V,u_PC_e,boundary_PC)

u_NC_e = Constant(bias)

bc_NC_e = DirichletBC(V,u_NC_e,boundary_NC)

bcs_e = [bc_PC_e,bc_NC_e]

# Define other variables
p = Constant(1.0E10)
n = Constant(0.9E10)
e = Constant(-1.602E-19)
eps = Constant(16.2*8.85E-14)
r = Expression('x[0]',degree=1)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V) 
f = Expression('((e/eps)*(n-p))*x[1]+((p*e)/eps)',degree=1, p=p, n=n, e=e, eps=eps)
# f = Constant((1.9E10*1.602E-19)/(16.2*8.85E-14))
a = dot(grad(u), grad(v))*r*dx
L = f*v*r*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs_e)

# Plot solution and mesh
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.set_title('miniPPC Electric Field')
plt.ylabel('z (cm)')
plt.xlabel('r (cm)')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)
E_abs = Function(V)
E_abs = project ( sqrt( dot(-1*grad(u), -1*grad(u)) ), V)
p2 = plot(E_abs)
#plot(mesh)
cbar2 = fig.colorbar(p2, ax=ax)
cbar2.ax.set_ylabel('Electric field (V/cm)')
cbar2.ax.yaxis.label.set_fontsize(15)
p2.set_cmap('plasma')
#plt.show()
plt.savefig('e_field.png')

# Find r and z components of the E_field
VS = VectorFunctionSpace(mesh,"CG",1)
E = Function(VS)
E = project(-grad(u),VS)
E_r, E_z = E.split(deepcopy=True)

#Save values in the plot to a file
effile = open('ef_data.dat','w')
i = 0
j = 0
while i<=1:
  while j<=1:
    #coolfile.write(f'{i:.1f}\t{j:.1f}\t{u(i,j):.1f}\n')
    #E_check = sqrt(E_r(i,j)**2 + E_z(i,j)**2)
    effile.write('%.2f %.2f %.6f %.6f %.6f %.6f\n' % (i*10,j*10,u(i,j),E_abs(i,j),E_r(i,j),E_z(i,j))) #r and z reported in mm for siggen
    j += step 
  i += step
  j = 0
  effile.write('\n')
effile.close()

#Save values of the potential to a file
#bpfile = open('bp_data.dat','w')
#i = 0
#j = 0
#while i<=1:
#  while j<=1:
#    #coolfile.write(f'{i:.1f}\t{j:.1f}\t{u(i,j):.1f}\n')
#    bpfile.write('%.3f %.3f %.6f\n' % (i,j,u(i,j)))
#    j += step 
#  i += step
#  j = 0
#bpfile.close()

#~~~~~~~~~Calculate the depletion voltage~~~~~~~~~~~~~#

#Read in the data!
import numpy as np
file1 = open('wp_data.dat')
dtype1 = np.dtype([('r','f2'),('z','f2'),('wp','f8')])
wp_array = np.loadtxt(file1, dtype=dtype1)
file1.close()

file2 = open('ef_data.dat')
dtype2 = np.dtype([('r','f2'),('z','f2'),('bp','f8'),('ef','f8'),('efr','f8'),('efz','f8')])
bp_array = np.loadtxt(file2, dtype=dtype2)
file2.close()

wp_values = wp_array['wp']
bp_values = bp_array['bp']

r_values = wp_array['r']
z_values = wp_array['z']

#Estimate depletion voltage
minv = bias
samples = 0
for ct in range(len(wp_array)):
  if(r_values[ct]<=1.2 and z_values[ct]<=0.2 and bp_values[ct]>0 and minv>bp_values[ct]/(1.0-wp_values[ct])):
    minv = bp_values[ct]/(1.0-wp_values[ct])
    samples += 1
  #if(r_values[ct]>=0.496 and r_values[ct]<=0.504 and z_values[ct]>=0.496 and z_values[ct]<=0.504):
  #  print("The e_field at the edge of the bulk (0.5,0.5) is: %.6f" % ef_values[ct])
depv = bias - minv
print("Estimated depletion voltage is %.8f volts." % depv)
print("Samples used in calculation: %i" % samples)

#Recreate ef plot from tabulated values
#import matplotlib.pyplot as pltt
#fig, ax = pltt.subplots()
#ax.set_title('miniPPC Electric Field')
#pltt.ylabel('z (cm)')
#pltt.xlabel('r (cm)')
#for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#             ax.get_xticklabels() + ax.get_yticklabels()):
#    item.set_fontsize(15)
#p3 = pltt.scatter(r_values, z_values, c=ef_values)
#cbar3 = fig.colorbar(p3, ax=ax)
#cbar3.ax.set_ylabel('Electric field (V/cm)')
#cbar3.ax.yaxis.label.set_fontsize(15)
#p3.set_cmap('plasma')
#pltt.show()


# Save solution to file in VTK format
#vtkfile = File('poisson/solution.pvd')
#vtkfile << u

# Compute error in L2 norm
#error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
#vertex_values_u_D = u_D.compute_vertex_values(mesh)
#vertex_values_u = u.compute_vertex_values(mesh)
#import numpy as np
#error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
#print('error_L2  =', error_L2)
#print('error_max =', error_max)

# Hold plot
# interactive()
