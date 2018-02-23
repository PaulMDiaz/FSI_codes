from fenics import *
from mshr import *

import numpy as np

# define meshes
# reference mesh
mesh_ref = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),10,10)

# first and second mesh to move. Without and with ALE
mesh_1 = RectangleMesh(Point(0.0,2.0),Point(1.0,3.0),10,10)
mesh_2 = RectangleMesh(Point(0.0,2.0),Point(1.0,3.0),10,10)

# Set up a vector function space on each mesh
V1 = VectorFunctionSpace(mesh_ref, "CG", 1)
V2 = VectorFunctionSpace(mesh_1, "CG", 1)
V3 = VectorFunctionSpace(mesh_2, "CG", 1)

# Create some non-trivial function on mesh 1. Define expr with higher
# degree, so that, when checking u2 vs expr, we have some non-zero error
# (that converges as the mesh is refined). (This mainly serves to check that
# we're not just getting zero by integrating over a wrong/empty domain.)
expr = Expression(("cos(x[0])","exp(x[0])"),degree=3)
u1 = interpolate(expr,V1)

# allow extrapolation on that function
u1.set_allow_extrapolation(True)

# Want to match values of function along the boundary with DirichletBC

# def boundary1(x, on_boundary):
#     return near(x[1],1.0)

def boundary2(x, on_boundary):
    return near(x[1],2.0)

def boundary3(x, on_boundary):
    return near(x[1],3.0)

bc0 = DirichletBC(V2, u1, boundary2)
bc1 = DirichletBC(V2, u1, boundary3)
bc2 = DirichletBC(V3, u1, boundary3)
bc3 = DirichletBC(V3, u1, boundary2)

# apply BC on mesh 2
u2 = Function(V2)
u2.vector().zero()
bc0.apply(u2.vector())

# # this should be zero, since u1 should equal u2 exactly on the shared
# # boundary

# find dofs
dofs_V1 = V1.tabulate_dof_coordinates().reshape((V1.dim(),-1))
dofs_V2 = V2.tabulate_dof_coordinates().reshape((V2.dim(),-1))
dofs_V3 = V3.tabulate_dof_coordinates().reshape((V3.dim(),-1))

# find indices corresponding to boundaries 1 and 2.
indices_1 = np.where((dofs_V1[:,1] <= 1.0 + DOLFIN_EPS)&(dofs_V1[:,1] >= 1.0 - DOLFIN_EPS))[0]
indices_2 = np.where((dofs_V2[:,1] <= 2.0 + DOLFIN_EPS)&(dofs_V2[:,1] >= 2.0 - DOLFIN_EPS))[0]
indices_3 = np.where((dofs_V3[:,1] <= 2.0 + DOLFIN_EPS)&(dofs_V3[:,1] >= 2.0 - DOLFIN_EPS))[0]

u1_border = u1.vector()[indices_1]
u2_border = u2.vector()[indices_2]

print "2 norm mesh 1 and mesh 2 velocities prior to mesh movement:", np.linalg.norm(u1_border - u2_border)/np.linalg.norm(u1_border)
# applies correct values somehow. Isn't expressly given top of first mesh as instruction.

# ### Move mesh_1
# # if I write this here, it still reads updated coordinates. The dofs are not read though.
# a = mesh_1.coordinates()[:]
# # # Adjust mesh coordinates.
for i_coords in range(len(mesh_1.coordinates())):
	mesh_1.coordinates()[i_coords] += 1.0

u3 = Function(V2)
u3.vector().zero()
# bc1.apply(u3.vector()) # poor results
bc0.apply(u3.vector())

u3_border = u3.vector()[indices_2]

print "2 norm mesh 1 and mesh 2 velocities post manual mesh movement:", np.linalg.norm(u1_border - u3_border)/np.linalg.norm(u1_border)
# Struggles if mesh coordinates are changed.

# move mesh with ALE
# apparently dofs and entities ordering will be preserved, and in consequence, the facetfunctions, function spaces etc...

class interface_1(SubDomain):
  def inside(self, x, on_boundary):
    return x[1] <= 2.0 + DOLFIN_EPS

class interface_2(SubDomain):
  def inside(self, x, on_boundary):
    return x[1] <= 3.0 + DOLFIN_EPS

interface_1 = interface_1()
interface_2 = interface_2()

mf = MeshFunction("size_t", mesh_2, 2)
mf.set_all(0)
interface_1.mark(mf,1)
interface_2.mark(mf,2)

facets = MeshFunction("size_t", mesh_2,1)
facets.set_all(0)
interface_1.mark(facets,1)
interface_2.mark(facets,2)

# Extract boundary mesh
bmesh = BoundaryMesh(mesh_2, "exterior", True)

for x in bmesh.coordinates():
    if interface.inside(x, True):
        x[0] += 1.0
        x[1] += 1.0

# update mesh...
ALE.move(mesh_2, bmesh)

# bc3 = DirichletBC(V3, u1, facets, 1) # puts out nothing
bc3 = DirichletBC(V3, u1, boundary3) # puts out something kinda interesting.
# only select values have elmements. Don't appear to be the right ones...

# bc3 = DirichletBC(V3, u1, facets, 2) # puts out wrong something kinda interesting.
# all values are filled though. So wrong.
# bc3 = DirichletBC(V3, u1, facets, 1) # puts out nothing

u4 = Function(V3)
u4.vector().zero()
# bc2.apply(u4.vector()) # achieves zeros no facets matching..
bc3.apply(u4.vector()) # achieves zeros no facets matching..

u4_border = u4.vector()[indices_3]

# I have some kind of dastadly error that needs a fix real bad.. really really bad. hmm.

print "2 norm mesh 1 and mesh 2 velocities post ALE mesh movement:", np.linalg.norm(u1_border - u4_border)/np.linalg.norm(u1_border)
print "border values", u4_border
