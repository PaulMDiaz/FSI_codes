from fenics import *
import numpy as np

# reference mesh that does not change in time.
mesh_ref = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),10,10)

# two test meshes, both offset from the initial mesh.
mesh_1 = RectangleMesh(Point(0.0,2.0),Point(1.0,3.0),10,10)
mesh_2 = RectangleMesh(Point(0.0,2.0),Point(1.0,3.0),10,10)

# Set up a vector function space on each mesh
V_ref = VectorFunctionSpace(mesh_ref, "CG", 1)
V1 = VectorFunctionSpace(mesh_1, "CG", 1)
V2 = VectorFunctionSpace(mesh_2, "CG", 1)

# simulate some solution on reference mesh
expr = Expression(("cos(x[0])","exp(x[0])"),degree=3)
u_ref = interpolate(expr,V_ref)

# allow extrapolation on that function
u_ref.set_allow_extrapolation(True)

# Want to match values of solution along the top boundary of the reference mesh to
# lower boundary of the separate mesh.

# itdentify boundary on mesh_1
def boundary2(x, on_boundary):
    return near(x[1],2.0)

bc0 = DirichletBC(V1, u_ref, boundary2)

################################################################################
########  Apply BC to separate mesh ##########
################################################################################

# Set up function on mesh_1 and apply BC.
u2 = Function(V1)
u2.vector().zero()
bc0.apply(u2.vector())

# Check error through identifiying nodes on boundary.

# find dofs
dofs_V1 = V_ref.tabulate_dof_coordinates().reshape((V_ref.dim(),-1))
dofs_V2 = V1.tabulate_dof_coordinates().reshape((V1.dim(),-1))
dofs_V3 = V2.tabulate_dof_coordinates().reshape((V2.dim(),-1))

# find indices corresponding to boundaries 1 and 2.
indices_1 = np.where((dofs_V1[:,1] <= 1.0 + DOLFIN_EPS)&(dofs_V1[:,1] >= 1.0 - DOLFIN_EPS))[0]
indices_2 = np.where((dofs_V2[:,1] <= 2.0 + DOLFIN_EPS)&(dofs_V2[:,1] >= 2.0 - DOLFIN_EPS))[0]
indices_3 = np.where((dofs_V3[:,1] <= 2.0 + DOLFIN_EPS)&(dofs_V3[:,1] >= 2.0 - DOLFIN_EPS))[0]

u1_border = u_ref.vector()[indices_1]
u2_border = u2.vector()[indices_2]

print "2 norm mesh 1 and mesh 2 velocities prior to mesh movement:", np.linalg.norm(u1_border - u2_border)/np.linalg.norm(u1_border)
# machine precision achieved.

################################################################################
########  Apply BC to separate mesh with coordinates moved ##########
################################################################################

# Move mesh_1 vertically by 1
for i_coords in range(len(mesh_1.coordinates())):
	mesh_1.coordinates()[i_coords] += [0.0, 1.0]

u3 = Function(V1)
u3.vector().zero()
# bc1.apply(u3.vector()) # poor results
bc0.apply(u3.vector())

u3_border = u3.vector()[indices_2]

print "2 norm mesh 1 and mesh 2 velocities post manual mesh movement:", np.linalg.norm(u1_border - u3_border)/np.linalg.norm(u1_border)

# Works perfectly when mesh is moved only in y.
# Fails when mesh moves in x.

################################################################################
########  Apply BC to separate mesh with coordinates moved in y using ALE ######
################################################################################
# move mesh with ALE

class Interface_1(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[1], 2.0)

interface_1 = Interface_1()

mf = MeshFunction("size_t", mesh_2, 2)
mf.set_all(0)
interface_1.mark(mf,1)

facets = MeshFunction("size_t", mesh_2,1)
facets.set_all(0)
interface_1.mark(facets,1)

# Extract boundary mesh
bmesh = BoundaryMesh(mesh_2, "exterior", True)
for x in bmesh.coordinates():
    # if interface.inside(x, True):
    x[0] += 1.0
    x[1] += 0.0

# update mesh...
ALE.move(mesh_2, bmesh)

bc3 = DirichletBC(V2, u_ref, facets, 1)

u4 = Function(V2)
u4.vector().zero()
bc3.apply(u4.vector()) # achieves zeros no facets matching..
u4_border = u4.vector()[indices_3]

print "2 norm mesh 1 and mesh 2 velocities post ALE y mesh movement:", np.linalg.norm(u1_border - u4_border)/np.linalg.norm(u1_border)

# works perfectly when mesh is only moved in y
# fails when there is an x component.
# further, moving with boundary does not fit FSI implementation so well. Want to be
# consistent with fluid mesh veloicty that is applied to fluid. 

################################################################################
########  ALE with x moved too ######
################################################################################
# y is already moved from previous step

for x in bmesh.coordinates():
    # if interface.inside(x, True):
    x[0] += 1.0
    x[1] += 0.0

ALE.move(mesh_2, bmesh)

bc3 = DirichletBC(V2, u_ref, facets, 1)

u4 = Function(V2)
u4.vector().zero()
bc3.apply(u4.vector()) # achieves zeros no facets matching..
u4_border = u4.vector()[indices_3]

print "2 norm mesh 1 and mesh 2 velocities post ALE x and y mesh movement:", np.linalg.norm(u1_border - u4_border)/np.linalg.norm(u1_border)
