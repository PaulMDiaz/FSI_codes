from dolfin import *

# class left(SubDomain):
#   def inside(self, x, on_boundary):
#     return x[0] > 1.0 - DOLFIN_EPS


# mesh_2 data array named "parent_vertex_indices" does not exist
mesh_1 = RectangleMesh(Point(0.0,0.0),Point(0.1,0.1),10,3)
mesh_2 = RectangleMesh(Point(0.0,0.2),Point(0.1,0.3),10,3)

V1 = VectorFunctionSpace(mesh_1, "CG", 1)
V2 = VectorFunctionSpace(mesh_2, "CG", 1)

dofs_V1 = V1.tabulate_dof_coordinates().reshape((V1.dim(),-1))
dofs_V2 = V2.tabulate_dof_coordinates().reshape((V2.dim(),-1))

# find indices corresponding to boundaries 1 and 2.
indices_1 = np.where((dofs_V1[:,1] <= 0.1 + DOLFIN_EPS)&(dofs_V1[:,1] >= 0.1 - DOLFIN_EPS))[0]
indices_2 = np.where((dofs_V2[:,1] <= 0.2 + DOLFIN_EPS)&(dofs_V2[:,1] >= 0.2 - DOLFIN_EPS))[0]

u1_border = u1.vector()[indices_1]
u2_border = u2.vector()[indices_2]

# Create some non-trivial function on mesh 1. Define expr with higher
# degree, so that, when checking u2 vs expr, we have some non-zero error
# (that converges as the mesh is refined). (This mainly serves to check that
# we're not just getting zero by integrating over a wrong/empty domain.)
expr = Expression(("cos(x[0])","exp(x[0])"),degree=3)
u1 = interpolate(expr,V1)

u1.set_allow_extrapolation(True)
# left = left()

def boundary1(x, on_boundary):
    return near(x[1],0.1)

def boundary2(x, on_boundary):
    return near(x[1],0.2)

def boundary3(x, on_boundary):
    return near(x[1],0.3)

bc3 = DirichletBC(V2, u1, boundary3)


# Subdomain marker
# mf = MeshFunction("size_t", mesh_2, 2)
# mf.set_all(1)
# left.mark(mf, 0)

# Define facet function over the new mesh_2
#ff = FacetFunction("size_t", mesh_2)
#Use MeshFunction<T>(mesh_2, mesh_2->topology().dim() - 1, value)
# maybe this:
# ff = MeshFunction("size_t", mesh_2, 1)
# ff.set_all(0)
# left.mark(ff, 1)

# Extract boundary mesh_2
bmesh = BoundaryMesh(mesh_2, "exterior", True)

# for x in bmesh.coordinates():
#     if left.inside(x, True):
#         # x[0] += 0.0
#         # x[1] *= 1.0 + 0.1*sin(0.1*t*x[1])
#         x[0] += 0.1
#         x[1] += 0.1

for x in bmesh.coordinates():
    # x[0] += 0.0
    # x[1] *= 1.0 + 0.1*sin(0.1*t*x[1])
    x[0] += 0.1
    x[1] += 0.1


# t = 0.0
# for i in range(10):
#   t += 0.1
#   for x in bmesh.coordinates():
#     if left.inside(x, True):
#       x[0] += 0.0
#       x[1] *= 1.0 + 0.1*sin(0.1*t*x[1])

ALE.move(mesh_2, bmesh)

u4 = Function(V2)
u4.vector().zero()
bc3.apply(u4.vector())

print "2 norm mesh 1 and mesh 2 velocities post ALE mesh movement:", np.linalg.norm(u1_border - u4_border)/np.linalg.norm(u1_border)
