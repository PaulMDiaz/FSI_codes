from fenics import *
from mshr import *
import numpy as np

# Two adjacent meshes that have coincident nodes on shared bouandary:
box_1 = Rectangle(Point(0, 0), Point(0.1, 0.1))
box_2 = Rectangle(Point(0, 0.1), Point(0.1, 0.2))

# Generate mesh
mesh_1 = generate_mesh(box_1, 4)
mesh_2 = generate_mesh(box_2, 4)

# Set up a vector function space on each mesh
Vector_space_1 = VectorFunctionSpace(mesh_1, "CG", 1)
Vector_space_2 = VectorFunctionSpace(mesh_2, "CG", 1)


# find the indices that correspond to matching nodes on the shared boundary.
# Independent meshes are designed so that these nodes do match.

# degree of freedom coordinates
dofs_1 = Vector_space_1.tabulate_dof_coordinates().reshape((Vector_space_1.dim(),-1))
dofs_2 = Vector_space_2.tabulate_dof_coordinates().reshape((Vector_space_2.dim(),-1))

# Locate coordinates on shared boundary.
indices_1 = np.where((dofs_1[:,1] <= 0.1 + DOLFIN_EPS) & (dofs_1[:,1] >= 0.1 - DOLFIN_EPS))[0]
indices_2 = np.where((dofs_2[:,1] <= 0.1 + DOLFIN_EPS) & (dofs_2[:,1] >= 0.1 - DOLFIN_EPS))[0]

# Treat x and y components seperately
indices_1_x = indices_1[::2]
indices_1_y = indices_1[1::2]
indices_2_x = indices_2[::2]
indices_2_y = indices_2[1::2]

# reorder indices
indices_1_x = indices_1_x[dofs_1[indices_1_x][:,0].argsort()]
indices_1_y = indices_1_y[dofs_1[indices_1_y][:,0].argsort()]
indices_2_x = indices_2_x[dofs_2[indices_2_x][:,0].argsort()]
indices_2_y = indices_2_y[dofs_2[indices_2_y][:,0].argsort()]

# concatenate
indices_1 = np.concatenate((indices_1_x,indices_1_y),axis = 0)
indices_2 = np.concatenate((indices_2_x,indices_2_y),axis = 0)

# Given some solution on mesh1/ vector_function_space_1

v_expression = Expression(("sin(pi*x[0]) + sin(pi*x[1])","sin(pi*x[0]) + sin(pi*x[1])"), degree = 1)
v_1 = interpolate(v_expression, Vector_space_1)

# allow extrapolation on that function
v_1.set_allow_extrapolation(True)
bc3.apply(fluid_solver.v_mesh.vector())

# use u1 as a BC on mesh 2
def boundary2(x, on_boundary):
    return near(x[1],0.1)
bc2 = DirichletBC(Vector_space_2, v_1, boundary2)

def boundary3(x, on_boundary):
    return near(x[1],0.19)
bc3 = DirichletBC(fluid_solver.V_space_mesh, v12, boundary3)

fluid_solver.v_mesh.vector().zero()

# apply BC on mesh 2
v_2 = Function(Vector_space_2)
v_2.vector().zero()
bc2.apply(v_2.vector())


# Project solution from mesh 1 onto mesh 2
# v_2 = project(v_1, Vector_space_1 , solver_type = "mumps",\
	# form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )

# Retrieve values from each mesh that lie on boundary.
#Indices go from left to right in x and then y. #dofs_1[indices_1]

v_1_boundary = v_1.vector()[indices_1]
v_2_boundary = v_2.vector()[indices_2]

# Compare values on boundary:
print "relative error boundary values v_1 and v_2", np.linalg.norm(v_1_boundary - v_2_boundary)/np.linalg.norm(v_1_boundary)
