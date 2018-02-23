# Cylinder mesh moving test case.

# Things to test.

#1. Mesh velocity on fluid interface.
#2. Fluid stress on stucture interface.

from fenics import *
import numpy as np

# Load meshes.
meshStructure = Mesh('cylinderbar_med.xml');
facetsStructure = MeshFunction("size_t", meshStructure, "cylinderbar_med_facet_region.xml")

meshFluid = Mesh('cylinderfluid_med.xml');
facetsFluid = MeshFunction("size_t", meshFluid, "cylinderfluid_med_facet_region.xml")

meshRef = Mesh('cylinderfluid_med.xml');
facetsRef = MeshFunction("size_t", meshRef, "cylinderfluid_med_facet_region.xml")

# set up function spaces.
# For now, address only mesh velocity and fluid velocity.

fluidVectorSpace = VectorFunctionSpace(meshFluid, 'CG', 2)
u0 = Function(fluidVectorSpace) 	# Previous time step
fluidScalarSpace = FunctionSpace(meshFluid, 'CG', 1)
p0 = Function(fluidScalarSpace)

meshVectorSpace = VectorFunctionSpace(meshRef, 'CG', 1)
u_mesh = Function(meshVectorSpace)

# Load solution from first time step of FSI problem:

# load data
u_temp = np.loadtxt('nodal_u_1')
u0.vector()[:] = u_temp
m_temp = np.loadtxt('nodal_m_1')
u_mesh.vector()[:] = m_temp

# set u_mesh to zero to check things. Invalid scalar encountered...
u_mesh.vector().zero()

p_temp = np.loadtxt('nodal_p_1')
p0.vector()[:] = p_temp
# U_temp = np.loadtxt('nodal_Ustr_1')
# U0.vector()[:] = U_temp

################################################################################
### Construct interface meshes.
################################################################################

# Boundary definition common to mesh and fluid in undeformed domain. Try this for now.
# will have to be more exact when applying to structure.
class FluidInterface(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] >= 0.24898979 - DOLFIN_EPS and x[0] < 2.5 and x[1] > 0.0 and x[1] < 0.41

fluidInterface = FluidInterface()

#################
# mesh
#################
boundariesMesh = MeshFunction("size_t", meshRef, meshRef.topology().dim() - 1)
boundariesMesh.set_all(0)
fluidInterface.mark(boundariesMesh, 1)

# restrict mesh to exterior boundaries only:
boundaryMesh = BoundaryMesh(meshRef, 'exterior')
interfaceMesh = SubMesh(boundaryMesh, fluidInterface)
interfaceMeshVectorSpace = VectorFunctionSpace(interfaceMesh, "CG", 1)

# thirdInterfaceFunction = Function(thirdInterfaceVectorSpace)
# thirdInterfaceFunction.vector().zero()
meshInterfaceSolution = interpolate(u_mesh,interfaceMeshVectorSpace)

interfaceMeshCoords = interfaceMeshVectorSpace.tabulate_dof_coordinates().reshape((interfaceMeshVectorSpace.dim(),-1))

# Interpolate solution of mesh onto the space defined only on boundary.

#################
# fluid
#################

boundariesFluid = MeshFunction("size_t", meshFluid, meshFluid.topology().dim() - 1)
boundariesFluid.set_all(0)
fluidInterface.mark(boundariesFluid, 1)

# restrict mesh to exterior boundaries only:
boundaryFluid = BoundaryMesh(meshFluid, 'exterior')
interfaceFluid = SubMesh(boundaryFluid, fluidInterface)
interfaceFluidVectorSpace = VectorFunctionSpace(interfaceFluid, "CG", 1)
interfaceFluidCoords = interfaceFluidVectorSpace.tabulate_dof_coordinates().reshape((interfaceFluidVectorSpace.dim(),-1))

fluidInterfaceSolution = Function(interfaceFluidVectorSpace)
fluidInterfaceSolution.vector().zero()

# does this get called at each time step??

class FluidInterfaceValues(Expression):
    def eval(self, values, x):
        try:
            values[:] = fluidInterfaceSolution(x)
        except:
            values[:] = 0
    def value_shape(self):
        return (2,)

fluidInterfaceValues = FluidInterfaceValues(degree=2)

# move to place with other bcs
# bcMove = DirichletBC(interfaceFluidVectorSpace,fluidInterfaceValues, boundariesFluid, 1)

##############################################################################
### Solve fluid problem with deformed mesh and right boundary conditions. ####
##############################################################################


##############################################################################
### Problem setup ####
##############################################################################

# delete p_s. and self.

nu_f = 0.001	# Fluid kinematic viscosity
nu_s = 0.4	# Structure Poisson coefficient
mu_s = 0.5e6 # structure first lame constant (very small)
#mu_s = 0.5e9 # structure first lame constant (very small)

rho_f = 1000.0	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)
#rho_f = 0.001
rho_s = 1000.0

# Numerical parameters
dt = 0.001# Time step
#T = 10.00	#  Set final time for iteration
T = 0.001#  Set final time for iteration

mu_f = rho_f*nu_f # dynamic visosity

# Geometric parameters
L = 2.5 	#Channel length
H = 0.41	#Channel height
l = 0.35 	#Bar length
h = 0.02	#Bar height

U_mean = 0.2

fluidTensorSpace = TensorFunctionSpace(meshFluid, 'CG', 1)

# u0 = Function(fluidVectorSpace) 	# Previous time step, alrady loaded.
u1  = Function(fluidVectorSpace)	# Current time step
p1 = Function(fluidScalarSpace)	# Current time step
# p0  = Function(fluidScalarSpace)	# Previous time step

# Define trial and test functions
u = TrialFunction(fluidVectorSpace)
p = TrialFunction(fluidScalarSpace)
v = TestFunction(fluidVectorSpace)
q = TestFunction(fluidScalarSpace)

# Define coefficients
k = Constant(dt)
U  = 0.5*(u0 + u)
n  = FacetNormal(meshFluid)
mu_f = Constant(mu_f)
rho_f = Constant(rho_f)
f = Constant((0.0,0.0))

cells = MeshFunction("size_t", meshFluid, meshFluid.topology().dim())
facets = MeshFunction("size_t", meshFluid, meshFluid.topology().dim() - 1)
dx = Measure('dx', domain = meshFluid, subdomain_data = cells)
ds = Measure('ds', domain = meshFluid, subdomain_data = facets)

##############################################################################
### Boundary Conditions
##############################################################################

# Create mesh function and mark all 0

inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')
bcuInlet = DirichletBC(fluidVectorSpace, Expression(inlet_profile, H = H, U_mean = U_mean, degree=2), facetsFluid, 16)
bcuWalls = DirichletBC(fluidVectorSpace, Constant((0, 0)), facetsFluid, 15)
bcuCylinder = DirichletBC(fluidVectorSpace, Constant((0, 0)), facetsFluid, 18)

# want to apply meshSolver.u_mesh, but coordinates differ.

# have to adjust to appropriate function space.

# bcuFSI = DirichletBC(fluidVectorSpace, u_mesh, facetsFluid, 20)
# bcuFSI = DirichletBC(fluidVectorSpace, fluidInterfaceValues, fluidInterface)
bcuFSI = DirichletBC(fluidVectorSpace,fluidInterfaceValues, boundariesFluid, 1)

bcu = [bcuInlet, bcuWalls, bcuCylinder, bcuFSI]

# DirichletBC(secondVectorSpace,secondInterfaceValues, secondInterface)
#bcuFSI = DirichletBC(fluidSolver.fluidVectorSpace, Constant((0, 0)), facetsFluid, 20)

# Pressure
bcpOutlet = DirichletBC(fluidScalarSpace, Constant(0), facetsFluid, 17)

# Set up the boundary conditions
# bcu = [bcuInlet, bcuWalls, bcuCylinder, bcuFSI]
bcp = [bcpOutlet]

##############################################################################
### Variational Form ####
##############################################################################

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu_f*epsilon(u) - p*Identity(len(u))

F1 = rho_f*dot((u - u0) / k, v)*dx \
		+ rho_f*dot(dot((u0-u_mesh), nabla_grad(u0)), v)*dx \
		+ inner(sigma(U, p0), epsilon(v))*dx \
		+ dot(p0*n, v)*ds \
		- dot(mu_f*nabla_grad(U)*n, v)*ds \
		- dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure correction
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p0), nabla_grad(q))*dx \
		- (1.0/k)*div(u1)*q*dx
	#-div(u1)*q*dx
	#-(1.0/k)*div(u1)*q*dx

# Velocity correction
a3 = dot(u, v)*dx
L3 = dot(u1, v)*dx\
		- k*dot(nabla_grad(p1-p0), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters["std_out_all_processes"] = False;
parameters['krylov_solver']['nonzero_initial_guess'] = True

##############################################################################
### Record indices prior to updating mesh.
# Index of boundary should be achieved in a similar way
##############################################################################

dofs_f_V = fluidVectorSpace.tabulate_dof_coordinates().reshape((fluidVectorSpace.dim(),-1))
dofs_m_V = meshVectorSpace.tabulate_dof_coordinates().reshape((meshVectorSpace.dim(),-1))

i_f_V_fsi1 = np.where((dofs_f_V[:,1] == 0.21) & (dofs_f_V[:,0] >= 0.24) & (dofs_f_V[:,0] <= 0.6))[0]
i_m_V_fsi1 = np.where((dofs_m_V[:,1] == 0.21) & (dofs_m_V[:,0] >= 0.24) & (dofs_m_V[:,0] <= 0.6))[0]

i_f_V_fsi2 = np.where((dofs_f_V[:,1] == 0.19) & (dofs_f_V[:,0] >= 0.24) & (dofs_f_V[:,0] <= 0.6))[0]
i_m_V_fsi2 = np.where((dofs_m_V[:,1] == 0.19) & (dofs_m_V[:,0] >= 0.24) & (dofs_m_V[:,0] <= 0.6))[0]

i_f_V_fsi3 = np.where((dofs_f_V[:,0] == 0.6) & (dofs_f_V[:,1] >= 0.19+DOLFIN_EPS) & (dofs_f_V[:,1] <= 0.21 -DOLFIN_EPS))[0]
i_m_V_fsi3 = np.where((dofs_m_V[:,0] == 0.6) & (dofs_m_V[:,1] >= 0.19+DOLFIN_EPS) & (dofs_m_V[:,1] <= 0.21 -DOLFIN_EPS))[0]

i_f_V_fsi = np.concatenate((i_f_V_fsi2,i_f_V_fsi3,i_f_V_fsi1),axis = 0)
i_m_V_fsi = np.concatenate((i_m_V_fsi2,i_m_V_fsi3,i_m_V_fsi1),axis = 0)


i_f_V_fsi_com = []
for i_s in range(dofs_m_V[i_m_V_fsi].shape[0]/2):
	for i_f in range(dofs_f_V[i_f_V_fsi].shape[0]/2):
		if dofs_m_V[i_m_V_fsi[2*i_s]][0] == dofs_f_V[i_f_V_fsi[2*i_f]][0] and dofs_m_V[i_m_V_fsi[2*i_s]][1] == dofs_f_V[i_f_V_fsi[2*i_f]][1]:
			i_f_V_fsi_com.append(i_f_V_fsi[2*i_f])
			i_f_V_fsi_com.append(i_f_V_fsi[2*i_f+1])


##############################################################################
### Update mesh and step forward in Time step ####
##############################################################################

# The interface is marked properly. Or at least marked. These factets are remembered.
# The challenge is applying the correct values. Use expressin to address this some way or another.
# New function is by definition on the coordiantes. how does that not work??

# move fluid mesh
# step through indices of coordinates. For each one update the fluid mesh coordinates
# as itself plus dt x the mesh velocity chosen at the correct coordinates.

# currently fails when mesh moves. I have not allowed extrapolation to be true.
# ie u.set_allow_extrapolation(True)
# Want values to work properly

for i_coords in range(len(meshFluid.coordinates())):
	meshFluid.coordinates()[i_coords] += dt*u_mesh(meshRef.coordinates()[i_coords])

for i_coords in range(len(interfaceFluid.coordinates())):
	interfaceFluid.coordinates()[i_coords] += dt*meshInterfaceSolution(interfaceMesh.coordinates()[i_coords])

### now look at deformed space.

# If i compute coordinates - nothing is changed.
interfaceFluidCoordsMoved = interfaceFluidVectorSpace.tabulate_dof_coordinates().reshape((interfaceFluidVectorSpace.dim(),-1))

# If I compute functionspace on moved mesh and look at coordinates they have moved
# This is quite alot of work...

# Match velocity of mesh to velocity of fluid. This will
# have to be done in each time step.
# dt = 0.001

# does this get called at each time step??
fluidInterfaceSolution.vector()[:] = meshInterfaceSolution.vector().get_local()

# all set to apply boundar condition and check out solution



print ""
print ""
print "ENTERING FLUID SOLVER"
print ''
print ''

# This should occure in initialization....  but requires defined boundaries ...

[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]
[bc.apply(A3) for bc in bcu]

### solver.
I = Identity(2) # could make more generic.

# Compute tentative velocity step
#begin("Computing tentative velocity")

################### BREAKS HERE ############################
b1 = assemble(L1)
[bc.apply(b1) for bc in bcu]
solve(A1, u1.vector(), b1, 'bicgstab', 'hypre_amg')
end()

# Pressure correction
#begin("Computing pressure correction")
b2 = assemble(L2)
[bc.apply(b2) for bc in bcp]
solve(A2, p1.vector(), b2, 'bicgstab', 'hypre_amg')
end()

# Velocity correction
#begin("Computing velocity correction")
b3 = assemble(L3)
[bc.apply(b3) for bc in bcu]
solve(A3, u1.vector(), b3, 'cg', 'sor')

#print('u test 2:',u1(0.1,0.2))
tau = mu_f*(grad(u1) + grad(u1).T)

# Sigma
sigma_FSI = project(-p1*I + tau, fluidTensorSpace, solver_type = "mumps", \
    form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )
#sigma_FSI = project(-p1*I + tau, tensorSpace, solver_type = "mumps", \
    #form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
print ""
print ""
print "EXITING FLUID SOLVER"

##############################################################################
### How did it go?
##############################################################################

# compare mesh, mesh_boundary, fluid, fluid boundary



# i_m_V_fsi_com = []
# for i_s in range(dofs_s_V[i_s_V_fsi].shape[0]/2):
# 	for i_m in range(dofs_m_V[i_m_V_fsi].shape[0]/2):
# 		if dofs_s_V[i_s_V_fsi[2*i_s]][0] == dofs_m_V[i_m_V_fsi[2*i_m]][0] and dofs_s_V[i_s_V_fsi[2*i_s]][1] == dofs_m_V[i_m_V_fsi[2*i_m]][1]:
# 			i_m_V_fsi_com.append(i_m_V_fsi[2*i_m])
# 			i_m_V_fsi_com.append(i_m_V_fsi[2*i_m+1])

fluidInterfaceSolution.vector()

u1.vector()[i_f_V_fsi]
u_mesh.vector()[i_m_V_fsi]

# print "2 norm mesh and fluid velocities :", np.linalg.norm(u1.vector()[i_f_V_fsi_com]  - u_mesh.vector()[i_m_V_fsi])/np.linalg.norm(u1.vector()[i_f_V_fsi])
print "2 norm mesh and fluid velocities :", np.linalg.norm(u1.vector()[i_f_V_fsi_com]  - u_mesh.vector()[i_m_V_fsi])

# it's kinda small, but not ideal.
