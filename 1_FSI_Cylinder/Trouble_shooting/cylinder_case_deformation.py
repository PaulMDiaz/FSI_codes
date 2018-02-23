"""
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for flow around a cylinder using the Incremental Pressure Correction
Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np
import scipy.io

# Start simulation from 8 s with medium mesh (projected onto desired mesh)
load_data = 0

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;
set_log_level(ERROR)

dt = 0.001
T = 0.002   # final time

nu = 0.001        # dynamic viscosity
rho = 1000.0         # density

#nu = mu/rho    # kinematic viscosity
mu = rho*nu     # dynamic viscosity

# Domain length and height
L = 2.5
H = 0.41

# mesh discretization
#N = 64

ufile = File("fluid_benchmark/velocity.pvd")
pfile = File("fluid_benchmark/pressure.pvd")

#U_max = 0.3
U_mean = 0.2

#mesh = Mesh('cylinderfluid.xml');
#facets_fluid = MeshFunction("size_t", mesh, "cylinderfluid_facet_region.xml")

mesh = Mesh('cylinderfluid_fine.xml');
facets_fluid = MeshFunction("size_t",mesh, "cylinderfluid_fine_facet_region.xml")
sub_domains_fluid = MeshFunction("size_t", mesh, "cylinderfluid_fine_physical_region.xml")


# mesh = Mesh('cylinderfluid_med.xml');
# facets_fluid = MeshFunction("size_t",mesh, "cylinderfluid_med_facet_region.xml")
# sub_domains_fluid = MeshFunction("size_t", mesh, "cylinderfluid_med_physical_region.xml")


# mesh = Mesh('cylinderfluid_course.xml');
# facets_fluid = MeshFunction("size_t",mesh, "cylinderfluid_course_facet_region.xml")
# sub_domains_fluid = MeshFunction("size_t", mesh, "cylinderfluid_course_physical_region.xml")

dx = dx(metadata={'quadrature_degree': 2})

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
T_space = TensorFunctionSpace(mesh, "Lagrange", 1)
S_space = FunctionSpace(mesh, "Lagrange", 1)
Q = FunctionSpace(mesh, "Lagrange", 1)

# space for mesh deformation
VFE = VectorElement("Lagrange",mesh.ufl_cell(),1)
VF = FunctionSpace(mesh,VFE)

# mess with mesh: move all nodes up one and one to right.
Fexpr = Expression(("1+x[0]","1+x[1]"),degree=2)
F = interpolate(Fexpr,VF)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u0 = Function(V)
u1  = Function(V)
p0 = Function(Q)
p1  = Function(Q)

meshVelocity = Function(VF)
meshVelocity.vector()[:] = F.vector().get_local()[:]

# dt = 1.0

meshDisplacement = Function(VF)
meshDisplacement.vector()[:] = meshVelocity.vector()[:]*dt

# Replace F with meshDisplacement (simulating what will happen in FSI problem)
F = meshDisplacement

# derivatives needed for change-of-variables
# deformation gradient tensor
DF = grad(F) + Identity(2)
J = det(DF)

if load_data == 1:
    # load saved data. Presently using a medium mesh with T = 8s.
    # Project onto mesh of this problem.
    # load data
    u_temp = np.loadtxt('nodal_u_med_t8')
    p_temp = np.loadtxt('nodal_p_med_t8')

    # load mesh
    mesh_load = Mesh('mesh_med.xml')

    # Define function spaces (P2-P1)
    V_load = VectorFunctionSpace(mesh_load, "Lagrange", 2)
    Q_load = FunctionSpace(mesh_load, "Lagrange", 1)

    # Define functions for solutions at previous and current time steps
    u0_load = Function(V_load)
    p0_load = Function(Q_load)

    u0_load.vector()[:] = u_temp
    p0_load.vector()[:] = p_temp

    u0_load.set_allow_extrapolation(True)
    #u0 = interpolate(u0_load, V)
    u0_load = project(u0_load, V)

    p0_load.set_allow_extrapolation(True)
    #p0 = interpolate(p0_load, Q, 'set_allow_extrapolation','true')
    p0_load = project(p0_load, Q)

    # Set to inital values. Needs to happen inside loop as well.
    u0 = u0_load
    p0 = p0_load


cells = CellFunction("size_t", mesh)

dx = Measure('dx', domain = mesh, subdomain_data = cells)
dx = dx(metadata={'quadrature_degree': 2})


ds = Measure('ds', domain = mesh, subdomain_data = facets_fluid)


def gradx(u):
    # (overkill for scalar-valued u)
    n = rank(u)
    ii = indices(n+2)
    # multivariate chain rule:
    # contract over last index of grad(u) and first index of
    # inv(DF) to change variables in derivative; should work for u of arbitrary
    # rank (scalar, vector, tensor)
    return as_tensor(grad(u)[ii[0:n+1]]*inv(DF)[ii[n],ii[n+1]],
                     ii[0:n]+(ii[n+1],))

def nabla_gradx(u):
    # (overkill for scalar-valued u)
    n = rank(u)
    ii = indices(n+2)
    # multivariate chain rule:
    # nabla_grad puts the new index First, in contrast to grad() which puts the
    # new index last. ie djdi as opposed to didj
    # not certain about this, need to test.
    # This idea is the same:
    # contract over last index of grad(u) and first index of
    # inv(DF) to change variables in derivative; should work for u of arbitrary
    # rank (scalar, vector, tensor)
    # effectively swapped what I think are the i and j from grad.
    return as_tensor(nabla_grad(u)[ii[1],ii[0]]*inv(DF)[ii[0],ii[n+1]],
                     (ii[1], ii[n+1]))

def divx(u):
    n = rank(u)
    ii = indices(n+2)
    # multivariate chain rule. yields a scalr so try without as_tensor
    return grad(u)[ii[0:n+1]]*inv(DF)[ii[1],ii[0]]

## load data
#u_temp = np.loadtxt('nodal_u_256_t3_8285')
#u0.vector()[:] = u_temp
#p_temp = np.loadtxt('nodal_p_256_t3_8285')
#p0.vector()[:] = p_temp

# Define inlet profile. 2 equivilent statements.
inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')

#inlet_profile = ('4*U_max*x[1]*(H-x[1])/pow(H, 2)', '0')

bcu_inlet = DirichletBC(V, Expression(inlet_profile, H = H, U_mean = U_mean, degree=2), facets_fluid, 16)
bcu_walls = DirichletBC(V, Constant((0, 0)), facets_fluid, 15)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), facets_fluid, 18)
bcu_fsi = DirichletBC(V, Constant((0, 0)), facets_fluid, 20)
bcp_outlet = DirichletBC(Q, Constant(0), facets_fluid, 17)


bcu = [bcu_inlet, bcu_walls, bcu_cylinder, bcu_fsi]
bcp = [bcp_outlet]

# may have to adjust n

# Define expressions used in variational forms
U  = 0.5*(u0 + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# adjust n, grad and ds...

def compute_forces(Mesh,mu,u,p, ds, J, DF):
	mesh = Mesh
	mu = mu
	# should include both pressure contribution and shear forces.
	# Face normals
	n = FacetNormal(mesh)
	# Stress tensor
	# Traction
	#sigma = nu*(grad(u)+grad(u).T) -p*Identity(2)
	sigma = mu*(gradx(u)+gradx(u).T) -p*Identity(2)

	# need to introduce the right normal
	# T = dot(sigma,n)
	# use first PK stress to transorm
	T = dot(J*sigma, inv(DF).T) * n

	# Traction is now in reference domain, so ds can be used as is.
	drag = -T[0]*ds(18)-T[0]*ds(20)
	lift = -T[1]*ds(18)-T[1]*ds(20)

	drag1 = assemble(drag); lift1 = assemble(lift);
	return drag1, lift1


# Define symmetric gradient
def epsilon(u):
    return sym(nabla_gradx(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u0) / k, v)*J*dx \
   + rho*dot(dot(u0, nabla_grad(u0)), v)*J*dx \
   + inner(sigma(U, p0), epsilon(v))*J*dx \
   + dot(p0*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*J*ds *inv(DF).T \
   - dot(f, v)*J*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_gradx(p), nabla_gradx(q))*J*dx
L2 = dot(nabla_gradx(p0), nabla_gradx(q))*J*dx - (1/k)*divx(u1)*q*J*dx

# Define variational problem for step 3
a3 = dot(u, v)*J*dx
L3 = dot(u1, v)*J*dx - k*dot(nabla_gradx(p1 - p0), v)*J*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)




# Create progress bar
#progress = Progress('Time-stepping')
#set_log_level(PROGRESS)


# initialize results matrix to store drag, lift and p_diff
# results = np.zeros((int(T/dt)+1,5))

# initialize results matrix to store drag, lift
results = np.zeros((int(T/dt)+1,2))


#results = np.loadtxt('cylinder_ipcs_nobc_256_0005_t3_8285')
# load results:
#cylinder_ipcs_nobc_256_0005_t3_829.mat'

# Time-stepping
t = 0
count = 0
# loaded data for this time and count
#t = 3.829
#count = 7658

# indices to compute traction.
dofs_f_S = S_space.tabulate_dof_coordinates().reshape((S_space.dim(),-1))

i_f_S_fsi1 = np.where((dofs_f_S[:,1] <= 0.21 + DOLFIN_EPS) & (dofs_f_S[:,1] >= 0.21 - DOLFIN_EPS))[0]
#i_f_S_fsi2 = np.where((dofs_f_S[:,1] == 0.19) & (dofs_f_S[:,0] >= 0.24) & (dofs_f_S[:,0] <= 0.6))[0]
i_f_S_fsi2 = np.where((dofs_f_S[:,1] <= 0.19 + DOLFIN_EPS) & (dofs_f_S[:,1] >= 0.19 - DOLFIN_EPS))[0]

#i_f_S_fsi3 = np.where((dofs_f_S[:,0] == 0.6) & (dofs_f_S[:,1] > 0.19) & (dofs_f_S[:,1] < 0.21))[0]
i_f_S_fsi3 = np.where((dofs_f_S[:,0] <= 0.6 + DOLFIN_EPS) & (dofs_f_S[:,0] >= 0.6 - DOLFIN_EPS)& (dofs_f_S[:,1] > 0.19) &(dofs_f_S[:,1] < 0.21 - DOLFIN_EPS))[0]

i_f_S_fsi1 = i_f_S_fsi1[dofs_f_S[i_f_S_fsi1][:,0].argsort()][::-1]
i_f_S_fsi2 = i_f_S_fsi2[dofs_f_S[i_f_S_fsi2][:,0].argsort()]
i_f_S_fsi3 = i_f_S_fsi3[dofs_f_S[i_f_S_fsi3][:,1].argsort()]

i_f_S_fsi = np.concatenate((i_f_S_fsi2,i_f_S_fsi3,i_f_S_fsi1),axis = 0)

traction_tensor = np.zeros((i_f_S_fsi.size,4,int(T/dt)+1))

f_normal_stresses = Function(S_space)
f_v=TestFunction(S_space)

[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]
[bc.apply(A3) for bc in bcu]

while t < T + DOLFIN_EPS:

    # Apply boundary conditions to matrices
    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u1.vector(), b1, 'bicgstab', 'hypre_amg')

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p1.vector(), b2, 'bicgstab', 'hypre_amg')

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    [bc.apply(b3) for bc in bcu]
    solve(A3, u1.vector(), b3, 'cg', 'sor')

    # deviatoric stress tnesor
    tau = mu*(gradx(u1) + gradx(u1).T)

    Dim = mesh.topology().dim()
    I = Identity(Dim)
    # Sigma
    sigma_FSI = project(-p1*I + tau, T_space, solver_type = "mumps", \
		form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )
    # Save to file
    ufile << u1
    pfile << p1

    # Update previous solution
    u0.assign(u1)
    p0.assign(p1)

    n_f = FacetNormal(mesh)


    # T_f = dot(sigma_FSI, n_f)
    # Apply first PK stress mapping
    T = dot(J*sigma_FSI, inv(DF).T) * n_f

    # Use reference ds now that mapping has been applied
    f_fx=(1/FacetArea(mesh))*f_v*T_f[0]*ds(20)
    f_fy=(1/FacetArea(mesh))*f_v*T_f[1]*ds(20)

    f_Traction_x=assemble(f_fx,tensor=f_normal_stresses.vector())[i_f_S_fsi]
    f_Traction_y=assemble(f_fy,tensor=f_normal_stresses.vector())[i_f_S_fsi]

    traction_tensor[:,:,count-1] = np.column_stack((f_Traction_x,f_Traction_y, f_Traction_x, f_Traction_x))

    drag, lift, = compute_forces(mesh, mu, u1, p1, ds, J, DF)

    # p_diff = p1(0.15, 0.2)-p1(0.25, 0.2)

    # C_D = 2/(pow(U_mean,2)*2*0.05)*drag
    # C_L = 2/(pow(U_mean,2)*2*0.05)*lift

    # results[count,:] = [p_diff, C_D, C_L, drag, lift]
    results[count,:] = [drag, lift]

    count += 1

    if count%500 == 0:
        scipy.io.savemat('cylinder_bar_fluid_med_1.mat', mdict={'results':results})

    print("t =", t)
    # print('u max:', u1.vector().array().max())
    # print('drag:',  drag)

scipy.io.savemat('cylinder_bar_fluid_med.mat', mdict={'results':results})
# scipy.io.savemat('traction_tensor_fluid.mat', mdict={'traction_tensor':traction_tensor})

# coordinates_fsi = dofs_f_S[i_f_S_fsi]
# scipy.io.savemat('coordinates_fluid.mat', mdict={'coordinates_fsi':coordinates_fsi})

#np.savetxt('cylinder_ipcs_bc_128_001_test', results)

nodal_values_u = u1.vector().get_local()
np.savetxt('nodal_u_med_t8', nodal_values_u)
nodal_values_p = p1.vector().get_local()
np.savetxt('nodal_p_med_t8', nodal_values_p)
File('mesh_med.xml')<< mesh
#File('mesh_256.xml')<<mesh
