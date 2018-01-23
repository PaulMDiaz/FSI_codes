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

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;
set_log_level(ERROR)

dt = 0.001
T = 15.00        # final time

mu = 0.001         # dynamic viscosity
rho = 1            # density

# Domain length and height
L = 2.5
H = 0.41


# mesh discretization
N = 64

ufile = File("fluid_benchmark/velocity.pvd")
pfile = File("fluid_benchmark/pressure.pvd")

# Create mesh
channel = Rectangle(Point(0, 0), Point(L, H))

cylinder = Circle(Point(0.2, 0.2), 0.05, N)
#U_max = 0.3
U_mean = 0.2

nu_f = 0.001    # kinematic viscosity
rho_f = 1

mu_f = nu_f*rho_f   # dynamic viscosity

geometry = channel - cylinder

mesh = generate_mesh(geometry, N)

h_cell = mesh.hmin()
refine_box = AutoSubDomain(lambda x: x[0] >= 0.1 - DOLFIN_EPS and x[0] <= 0.8 + DOLFIN_EPS and x[1] >= 0.1 -DOLFIN_EPS and x[1] <= 0.3+DOLFIN_EPS)

for i_mark in range(1):
	cell_markers_1 = CellFunction("bool", mesh, False)
	refine_box.mark(cell_markers_1,True)
	mesh = refine(mesh, cell_markers_1)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)


class Inlet(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

class Outlet(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > (L - DOLFIN_EPS) and on_boundary

class Walls(SubDomain):
    def inside(self,x,on_boundary):
        return 'near(x[1], 0) || near(x[1], 0.41)' and on_boundary

class Cylinder(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > DOLFIN_EPS and x[0] < (L - DOLFIN_EPS) and \
        x[1] > DOLFIN_EPS and x[1] < (H - DOLFIN_EPS) and on_boundary

# Compile subdomains
walls = Walls()
inlet = Inlet()
outlet = Outlet()
cylinder = Cylinder()

# Mark subdomains
facets = FacetFunction("size_t", mesh)
facets.set_all(0)
cells = CellFunction("size_t", mesh)

walls.mark(facets, 4)
inlet.mark(facets, 1)
outlet.mark(facets, 2)
cylinder.mark(facets, 3)

walls    = 'near(x[1], 0) || near(x[1], 0.41)'

dx = Measure('dx', domain = mesh, subdomain_data = cells)
ds = Measure('ds', domain = mesh, subdomain_data = facets)

# Define inlet profile. 2 equivilent statements.
inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')
#inlet_profile = ('4*U_max*x[1]*(H-x[1])/pow(H, 2)', '0')

bcu_inlet = DirichletBC(V, Expression(inlet_profile, H = H, U_mean = U_mean, degree=2), inlet)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)

bcp_outlet = DirichletBC(Q, Constant(0), outlet)

bcu = [bcu_inlet, bcu_walls, bcu_cylinder]
bcp = [bcp_outlet]

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u0 = Function(V)
u1  = Function(V)
p0 = Function(Q)
p1  = Function(Q)

## load data
#u_temp = np.loadtxt('nodal_u_256_t3_8285')
#u0.vector()[:] = u_temp
#p_temp = np.loadtxt('nodal_p_256_t3_8285')
#p0.vector()[:] = p_temp

# Define expressions used in variational forms
U  = 0.5*(u0 + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

def compute_forces(mesh,nu, u, p, ds):
    # should include both pressure contribution and shear forces.
    # Face normals
    n = FacetNormal(mesh)
    # Stress tensor
    # Traction
    sigma = nu*(grad(u)+grad(u).T) -p*Identity(2)
    T = dot(sigma,n)
    drag = -T[0]*ds(3)
    lift = -T[1]*ds(3)
    # if accounting for pressure only
    #drag = -p*n[0]*ds(1)
    #lift = p*n[1]*ds(1)
    drag = assemble(drag); lift = assemble(lift);
    return drag, lift

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
# Tentative velocity step
F1 = (1/k)*rho*inner(u - u0, v)*dx + rho*inner(grad(u0)*u0, v)*dx + \
     nu_f*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
#[bc.apply(A1) for bc in bcu]
#[bc.apply(A2) for bc in bcp]

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"


# initialize results matrix to store drag, lift and p_diff
results = np.zeros((int(T/dt)+1,5))

#results = np.loadtxt('cylinder_ipcs_nobc_256_0005_t3_8285')
# load results:
#cylinder_ipcs_nobc_256_0005_t3_829.mat'

# Time-stepping
t = 0
count = 0
# loaded data for this time and count
#t = 3.829
#count = 7658

while t < T + DOLFIN_EPS:

    # Update current time
    t += dt

    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "cg", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")

    # Save to file
    #ufile << u1
    #pfile << p1

    # Update previous solution
    u0.assign(u1)
    p0.assign(p1)

    drag, lift, = compute_forces(mesh, nu_f, u1, p1, ds)

    p_diff = p1(0.15, 0.2)-p1(0.25, 0.2)

    C_D = 2/(pow(U_mean,2)*2*0.05)*drag
    C_L = 2/(pow(U_mean,2)*2*0.05)*lift

    results[count,:] = [p_diff, C_D, C_L, drag, lift]

    count += 1

    if count%500 == 0:
        scipy.io.savemat('cylinder_chorin_64_001_check.mat', mdict={'results':results})

    print("t =", t)
    #print('u max:', u1.vector().array().max())

scipy.io.savemat('cylinder_chorin_64_001.mat', mdict={'results':results})

#np.savetxt('cylinder_chorin_32_01', results)

#nodal_values_u = u1.vector().get_local()
#np.savetxt('nodal_u_256_t3_8285', nodal_values_u)
#nodal_values_p = p1.vector().get_local()
#np.savetxt('nodal_p_256_t3_8285', nodal_values_p)
