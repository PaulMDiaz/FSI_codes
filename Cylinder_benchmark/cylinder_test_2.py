"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

from dolfin import *
import numpy as np
from mshr import *
import scipy.io

# Validate flow past Cylinder
# Compare drag and lift on cylinder and the pressure difference between the front
# and back to Turek et al.
# Drag lower and upper bound [5.57 5.59]
# Lift lower and upper bound [0.0104 0.0110]
# p_diff lower and upper bound [0.1172 0.1176]
from dolfin import *
from dolfin_utils import meshconvert
##dolfin-convert cylinder.msh cylinder.xml
import os
os.system('dolfin-convert cylinder_2.msh cylinder_2.xml')


# Variables to generate files


pwd = './Results_cylinder/'
file_v_comp = File(pwd + 'u_f.pvd')
file_p_comp = File(pwd + 'p_f.pvd')

# Domain length and height
L = 2.2
H = 0.41

# mesh discretization
#N = 32

#channel = Rectangle(Point(0, 0), Point(L, H))

#cylinder = Circle(Point(0.2, 0.2), 0.05, N)

# Two equivilent ways to define velocity inlet parabola
# mean inlet velocity
U_mean = 0.2
# max velocity at inlet
U_m = 0.3


# Set parameter values
dt = 0.001
T = 8.000

nu_f = 0.001    # kinematic viscosity
rho_f = 1

mu_f = nu_f*rho_f   #dynamic viscosity

#geometry = channel - cylinder

#mesh = generate_mesh(geometry, N)

mesh = Mesh("cylinder.xml")
# refine mesh on rectangle that captures cylinder and recirculation region
#h_cell = mesh.hmin()
#h_cell = 0.1
#center1 = Point(0.2, 0.2)
#center2 = Point(0.25, 0.2)
#center3 = Point(0.30, 0.2)
#center4 = Point(0.35, 0.2)
#center5 = Point(0.40, 0.2)

#refine_line = AutoSubDomain(lambda x: x[0] >= 0.2 and x[0] <= 0.5 and x[1] == 0.2)

#cell_f = CellFunction("bool", mesh, False)

#h = mesh.hmin()
#center = Point(c_x, c_y)
#cell_f = CellFunction('bool', mesh, False)
#for cell in cells(mesh):
#    if cell.midpoint().distance(center1) < 0.07:
#        cell_f[cell] = True
#    if cell.midpoint().distance(center2) < 0.07:
#        cell_f[cell] = True
#    if cell.midpoint().distance(center3) < 0.07:
#        cell_f[cell] = True
#    if cell.midpoint().distance(center4) < 0.07:
#        cell_f[cell] = True
#    if cell.midpoint().distance(center5) < 0.07:
#        cell_f[cell] = True

#mesh = refine(mesh, cell_f)

## Load saved results
#nodal_values_u = np.loadtxt('nodal_u', dtype = float)
#u_n.vector()[:] = nodal_values_u
#nodal_values_p = np.loadtxt('nodal_p', dtype = float)
#p_n.vector()[:] = nodal_values_p

# refine mesh about cylinder
#center = Point(0.2, 0.2)

# minimum cell diameter
#h = mesh.hmin()
#cell_f = CellFunction("bool", mesh, False)
#for cell in cells(mesh):
#    if cell.midpoint().distance(center) < 0.05 + h:#
#        cell_f[cell] = True
#mesh = refine(mesh, cell_f)

#plot(mesh, interactive = True)

# Variables to generate files
pwd = './Results_Cylinder_benchmark/'
file_v_f = File(pwd + 'v_f.pvd')
file_p_f = File(pwd + 'p_f.pvd')

# Define boundaries

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




## Refine mesh on cylinder facets
#cell_markers = FacetFunction("bool", mesh, False)
#cylinder.mark(cell_markers, True)
#mesh = refine(mesh, cell_markers)


# Redefine cylinder...
#class Cylinder(SubDomain):
#    def inside(self,x,on_boundary):
#        return x[0] > DOLFIN_EPS and x[0] < (L - DOLFIN_EPS) and \
#        x[1] > DOLFIN_EPS and x[1] < (H - DOLFIN_EPS) and on_boundary
#cylinder = Cylinder()

# expect difficulties here if not earlier...

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


# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)



# Define inlet profile
#inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')
inlet_profile = ('4*U_m*x[1]*(H-x[1])/pow(H, 2)', '0')

# Define boundary conditions
bcu_inlet = DirichletBC(V, Expression(inlet_profile, H = H, U_m = U_m, degree=2), inlet)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)

#bcu_outlet = DirichletBC(V, Expression(inlet_profile, H = H, U_m = U_m, degree=2), outlet)

bcp_outlet = DirichletBC(Q, Constant(0), outlet)

bcu = [bcu_inlet, bcu_walls, bcu_cylinder]
#bcu = [bcu_inlet, bcu_walls, bcu_cylinder, bcu_outlet]
bcp = [bcp_outlet]






## Functions for results
u_res = Function(V, name = 'u')
p_res = Function(Q, name = 'p')


def compute_forces(mesh,nu, u, p, ds):
    # should include both pressure contribution and shear forces.
    # Face normals
    n = FacetNormal(mesh)
    # Stress tensor
    # Traction
    sigma = nu*(grad(u)+grad(u).T) -p*Identity(2)
    T = dot(sigma,n)
    drag = -T[0]*ds(3)
    lift = T[1]*ds(3)
    # if accounting for pressure only
    #drag = -p*n[0]*ds(1)
    #lift = p*n[1]*ds(1)
    drag = assemble(drag); lift = assemble(lift);
    return drag, lift

#(1/k)*inner(u - u0, v)*dx + inner(grad(u0)*(u0-u_mesh), v)*dx + \
# Tentative velocity step

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu_f = Constant(mu_f)
nu_f = Constant(nu_f)
rho_f = Constant(rho_f)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu_f*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho_f*dot((u - u_n) / k, v)*dx \
   + rho_f*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu_f*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx

 # doesn't have rho_f in it. 1 in this problem so okay.
#F1 = (1/k)*inner(u - u_n, v)*dx \
#    + inner(grad(u_n)*u_n, v)*dx \
#    + nu_f*inner(grad(u), grad(v))*dx \
#    - inner(f, v)*dx

a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx
#a2 = inner(grad(p), grad(q))*dx
#L2 = -(1/k)*div(u_)*q*dx

# Velocity update
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx
#a3 = inner(u, v)*dx
#L3 = inner(u_, v)*dx - k*inner(grad(p_), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
#[bc.apply(A1) for bc in bcu]
#[bc.apply(A2) for bc in bcp]

# Time-stepping
t = 0.0

# Tabulate dofs for velocity and pressure
dofs_f_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
dofs_f_Q = Q.tabulate_dof_coordinates().reshape((Q.dim(),-1))

# identify nodes for veloicty inlet
i_f_V_in = np.where((dofs_f_V[:,0] == 0.0))[0]

count = 0;
# initialize results matrix to store drag, lift and p_diff
results = np.zeros((int(T/dt)+1,3))

while t < T + DOLFIN_EPS:
    # change from ilu, amg_hypre, ilu to "default"

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    #[bc.apply(b1) for bc in bcu]
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'ilu')
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    #[bc.apply(b2) for bc in bcp]
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p_.vector()) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'ilu')
    end()

    # applying A3 BC3 is new. Not sure will work.

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u_.vector(), b3, 'cg')
    end()

    # Update previous solution
    u_n.assign(u_)
    file_v_comp << (u_n, t)
    p_n.assign(p_)
    file_p_comp << (p_n, t)


    print('time:', t)
    #print u_.vector().array()[i_f_V_in]

    drag, lift, = compute_forces(mesh, nu_f, u_, p_, ds)

    p_diff = p_(0.15, 0.2)-p_(0.25, 0.2)

    C_D = 2/(pow(U_mean,2)*2*0.05)*drag
    C_L = 2/(pow(U_mean,2)*2*0.05)*lift

    print "cylinder drag and lift coefficients = ", C_D, C_L

    u_res.assign(u_)
    file_v_f << (u_res, t)
    p_res.assign(p_)
    file_p_f << (p_res, t)

    results[count,:] = [p_diff, C_D, C_L]

    count += 1
    t += dt

#nodal_values_u = u_.vector().array()
#np.savetxt('nodal_u', nodal_values_u)
#nodal_values_p = p_.vector().array()
#np.savetxt('nodal_p', nodal_values_p)

scipy.io.savemat('results_32_1.mat', mdict={'results':results})
