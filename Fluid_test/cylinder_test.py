"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

# make this work for mesh velocity. Work out how.
# same thing for structure stres. Work out how.

from dolfin import *
import numpy as np
from mshr import *
import scipy.io


# Variables to generate files
pwd = './Results_cylinder/'
file_v_comp = File(pwd + 'u_f.pvd')
file_p_comp = File(pwd + 'p_f.pvd')

# Domain length and height
L = 2.5
H = 0.41

channel = Rectangle(Point(0, 0), Point(L, H))

cylinder = Circle(Point(0.2, 0.2), 0.05,50)

l = 0.35
h = 0.02

# mean inlet velocity
U_mean = 0.2

# Set parameter values
dt = 0.001
T = 0.005

nu_f = 0.001    # kinematic viscosity
rho_f = 1000

mu_f = nu_f*rho_f   #dynamic viscosity

# left end is fully attatched to the fixed cylinder...
# have to work out how to do that well.
structure = Rectangle(Point(0.6-l-0.01,0.19), Point(0.6, 0.19+h))

#structure = Rectangle(Point(0.6-l,0.19), Point(0.6, 0.19+h))

A_point = Point(0.6, 0.2)

# this will have to change some...
domain = channel - cylinder - structure

mesh = generate_mesh(domain, 64)

cell_markers = CellFunction("bool", mesh)
cell_markers.set_all(False)

# This maybe does something
for cell in cells(mesh):
    def set_cell(self, x, cell):
        # need a better if statement...
        if x[0] >= 0.2 + 0.05 and x[0] <= 0.6 and x[1] >= 0.19 and x[1] <= 0.19 + 0.02:
            cell_markers[cell] = True

    #if the cell's error indicator is greater than some criteria.
    # in this case just identify cells on structure and maybe sphere.
mesh = refine(mesh, cell_markers)


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

#v_ = Function(V, name = 'u')
#p_ = Function(Q, name = 'p')

## Create functions
#u0 = Function(V)
#u1 = Function(V)
#p1 = Function(Q)

# reintroduce as first line when mesh is moving.

# Define boundaries



inlet   = CompiledSubDomain('near(x[0], 0)')
outlet  = CompiledSubDomain('near(x[0], L)', L = L)
walls   = CompiledSubDomain('near(x[1], 0) || near(x[1], H)', H = H)

# structure is any object in interior of domain
fsi = CompiledSubDomain('on_boundary && x[0] > DOLFIN_EPS && x[0]< L - DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H)

# Define inflow profile
inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')

# Define boundary conditions
bcu_inlet = DirichletBC(V, Expression(inlet_profile, H = H, U_mean = U_mean, degree=2), inlet)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_fsi = DirichletBC(V, Constant((0, 0)), fsi)

# set the reference pressure of the outflow to have zero mean value.
bcp_outlet = DirichletBC(Q, Constant(0), outlet)

bcu = [bcu_inlet, bcu_walls, bcu_fsi]
bcp = [bcp_outlet]

#(1/k)*inner(u - u0, v)*dx + inner(grad(u0)*(u0-u_mesh), v)*dx + \
# Tentative velocity step

#Load previously saved solution
#nodal_values_u = np.loadtxt('nodal_u', dtype = float)
#u_n.vector()[:] = nodal_values_u
#nodal_values_p = np.loadtxt('nodal_p', dtype = float)
#p_n.vector()[:] = nodal_values_p

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)

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
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Velocity update
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Time-stepping
t = 0.0

# Create progress bar
#progress = Progress('Time-stepping')
#set_log_level(PROGRESS)

#dofs_f_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
#i_f_V_top = np.where((dofs_f_V[:,1] == H))[0]
#i_f_V_base = np.where((dofs_f_V[:,1] == h))[0]


dofs_f_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
# : -all rows. Column 0 is x, colum 1 is y
# indices for x and y
i_f_V_in = np.where((dofs_f_V[:,0] == 0.0))[0]

while t < T + DOLFIN_EPS:
    # change from ilu, amg_hypre, ilu to "default"

    print('time:', t)

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'ilu')
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'ilu')
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg')
    end()

    # Plot solution
    #plot(p1, title="Pressure", rescale=True)
    #plot(u1, title="Velocity", rescale=True)

    # Save to file
    #ufile << u1
    u_inlet = u_.vector()[i_f_V_in]
    print "fluid velocity at inlet = ", u_inlet

    #u_top = u1.vector()[i_f_V_top]
    #print "fluid velocity on top plate = ", u_top

    #u_base = u1.vector()[i_f_V_base]
    #np.savetxt('u_base', u_base)
    #u_base_coords = dofs_f_V[i_f_V_base]
    #np.savetxt('u_base_coords', u_base_coords)

    #print "fluid velocity at base = ", u_base

    # Move to next time step

    # Update previous solution
    u_n.assign(u_)
    file_v_comp << (u_n, t)
    p_n.assign(p_)
    file_p_comp << (p_n, t)


    t += dt

# Hold plot
#interactive()
# Record y velocity along mid y line. Record x velocity for x  = 0

#u_v_chorin = np.column_stack((np.linspace(0.0, 1.0, 101), np.zeros(101)))
#u_u_chorin = np.column_stack((np.linspace(0.0, 1.0, 101), np.zeros(101)))

#for i_x in range(len(u_v_chorin)):
	#u_v_chorin[i_x,1] = u1(u_v_chorin[i_x,0],0.5)[1]
	#u_u_chorin[i_x,1] = u1(0.5,u_u_chorin[i_x,0])[0]

#scipy.io.savemat('u_v_chorin.mat', mdict={'u_v_chorin':u_v_chorin})
#scipy.io.savemat('u_u_chorin.mat', mdict={'u_u_chorin':u_u_chorin})
#File('./Results_cylinder/' + 'p1.pvd') << (p_)
#File('./Results_cylinder/' + 'u1.pvd') << (u_)
