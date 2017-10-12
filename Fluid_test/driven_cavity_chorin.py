"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

from dolfin import *
import numpy as np
import scipy.io


# Variables to generate files
pwd = './Results_Driven_Cavity_Ghia/'
file_v_comp = File(pwd + 'v_f.pvd')
file_p_comp = File(pwd + 'p_f.pvd')


# Load mesh from file
#mesh = UnitSquare(20,20)
mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0,1.0), 60, 48)
# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.5
T = 6
nu = 0.01

# Define boundary conditions
noslip  = DirichletBC(V, (0, 0), "x[0] < DOLFIN_EPS || x[0] > 1.0 - DOLFIN_EPS || x[1] < DOLFIN_EPS")
lid  = DirichletBC(V, (1,0), "x[1] > 1.0 - DOLFIN_EPS")
bcu = [noslip, lid]

# Functions for results
v_ = Function(V, name = 'u')
p_ = Function(Q, name = 'p')

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
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

# Create files for storing solution
#ufile = File("velocity.pvd")

# Time-stepping
t = 0.0
#p = Progress("Time-stepping")

dofs_f_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
i_f_V_top = np.where((dofs_f_V[:,1] == 1.0))[0]

while t < T + DOLFIN_EPS:
    # change from ilu, amg_hypre, ilu to "default"
    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "ilu")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    solve(A2, p1.vector(), b2, "gmres", "ilu")
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "ilu")
    end()

    # Plot solution
    #plot(p1, title="Pressure", rescale=True)
    #plot(u1, title="Velocity", rescale=True)

    # Save to file
    #ufile << u1

    u_top = u1.vector()[i_f_V_top]
    print "fluid velocity on top plate = ", u_top

    # Move to next time step
    u0.assign(u1)
    #p.update(t / T)

    v_.assign(u1)
    file_v_comp << (v_, t)
    p_.assign(p1)
    file_p_comp << (p_, t)

    t += dt

# Hold plot
#interactive()
# Record y velocity along mid y line. Record x velocity for x  = 0

u_v_chorin = np.column_stack((np.linspace(0.0, 1.0, 101), np.zeros(101)))
u_u_chorin = np.column_stack((np.linspace(0.0, 1.0, 101), np.zeros(101)))

for i_x in range(len(u_v_chorin)):
	u_v_chorin[i_x,1] = u1(u_v_chorin[i_x,0],0.5)[1]
	u_u_chorin[i_x,1] = u1(0.5,u_u_chorin[i_x,0])[0]

scipy.io.savemat('u_v_chorin.mat', mdict={'u_v_chorin':u_v_chorin})
scipy.io.savemat('u_u_chorin.mat', mdict={'u_u_chorin':u_u_chorin})
