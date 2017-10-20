"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

# make this work for mesh velocity. Work out how.
# same thing for structure stres. Work out how.

from dolfin import *
import numpy as np
import scipy.io


# Variables to generate files
pwd = './Results_Driven_Cavity_Ghia/'
file_v_comp = File(pwd + 'v_f.pvd')
file_p_comp = File(pwd + 'p_f.pvd')

H = 2.0
h = 0.5
W = 2.0

N_w = 64
N_h = 48
# Load mesh from file
#mesh = UnitSquare(20,20)
mesh = RectangleMesh(Point(0.0, h), Point(W,H), N_w, N_h)
# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.1
T = 0.0
nu = 0.2

# top boundary conditions
class FSS(Expression):
    def eval(self,values, x):
        if between(x[0], (0.0, 0.3)) and near(x[1], 2):
            values[0] = 0.5*sin(pi*x[0]/0.6)**2
            values[1] = 0.0
        elif between(x[0], (0.3, 1.7)) and near(x[1], 2):
            values[0] = 0.5
            values[1] = 0.0
        else:
            values[0] = 0.5*sin(pi*(x[0]-2)/0.6)**2
            values[1] = 0.0

    def value_shape(self):
        return (2,)

U_top = interpolate(FSS(degree = 0), V)

#U_top = Expression(("1" , "0.0"), degree = 0)	# Top veloc

top = CompiledSubDomain('near(x[1], H1) && on_boundary', H1 = H)
left = CompiledSubDomain('near(x[0], 0.0) && on_boundary')
right = CompiledSubDomain('near(x[0],W1) && on_boundary', W1 = W)
fsi = CompiledSubDomain('near(x[1], h1)  && on_boundary', h1 = h)

# Define boundary conditions
noSlipLeft = DirichletBC(V, Constant((0.0, 0.0)), left)
noSlipRight = DirichletBC(V, Constant((0.0, 0.0)), right)

freestreamV = DirichletBC(V, U_top, top)

dofs_f_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
i_f_V = np.where((dofs_f_V[:,1] == h))[0] #  & (x <= 0.5)
 #  & (x <= 0.5)

#mesh is 25026.
#nodal_M_u = np.loadtxt('nodal_M_u', dtype = float)
#nodal_M_u = np.loadtxt('nodal_M_u_1', dtype = float)
nodal_M_u = np.loadtxt('nodal_M_a_incomp', dtype = float)

u_mesh = Function(V)

u_mesh.vector()[:] = nodal_M_u
u_mesh_FSI = u_mesh.vector()[i_f_V]
u_mesh_FSI1 = u_mesh_FSI
np.savetxt('u_mesh_a_incomp', u_mesh_FSI1)

dofs_f_V[i_f_V]

u_mesh_print = u_mesh.vector()[np.where((dofs_f_V[:,0] == W))[0]]

#class u_mesh_2(Expression):
#    def eval(self,values, x):
#        #between(x[0], (0.0, H)) and near(x[1], 2):
#        for i_c in range(len(u_mesh_FSI)/2)
#            values[0] = u_mesh_FSI[i_c*2]
#            values[1] = u_mesh_FSI[i_c*2]
#    def value_shape(self):
#        return (2,)

#failed to converge using u_mesh.
# try sinusoidal expression.
# works with sinusoid. u_mesh doesn't identify component on boundary.. perhaps.


class u_mesh_1(Expression):
    def eval(self,values, x):
        #between(x[0], (0.0, H)) and near(x[1], 2):
        # one values for symmetric sin, one for asymettric
        #values[1] = 0.2*sin(2*pi*x[0]/W)
        values[1] = (1+0.2*x[0])*0.05*sin(pi*x[0])
        values[0] = 0.0
    def value_shape(self):
        return (2,)

u_mesh_boundary = interpolate(u_mesh_1(degree = 0), V)

#u_mesh_boundary = interpolate(u_mesh_FSI, V)

fluidFSI = DirichletBC(V, u_mesh_boundary, fsi)

#fluidFSI = DirichletBC(V, u_mesh, fsi)

# I think that u_mesh is not physical.
#fluidFSI = DirichletBC(V, Constant((0,-0.1)), fsi)

#fluidFSI = DirichletBC(V, u_mesh, fsi)
#fluidFSI = DirichletBC(V, Constant((0.0, 0.0)), fsi)

bcu = [noSlipLeft, noSlipRight, freestreamV, fluidFSI]

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
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*(u0-u_mesh), v)*dx + \
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
i_f_V_top = np.where((dofs_f_V[:,1] == H))[0]
i_f_V_base = np.where((dofs_f_V[:,1] == h))[0]

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

    #u_top = u1.vector()[i_f_V_top]
    #print "fluid velocity on top plate = ", u_top

    u_base = u1.vector()[i_f_V_base]
    np.savetxt('u_base', u_base)
    u_base_coords = dofs_f_V[i_f_V_base]
    np.savetxt('u_base_coords', u_base_coords)

    print "fluid velocity at base = ", u_base

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

#u_v_chorin = np.column_stack((np.linspace(0.0, 1.0, 101), np.zeros(101)))
#u_u_chorin = np.column_stack((np.linspace(0.0, 1.0, 101), np.zeros(101)))

#for i_x in range(len(u_v_chorin)):
	#u_v_chorin[i_x,1] = u1(u_v_chorin[i_x,0],0.5)[1]
	#u_u_chorin[i_x,1] = u1(0.5,u_u_chorin[i_x,0])[0]

#scipy.io.savemat('u_v_chorin.mat', mdict={'u_v_chorin':u_v_chorin})
#scipy.io.savemat('u_u_chorin.mat', mdict={'u_u_chorin':u_u_chorin})

File('./Results_fluid_test/' + 'u1.pvd') << (u1)
