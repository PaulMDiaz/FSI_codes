"""This program solves the incompressible Navier-Stokes equations
on flow past a cyilnder using Chorin's splitting method."""

# Copyright (C) 2010-2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Mikael Mortensen 2011
#
# First added:  2010-08-30
# Last changed: 2011-06-30

# Begin demo

from __future__ import print_function
from dolfin import *
import numpy as np
from mshr import *
import scipy.io

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
#mesh = Mesh("../lshape.xml.gz")
#mesh = Mesh("lshape.xml")

# Set parameter values
dt = 0.01
T = 8.00

# Domain length and height
L = 2.5
H = 0.41

# mesh discretization
N = 64

channel = Rectangle(Point(0, 0), Point(L, H))

cylinder = Circle(Point(0.2, 0.2), 0.05, 2*N)

#U_max = 0.3
U_mean = 0.2

nu_f = 0.001    # kinematic viscosity
rho_f = 1

mu_f = nu_f*rho_f   # dynamic viscosity

geometry = channel - cylinder

mesh = generate_mesh(geometry, N)



# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

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

# Create functions
u0 = Function(V)
u1 = Function(V)
p0 = Function(Q)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

U = 0.5*(u0+u)
n = FacetNormal(mesh)

# beta is a flag. set to zero when periodic boundary conditions are used.
beta = 1

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu_f*epsilon(u) - p*Identity(len(u))


# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx \
    + inner(grad(u0)*u0, v)*dx \
    + inner(sigma(U, p0), epsilon(v))*dx \
    + inner(p0*n, v)*ds \
    - beta*nu_f*inner(grad(U).T*n, v)*ds \
    - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure correction
a2 = inner(grad(p), grad(q))*dx
L2 = inner(grad(p0), grad(q))*dx \
    -(1.0/k)*div(u1)*q*dx

# Velocity correction
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1-p0), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("results_cylinder/velocity.pvd")
pfile = File("results_cylinder/pressure.pvd")


count = 0;
# initialize results matrix to store drag, lift and p_diff
results = np.zeros((int(T/dt)+1,5))

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    #p_in.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "bicgstab", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p1.vector()) for bc in bcp]
    solve(A2, p1.vector(), b2, "bicgstab", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "bicgstab", "default")
    end()

    # Plot solution
    #plot(p1, title="Pressure", rescale=True)
    #plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    p0.assign(p1)
    t += dt

    drag, lift, = compute_forces(mesh, nu_f, u1, p1, ds)

    p_diff = p1(0.15, 0.2)-p1(0.25, 0.2)

    C_D = 2/(pow(U_mean,2)*2*0.05)*drag
    C_L = 2/(pow(U_mean,2)*2*0.05)*lift

    results[count,:] = [p_diff, C_D, C_L, drag, lift]

    count += 1

    print("t =", t)

# Hold plot
#interactive()
scipy.io.savemat('results_cylinder_64.mat', mdict={'results':results})
