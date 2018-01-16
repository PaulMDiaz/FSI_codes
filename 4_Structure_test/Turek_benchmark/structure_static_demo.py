from dolfin import *
from mshr import *
from numpy import array, loadtxt
import numpy as np
import scipy.io

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True

# Create mesh and define function space
N = 612
cylinder = Circle(Point(0.2, 0.2), 0.05, N)
# 32 corresponds to 546
l = 0.35
h = 0.02
L = 2.5
H = 0.41

Bar = Rectangle(Point(0.0,0.0), Point(l, h))
# Too long bar, cut off with Circle
bar_2 = Rectangle(Point(0.2,0.19), Point(0.6, 0.21))

bar = bar_2 - cylinder

mesh = generate_mesh(bar, N)

rho0 = Constant(1000.0)
nu_s = Constant(0.4)
mu_s = Constant(0.5e6)

E_s = mu_s*2*(1+nu_s)	# structure elastic modulus

lmbda = nu_s*E_s/((1+nu_s)*(1-2*nu_s)) # 2nd lame constant

V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Mark boundary subdomians
class Fsi(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > l - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < (0.19 + DOLFIN_EPS) and on_boundary or \
		x[0] > l - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > (0.21 - DOLFIN_EPS) and x[1] < H - DOLFIN_EPS and on_boundary or \
		x[0] > 0.6 - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < H - DOLFIN_EPS and on_boundary

# find boundary with cylinder by using cylinder geometry..
# For now do interior nodes, not the one in contact with the cylinder...

#class Left(SubDomain):
#    def inside(self,x,on_boundary):
#        return x[0] < l + DOLFIN_EPS and \
#        x[1] > 0.19+DOLFIN_EPS and x[1] < 0.21 - DOLFIN_EPS and on_boundary #or\

class Left(SubDomain):
    def inside(self,x,on_boundary):
        return (x[0]-0.2)**2 + (x[1]-0.2)**2 < 0.05**2 and on_boundary #or\

fsi = Fsi()
left = Left()

facets = FacetFunction("size_t", mesh)
facets.set_all(0)
fsi.mark(facets, 2)

#left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)

Fixed = Constant((0.0, 0.0))
bcl = DirichletBC(V, Fixed, left)
bcl = [bcl]


# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0, -2000.0))      # Body force per unit volume
#B  = Constant((0.0, -2000.0))      # Body force per unit volume

T_hat  = Constant((0.0,  0.0))  # Traction force on the boundary

# Kinematics
d = u.geometric_dimension()
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

E = (C - I)/2
E = variable(E)

# Stored strain energy density (compressible neo-Hookean model)
#psi = (mu_s/2)*(Ic - 3) - mu_s*ln(J) + (lmbda/2)*(ln(J))**2



# St venant model for stored strain energy
psi = lmbda/2*(tr(E)**2) + mu_s*tr(E*E)

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T_hat, u)*ds(2)

F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

#solve(F == 0, u, bcl, J=J)

problem = NonlinearVariationalProblem(F, u, bcl, J=J)
solver = NonlinearVariationalSolver(problem)

solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-9
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-9
solver.parameters["newton_solver"]["maximum_iterations"] = 50

solver.solve()

#problem = LinearVariationalProblem(F, u, bcl, J)
#solver = LinearVariationalSolver(problem)

#solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
#solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
#solver.parameters["newton_solver"]["maximum_iterations"] = 1000

#solver.solve()



#solver.solve()
# Save solution in VTK format
file = File("Results_static_demo/displacement.pvd");
file << u;

results = [u(0.6,0.2)[0], u(0.6,0.2)[1]]
# Plot and hold solution
#plot(u, mode = "displacement", interactive = True)
#file << u;

# save displacement field
#u_txt = u.vector().array()
#u_txt = u.vector().get_local()

#np.savetxt('twisty.txt', u_txt)

# Plot and hold solution
#plot(u, mode = "displacement", interactive = True)
