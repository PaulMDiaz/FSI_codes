from dolfin import *
from mshr import *
from numpy import array, loadtxt
import numpy as np
import scipy.io
from nonlinear_solver import *

# Optimization options for the form compiler
parameters["form_compiler"]["optimize"] = True
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4

def default_parameters():
   "Return default solver parameters."
   p = Parameters("solver_parameters")
   p.add("plot_solution", True)
   p.add("save_solution", False)
   p.add("store_solution_data", False)
   p.add("element_degree",2)
   p.add("problem_formulation",'displacement')
   new = Parameters("newton_solver")
   new.add("value", 1.0)
   new.add("adaptive", True)
   new.add("loading_number_of_steps", 1)
   p.add(new)

   return p
# Create mesh and define function space
N = 64
# 32 corresponds to 546
Bar = Rectangle(Point(0.0,0.0), Point(0.35, 0.02))

mesh = generate_mesh(Bar, N)

rho0 = Constant(1000.0)
nu_s = Constant(0.4)
mu_s = Constant(0.5e6)

E_s = mu_s*2*(1+nu_s)	# structure elastic modulus

lmbda = nu_s*E_s/((1+nu_s)*(1-2*nu_s)) # 2nd lame constant

V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Mark boundary subdomians
class Fsi(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > DOLFIN_EPS and on_boundary

fsi = Fsi()

facets = FacetFunction("size_t", mesh)
facets.set_all(0)
fsi.mark(facets, 2)

left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)

Fixed = Constant((0.0, 0.0))
bcl = DirichletBC(V, Fixed, left)
bcl = [bcl]


# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0, -0.1))      # Body force per unit volume
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


# St venant model for stored strain energy
psi = lmbda/2*(tr(E)**2) + mu_s*tr(E*E)
# 2nd PK stress
S = diff(psi, E)
# 1st PK stress
P = F*S

theta = Constant(1.0)

L = inner(P, grad(v))*dx -inner(B, v)*dx

L = L - inner(T_hat, v)*ds(2)

a = derivative(L, u, du)

#solver = AugmentedNewtonSolver(L, u, a, bcl,\
#                               load_increment = theta)

problem = NonlinearVariationalProblem(L, u, bcl, a)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
solver.parameters["newton_solver"]["maximum_iterations"] = 1000
#prm["newton_solver"]["relaxation_parameter"] = theta

#newton_parameters = parameters["newton_solver"]
#solver.parameters['newton_solver'].update(newton_parameters)

solver.solve()

# Save solution in VTK format
file = File("Results_static/displacement.pvd");
file << u;

# Plot and hold solution
#plot(u, mode = "displacement", interactive = True)
#file << u;

# save displacement field
#u_txt = u.vector().array()
#u_txt = u.vector().get_local()

#np.savetxt('twisty.txt', u_txt)

# Plot and hold solution
#plot(u, mode = "displacement", interactive = True)
