from dolfin import *
from mshr import *
from numpy import array, loadtxt
import numpy as np
import scipy.io

set_log_level(ERROR)

parameters["form_compiler"]["cpp_optimize"] = True
# parameters["form_compiler"]["optimize"] = True

# Three tests to complete
# CSM 3 is dynamic

N = 64
# 32 corresponds to 546
Bar = Rectangle(Point(0.0,0.0), Point(0.35, 0.02))

mesh = generate_mesh(Bar, N)

vector = VectorFunctionSpace(mesh, "CG", 1)

# Inserted degree for expression
# gravity acting in body force.
B  = Expression(("0.0", "-2"), degree=1)
T_hat  = Expression(("0.0", "0.0"), degree=1)


rho0 = Constant(1000.0)
nu_s = Constant(0.4)
mu_s = Constant(0.5e6)
#mu_s = Constant(2.0e6)

E_s = mu_s*2*(1+nu_s)	# structure elastic modulus

lmbda = nu_s*E_s/((1+nu_s)*(1-2*nu_s)) # 2nd lame constant

# Mixed space for computing both displacement and rate of change of displacement
V_element = VectorElement("CG", mesh.ufl_cell() , 2)
mixed_space = FunctionSpace(mesh, MixedElement([V_element, V_element]))

#mixed_space = MixedFunctionSpace([vector, vector])
V = TestFunction(mixed_space)
dU = TrialFunction(mixed_space)
U = Function(mixed_space)
U0 = Function(mixed_space)

xi, eta = split(V)
u, v = split(U)

v0 = Constant((0,)*vector.mesh().geometry().dim())
u0 = Constant((0,)*vector.mesh().geometry().dim())

cells = CellFunction("size_t", mesh)
dx = Measure('dx', domain = mesh, subdomain_data = cells)

# Project u0 and v0 into U0
a_proj = inner(dU, V)*dx
L_proj = inner(u0, xi)*dx + inner(v0, eta)*dx
solve(a_proj == L_proj, U0)

u0, v0 = split(U0)

left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)

Fixed = Constant((0.0, 0.0))
bcl = DirichletBC(mixed_space.sub(0), Fixed, left)
bcs = [bcl]

# Kinematics
u_mid = 0.5*(u0 + u)
v_mid = 0.5*(v0 + v)

d = u.geometric_dimension()
I = Identity(d)         #identity
F = I + grad(u_mid)     #deformation gradient
C = F.T*F               #Right Cauchy-Green tensor
E = (C - I)/2           #Green-Lagrange strain tensor
E = variable(E)

# Elasticity parameters - defined earlier
#mu    = Constant(3.85)
#lmbda = Constant(5.77)

# St venant model for stored strain energy
psi = lmbda/2*(tr(E)**2) + mu_s*tr(E*E)

# Second Piola-Kirchhoff stress
S = diff(psi, E)
# First Piola-Kirchhoff stress
P = F*S

# replace with expression from Turek
# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)
P = F *( lmbda*tr(E)*I + 2* mu_s*E)*F.T*inv(F.T)


dt = 0.01
k = Constant(dt)

# Variational form dynamic hyperelasticity
L = rho0*inner(v - v0, xi)*dx \
    + k*inner(P, grad(xi))*dx \
    - k*inner(B, xi)*dx \
    + inner(u - u0, eta)*dx \
    - k*inner(v_mid, eta)*dx

L = L - k*inner(T_hat, xi)*ds

a = derivative(L, U, dU)

t = 0.0
T = 10.00

#problem = VariationalProblem(a, L, bcl, nonlinear = True)
problem = NonlinearVariationalProblem(L, U, bcl, a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
solver.parameters["newton_solver"]["maximum_iterations"] = 100

file = File("./Results_twisty/displacement.pvd")

results = np.zeros((int(T/dt)+1,3))
# is the +1 necessary?
count = 0

while t < T:


    t = t + float(k)
    print 'Time ', t
    #problem.solve(U)
    solver.solve()
    u, v = U.split()
    file << u

    U0.assign(U)

    u11, v11 = U.split(deepcopy = True)

    results[count,:] = [t, u11(0.35,0.01)[0], u11(0.35,0.01)[1]]
    count += 1



scipy.io.savemat('results.mat', mdict={'results':results})
## Check to see intial displacement is loaded
#U.assign(U0)
#u, v = U.split()
#file << u;
