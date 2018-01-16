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

N = 240
cylinder = Circle(Point(0.2, 0.2), 0.05, N)
# 32 corresponds to 546

l = 0.35
h = 0.02
L = 2.5
H = 0.41

bar_2 = Rectangle(Point(0.2,0.19), Point(0.6, 0.21))
bar = bar_2 - cylinder

mesh = generate_mesh(bar, N)

vector = VectorFunctionSpace(mesh, "CG", 1)

# Inserted degree for expression
# gravity acting in body force.
#B  = Expression(("0.0", "-2.0"), degree=2)
rho0 = Constant(1000.0)
B  = Constant((0.0, -2000.0))
# B is -2000: - rho * g with g = 2.

T_hat  = Constant((0.0, 0.0))

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

class Fsi(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > l - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < (0.19 + DOLFIN_EPS) and on_boundary or \
		x[0] > l - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > (0.21 - DOLFIN_EPS) and x[1] < H - DOLFIN_EPS and on_boundary or \
		x[0] > 0.6 - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < H - DOLFIN_EPS and on_boundary

class Left(SubDomain):
    def inside(self,x,on_boundary):
        return (x[0]-0.2)**2 + (x[1]-0.2)**2 < 0.05**2 and on_boundary #or\

fsi = Fsi()
left = Left()

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

dt = 0.02
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
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-8
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-8
solver.parameters["newton_solver"]["maximum_iterations"] = 1000


file = File("./Results_dyanmic/displacement.pvd")

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

    results[count,:] = [t, u11(0.6,0.2)[0], u11(0.6,0.2)[1]]
    count += 1

#N = 64
# scipy.io.savemat('results_7906_02.mat', mdict={'results':results})
# N = 128
#scipy.io.savemat('results_29850_005.mat', mdict={'results':results})
# N = 240
scipy.io.savemat('results_105134_02.mat', mdict={'results':results})

## Check to see intial displacement is loaded
#U.assign(U0)
#u, v = U.split()
#file << u;
