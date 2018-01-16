from dolfin import *
from mshr import *
from numpy import array, loadtxt
import numpy as np
parameters["form_compiler"]["cpp_optimize"] = True
# parameters["form_compiler"]["optimize"] = True

mesh = UnitCubeMesh(8, 8, 8)
vector = VectorFunctionSpace(mesh, "CG", 1)

# Inserted degree for expression
B  = Expression(("0.0", "0.0", "0.0"), degree=2)
T  = Expression(("0.0", "0.0", "0.0"), degree=2)

# Mixed space for computing both displacement and rate of change of displacement
V_element = VectorElement("CG", mesh.ufl_cell() , 1)
mixed_element = FunctionSpace(mesh, MixedElement([V_element, V_element]))

#mixed_element = MixedFunctionSpace([vector, vector])
V = TestFunction(mixed_element)
dU = TrialFunction(mixed_element)
U = Function(mixed_element)
U0 = Function(mixed_element)

xi, eta = split(V)
u, v = split(U)

u0 = Function(vector)
# This txt file, saved from demo_hyperelasticity works with projection method.
_u0 = np.loadtxt('twisty.txt', dtype = float)
u0.vector()[:] = _u0[:]

v0 = Constant((0,)*vector.mesh().geometry().dim())

cells = CellFunction("size_t", mesh)
dx = Measure('dx', domain = mesh, subdomain_data = cells)

# Project u0 and v0 into U0
a_proj = inner(dU, V)*dx
L_proj = inner(u0, xi)*dx + inner(v0, eta)*dx
solve(a_proj == L_proj, U0)

## Try alternative for twist data from web
#_u0 = np.loadtxt('twisty_downloaded.txt', dtype = float)
#_u02 = np.concatenate((_u0 ,0.0*_u0),axis = 0)
#U0.vector()[:] = _u02


u0, v0 = split(U0)

# deecopy appears to make it a permanant command, ie solution doesn't change.
#u0, v0 = U0.split(deepcopy = True)


clamp = Constant((0.0, 0.0, 0.0))
#left = compile_subdomains("x[0] == 0")
#class Left(SubDomain):#
#    def inside(self,x,on_boundary):
#        return x[0] < DOLFIN_EPS and on_boundary
## Compile boundary
#left = Left()

left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)



bcl = DirichletBC(mixed_element.sub(0), clamp, left)

u_mid = 0.5*(u0 + u)
v_mid = 0.5*(v0 + v)

#I = Identity(v.cell().d)
I = Identity(v.__len__())
F = I + grad(u_mid)
C = F.T*F
E = (C - I)/2
E = variable(E)

mu    = Constant(3.85)
lmbda = Constant(5.77)

# St venant model for stored strain energy
psi = lmbda/2*(tr(E)**2) + mu*tr(E*E)


S = diff(psi, E)
P = F*S

rho0 = Constant(1.0)
dt = Constant(0.1)

L = rho0*inner(v - v0, xi)*dx + dt*inner(P, grad(xi))*dx \
    - dt*inner(B, xi)*dx - dt*inner(T, xi)*ds \
    + inner(u - u0, eta)*dx - dt*inner(v_mid, eta)*dx
a = derivative(L, U, dU)

t = 0.0
T = 2.0

#problem = VariationalProblem(a, L, bcl, nonlinear = True)
problem = NonlinearVariationalProblem(L, U, bcl, a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
solver.parameters["newton_solver"]["maximum_iterations"] = 100

file = File("./Results_twisty/displacement.pvd")

while t < T:

    t = t + float(dt)

    #problem.solve(U)
    solver.solve()
    u, v = U.split()
    file << u

    U0.assign(U)

## Check to see intial displacement is loaded
#U.assign(U0)
#u, v = U.split()
#file << u;
