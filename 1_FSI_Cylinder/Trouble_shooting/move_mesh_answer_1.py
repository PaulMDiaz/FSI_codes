from dolfin import *

mesh = UnitSquareMesh(10,10)

# (including the change-of-variables in the form can cause the automatically-
# computed quadrature degree to blow up unnecessarily)
dx = dx(metadata={'quadrature_degree': 2})

# space for solution
VuE = FiniteElement("Lagrange",mesh.ufl_cell(),1)
Vu = FunctionSpace(mesh,VuE)

# space for mesh deformation
VFE = VectorElement("Lagrange",mesh.ufl_cell(),1)
VF = FunctionSpace(mesh,VFE)

# set some mapping, F, from the mesh to physical space
# (deforms square into parallelogram in this case)
Fexpr = Expression(("x[1]+x[0]","x[1]"),degree=2)
F = interpolate(Fexpr,VF)

# This could be equivilent to applied mesh velocities. ie:
# Some mesh velocity is solved for
meshVelocity = Function(VF)
# meshVelocity.vector().get_local()[:] = F.vector().get_local()[:]
# Test case, use F as values in mesh velocity, set dt as 1.

meshVelocity.vector()[:] = F.vector().get_local()[:]
# This is then applied to mapping from reference to deformed.
# More straightforward in uniform, but should be okay for unstructured too.
# constructed on same space as F.
dt = 1.0

meshDisplacement = Function(VF)
meshDisplacement.vector()[:] = meshVelocity.vector()[:]*dt

# check this worked:
meshDisplacement.vector().get_local()

# Replace F with meshDisplacement (simulating what will happen in FSI problem)
F = meshDisplacement

# derivatives needed for change-of-variables
DF = grad(F) + Identity(2)

J = det(DF)
def gradx(u):
    # (overkill for scalar-valued u)
    n = rank(u)
    ii = indices(n+2)
    # multivariate chain rule:
    # contract over last index of grad(u) and first index of
    # inv(DF) to change variables in derivative; should work for u of arbitrary
    # rank (scalar, vector, tensor)
    return as_tensor(grad(u)[ii[0:n+1]]*inv(DF)[ii[n],ii[n+1]],
                     ii[0:n]+(ii[n+1],))

# pose Poisson problem in physical space
u = TrialFunction(Vu)
v = TestFunction(Vu)

# use gradients w.r.t. physical space (gradx), and include Jacobian determinant
# (J) in integration measure
a = inner(gradx(u),gradx(v))*J*dx

# NOTE: some care must be taken with more complicated Expressions, as "x[0]"
# and "x[1]" refer to coordinates in the UN-deformed mesh, not physical space
L = inner(Constant(1.0),v)*J*dx

# boundary is defined on the UN-deformed mesh
def boundary(x,on_boundary):
    return on_boundary

u = Function(Vu)
solve(a==L,u,[DirichletBC(Vu,Constant(0.0),boundary),])

# Can plot solution on deformed mesh in Paraview by outputting u and d to
# separate files, loading both of those files, applying the "Calculator"
# filter to compute the identity map for each one, selecting both of the
# Calculator results, applying the "Append Attributes" filter, then applying
# "Warp by Vector" using d, and coloring based on u.
