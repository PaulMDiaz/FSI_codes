from fenics import *
from mshr import *

# mshr generates non-regular meshes, so I've just used RectangleMesh()
# to make sure the nodes line up as expected.

# # Two adjacent meshes that have coincident nodes on shared bouandary:
# box_1 = Rectangle(Point(0, 0), Point(0.1, 0.1))
# box_2 = Rectangle(Point(0, 0.1), Point(0.1, 0.2))
#
# # Generate mesh
# mesh_1 = generate_mesh(box_1, 4)
# mesh_2 = generate_mesh(box_2, 4)

mesh_1 = RectangleMesh(Point(0.0,0.0),Point(0.1,0.1),3,3)
mesh_2 = RectangleMesh(Point(0.0,0.1),Point(0.1,0.2),3,3)

# Set up a vector function space on each mesh
V1 = VectorFunctionSpace(mesh_1, "CG", 1)
V2 = VectorFunctionSpace(mesh_2, "CG", 1)

# Create some non-trivial function on mesh 1. Define expr with higher
# degree, so that, when checking u2 vs expr, we have some non-zero error
# (that converges as the mesh is refined). (This mainly serves to check that
# we're not just getting zero by integrating over a wrong/empty domain.)
expr = Expression(("cos(x[0])","exp(x[0])"),degree=3)
u1 = interpolate(expr,V1)

# allow extrapolation on that function
u1.set_allow_extrapolation(True)

# use u1 as a BC on mesh 2
def boundary2(x, on_boundary):
    return near(x[1],0.1)
bc2 = DirichletBC(V2, u1, boundary2)

# apply BC on mesh 2
u2 = Function(V2)
u2.vector().zero()
bc2.apply(u2.vector())


# now check the accuracy:

x1 = Expression("x[1]",degree=1)

# Look at L2 error on the shared boundary

# this should be zero, since u1 should equal u2 exactly on the shared
# boundary
print("Error between u2 and u1 (should be machine zero): "
      # (conditional isolates bottom boundary of mesh_2)
      # conditional(condition, true_value, false_value)
      # gt(a,b) = (a > b)
      # this therefore appears to apply true value above boundary...
      +str(assemble(conditional
                    (gt(0.1-DOLFIN_EPS,x1),
                     inner(u2-u1, u2-u1),
                     Constant(0.0))*ds(domain=mesh_2))))

# this should be nonzero, but converging with mesh refinement, since u1
# (and hence u2) is not exactly equal to the higher-degree expression expr
print("Error between u2 and expr (should be nonzero but convergent): "
      +str(assemble(conditional
                    (gt(0.1-DOLFIN_EPS,x1),
                     inner(u2-expr, u2-expr),
                     Constant(0.0))*ds(domain=mesh_2))))
