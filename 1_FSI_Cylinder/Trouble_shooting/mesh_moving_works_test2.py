from fenics import *
from mshr import *

import numpy as np

# from dolfin import *
# from fenicstools import interpolate_nonmatching_mesh
################################################################################
########  Reference mesh ##########
################################################################################

# perhaps have this line.
# parameters["reorder_dofs_serial"] = False

refMesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),10,10)
refVectorSpace = VectorFunctionSpace(refMesh, "CG", 1)

# some test solution
expr = Expression(("cos(x[0])","exp(x[0])"),degree=3)
refSolution = interpolate(expr,refVectorSpace)
refSolution.set_allow_extrapolation(True)

# reference interface
class refInterface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)

# mark interface
refInterface = refInterface()
boundariesRef = MeshFunction("size_t", refMesh, refMesh.topology().dim() - 1)
boundariesRef.set_all(0)
refInterface.mark(boundariesRef, 1)

# Create subMesh on interface
boundaryMeshRef = BoundaryMesh(refMesh, 'exterior')
refInterfaceMesh = SubMesh(boundaryMeshRef, refInterface)
refInterfaceVectorSpace = VectorFunctionSpace(refInterfaceMesh, "CG", 1)

# Introduce interface solution
refInterfaceSolution = interpolate(refSolution,refInterfaceVectorSpace)

# Check refValues coordinates
refInterfaceCoords = refInterfaceVectorSpace.tabulate_dof_coordinates().reshape((refInterfaceVectorSpace.dim(),-1))

dofsRef = refVectorSpace.tabulate_dof_coordinates().reshape((refVectorSpace.dim(),-1))
indicesRef = np.where((dofsRef[:,1] <= 1.0 + DOLFIN_EPS)&(dofsRef[:,1] >= 1.0 - DOLFIN_EPS))[0]

# refInterfaceSolution has a funny order, but correct values.

# compare solution:
# indicesRef appear to be in opposite direction to
refInterfaceSolution.vector().get_local()
refSolution.vector()[indicesRef]


################################################################################
######## Mesh_2 : Different coords - same degree space  ##########
################################################################################

mesh_2 = RectangleMesh(Point(1.0,2.0),Point(2.0,3.0),10,10)

# Set up a vector function space on each mesh
secondVectorSpace = VectorFunctionSpace(mesh_2, "CG", 1)

class SecondInterface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 2.0)

#The simplest approach is to specify the boundary by a SubDomain object,
# using the inside() function to specify on which facets the boundary conditions should be applied.
# The boundary facets will then be searched for and marked only on the first call to apply.
#This means that the mesh could be moved after the first apply and the boundary markers would still remain intact.

#Note that there may be caching employed in BC computation for performance reasons.
# In particular, applicable DOFs are cached by some methods on a first apply().
# This means that changing a supplied object (defining boundary subdomain) after
# first use may have no effect. But this is implementation and method specific.
secondInterface = SecondInterface()

boundariesMesh2 = MeshFunction("size_t", mesh_2, mesh_2.topology().dim() - 1)
boundariesMesh2.set_all(0)
secondInterface.mark(boundariesMesh2, 1)
# ds = Measure('ds')[boundariesRef]

		# structureSolver.dx = Measure('dx', domain = structureSolver.mesh, subdomain_data = structureSolver.cells)
# ds = Measure('ds', domain = mesh_2, subdomain_data = boundariesMesh2)

boundaryMesh2 = BoundaryMesh(mesh_2, 'exterior')
interfaceMesh2 = SubMesh(boundaryMesh2, secondInterface)
secondInterfaceVectorSpace = VectorFunctionSpace(interfaceMesh2, "CG", 1)
secondInterfaceFunction = Function(secondInterfaceVectorSpace)


# populate secondInterfaceFunction with values from first boundary solution:
secondInterfaceFunction.vector()[:] = refInterfaceSolution.vector()[:]
# ^ this copies exactly.

# secondInterfaceFunction.vector()[:] = refInterfaceSolution.vector()[:]
secondInterfaceFunction.vector().get_local()
secondInterfaceFunction(2.0,2.0)


# is in reverse order to values on original mesh.
# Use this solution in application of boundary condition this mesh:
# what if mesh moves? First mesh, second? Figure out. Need flexible and wise.

# use in definition of the variational problem, ie apply as BC.
# secondInterfaceFunction.set_allow_extrapolation(True)

# Set up expression to use in BC application
class SecondInterfaceValues(Expression):
    def eval(self, values, x):
        try:
            values[:] = secondInterfaceFunction(x)
        except:
            values[:] = 0
    def value_shape(self):
        return (2,)

secondInterfaceValues = SecondInterfaceValues(degree=2)
#
bcNew = DirichletBC(secondVectorSpace,secondInterfaceValues, secondInterface)
# bcNew = DirichletBC(secondVectorSpace,secondInterfaceFunction, secondInterface)

secondSolution = Function(secondVectorSpace)
secondSolution.vector().zero()
#
bcNew.apply(secondSolution.vector())

secondInterfaceCoords = secondInterfaceVectorSpace.tabulate_dof_coordinates().reshape((secondInterfaceVectorSpace.dim(),-1))

dofsSecond = secondVectorSpace.tabulate_dof_coordinates().reshape((secondVectorSpace.dim(),-1))
indicesSecond = np.where((dofsSecond[:,1] <= 2.0 + DOLFIN_EPS)&(dofsSecond[:,1] >= 2.0 - DOLFIN_EPS))[0]
 # opposite order to secondInterfaceCoords
#
secondSolution.vector()[indicesSecond]
secondInterfaceFunction.vector().get_local()

dofsSecond[indicesSecond]
secondSolution.vector()[indicesSecond]

dofsRef[indicesRef]
refSolution.vector()[indicesRef]

print "2 norm refMesh and mesh_2 :", np.linalg.norm(secondSolution.vector()[indicesSecond]  - refSolution.vector()[indicesRef])/np.linalg.norm(refSolution.vector()[indicesRef])

################################################################################
######## Move mesh, start with same coordinates then move, same degree ##########
################################################################################

mesh_3 = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),10,10)

thirdVectorSpace = VectorFunctionSpace(mesh_3, "CG", 1)

thirdSolution = Function(thirdVectorSpace)
thirdSolution.vector().zero()

class ThirdInterface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)

# The boundary facets will then be searched for and marked only on the first
# call to apply. This means that the mesh could be moved after the first apply
# and the boundary markers would still remain intact.

thirdInterface = ThirdInterface()

boundariesMesh3 = MeshFunction("size_t", mesh_3, mesh_3.topology().dim() - 1)
boundariesMesh3.set_all(0)
thirdInterface.mark(boundariesMesh3, 1)

# created a functionspace on the third mesh prior to movement.
boundaryMesh3 = BoundaryMesh(mesh_3, 'exterior')
interfaceMesh3 = SubMesh(boundaryMesh3, thirdInterface)
thirdInterfaceVectorSpace = VectorFunctionSpace(interfaceMesh3, "CG", 1)
thirdInterfaceFunction = Function(thirdInterfaceVectorSpace)
thirdInterfaceFunction.vector().zero()

# thirdInterfaceCoords1 = thirdInterfaceVectorSpace.tabulate_dof_coordinates().reshape((thirdInterfaceVectorSpace.dim(),-1))

class ThirdInterfaceValues(Expression):
    def eval(self, values, x):
        try:
            values[:] = thirdInterfaceFunction(x)
        except:
            values[:] = 0
    def value_shape(self):
        return (2,)

thirdInterfaceValues = ThirdInterfaceValues(degree = 2)


bcMove = DirichletBC(thirdVectorSpace,thirdInterfaceValues, boundariesMesh3, 1)


# move mesh
for i_coords in range(len(mesh_3.coordinates())):
	mesh_3.coordinates()[i_coords] += [1.0, 1.0]

for i_coords in range(len(interfaceMesh3.coordinates())):
	interfaceMesh3.coordinates()[i_coords] += [1.0, 1.0]

# current things that happen post movement:
# assigning values. May have to stay
# applying bc - should probably happen here. May have to happen before movement too..

thirdInterfaceFunction.vector()[:] = refInterfaceSolution.vector().get_local()

bcMove.apply(thirdSolution.vector())

indicesThird = indicesRef

thirdInterfaceCoords = thirdInterfaceVectorSpace.tabulate_dof_coordinates().reshape((thirdInterfaceVectorSpace.dim(),-1))
thirdInterfaceCoords

thirdSolution.vector()[indicesThird]
thirdInterfaceFunction.vector().get_local()
thirdSolution.vector().get_local()
# dofsSecond[indicesSecond]
# secondSolution.vector()[indicesSecond]

dofsRef[indicesRef]
refSolution.vector()[indicesRef]

print "2 norm refMesh and mesh_moved :", np.linalg.norm(thirdSolution.vector()[indicesThird]  - refSolution.vector()[indicesRef])/np.linalg.norm(refSolution.vector()[indicesRef])

# Achieves the right solution in the right place...
# Cost: setting up new subdomain, boundarymesh, interfacemesh, functionspace
# and function. Also expression necessary in order to apply to BC.

# How to make this more elegant? Magic toward this expression?

# make new file
# move the boundary mesh too.
