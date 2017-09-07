# Test that structure solver gives physical results
# Import structure solver
from StructureSolver_1 import Structure_Solver
# Import dolfin and numpy
from fenics import *
import numpy as np
import pylab as plt
import math

# class for structure parameters
class Structure_params:
    def __init__(self): #Initialize
        # Lame' constants
        self.mu_s = 2 #self.E_s/(2.0*(1.0 + self.nu_s))
        self.lambda_s = 2.7*10**2  # non-zero for compressible structures?
        # Number of elements in each dimension
        self.N_x = 32
        self.N_y = 8
        # height of structure
        self.h = 0.5
        # Numerical parameters
        self.dt = 0.05	# Time step
        self.T = 0.00   # Set final time for iteration ie 0.05
        # Set up a variable for time
        self.t = 0
        self.iter_t = 1
        # Mesh definition
        self.mesh = RectangleMesh(Point(0.0, 0.0), Point(2.0, self.h), self.N_x, self.N_y)
        self.Dim = self.mesh.topology().dim() #Not sure what this does?

    # Boundary conditions function for Structure_params object-1
    def Define_Boundary_Conditions(self, S):
        h = self.h
        # Left boundary for structure and fluid
        class Left(SubDomain):
            def inside(self, x, on_boundary):
				return near(x[0], 0.0)
	# Right boundary for structure and fluid
        class Right(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[0], 2.0)
	# FSI boundary for fluid and structure
        class Top(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[1], h)
	# Bottom boundary for structure
        class Bottom(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], 0.0)

        # Initialize structure boundary parts
        S.boundary_parts = MeshFunction("size_t", S.mesh, S.mesh.topology().dim()-1)

        # Initialize boundary objects
        S.left = Left()
        S.right = Right()
        S.top = Top()
        S.bottom = Bottom()

        # Mark boundaries
        S.left.mark(S.boundary_parts, 1)
        S.right.mark(S.boundary_parts, 2)
        S.top.mark(S.boundary_parts, 3)
        S.bottom.mark(S.boundary_parts, 0)

        # Define the boundary integral variables
        S.dA = ds(subdomain_data=S.boundary_parts)
        S.dV = dx

        S.N = FacetNormal(S.mesh)

        #  Set up the boundary conditions (Should we have a top BC?)
    	#  BCs for the left side (no displacement)
        LeftBC = DirichletBC(S.V_space, Constant((0.0, 0.0)), S.left)

		#  BCs for the right side (no displacement)
        RightBC = DirichletBC(S.V_space, Constant((0.0, 0.0)), S.right)

        S.bcs = [LeftBC, RightBC]

#  Create Structure_params object
SP = Structure_params()


# an artifical 'Fluid_Solver' object is created to provide the pressure field at the
# 'interface' a.k.a. top boundary of the structure
class Fluid_Solver(object):
    def __init__(self, SP):

        self.mesh = SP.mesh #This is the structure mesh
        self.T_space = TensorFunctionSpace(self.mesh, "Lagrange", 2)
        I = Identity(SP.Dim)
        #self.S_space = FunctionSpace(self.mesh, "Lagrange", 1)
        #self.p1 = Function(self.S_space)
        #I = np.zeros((2,2))
        #zero_function = Expression("value", value = 0.0, degree = 1)
        # Here is where we attempt to implement a zero pressure field, which should correspond to zero stress on the top
        # boundary of the structure
        self.sigma_FSI = project(I-I, self.T_space, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )



F = Fluid_Solver(SP)


## Structure Solver Setup

#  Set the Structure Element Type
StructureElementType = "CG" # CG is equivilent to lagrange.

#  Set the Structure Element Degree
StructureElementDegree = 1 # Standard linear lagrange element.

# Set the solver used for the structure problem
# The options are "Linear" or "NeoHookean" for the structure solver
# Note that the Linear solver does not work and didn't work in the original code
#StructureSolverMethod = "NeoHookean"
StructureSolverMethod = "Linear"

# Body forces on the structure
StructureBodyForce = Constant((0.0, 0.0))

# Set the equations for the variational solver and boundary conditions and creat structure solver object
S = Structure_Solver(SP.mesh, StructureElementType,
	StructureElementDegree, StructureSolverMethod, StructureBodyForce)

# Define boundary conditions on the structure
SP.Define_Boundary_Conditions(S)

# Structure array of coordinates
dofs_s_S = S.S_space.tabulate_dof_coordinates().reshape((S.S_space.dim(),-1))
dofs_s_V = S.V_space.tabulate_dof_coordinates().reshape((S.V_space.dim(),-1))
dofs_s_T = S.T_space.tabulate_dof_coordinates().reshape((S.T_space.dim(),-1))

#Extract dof indices for values on boundary.
# y = 0.5 if mesh is not deformed.

i_s_S = np.where((dofs_s_S[:,1] == 0.5))[0] #  & (x <= 0.5)
i_s_V = np.where((dofs_s_V[:,1] == 0.5))[0] #  & (x <= 0.5)
i_s_T = np.where((dofs_s_T[:,1] == 0.5))[0] #  & (x <= 0.5)


# Begin solving over time interval
while SP.t < SP.T + DOLFIN_EPS:
    print ' STARTING TIME LOOP ITERATION ', SP.iter_t
    SP.t += SP.dt
    for ii in range(1):
        print ''
        print ''
        print 'Time loop iteration number = ', SP.iter_t
        print 'Loop iteration time = ' , SP.t

        d_ = S.d.vector()[i_s_S]
        print 'structure deflection on interface prior to structure solver = ', d_

        #sigma_FSI = F.sigma_FSI.vector()[i_s_T]
        #print 'stress at interface = ', sigma_FSI

        # have to edit this line a fair bit.
        S.Structure_Problem_Solver(SP, F)

        d_  = S.d.vector()[i_s_S]
        print 'structure deflection on interface after structure solver = ', d_

    S.d0.assign(S.d) # set current time to previous for displacement
