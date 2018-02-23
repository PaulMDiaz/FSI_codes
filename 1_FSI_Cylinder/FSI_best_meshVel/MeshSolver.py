#  This class contains the functions necessary to setup the FSI problem
#  and to handle the transfer of information between the fluid and structure
#  solvers during the time-stepping iteration.

#  Created: 28 March 2014
#  Modified: 23 August 2015
#from dolfin import *
from fenics import *
import numpy as np
import pylab as plt

class MeshSolver(object):

	def __init__(self, meshElementType, meshElementDegree, p_s):

		self.mesh = p_s.meshRef	# reference mesh
		self.elementType = meshElementType
		self.elementDegree = meshElementDegree

		self.tensorSpace = TensorFunctionSpace(self.mesh, self.elementType, self.elementDegree)
		self.vectorSpace = VectorFunctionSpace(self.mesh, self.elementType, self.elementDegree)
		self.scalarSpace = FunctionSpace(self.mesh, self.elementType, self.elementDegree)

		self.meshVelocity = Function(self.vectorSpace) # Current time step
		self.meshDisplacement1 = Function(self.vectorSpace) # Current time step
		self.meshDisplacement1.vector().zero()
		self.meshDisplacement0 = Function(self.vectorSpace) # previos time step
		self.meshDisplacement0.vector().zero()
		self.meshDisplacement = Function(self.vectorSpace) # previos time step
		self.meshDisplacement.vector().zero()

		# self.v_mesh = Function(self.vectorSpace) # for projection of structure velocity

		# Trial and test functions for mesh motion
		self.del_u_mesh = TestFunction(self.vectorSpace)
		self.du_mesh = TrialFunction(self.vectorSpace)

		self.cells = CellFunction("size_t", self.mesh)
		self.facets = FacetFunction("size_t", self.mesh)
		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)
		# self.ds = Measure('ds', domain = self.mesh, subdomain_data = self.facets)

		self.interfaceMeshVectorSpace = VectorFunctionSpace(p_s.interfaceMesh, "CG", 1)
		self.meshInterfaceVelocity = interpolate(self.meshVelocity,self.interfaceMeshVectorSpace)
		self.interfaceMeshCoords = self.interfaceMeshVectorSpace.tabulate_dof_coordinates().reshape((self.interfaceMeshVectorSpace.dim(),-1))



	def meshProblemSolver(self, structureSolver, fluidSolver, p_s):
		# def Move_Mesh(self, Structure, Fluid):
		# Description: Moves structure and fluid meshes given structure displacement
		# Inputs: Structure - structure solver object
		#		  Fluid - fluid solver object

		print ''
		print ''
		print 'ENTERING MESH MOVEMENT'
		print ''
		print ''


		# artificial viscosity
		a = 1.0e-11		# MPa/s
		# a = 1.0e-9		# MPa/s

		self.Form_m = a*sym(grad(self.meshVelocity))[i, j]*self.del_u_mesh[i].dx(j)*self.dx
		self.Gain_m = derivative(self.Form_m, self.meshVelocity, self.du_mesh)

		solve(self.Form_m == 0, self.meshVelocity, self.bcMesh, J = self.Gain_m,
			solver_parameters = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-10} }, \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )

		# try with no mesh movement to trouble shoot fluid solver.
		# self.meshVelocity.vector().zero()

		# interpolate mesh solution to FSI interface
		self.meshInterfaceVelocity = interpolate(self.meshVelocity,self.interfaceMeshVectorSpace)

		# # update mesh displacement using mesh velocity.
		# # first, set previos mesh displacement to present
		# self.meshDisplacement0.vector()[:] = self.meshDisplacement1.vector().get_local()

		# this should only occur when time step changes!!! not within time step.. ie displacement remains at 0 until mesh is moved.
		# Update present
		# self.meshDisplacement1.vector()[:] += p_s.dt*self.meshVelocity.vector().get_local()
        #
		# # take average over time step for use in piola transform of fluid stress
		# # may have to initialize meshDisplacment.
		# self.meshDisplacement.vector()[:] = 0.5*(self.meshDisplacement0.vector().get_local()+ self.meshDisplacement1.vector().get_local())

	def Generate_Files(Structure, Fluid):
		File("fluid_pressure.pvd") << Fluid.p1
		File("fluid_velocity.pvd") << Fluid.u1
		File("fluid_mesh.pvd") << Fluid.mesh
		File("mesh.pvd") << mesh
		File("structure_mesh.pvd") << Structure.mesh
		File("structure_displacement.pvd") << Structure.d
		print ""
		print "Simulation Finished"
		print "Close plots to exit"
