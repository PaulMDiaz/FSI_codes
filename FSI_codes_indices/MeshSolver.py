#  This class contains the functions necessary to setup the FSI problem
#  and to handle the transfer of information between the fluid and structure
#  solvers during the time-stepping iteration.

#  Created: 28 March 2014
#  Modified: 23 August 2015
#from dolfin import *
from fenics import *
import numpy as np
import pylab as plt

class Mesh_Solver:

	def Move_Mesh(self, S, F):
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
		bc_m = [DirichletBC(F.V_space, S.d_dot, F.facets, 3)]
		Form_m = a*sym(grad(F.u_mesh))[i, j]*F.del_u_mesh[i].dx(j)*F.dx
		Gain_m = derivative(Form_m, F.u_mesh, F.du_mesh)

		solve(Form_m == 0, F.u_mesh, bc_m, J = Gain_m, 
			solver_parameters = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )


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
