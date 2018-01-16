def s_solve(u0, u00, sigma_f, t):
	rho_s = 8.3E-9		# tonne/mm^3
	nu_s = 0.3
	E_s = 200E3		# in MPa
	mu_s = E_s/(2.*(1.+nu_s))
	lam_s = 2.*mu_s*nu_s/(1.-2.*nu_s)

        print 'before sigma_s projection'
        # mesh.init(1)
        # print mesh.num_cells(),mesh.num_edges(),mesh.num_faces(),mesh.num_facets()

	sigma_s = project(sigma_f, T_s_space, solver_type = "mumps",\
		form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
        print 'passed sigma_s projection'

	# plot ( facets_s, interactive = True)
	displ_s_bottom = Expression(('A*sin(2.*pi*nu*t)', '0.0'), A = 0.50, nu = 10., t = t, degree = 2)
	bc_s = [DirichletBC(V_s_space, displ_s_bottom, facets_s, 1)]
	u_s = Function(V_s_space)
	del_u = TestFunction(V_s_space)
	du = TrialFunction(V_s_space)

	F_s = as_tensor( u_s[k].dx(i) + delta[k, i], (k, i) )
	J_s = det(F_s)

	C_s = as_tensor( F_s[k, i]*F_s[k, j], (i, j) )
	E_s = as_tensor(1./2.*(C_s[i, j] - delta[i, j]), (i, j) )
	S_s = as_tensor( lam_s*E_s[k, k]*delta[i, j] + 2.*mu_s*E_s[i, j], (i, j))
	P_s = as_tensor( F_s[i, j]*S_s[k, j], (k, i) )

	t_hat = as_tensor( J_s*inv(F_s)[k, j]*sigma_s[j, i]*N[k] , (i, ) )

	Form_s = ( rho_s*(u_s-2.*u0_s+u00_s)[i]/(dt*dt)*del_u[i] + P_s[k, i]*del_u[i].dx(k) - rho_s*f[i]*del_u[i] )*dV - \
			 t_hat[i]*del_u[i]*dA(2)

	Gain_s = derivative(Form_s, u_s, du)

	solve(Form_s == 0, u_s, bc_s, J = Gain_s, \
		solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
		form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

	v_s = project( (u_s - u0_s)/dt, V_s_space, solver_type = "mumps", \
		form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

	# plot(v_s, interactive = True)

	return u_s, v_s
