#__author__ = "Harish Narayanan"
#__copyright__ = "Copyright (C) 2009 Simula Research Laboratory and %s" % __author__
#__license__  = "GNU GPL Version 3 or any later version"

from dolfin import *

# Infinitesimal strain tensor
def infinitesimalStrain(u):
	return variable(0.5*(grad(u) + grad(u).T))
	# Second order identity tensor
def secondOrderIdentity(u):
	return variable(Identity(u.geometric_dimension()))

# Deformation gradient
def deformationGradient(u):
	I = secondOrderIdentity(u)
	return variable(I + grad(u))

# Determinant of the deformation gradient
def jacobian(u):
	F = deformationGradient(u)
	return variable(det(F))

# Right Cauchy-Green tensor
def rightCauchyGreen(u):
	F = deformationGradient(u)
	return variable(F.T*F)

# Green-Lagrange strain tensor
def greenLagrangeStrain(u):
    I = secondOrderIdentity(u)
    C = rightCauchyGreen(u)
    return variable(0.5*(C - I))

# Left Cauchy-Green tensor
def leftCauchyGreen(u):
	F = deformationGradient(u)
	return variable(F*F.T)

# Euler-Almansi strain tensor
def eulerAlmansiStrain(u):
	I = secondOrderIdentity(u)
	b = leftCauchyGreen(u)
	return variable(0.5*(I - inv(b)))

# Invariants of an arbitrary tensor, A
def invariants(A):
    I1 = tr(A)
    I2 = 0.5*(tr(A)**2 - tr(A*A))
    I3 = det(A)
    return [I1, I2, I3]

# Invariants of the (right/left) Cauchy-Green tensor
def cauchyGreenInvariants(u):
	C = rightCauchyGreen(u)
	[I1, I2, I3] = invariants(C)
	return [variable(I1), variable(I2), variable(I3)]

	# Isochoric part of the deformation gradient
def isochoricDeformationGradient(u):
	F = deformationGradient(u)
	J = jacobian(u)
	return variable(J**(-1.0/3.0)*F)

# Isochoric part of the right Cauchy-Green tensor
def isochoricRightCauchyGreen(u):
	C = rightCauchyGreen(u)
	J = jacobian(u)
	return variable(J**(-2.0/3.0)*C)

# Invariants of the ischoric part of the (right/left) Cauchy-Green
# tensor. Note that I3bar = 1 by definition.
def isochoricCauchyGreenInvariants(u):
	Cbar = isochoricRightCauchyGreen(u)
	[I1bar, I2bar, I3bar] = invariants(Cbar)
	return [variable(I1bar), variable(I2bar)]



# Principal stretches
def principalStretches(u):
    C = rightCauchyGreen(u)
    S = FunctionSpace(u.function_space().mesh(), "CG", 1)
    if (u.cell().d == 2):
        D = sqrt(tr(C)*tr(C) - 4.0*det(C))
	eig1 = sqrt(0.5*(tr(C) + D))
	eig2 = sqrt(0.5*(tr(C) - D))
	return [variable(eig1), variable(eig2)]
    if (u.cell().d == 3):
	c = (1.0/3.0)*tr(C)
	D = C - c*secondOrderIdentity(u)
	q = (1.0/2.0)*det(D)
	p = (1.0/6.0)*inner(D, D)
	ph = project(p, S)
	if (norm(ph) < DOLFIN_EPS):
            eig1 = sqrt(c)
	    eig2 = sqrt(c)
	    eig3 = sqrt(c)
        else:
	    phi = (1.0/3.0)*atan(sqrt(p**3.0 - q**2.0)/q)
	    if (phi < 0.0):
                phi = phi + DOLFIN_PI/3.0
	    end
	    eig1 = sqrt(c + 2*sqrt(p)*cos(phi))
	    eig2 = sqrt(c - sqrt(p)*(cos(phi) + sqrt(3)*sin(phi)))
	    eig3 = sqrt(c - sqrt(p)*(cos(phi) - sqrt(3)*sin(phi)))
        return [variable(eig1), variable(eig2), variable(eig3)]

# Pull-back of a two-tensor from the current to the reference
# configuration
def piolaTransform(A, u):
    J = jacobian(u)
    F = deformationGradient(u)
    B = J*A*inv(F).T
    return B

# Push-forward of a two-tensor from the reference to the current
# configuration
def inversePiolaTransform(A, u):
    J = jacobian(u)
    F = deformationGradient(u)
    B = (1/J)*A*F.T
    return B
