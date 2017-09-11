import numpy as np
from numpy.linalg import det
import sys
from sympy.abc import *
from sympy import *
from sympy.matrices import *

class Atom :
	def __init__(self, Z,x,y,z) :
		self.Z,self.x,self.y,self.z = Z,x,y,z
	def __str__(self) :
		return self.Z + " (" + self.x + "," + self.y + "," + self.z + ")"

class Primitive :
	def __init__(self, C,i,j,k,a) :
		self.C,self.i,self.j,self.k,self.a = C,i,j,k,a
	def __call__(self, r) :
		return self.C * r[0]**self.i * r[1]**self.j * r[2]**self.k \
				* exp(-self.a * np.dot(r,r))
	def __str__(self) :
		return  "%10.6g (%1d,%1d,%1d) exp(-%8.5g r²)" % \
				(self.C,self.i,self.j,self.k,self.a)

class Contracted :
	def __init__(self, n,Ax,Ay,Az) :
		self.numberOfPrimitives = n
		self.Ax,self.Ay,self.Az = Ax,Ay,Az
		self.primitives = []
	def addPrimitive(self, primitive) :
		self.primitives.append(primitive)
	def __call__(self, x,y,z) :
		x,y,z = x-self.Ax,y-self.Ay,z-self.Az
		v = self.primitives[0]([x,y,z])
		for i in range(1,self.numberOfPrimitives) :
			v += self.primitives[i]([x,y,z])
		return v
	def __str__(self) :
		string = "Contracted @(%3.1g,%3.1g,%3.1g) :\n" % \
				(self.Ax,self.Ay,self.Az)
		for i in range(self.numberOfPrimitives) :
			string += "   " + self.primitives[i].__str__() + "\n"
		return string

class Basis :
	def __init__(self, Cu, Cd) :
		self.Cu = Cu
		self.Cd = Cd
		self.basis = []
	def addBasisFunction(self, f) :
		self.basis.append(f)
	def __call__(self, x,y,z,index,spin) :
		C = self.Cu
		if spin==0 :
			C = self.Cd
		v = C[index][0] * self.basis[0](x,y,z)
		for i in range(1,len(self.basis)) :
			v += C[index][i] * self.basis[i](x,y,z)
		return v


fileName = '/Users/morten/Documents/Master/HartreeFock/HartreeFockBases/Be-3-21G'
inFile   = open(fileName,'r') 

l = inFile.readline().split()
numberOfElectrons 			= int(l[0])
numberOfSpinUpElectrons 	= int(l[1])
numberOfSpinDownElectrons 	= int(l[2])
basisSize 					= int(l[3])
numberOfAtoms 				= int(l[4])

atoms = []
for i in range(numberOfAtoms) :
	l = inFile.readline().split()
	Z,x,y,z = int(l[0]),float(l[1]), float(l[2]), float(l[3])
	atoms.append(Atom(Z,x,y,z))

b = []
for i in range(basisSize) :
	numberOfPrimitives = int(inFile.readline())
	l = inFile.readline().split()
	x,y,z = float(l[0]), float(l[1]), float(l[2])
	c = Contracted(numberOfPrimitives,x,y,z)

	for i in range(numberOfPrimitives) :
		l = inFile.readline().split()
		i,j,k = int(l[0]),int(l[1]),int(l[2])
		alpha,C = float(l[3]), float(l[4])
		p = Primitive(C,i,j,k,alpha)
		c.addPrimitive(p)	
	b.append(c)

coefficientsUp   = np.zeros((numberOfSpinUpElectrons,   basisSize))
coefficientsDown = np.zeros((numberOfSpinDownElectrons, basisSize))
basis = Basis(coefficientsUp,coefficientsDown)
for i in range(len(b)) :
	basis.addBasisFunction(b[i]) 

for i in range(numberOfSpinUpElectrons) :
	l = inFile.readline().split()
	for j in range(basisSize) :
		coefficientsUp[i][j] = float(l[j])

for i in range(numberOfSpinDownElectrons) :
	l = inFile.readline().split()
	for j in range(basisSize) :
		coefficientsDown[i][j] = float(l[j])

inFile.close()

#===========================================================================
R = [[-2.232577556201, -2.027172478823, -1.847618479869],
     [-0.09484430410055, 0.6015963661268, -0.4881592969004],
     [-0.5670967224069, -3.856472063511, -0.4809414284783],
     [0.5196217494483, 0.621255038224, 1.99644060499]]
Ru = [	[R[0][0], R[0][1], R[0][2]],
		[R[1][0], R[1][1], R[1][2]]]
Rd = [	[R[2][0], R[2][1], R[2][2]],
		[R[3][0], R[3][1], R[3][2]]]
#===========================================================================
x1,y1,z1 = symbols('x_1 y_1 z_1')
x2,y2,z2 = symbols('x_2 y_2 z_2')
x3,y3,z3 = symbols('x_3 y_3 z_3')
x4,y4,z4 = symbols('x_4 y_4 z_4')
x,y,z    = symbols('x y z')
r 		 = [x,y,z]
r1 		 = [x1,y1,z1]
r2 		 = [x2,y2,z2]
r3 		 = [x3,y3,z3]
r4 		 = [x4,y4,z4]
rr 		 = [r1,r2,r3,r4]
ru 		 = [r1,r2]
rd 		 = [r3,r4]
#===========================================================================

class Det :
	def __init__(self, spin=None) :
		self.spin = spin
		self.du = [[None for i in range(numberOfSpinUpElectrons)] \
						 for j in range(numberOfSpinUpElectrons)]
		self.dd = [[None for i in range(numberOfSpinDownElectrons)] \
						 for j in range(numberOfSpinDownElectrons)]
		if spin==None :
			self.setupDu()
			self.setupDd()
		else :
			if spin==1 :
				self.setupDu()
			else :
				self.setupDd()
	def setupDu(self) :
		for i in range(numberOfSpinUpElectrons) :
			for j in range(numberOfSpinUpElectrons) :
				self.du[i][j] = basis(*ru[i],j,1)
		self.du = Matrix(self.du)
	def setupDd(self) :
		for i in range(numberOfSpinDownElectrons) :
			for j in range(numberOfSpinDownElectrons) :
				self.dd[i][j] = basis(*rd[i],j,0)
		self.dd = Matrix(self.dd)
	def __call__(self, spin=None) :
		if spin==None :
			if self.spin != None :
				spin = self.spin
			else :
				spin = 1
		return self.du if (spin==1) else self.dd
	def value(self, spin=None) :
		if spin==None :
			if self.spin != None :
				spin = self.spin
			else :
				spin = 1
		if spin==1 :
			return self.du.subs({x1:R[0][0],y1:R[0][1],z1:R[0][2],
								 x2:R[1][0],y2:R[1][1],z2:R[1][2],
								 x3:R[2][0],y3:R[2][1],z3:R[2][2],
								 x4:R[3][0],y4:R[3][1],z4:R[3][2]})
		else :
			return self.dd.subs({x1:R[0][0],y1:R[0][1],z1:R[0][2],
								 x2:R[1][0],y2:R[1][1],z2:R[1][2],
								 x3:R[2][0],y3:R[2][1],z3:R[2][2],
								 x4:R[3][0],y4:R[3][1],z4:R[3][2]})

def laplacian(f) :
	return 	diff(f,x1,x1) + diff(f,y1,y1) + diff(f,z1,z1) + \
			diff(f,x2,x2) + diff(f,y2,y2) + diff(f,z2,z2) + \
			diff(f,x3,x3) + diff(f,y3,y3) + diff(f,z3,z3) + \
			diff(f,x4,x4) + diff(f,y4,y4) + diff(f,z4,z4)

def value(f) :
	return f.subs({x1:R[0][0],y1:R[0][1],z1:R[0][2],
				   x2:R[1][0],y2:R[1][1],z2:R[1][2],
				   x3:R[2][0],y3:R[2][1],z3:R[2][2],
				   x4:R[3][0],y4:R[3][1],z4:R[3][2]})

#Du = Det(1)
#Dd = Det(0)
#
#Su = det(Du.du)
#Sd = det(Dd.dd)
#
#print(value(laplacian(Su)/Su))
#print(value(laplacian(Sd)/Sd))
#
#sys.exit(1)


x     = Symbol('x',     real=True)
y     = Symbol('y',     real=True)
z     = Symbol('z',     real=True)
a     = Symbol('alpha', real=True, positive=True)
r     = Symbol('r',     real=True, positive=True)
phi   = Symbol('phi',   real=True)
theta = Symbol('theta', real=True)


psi_1s  = exp(-a*sqrt(x**2+y**2+z**2))
psi_2s  = (1-a*r/2)*exp(-a*r/2)
psi_2px = x*exp(-a*r/2)
psi_2py = y*exp(-a*r/2)
psi_2pz = z*exp(-a*r/2)

def integrand(f) :
	return f*r**2*sin(theta)

def sphericalIntegral(f) :
	f = f.subs(sqrt(x**2+y**2+z**2), r)
	f = f.subs(x, r*sin(theta)*cos(phi))
	f = f.subs(y, r*sin(theta)*sin(phi))
	f = f.subs(z, r*cos(theta))
	return integrate(integrate(integrate(integrand(f),(theta,0,pi)),(phi,0,2*pi)),(r,0,oo))

N_1s  = simplify(sqrt(1/sphericalIntegral(psi_1s**2)))
N_2s  = simplify(sqrt(1/sphericalIntegral(psi_2s**2)))
N_2px = simplify(sqrt(1/sphericalIntegral(psi_2px**2)))
N_2py = simplify(sqrt(1/sphericalIntegral(psi_2py**2)))
N_2pz = simplify(sqrt(1/sphericalIntegral(psi_2pz**2)))

print("", end="\n       ")
print("N_1s:");   pprint(N_1s);  print("",end="\n       ")
print("N_2s:");   pprint(N_2s);  print("",end="\n       ")
print("N_2px:");  pprint(N_2px); print("",end="\n       ")
print("N_2py:");  pprint(N_2py); print("",end="\n       ")
print("N_2pz:");  pprint(N_2pz); print("",end="\n       ")




pprint(diff(N_1s*psi_1s,x).subs(sqrt(x**2+y**2+z**2), r))




