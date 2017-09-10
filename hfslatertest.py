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
		v = self.primitives[0](x,y,z)
		for i in range(1,self.numberOfPrimitives) :
			v += self.primitives[i](x,y,z)
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
	def addBasisFunction(f) :
		self.basis.append(f)
	def __call__(self, x,y,z,index,spin) :
		C = self.Cu
		if spin==0 :
			C = self.Cd
		v = self.Cu[index][0] * self.basis[0](x,y,z)
		for i in range(1,len(self.basis)) :
			v += self.Cu[index][i] * self.basis[i](x,y,z)
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
R = [[0.278167628856, -4.238416802148, -0.5432147532112],
     [-1.88033748609, 1.867496255763, 2.707593002034],
     [-2.993552651848, -1.081825109678, -1.583887689758],
     [-1.880942784354, 2.157179898052, -1.021396690302]]
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
#===========================================================================

class Det :
	def __init__(self) :
		self.du = [[None for i in range(numberOfSpinUpElectrons)] \
					for j in range(numberOfSpinUpElectrons)]
		self.dd = [[None for i in range(numberOfSpinDownElectrons)] \
					for j in range(numberOfSpinDownElectrons)]
	def setupDu(self) :
		




sys.exit(1)






i,j,k 		= symbols('i j k', integer=True)
C, alpha 	= symbols('C alpha')
x,y,z 		= symbols('x y z')
r 			= [x,y,z]

def Gp(C,i,j,k,a,R) :
	r2 = np.dot(R,R)
	return C*r[0]**i * r[1]**j * r[2]**k * exp(-a * r2)

def Gc(p,R) :
	if isinstance(p[0],list) :
		s = Gp(*p[0],R)
		for i in range(1,len(p)) :
			s += Gp(*p[i],R)
		return s
	else :
		return Gp(*p,R)

#===========================================================================
electrons = 4
eUp       = 2
eDown     = 2
#===========================================================================
HF = [[-0.99281, -0.076425,  0.028727], 
	  [ 0.21571, -0.229340, -0.822350],
	  [-0.99281, -0.076425,  0.028727], 
	  [-0.21571,  0.229340,  0.822350]]
#===========================================================================
contracted1 = 	[[ 0.06442630000000,0,0,0,71.887600000000006], 
				 [ 0.36609600000000,0,0,0,10.728899999999999],
				 [ 0.69593400000000,0,0,0, 2.222050000000000]]
contracted2 = 	[[-0.42106400000000,0,0,0, 1.295480000000000],
				 [ 1.22407000000000,0,0,0, 0.268881000000000]]
contracted3 = 	[[ 1.00000000000000,0,0,0, 0.077350000000000]]
#===========================================================================
basis = [contracted1,contracted2,contracted3]
#===========================================================================



x1,y1,z1 = symbols('x_1 y_1 z_1')
x2,y2,z2 = symbols('x_2 y_2 z_2')
x3,y3,z3 = symbols('x_3 y_3 z_3')
x4,y4,z4 = symbols('x_4 y_4 z_4')

Rs = [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],[x4,y4,z4]]

def value(f,X,Y,Z) :
	return f.subs({x:X, y:Y, z:Z})

def Du(orbital,R) :
	d = HF[orbital][0]*Gc(basis[0],R)
	for i in range(1,3) :
		d += HF[orbital][i]*Gc(basis[i],R)
	return d


SlaterUp      = [[None for i in range(2)] for j in range(2)]
SlaterDown    = [[None for i in range(2)] for j in range(2)]

for i in range(eUp) : # electron
	for j in range(eUp) : # orbital
		SlaterUp[i][j] = Du(j+2,Rs[i][:])

for i in range(eUp) : # electron
	for j in range(eUp) : # orbital
		SlaterDown[i][j] = Du(j, Rs[2+i][:])


def subsVal(a) :
	return a.subs({x1:R[0][0],y1:R[0][1],z1:R[0][2],x2:R[1][0],y2:R[1][1],z2:R[1][2],x3:R[2][0],y3:R[2][1],z3:R[2][2],x4:R[3][0],y4:R[3][1],z4:R[3][2]})


def laplacian(D) :
	return 	diff(D,x1,x1)+diff(D,x2,x2)+diff(D,x3,x3) + \
			diff(D,y1,y1)+diff(D,y2,y2)+diff(D,y3,y3) + \
			diff(D,z1,z1)+diff(D,z2,z2)+diff(D,z3,z3)

matUp 	= Matrix(SlaterUp)
matDown = Matrix(SlaterDown)

pprint(subsVal(matDown))
pprint(subsVal(matUp))

laplacianUp   = laplacian(det(matUp))
laplacianDown = laplacian(det(matDown))

Ku = subsVal(laplacianUp/det(matUp))
Kd = subsVal(laplacianDown/det(matDown))
print(Ku,Kd,Ku+Kd)














