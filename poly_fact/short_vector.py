import numpy as np 
import copy

"""
FUNCTIONS IMPLEMENTED: BASIS_REDUCTION

	- INPUT: 
		- b 		- a basis, with integer coefficients
					- (refer to helper class for declaration)

	- OUTPUT:
		- g 		- a reduced basis of the lattice spanned by b
			- reduced basis has a corresponding Gram-Schmidt orthogonal basis
			  (f1*,f2*,f3*, ... , fn*) such that ||fi*||^2 <= 2 ||f(i+1)*||^2




HELPER CLASS: BASIS
		- stores a list of basis vectors 
		- __init__(self,N=23,mat=None)
			- does not check linear independence
			- mat can be specified as a 2D square np array

	- ATTRIBUTES
		- N 			- number of dimensions spanned (assuming linearly independent)
	
	- METHODS
		- get_vec(i)	
		- set_vec(i)	
		- get_val(i,j)	
		- set_val(i,j)	
		- GSO()			- computes Gram-Schmidt orthogonalisation, returns as a basis
		- swap_vec(i,j)

HELPER CLASS: VECTOR
		- __init__(self,val)
			- val is a 1D np array

	- METHODS
		- norm()		- computes ||f||^2

HELPER FUNCTIONS: 
	- inner(v1,v2) 		- computes inner product of 2 vectors
	
"""

class basis:
	def __init__(self,N=23,mat=None):
		assert float(N).is_integer()
		if mat is None:
			self.mat = np.zeros((N,N))
		else:
			assert mat.shape[0] == mat.shape[1]
			self.mat = mat
			self.N = mat.shape[0]
			return

		self.N = N

	def __str__(self):
		return self.mat.__str__()

	def get_vec(self,i):
		return copy.deepcopy(vector(self.mat[i,:]))

	def set_vec(self,i,v):
		assert v.shape == self.mat[i,:].shape
		assert isinstance(v,vector)
		self.mat[i,:] = v.val

	def get_val(self,i,j):
		return self.mat[i,j]

	def set_val(self,i,j,val):
		self.mat[i,j] = val

	def GSO(self):
		gso = basis(N=self.N)
		M = basis(mat=np.eye(self.N))

		for i in range(self.N):
			v = self.get_vec(i)
			vi = self.get_vec(i)
			for j in range(0,i):
				vj = self.get_vec(j)
				mu = inner(vi,vj)/inner(vj,vj)
				v = v - vj*mu

				M.set_val(i,j,mu)
			gso.set_vec(i,v)
		return gso,M

	def swap_vec(self,i,j):
		t_v = self.get_vec(i)
		self.set_vec(i,self.get_vec(j))
		self.set_vec(j,t_v)


def inner(v,w):
	assert isinstance(v,vector) and isinstance(w,vector)

	v_ = np.inner(v.val,w.val)
	return v_

class vector:
	def __init__(self,val):
		self.val = val
		self.shape = val.shape

	def __str__(self):
		return self.val.__str__()

	def __add__(self,v):
		assert isinstance(v,vector)
		assert self.shape == v.shape
		return vector(self.val+v.val)

	def __sub__(self,v):
		assert isinstance(v,vector)
		assert self.shape == v.shape
		return vector(self.val-v.val)

	def __mul__(self,s):
		assert not isinstance(s,vector)
		return vector(self.val*s)

	def norm(self):
		return np.inner(self.val,self.val)


def basis_reduction(b):
	assert np.all(np.equal(np.mod(b.mat, 1), 0))

	g = copy.deepcopy(b)
	i = 1
	while i<b.N:
		gso,M = g.GSO()
		for j in range(i):
			v = g.get_vec(i) - g.get_vec(j)*int(M.get_val(i,j))
			g.set_vec(i,v)
		gso,M = g.GSO()
		if i>0 and gso.get_vec(i-1).norm()>2*gso.get_vec(i).norm():
			g.swap_vec(i,i-1)
			i-=1
		else:
			i+=1
	return g



if __name__ == "__main__":
	b = basis(N=3)
	v1 = vector(np.array([1,1,1]))
	v2 = vector(np.array([-1,0,2]))
	v3 = vector(np.array([3,5,6]))
	b.set_vec(0,v1)
	b.set_vec(1,v2)
	b.set_vec(2,v3)
	print(b)

	print(basis_reduction(b))

	v4 = vector(np.array([0,1,0]))



