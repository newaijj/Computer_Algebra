from polynomial import *
from random import randint

def rand_polynomial(max_degree):
	p = polynomial()
	for i in range(0,max_degree+1):
		p.set_coef(i,randint(0,N))
	return p

def filter_degree(p,min_degree=1):
	return p.degree>=min_degree

"""
Inputs:
- p : squarefree polynomial of degree n
- d : target degree of factors (integer divisor of n)

Output:
- a factor of p of degree d or failure exception (happens probabilistically)
"""

class UnluckyStartError(Exception):
	pass

def equal_degree_splitting(p,d):
	a = rand_polynomial(p.degree-1)
	if not filter_degree(a):
		raise UnluckyStartError

	g1 = gcd_polynomial(a,p)
	p_one = polynomial(size=g1.size)
	p_one.set_coef(0,1)
	if(g1.degree>0):#not np.all(np.equal(p_one.coef, g1.coef))):
		return g1

	b = exp_polynomial_rem(a,int((N**d-1)/2),p)
	b.set_coef(0,b.get_coef(0)-1)

	g2 = gcd_polynomial(b,p)
	p_one = polynomial(size=g2.size)
	p_one.set_coef(0,1)
	if(g2.degree>0 and not np.all(np.equal(p.coef, g2.coef))): #not np.all(np.equal(p_one.coef, g2.coef)) and not np.all(np.equal(p.coef, g2.coef))):
		return g2
	else:
		raise UnluckyStartError

def equal_degree_factorisation(p,d):
	if(p.degree==d):
		return [p]
	while(True):
		try:
			fac = equal_degree_splitting(p,1)
		except UnluckyStartError:
			continue
		break
	return [fac,*equal_degree_factorisation(mod_polynomial(p,fac)[0],d)]

"""
Implements polynomial factorisation over a polynomial with coefficients in finite field Fn
Input: polynomial p (integer coefficients with finite field global N)
Output: list of polynomial factors of p
"""
def fact_polynomial(p):
	h = polynomial()
	h.set_coef(1,1)

	#keeps track of remainder of polynomial
	#lc_inv = polynomial()
	#lc_inv.set_coef(0,mul_inv(p.get_coef(p.degree)))
	v = p#mul_polynomial(p,lc_inv)

	i=0

	#set of factors found
	U = []
	#const = polynomial()
	#const.set_coef(0,p.get_coef(p.degree))

	v_one = polynomial()
	v_one.set_coef(0,1)
	v_zero = polynomial()

	while(v.degree>0):#not equal_polynomial(v_one,v)):
		i+=1

		"""
		one distinct degree factorisation step
		(removes all factors of degree i)
		"""
		h = exp_polynomial_rem(h,N,p)
		x = polynomial()
		x.set_coef(1,1)
		g = gcd_polynomial(sub_polynomial(h,x),v)

		if(g.degree>0):
			"""
			call equal degree factorisation algorithm to separate factors fo degree i
			"""
			facs = equal_degree_factorisation(g,i)

			"""
			determine multiplicites of factors found
			"""
			for fac in facs:
				while(equal_polynomial(mod_polynomial(v,fac)[1],v_zero)):
					U.append(fac)
					v = mod_polynomial(v,fac)[0]
	#v = mod_polynomial(v,g)[0]
	#lc_inv = mul_polynomial(lc_inv,g)


	#lc_inv = mul_polynomial(lc_inv,v)
	#U.append(lc_inv)
		
	return U



p0 = polynomial()
p3 = polynomial()
p3.set_coef(0, 2)
p3.set_coef(1, 3)
p3.set_coef(2, 1)
p4 = polynomial()
p4.set_coef(0, 1)
p4.set_coef(1, 3)
p4.set_coef(2, 3)
p4.set_coef(3, 1)
p_prob  = polynomial()
p_prob.set_coef(0,8)
p_prob.set_coef(1,12)
p_prob.set_coef(2,19)
p_prob.set_coef(3,20)
p_prob.set_coef(4,6)
p_prob.set_coef(5,18)
#print("p3:", p3)
p_rand = rand_polynomial(6)
print("random_polynomial: ", p_rand)
#print("p4:", p4)

#print("equal_degree_factorisation",equal_degree_factorisation(p4,1))
print("factorisation: ")
facs = fact_polynomial(p_rand)

for f in facs:
	print(f)
print("verify: ",mul_polynomial(*facs))

#print("split: ", gcd_polynomial(p3, p4))

#p5 = rand_polynomial(4)
#print("p_rand: ", filter_degree(p0))


