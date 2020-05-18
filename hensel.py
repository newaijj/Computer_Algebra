from polynomial import *
from polynomial_fact import *

"""
FUNCTIONS IMPLEMENTED: HENSEL_STEP
	- raises g,h,s,t from mod m to mod m^2
	- f = g*h mod m
	- s*g + t*h = 1 mod m

	- INPUT: 
		- f,g,h,s,t (all elements of R[x])

	- f_ = g_*h_ mod m^2
	- s_*g_ + t_*h_ = 1 mod m^2

	- OUTPUT:
		- g_,h_,s_,t_ (all elements of R[x])
"""


def Hensel_step(f, g, h, s, t, N=N):
    assert is_instance(f, polynomial)
    assert is_instance(g, polynomial)
    assert is_instance(h, polynomial)
    assert is_instance(s, polynomial)
    assert is_instance(t, polynomial)
    assert equal_polynomial(f, mul_polynomial(g, h, N=N))
    assert equal_polynomial
    return
