from polynomial import *
from random import randint

"""
FUNCTIONS IMPLEMENTED: FACTORISE_POLYNOMIAL_INT_FINITE
	- factorises integer polynomial with coefficients in finite field fully 
	  into irreducible factors (irreducible in Fn[x])

	- INPUT 	: polynomial with integer coefficients in finite field
	- OUTPUT 	: list of factors

- HELPER FUNCTIONS
	- equal_degree_factorisation(p, d) 	- factorises squarefree polynomial (p) with factors only of equal degree (d)
	- equal_degree_splitting(p, d)		- returns one factor of squarefree polynomial (p) with factors only of equal degree (d)

- GENERAL PURPOSE FUNCTIONS
	- rand_polynomial(max_degree)		- random polynomial in finite field FN of maximum degree (max_degree)
	- filter_degree(p, min_degree=1)	- returns true if degree of p >= min_degree

"""


def rand_polynomial(max_degree):
    p = polynomial()
    for i in range(0, max_degree + 1):
        p.set_coef(i, randint(0, N))
    return p


def filter_degree(p, min_degree=1):
    return p.degree >= min_degree


class UnluckyStartError(Exception):
    pass


"""
FUNCTION: EQUAL_DEGREE_SPLITTING
Inputs:
- p : squarefree polynomial of degree n
- d : target degree of factors (integer divisor of n)

Output:
- a factor of p of degree d or failure exception (happens probabilistically)
"""


def equal_degree_splitting(p, d):
    a = rand_polynomial(p.degree - 1)
    if not filter_degree(a):
        raise UnluckyStartError

    g1 = gcd_polynomial(a, p)
    p_one = polynomial(size=g1.size)
    p_one.set_coef(0, 1)
    if g1.degree > 0:
        return g1

    b = exp_polynomial_rem(a, int((N ** d - 1) / 2), p)
    b.set_coef(0, b.get_coef(0) - 1)

    g2 = gcd_polynomial(b, p)
    p_one = polynomial(size=g2.size)
    p_one.set_coef(0, 1)
    if g2.degree > 0 and not np.all(np.equal(p.coef, g2.coef)):
        return g2
    else:
        raise UnluckyStartError


def equal_degree_factorisation(p, d):
    if p.degree == d:
        return [p]
    while True:
        try:
            fac = equal_degree_splitting(p, 1)
        except UnluckyStartError:
            continue
        break
    return [fac, *equal_degree_factorisation(mod_polynomial(p, fac)[0], d)]


"""
Implements polynomial factorisation over a polynomial with coefficients in finite field Fn
Input: polynomial p (integer coefficients with finite field global N)
Output: list of polynomial factors of p
"""


def factorise_polynomial_int_finite(p):
    h = polynomial()
    h.set_coef(1, 1)

    v = p

    i = 0

    # set of factors found
    U = []

    v_one = polynomial()
    v_one.set_coef(0, 1)
    v_zero = polynomial()

    while v.degree > 0:
        i += 1

        """
		one distinct degree factorisation step
		(removes all factors of degree i)
		"""
        h = exp_polynomial_rem(h, N, p)
        x = polynomial()
        x.set_coef(1, 1)
        g = gcd_polynomial(sub_polynomial(h, x), v)

        if g.degree > 0:
            """
			call equal degree factorisation algorithm to separate factors fo degree i
			"""
            facs = equal_degree_factorisation(g, i)

            """
			determine multiplicites of factors found
			"""
            for fac in facs:
                while equal_polynomial(mod_polynomial(v, fac)[1], v_zero):
                    U.append(fac)
                    v = mod_polynomial(v, fac)[0]

    return U


if __name__ == "__main__":

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
    p_prob = polynomial()
    p_prob.set_coef(0, 8)
    p_prob.set_coef(1, 12)
    p_prob.set_coef(2, 19)
    p_prob.set_coef(3, 20)
    p_prob.set_coef(4, 6)
    p_prob.set_coef(5, 18)
    # print("p3:", p3)
    p_rand = rand_polynomial(6)
    print("random_polynomial: ", p_rand)
    print("factorisation: ")
    facs = factorise_polynomial_int_finite(p_rand)

    for f in facs:
        print(f)
    print("verify: ", mul_polynomial(*facs))
