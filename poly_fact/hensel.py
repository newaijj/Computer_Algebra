from polynomial import *
from polynomial_fact_finite import *
from factor_tree import *

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
    assert isinstance(f, polynomial)
    assert isinstance(g, polynomial)
    assert isinstance(h, polynomial)
    assert isinstance(s, polynomial)
    assert isinstance(t, polynomial)

    assert equal_polynomial(f, mul_polynomial(g, h, N=N), N=N)

    print(
        add_polynomial(
            mul_polynomial(s, g, N=N), mul_polynomial(h, t, N=N), N=N
        )
    )
    print(N)

    assert equal_polynomial(
        polynomial((0, 1), N=N),
        add_polynomial(
            mul_polynomial(s, g, N=N), mul_polynomial(h, t, N=N), N=N
        ),
        N=N,
    )

    N_ = N ** 2

    e = sub_polynomial(f, mul_polynomial(g, h, N=N_), N=N_)
    se = mul_polynomial(s, e, N=N_)
    q, r = mod_polynomial(se, h, N=N_)
    g_ = add_polynomial(
        g, mul_polynomial(t, e, N=N_), mul_polynomial(q, g, N=N_), N=N_
    )
    h_ = add_polynomial(h, r, N=N_)

    b = add_polynomial(
        mul_polynomial(s, g_, N=N_),
        mul_polynomial(t, h_, N=N_),
        polynomial((0, -1), N=N_),
        N=N_,
    )
    sb = mul_polynomial(s, b, N=N_)
    c, d = mod_polynomial(sb, h_, N=N_)
    s_ = sub_polynomial(s, d, N=N_)
    t_ = sub_polynomial(
        t,
        add_polynomial(
            mul_polynomial(t, b, N=N_), mul_polynomial(c, g_, N=N_), N=N_
        ),
        N=N_,
    )

    assert equal_polynomial(f, mul_polynomial(g_, h_, N=N_), N=N_)
    assert equal_polynomial(
        polynomial((0, 1), N=N_),
        add_polynomial(
            mul_polynomial(s_, g_, N=N_), mul_polynomial(h_, t_, N=N_), N=N_
        ),
        N=N_,
    )

    return g_, h_, s_, t_


"""
IMPLEMENTED FUNCTION: MULTIFACTOR_HENSEL_LIFTING

	-INPUT:
		- N 			- mod
		- f 			- R[x] of degree n (to be factorised)
		- a  			- multiplicative inverse of lc(f) (mod N)
		- l 			- factor to Hensel Lift by (i.e. N**l)
							-HAS TO BE A FACTOR OF 2 OTHERWISE IT BEHAVES WEIRDLY
		- T 			- factor tree of normal(f) (mod N)

	-OUTPUT:
		- a 			- multiplicative inverse of lc(f) (mod N**l)
		- T 			- factor tree of normal(f) (mod N**l)

"""


def Multifactor_Hensel_Lifting(N, f, a, l, T):
    d = int(np.ceil(np.log(l) / np.log(2)))

    for j in range(1, d + 1, 1):
        N_ = N ** (2 ** j)
        a = (2 * a - lc(f) * a * a) % N_
        T.value = mul_polynomial(polynomial((0, a), N=N_), f, N=N_)

        for node in T.get_walk():
            if not node.isleaf():
                node.L.value, node.R.value, node.s, node.t = Hensel_step(
                    node.value,
                    node.L.value,
                    node.R.value,
                    node.s,
                    node.t,
                    N=N ** (2 ** (j - 1)),
                )

    return a, T


if __name__ == "__main__":

    """
    HENSEL STEP VERIFICATION
    p1 = polynomial()
    p1.set_coef(0, 1)
    p1.set_coef(1, 2)
    p1.set_coef(2, 1)
    print("g: ", p1)
    p2 = polynomial()
    p2.set_coef(0, 1)
    p2.set_coef(1, 2)
    print("h: ", p2)
    print("gcd: ", gcd_polynomial(p1, p2))
    p3 = mul_polynomial(p1, p2)
    print("f (g*h): ", p3)

    _, s, t, _ = extended_euclidean(p1, p2)
    print("Hensel result: ", *Hensel_step(p3, p1, p2, s, t))

    print(
        "verification f: ",
        mul_polynomial(*Hensel_step(p3, p1, p2, s, t)[0:2], N=N ** 2),
    )
    g, h, s_, t_ = Hensel_step(p3, p1, p2, s, t)
    print(
        "verification gs*ht mod N**2: ",
        add_polynomial(
            mul_polynomial(g, s_, N=N ** 2),
            mul_polynomial(h, t_, N=N ** 2),
            N=N ** 2,
        ),
    )
    """
    """
    MULTIFACTOR HENSEL LIFTING VERIFICATION
    p_prob = polynomial()
    p_prob.set_coef(0, 21)
    p_prob.set_coef(1, 20)
    p_prob.set_coef(2, 21)
    p_prob.set_coef(3, 9)
    p_prob.set_coef(4, 15)
    p_prob.set_coef(5, 19)
    p_prob.set_coef(6, 8)
    p_prob.set_coef(7, 9)
    p_prob.set_coef(8, 19)
    p_prob.set_coef(9, 4)
    p_rand = rand_polynomial(10)
    print("random_polynomial: ", p_prob)

    a = mul_inv(lc(p_prob))
    p_norm = mul_polynomial(p_prob,polynomial((0,a)))
    print("normalised: ",p_norm)

    print("factorisation: ")
    facs = factorise_polynomial_int_finite(p_norm)
    for f in facs:
        print(f)
    print("verify: ", mul_polynomial(*facs))


    T = make_tree(facs)
    a,t = Multifactor_Hensel_Lifting(N,p_prob,a,4,T)

    print("a: ", a)
    print("T.value: ",mul_polynomial(T.value,polynomial((0,mul_inv(a,N=N**4)),N=N**4)))

    verification=mul_polynomial(T.L.value,T.R.value,N=N**4)
    verification2 = mul_polynomial(T.L.L.value,T.L.R.value,T.R.L.value,T.R.R.value,N=N**4)
    print("verify: ",verification)
    print("verify2: ", verification2)
    """
