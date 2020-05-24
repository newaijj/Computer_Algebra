import numpy as np
from logging import warning
from itertools import takewhile
import copy

"""
CLASS IMPLEMENTED: POLYNOMIAL
    - univariable only
    - coefficients mod N
    - __init__(self, initial=None, size=64, N=N)

- ATTRIBUTES
    - size      - total number of coefficients stored (defaults to 64)
    - degree    - degree of polynomial (maintained on insertions, but not for direct access to p.coef)
    - N         - current N of this polynomial

- METHODS
    - set_coef(power, coefficient)
    - get_coef(power)
    - update_degree(power=None) - checks if degree needs to be updated. If power is specified, checks if that power is larger than current degree and sets it if it is
    - copy_meta()               - returns polynomial with copy of all meta variables but empty coefficients
    - update_N(N)               - updates self.N if N is bigger than self.N
    - max_norm                  - return max_norm (max coefficient)
    - lc                        - wrapper for lc function

- FUNCTIONS - POLYNOMIAL
    - add_polynomial()                      - accepts any number of polynomial objects as input
    - mul_polynomial()                      - accepts any number of polynomial objects as input
    - sub_polynomial(p1,p2,N=N)             - (p1 - p2)
    - mod_polynomial(p1,p2,N=N)             - (p1 % p2); returns (p_quotient , p_remainder) given global N
    - exp_polynomial(p,e,N=N)               - returns (p^e), e must be integer
    - exp_polynomial_rem(p, e, f,N=N)       - returns (p^e mod f), e must be integer, p and f must be polynomials
    - gcd_polynomial(p1, p2, N=N)           - return GCD of the two polynomials given global N
    - equal_polynomial(p1, p2, N=N)         - returns if polynomials are equal
    - lc(p)                                 - returns leading coefficient of polynomial
    - extended_euclidean(p1, p2, N=N)       - return p,s,t,r (Algo 3.6), s*p1+t*p2=r=gcd(p1,p2)     

- FUNCTION - MISC
    - mul_inv(a)    -   calculates multiplicative inverse given global N
"""

N = 23


def add_polynomial(*args, N=N):
    max_size = 0
    for p in args:
        assert isinstance(p, polynomial)
        max_size = max(max_size, p.size)
    p_out = polynomial(size=max_size, N=N)

    for p in args:
        out = np.zeros(max_size)
        out[0 : (p.degree + 1)] = p.coef[0 : (p.degree + 1)]
        p_out.coef += out
        p_out.coef = p_out.coef % N

    p_out.update_degree()
    return p_out


def add_inv_polynomial(p1, N=N):
    assert isinstance(p1, polynomial)
    p_out = copy.deepcopy(p1)
    p_out.coef = -p_out.coef % N
    return p_out


def sub_polynomial(p1, p2, N=N):
    return add_polynomial(p1, add_inv_polynomial(p2, N=N), N=N)


def mul_2_polynomial(p1, p2, N=N):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)

    np_out = np.convolve(p1.coef[0 : p1.degree + 1], p2.coef[0 : p2.degree + 1])
    p_out = polynomial(size=max(p1.size, p2.size, np_out.shape[0]), N=N)
    p_out.degree = np_out.shape[0] - 1
    p_out.coef[0 : np_out.shape[0]] = np_out
    p_out.coef = p_out.coef % N

    return p_out


def mul_polynomial(*args, N=N):
    max_size = 0
    for p in args:
        assert isinstance(p, polynomial)
        max_size = max(max_size, p.size)
    p_out = polynomial(size=max_size, N=N)
    p_out.coef[0 : args[0].degree + 1] = args[0].coef[0 : args[0].degree + 1]
    p_out.degree = args[0].degree

    for p in args[1:]:
        p_out = mul_2_polynomial(p_out, p, N=N)

    return p_out


# calc multiplicative inverse of a number mod N based on solving Bézout's identity
def mul_inv(a, N=N):
    a = a % N
    assert a >= 0
    other_coef = np.array([0, 1])
    cur_coef = np.array([1, 0])
    other = N
    cur = a
    while cur != 0:
        while other >= cur:
            other -= cur
            other_coef -= cur_coef
        other, cur = cur, other
        other_coef, cur_coef = cur_coef, other_coef
    return other_coef[0] % N


def mod_polynomial(p1, p2, N=N):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)

    p_q = polynomial(size=max(p1.size, p2.size), N=N)
    p_r = copy.copy(p1)
    for p in range(p1.degree, p2.degree - 1, -1):
        fac = (p_r.get_coef(p) * mul_inv(p2.get_coef(p2.degree), N=N)) % N
        p_q.set_coef(p - p2.degree, fac)
        if fac == 0:
            continue
        p_fac = polynomial(size=p2.size, N=N)
        p_fac.set_coef(p - p2.degree, fac)
        p_r = sub_polynomial(p_r, mul_polynomial(p2, p_fac, N=N), N=N)

    return p_q, p_r


# exponent by repeated squaring
def exp_polynomial(p, e, N=N):
    assert isinstance(p, polynomial) and float(e).is_integer()

    p_out = p.copy_meta()
    p_out.set_coef(0, 1)
    if e == 0:
        return p_out
    else:
        for a in bin(e).replace("0b", ""):
            p_out = mul_polynomial(p_out, p_out, N=N)
            if a == "1":
                p_out = mul_polynomial(p_out, p, N=N)
    return p_out


# exponent by repeated squaring, modded some polynomial f
def exp_polynomial_rem(p, e, f, N=N):
    assert (
        isinstance(p, polynomial)
        and isinstance(f, polynomial)
        and float(e).is_integer()
    )

    p_out = p.copy_meta()
    p_out.set_coef(0, 1)
    if e == 0:
        return p_out
    else:
        for a in bin(e).replace("0b", ""):
            p_out = mod_polynomial(mul_polynomial(p_out, p_out, N=N), f, N=N)[1]
            if a == "1":
                p_out = mod_polynomial(mul_polynomial(p_out, p, N=N), f, N=N)[1]
    return p_out


def gcd_polynomial(p1, p2, N=N):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    if p2.degree > p1.degree:
        p1, p2 = p2, p1
    p_zero = p2.copy_meta()

    while not np.all(np.equal(p2.coef, p_zero.coef)):
        p_t = mod_polynomial(p1, p2, N=N)[1]
        p1 = p2
        p2 = p_t
    return p1


def equal_polynomial(p1, p2, N=N):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    return np.all(np.equal(p1.coef % N, p2.coef % N))


def lc(p):
    assert isinstance(p, polynomial)
    return p.get_coef(p.degree)


def extended_euclidean(p1, p2, N=N):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    p = [polynomial((0, lc(p1)), N=N), polynomial((0, lc(p2)), N=N)]
    s = [polynomial((0, mul_inv(lc(p1), N=N)), N=N), polynomial(N=N)]
    t = [polynomial(N=N), polynomial((0, mul_inv(lc(p2), N=N)), N=N)]
    r1 = mul_polynomial(p1, s[0], N=N)
    r2 = mul_polynomial(p2, t[1], N=N)
    r = [r1, r2]
    q = [0, 0]

    cur = 1
    other = 0
    p_zero = polynomial(N=N)
    while not equal_polynomial(r[other], p_zero):
        temp = cur
        cur = other
        other = temp

        q[cur] = mod_polynomial(r[other], r[cur], N=N)[0]
        rem = mod_polynomial(r[other], r[cur], N=N)[1]
        p[other] = polynomial((0, lc(rem)), N=N)
        r[other] = mul_polynomial(
            polynomial((0, mul_inv(lc(rem), N=N)), N=N), rem
        )
        s[other] = mul_polynomial(
            sub_polynomial(s[other], mul_polynomial(s[cur], q[cur], N=N), N=N),
            polynomial((0, mul_inv(p[other].get_coef(0), N=N))),
            N=N,
        )
        t[other] = mul_polynomial(
            sub_polynomial(t[other], mul_polynomial(t[cur], q[cur], N=N), N=N),
            polynomial((0, mul_inv(p[other].get_coef(0), N=N))),
            N=N,
        )

    return p[cur], s[cur], t[cur], r[cur]


class polynomial:
    # polynomial up to 64 terms
    def __init__(self, initial=None, size=64, N=N):
        self.coef = np.zeros(size)
        self.size = size
        self.degree = 0
        self.N = N
        if initial != None:
            assert isinstance(initial, tuple)
            self.set_coef(initial[0], initial[1])

    def copy_meta(self):
        p = polynomial(size=self.size)
        p.size = self.size
        p.degree = self.degree
        p.N = self.N
        return p

    def update_N(N):
        assert float(N).is_integer()
        if N < self.N:
            self.coef %= N
        update_degree()
        self.N = N

    def update_degree(self, power=None):
        if power == None:
            self.degree = self.max_nonzero_pow()
        else:
            assert float(power).is_integer()
            if self.coef[power] != 0 and self.degree < power:
                self.degree = power

    def set_coef(self, power, coef):
        assert power >= 0 and power < self.size
        self.coef[power] = coef % self.N
        self.update_degree(power=power)

    def get_coef(self, power):
        assert power >= 0 and power < self.size
        return self.coef[power]

    def max_nonzero_pow(self):
        try:
            return np.max(np.nonzero(self.coef))
        except ValueError:
            return 0

    def max_norm(self):
        return np.max(self.coef)

    def lc(self):
        return lc(self)

    def __str__(self):
        nzero = np.nonzero(self.coef)
        if nzero[0].size == 0:
            r = "0"
        else:
            r = ""
            it = np.nditer(nzero)
            nxt = next(it)
            r += str(self.get_coef(nxt)) + "x^" + str(nxt) + " "

            while True:
                try:
                    nxt = next(it)
                    r += "+ " + str(self.get_coef(nxt)) + "x^" + str(nxt) + " "
                except StopIteration:
                    break
        return r


if __name__ == "__main__":
    """
    p3 = polynomial()
    p3.set_coef(0, 2)
    p3.set_coef(1, 3)
    p3.set_coef(2, 1)
    p4 = polynomial()
    p4.set_coef(0, 2)
    p4.set_coef(1, 6)
    p4.set_coef(2, 6)
    p4.set_coef(3, 2)
    p, s, t, r = extended_euclidean(p3, p4)
    jkj = p3, p4
    for i in range(1):
        print(p)
        print(s)
        print(t)
        print(r)

    p1 = polynomial()
    p1.set_coef(0, 3)
    p1.set_coef(1, 4)
    p1.set_coef(2, 7)
    p1 = mul_polynomial(p1, polynomial((0, mul_inv(5))))
    print(p1)
    p2 = polynomial()
    p2.set_coef(0, 4)
    p2.set_coef(1, 3)
    print(p2)
    print(*extended_euclidean(p1, p2))
    _, s, t, _ = extended_euclidean(p1, p2)
    print(add_polynomial(mul_polynomial(p1, s), mul_polynomial(p2, t)))

    print(polynomial((0, -1)))
    """

    # 18.0x^0 + 2.0x^1  14.0x^0
    # print("p2-p1", *mod_polynomial(p2, p1))
    # print(gcd_polynomial(p2,p1))
