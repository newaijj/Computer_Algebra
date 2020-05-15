import numpy as np
from logging import warning
from itertools import takewhile
import copy

"""
CLASS IMPLEMENTED: POLYNOMIAL
    -univariable only
    -coefficients mod N

- ATTRIBUTES
    - size  - total number of coefficients stored (defaults to 64)
    - degree - degree of polynomial (maintained on insertions, but not for direct access to p.coef)

- METHODS
    - set_coef(power, coefficient)
    - get_coef(power)
    - update_degree(power=None) - checks if degree needs to be updated. If power is specified, checks if that power is larger than current degree and sets it if it is
    - copy_meta()   - returns polynomial with copy of all meta variables but empty coefficients

- FUNCTIONS - POLYNOMIAL
    - add_polynomial() -                    accepts any number of polynomial objects as input
    - mul_polynomial() -                    accepts any number of polynomial objects as input
    - sub_polynomial(p1,p2) -               (p1 - p2)
    - mod_polynomial(p1,p2) -               (p1 % p2); returns (p_quotient , p_remainder) given global N
    - exp_polynomial(p,e)   -               returns (p^e), e must be integer
    - exp_polynomial_rem(p, e, f) -         returns (p^e mod f), e must be integer, p and f must be polynomials
    - gcd_polynomial(p1,p2) -               return GCD of the two polynomials given global N
    - equal_polynomial(p1,p2) -             returns if polynomials are equal

- FUNCTION - MISC
    - mul_inv(a)    -   calculates multiplicative inverse given global N
"""

N = 23


def add_polynomial(*args):
    max_size = 0
    for p in args:
        assert isinstance(p, polynomial)
        max_size = max(max_size, p.size)
    p_out = polynomial(size=max_size)

    for p in args:
        out = np.zeros(max_size)
        out[0 : (p.degree + 1)] = p.coef[0 : (p.degree + 1)]
        p_out.coef += out
        p_out.coef = p_out.coef % N

    p_out.update_degree()
    return p_out


def add_inv_polynomial(p1):
    assert isinstance(p1, polynomial)
    p_out = copy.deepcopy(p1)
    p_out.coef = -p_out.coef
    return p_out


def sub_polynomial(p1, p2):
    return add_polynomial(p1, add_inv_polynomial(p2))


def mul_2_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)

    np_out = np.convolve(p1.coef[0 : p1.degree + 1], p2.coef[0 : p2.degree + 1])
    p_out = polynomial(size=max(p1.size, p2.size, np_out.shape[0]))
    p_out.degree = np_out.shape[0] - 1
    p_out.coef[0 : np_out.shape[0]] = np_out
    p_out.coef = p_out.coef % N

    return p_out


def mul_polynomial(*args):
    max_size = 0
    for p in args:
        assert isinstance(p, polynomial)
        max_size = max(max_size, p.size)
    p_out = polynomial(size=max_size)
    p_out.coef[0 : args[0].degree + 1] = args[0].coef[0 : args[0].degree + 1]
    p_out.degree = args[0].degree

    for p in args[1:]:
        p_out = mul_2_polynomial(p_out, p)

    return p_out


# calc multiplicative inverse of a number mod N based on solving Bézout's identity
def mul_inv(a):
    assert a < N
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


def mod_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)

    p_q = polynomial(size=max(p1.size, p2.size))
    p_r = copy.copy(p1)
    for p in range(p1.degree, p2.degree - 1, -1):
        fac = (p_r.get_coef(p) * mul_inv(p2.get_coef(p2.degree))) % N
        p_q.set_coef(p - p2.degree, fac)
        if fac == 0:
            continue
        p_fac = polynomial(size=p2.size)
        p_fac.set_coef(p - p2.degree, fac)
        p_r = sub_polynomial(p_r, mul_polynomial(p2, p_fac))

    return p_q, p_r


# exponent by repeated squaring
def exp_polynomial(p, e):
    assert isinstance(p, polynomial) and float(e).is_integer()

    p_out = p.copy_meta()
    p_out.set_coef(0, 1)
    if e == 0:
        return p_out
    else:
        for a in bin(e).replace("0b", ""):
            p_out = mul_polynomial(p_out, p_out)
            if a == "1":
                p_out = mul_polynomial(p_out, p)
    return p_out


# exponent by repeated squaring, modded some polynomial f
def exp_polynomial_rem(p, e, f):
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
            p_out = mod_polynomial(mul_polynomial(p_out, p_out), f)[1]
            if a == "1":
                p_out = mod_polynomial(mul_polynomial(p_out, p), f)[1]
    return p_out


def gcd_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    if p2.degree > p1.degree:
        p1, p2 = p2, p1
    p_zero = p2.copy_meta()
    while not np.all(np.equal(p2.coef, p_zero.coef)):
        p_t = mod_polynomial(p1, p2)[1]
        p1 = p2
        p2 = p_t
    return p1


def equal_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    return np.all(np.equal(p1.coef, p2.coef))


class polynomial:
    # polynomial up to 64 terms
    def __init__(self, size=64):
        self.coef = np.zeros(size)
        self.size = size
        self.degree = 0

    def copy_meta(self):
        p = polynomial(size=self.size)
        p.size = self.size
        p.degree = self.degree
        return p

    def update_degree(self, power=None):
        if power == None:
            self.degree = self.max_nonzero_pow()
        else:
            assert float(power).is_integer()
            if self.coef[power] != 0 and self.degree < power:
                self.degree = power

    def set_coef(self, power, coef):
        assert power >= 0 and power < self.size
        self.coef[power] = coef % N
        self.update_degree(power=power)

    def get_coef(self, power):
        assert power >= 0 and power < self.size
        return self.coef[power]

    def max_nonzero_pow(self):
        try:
            return np.max(np.nonzero(self.coef))
        except ValueError:
            return 0

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
    p1 = polynomial()
    p1.set_coef(0, 3)
    p1.set_coef(2, 7)
    p2 = add_inv_polynomial(p1)
    p2 = mul_polynomial(p2, p1)
    p2.set_coef(3, 1)
    p3 = polynomial()
    p3.set_coef(0, 2)
    p3.set_coef(1, 3)
    p3.set_coef(2, 1)
    p4 = polynomial()
    p4.set_coef(0, 1)
    p4.set_coef(1, 3)
    p4.set_coef(2, 3)
    p4.set_coef(3, 1)
    print("p3:", p3)
    print("p4:", p4)
    print("gcd: ", gcd_polynomial(p3, p4))

    p5 = polynomial()
    p5.set_coef(0, 18)
    p5.set_coef(1, 2)
    p6 = polynomial()
    p6.set_coef(0, 14)

    print(equal_polynomial(p5, p5))
    print(equal_polynomial(p6, p6))
    print(6 * mul_inv(9) % N)
    # 18.0x^0 + 2.0x^1  14.0x^0
    # print("p2-p1", *mod_polynomial(p2, p1))
    # print(gcd_polynomial(p2,p1))
