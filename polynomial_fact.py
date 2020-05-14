import numpy as np
from logging import warning
from itertools import takewhile
import copy

"""
CLASS IMPLEMENTED: POLYNOMIAL
    -univariable only

- ATTRIBUTES
    - size  - total number of coefficients stored (defaults to 64)
    - min_power     - minimum power coefficient that can be stored (defaults to -32)
    - max_power     - maximum power coefficient that can be stored (defaults to 31)

- METHODS
    - set_coef(power, coefficient)
    - get_coef(power)
    - max_nonzero_pow()     - maximum power with non-zero coefficient stored

- FUNCTIONS
    - add_polynomial() -        accepts any number of polynomial objects as input
    - mul_polynomial() -        accepts any number of polynomial objects as input
    - sub_polynomial(p1,p2) -   (p1 - p2)
    - mod_polynomial(p1,p2) -   (p1 % p2); returns (p_quotient , p_remainder) 


"""
"""

def np_shift(np_arr, val):
    assert np_arr.ndim == 1
    assert np_arr.shape[0] > val

    shift = np.array(np_arr)
    np_arr[0 : -val - 1] = shift[val:-1]
    np_arr[-val:-1] = shift[0 : val - 1]
    return np_arr
"""


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


def mod_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)

    p_q = polynomial(size=max(p1.size, p2.size))
    p_r = copy.copy(p1)
    for p in range(p1.degree, p2.degree - 1, -1):
        fac = p_r.get_coef(p) / p2.get_coef(p2.degree)
        p_q.set_coef(p - p2.degree, fac)
        if fac == 0:
            continue
        p_fac = polynomial(size=p2.size)
        p_fac.set_coef(p - p2.degree, fac)
        p_r = sub_polynomial(p_r, mul_polynomial(p2, p_fac))

    return p_q, p_r


class polynomial:
    # polynomial up to 64 terms
    def __init__(self, size=64):
        self.coef = np.zeros(size)
        self.size = size
        self.degree = 0

    def update_degree(self, power=None):
        if power == None:
            self.degree = self.max_nonzero_pow()
        else:
            assert float(power).is_integer()
            if self.degree < power:
                self.degree = power

    def set_coef(self, power, coef):
        assert power >= 0 and power < self.size
        self.coef[power] = coef
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


p1 = polynomial()
p1.set_coef(0, 3)
p1.set_coef(2, 7)
p2 = add_inv_polynomial(p1)
p2 = mul_polynomial(p2, p1)
p2.set_coef(3, 1)
print(p1)
print(p2)
print("p2-p1", *mod_polynomial(p2, p1))
