import numpy as np
from logging import warning
from itertools import takewhile

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


def np_shift(np_arr, val):
    assert np_arr.ndim == 1
    assert np_arr.shape[0] > val

    shift = np.array(np_arr)
    np_arr[0 : -val - 1] = shift[val:-1]
    np_arr[-val:-1] = shift[0 : val - 1]
    return np_arr


def add_polynomial(*args):
    for p in args:
        assert isinstance(p, polynomial)
        assert p.size == args[0].size
    p_out = polynomial(size=args[0].size)

    for p in args:
        p_out.coef += p.coef

    return p_out


def add_inv_polynomial(p1):
    assert isinstance(p1, polynomial)
    p_out = polynomial(size=p1.size)
    p_out.coef = -p1.coef
    return p_out


def sub_polynomial(p1, p2):
    return add_polynomial(p1, add_inv_polynomial(p2))


def mul_2_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    assert p1.size == p2.size

    p_out = polynomial(size=p1.size)
    p_out.coef = np.convolve(p1.coef, p2.coef, "same")
    """
	shift = np.array(p_out.coef)
	p_out.coef[0:-2] = shift[1:-1]
	p_out.coef[-1] = shift[0]

	"""
    p_out.coef = np_shift(p_out.coef, 1)
    return p_out


def mul_polynomial(*args):
    for p in args:
        assert isinstance(p, polynomial)
        assert p.size == args[0].size
    p_out = args[0]

    for p in args[1:]:
        p_out = mul_2_polynomial(p_out, p)

    return p_out


"""
def mul_inv_polynomial(p1):
    assert isinstance(p1, polynomial)
    p_out = polynomial(size=p1.size)
    p_out.coef = np.divide(
        np.ones(p1.coef.size, dtype=p1.coef.dtype),
        p1.coef,
        out=np.zeros_like(p1.coef),
        where=p1.coef != 0.0,
    )
    p_out.coef = np.flip(p_out.coef)
    p_out.coef = np_shift(p_out.coef, -1)
    return p_out


def div_polynomial(p1, p2):
    assert isinstance(p1, polynomial) and isinstance(p2, polynomial)
    assert p1.size == p2.size

    p_out = mul_polynomial(p1, mul_inv_polynomial(p2))
    return p_out
"""


def mod_polynomial(p1, p2):
    for p in range(p1.min_pow, -1):
        assert p1.get_coef(p) == 0
    for p in range(p2.min_pow, -1):
        assert p2.get_coef(p) == 0
    assert p1.max_nonzero_pow() >= p2.max_nonzero_pow()

    p1_max = p1.max_nonzero_pow()
    p2_max = p2.max_nonzero_pow()
    p_q = polynomial()
    for p in range(p1_max, p2_max - 1, -1):
        fac = p1.get_coef(p) / p2.get_coef(p2_max)
        p_q.set_coef(p - p2_max, fac)
        if fac == 0:
            continue
        p_fac = polynomial(size=p2.size)
        p_fac.set_coef(p - p2_max, fac)
        p1 = sub_polynomial(p1, mul_polynomial(p2, p_fac))

    return p_q, p1


class polynomial:
    # polynomial up to 64 terms
    def __init__(self, size=64):
        self.coef = np.zeros(size)
        self.size = size
        self.min_pow = int(-size / 2)
        self.max_pow = int(size / 2 - 1)

    def set_coef(self, power, coef):
        assert abs(power < self.size / 2)
        """
      indexing of coef: [-(size/2)+1, -(size/2)+2, ..., (size/2)-1]
      index: exponent
      value at index: coefficient
      """
        self.coef[int(power + self.size / 2)] = coef

    def get_coef_index(self, index):
        return self.coef[index]

    def get_coef(self, power):
        return self.coef[self.get_index(power)]

    def get_pow(self, index):
        return int(
            index - self.size / 2
            if index >= (self.size / 2)
            else index - self.size / 2
        )

    def get_index(self, power):
        if power >= 0:
            assert power <= (self.size / 2) - 1
        else:
            assert power >= -(self.size / 2)
        return int(self.size / 2) + power

    def max_nonzero_pow(self):
        return self.get_pow(np.max(np.nonzero(self.coef)))

    def __str__(self):
        nzero = np.nonzero(self.coef)
        if nzero[0].size == 0:
            r = "0"
        else:
            r = ""
            it = np.nditer(nzero)
            nxt = next(it)
            r += (
                str(self.get_coef_index(nxt))
                + "x^"
                + str(self.get_pow(nxt))
                + " "
            )

            while True:
                try:
                    nxt = next(it)
                    r += (
                        "+ "
                        + str(self.get_coef_index(nxt))
                        + "x^"
                        + str(self.get_pow(nxt))
                        + " "
                    )
                except StopIteration:
                    break

        return r


import time

p1 = polynomial()
for i in range(-32, 32):
    p1.set_coef(i, i)
t0 = time.time()
ans = mul_polynomial(p1, p1)
t1 = time.time()

print("p1*p1: ", ans)
print(t1 - t0)
