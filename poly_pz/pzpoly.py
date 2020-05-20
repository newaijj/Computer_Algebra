"""
structure of polynomial: {exp: coef}
"""
from collections import defaultdict


class polynomial:
    def __init__(self):
        self.eqn = defaultdict(lambda: 0)

    def set_coef(self, coef, power):
        if coef == 0:
            return
        self.eqn[power] = coef

    def __str__(self):
        # tested to be 2x faster for n=100 than loop
        return " + ".join(
            [
                str(coef) + "x^" + str(exp)
                for exp, coef in sorted(self.eqn.items())
            ]
        )

    def __add__(self, other):
        result = polynomial()
        result.eqn = self.eqn
        for exp, coef in other.eqn.items():
            result.eqn[exp] += coef
        return result

    def __mul__(self, other):
        result = polynomial()
        for e1, co1 in self.eqn.items():
            for e2, co2 in other.eqn.items():
                result.eqn[e1 + e2] += co1 * co2
        return result

    def __lshift__(self, other):
        if other[0]:
            self.eqn[other[1]] = other[0]
        return self


p1 = polynomial()
p1 << [1, 1] << [0, 2] << [-3, 4]
print("p1", p1)
print("p1*p1", p1 * p1)
