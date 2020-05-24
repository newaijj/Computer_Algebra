from polynomial import *

"""
CLASS IMPLEMENTED: FACTOR_TREE
	- binary tree of factors of a polynomial (mod N)
		- polynomial at each node
		- 2 branches, factors
		- s,t at each node, solutions to p1*s + p2*t = 1 (mod N)
	- used for 15.17 (Multifactor Hensel Lifting)
	- __init__(self, value=None, parent=None)

ATTRIBUTES
	- L 	- Left branch
	- R 	- Right branch
	- P 	- Parent
	- value - Polynomial at current node
	- s
	- t     - solutions to p1*s + p2*t = 1 (mod N)

METHODS
	- get_walk(self)	- returns list of nodes at that level and below (like os.walk)

	- set_value(self, t1, t2, N=N) 	- sets node based on two children branches
	- isleaf(self)					- returns bool based on whether L and R are None
	- get_L, get_R
	- set_L, set_R

FUNCTIONS
	- make_tree(factors, N=N)
		- INPUT: 	list of factors
			- must be coprime (product is squarefree)
		- OUTPUT: 	root node of factor_tree 


"""


class factor_tree:
    def __init__(self, value=None, parent=None):
        # assert isinstance(value,polynomial)
        self.L = None
        self.R = None
        self.P = parent
        self.value = value
        self.s = None
        self.t = None

    def isleaf(self):
        return self.L == None and self.R == None

    def set_L(self, L):
        assert isinstance(L, polynomial)
        self.L = L

    def set_R(self, R):
        assert isinstance(R, polynomial)
        self.R = R

    def get_L(self):
        return self.L

    def get_R(self):
        return self.R

    def set_value(self, t1, t2, N=N):
        assert isinstance(t1, factor_tree) and isinstance(t2, factor_tree)
        self.value = mul_polynomial(t1.value, t2.value, N=N)
        self.L = t1
        self.R = t2
        _, self.s, self.t, _ = extended_euclidean(t1.value, t2.value, N=N)
        t1.parent = self
        t2.parent = self

    def get_walk(self):
        out = [self]
        if self.R != None:
            out += self.R.get_walk()
        if self.L != None:
            out += self.L.get_walk()
        return out


def make_tree(factors, N=N):
    assert isinstance(factors, list)
    for fac in factors:
        assert isinstance(fac, polynomial)

    others = []
    f = factor_tree(factors[0])
    others.append((f, 1))
    for fac in factors[1:]:
        node_c = factor_tree(fac)

        if others[-1][1] == 1:
            i = 0
            while len(others) > 0 and others[-1][1] == i + 1:
                node_new = factor_tree()
                node_new.set_value(others[-1][0], node_c)
                others = others[:-1]
                node_c = node_new
                i += 1
            i += 1
            others.append((node_c, i))
        else:
            others.append((node_c, 1))

    for i in range(len(others) - 1, 0, -1):
        node_new = factor_tree()
        node_new.set_value(others[i - 1][0], others[i][0])
        others = others[:-2]
        others.append((node_new, 0))

    return others[0][0]


if __name__ == "__main__":
    from polynomial_fact_finite import *

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
    p_prob.set_coef(9, 1)
    p_rand = rand_polynomial(10)
    print("random_polynomial: ", p_prob)
    print("factorisation: ")
    facs = factorise_polynomial_int_finite(p_prob)
    for f in facs:
        print(f)
    print("verify: ", mul_polynomial(*facs))

    t = make_tree(facs)
    print("t.value: ", t.value)
    print("t.L.value: ", t.L.value)
    print("t.R.value: ", t.R.value)
    print("t.L.L.value: ", t.L.L.value)
    print("t.L.R.value: ", t.L.R.value)
    print("t.R.L.value: ", t.R.L.value)
    print("t.R.R.value: ", t.R.R.value)

    print("t.L.s: ", t.L.s)
    print("t.L.t: ", t.L.t)
    print(*extended_euclidean(t.L.L.value, t.L.R.value))
    print(*t.get_walk())
