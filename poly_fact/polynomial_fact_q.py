from polynomial import *
from gen_primes import *
from polynomial_fact_finite import *
from factor_tree import *
from hensel import Multifactor_Hensel_Lifting
import logging


def deriv_polynomial(p):
    assert isinstance(p, polynomial)
    p_out = p.copy_meta()

    p_out.set_coef(p.degree, 0)
    for i in range(p.degree):
        p_out.set_coef(i, p.get_coef(i + 1) * (i + 1))
    p_out.update_degree()
    return p_out


def mod_const(p, n):
    assert isinstance(p, polynomial) and float(n).is_integer()
    p_ = p.copy_meta()

    for i in range(p.degree + 1):
        p_.set_coef(i, p.get_coef(i) % n)
    return p_


def poly_fact_z(f):
    assert isinstance(f, polynomial)
    assert f.degree > 0

    output = set()
    logger.debug("poly_fact_z called")

    # 1
    if f.degree == 1:
        output.add(f)
        return output
    b = f.lc()
    A = f.max_norm()
    n = f.degree
    B = ((n + 1) ** 0.5) * (2 ** n) * A
    C = ((n + 1) ** (2 * n)) * (A ** (2 * n - 1))
    y = np.ceil(2 * (np.log(C) / np.log(2)))

    logger.debug("#1 completed")

    # 2
    prime = gen_primes()

    while True:
        p = next(prime)
        f_bar = mod_const(f, p)
        while not (
            b % p != 0
            and gcd_polynomial(
                f_bar, mod_const(deriv_polynomial(f_bar), p), N=p
            ).degree
            == 0
        ):  # equal_polynomial(gcd_polynomial(f_bar,mod_const(deriv_polynomial(f_bar),p),N=p),polynomial((0,1)),N=p))):
            p = next(prime)
            f_bar = mod_const(f, p)
            assert p <= 2 * y * np.log(y)
        l = np.ceil(np.log(2 ** (n ** 2) * B ** (2 * n)) / np.log(p))
        i = 0
        while l >= 2 ** i:
            i += 1
        l = 2 ** i

        logger.debug("#2 completed: possible p={}".format(p))

        # 3
        b_inv_poly = polynomial((0, mul_inv(b, N=p)), N=p)
        f_lc = mul_polynomial(f, b_inv_poly, N=p)

        print("here")
        print(f_lc, p)
        h_list = factorise_polynomial_int_finite(f_lc, N=p)
        cond = True
        for h in h_list:
            if h.max_norm() > np.ceil(p / 2):
                cond = False
        if cond:
            break

    b_poly = polynomial((0, b), N=p)
    facs = [b_poly] + h_list
    assert equal_polynomial(mul_polynomial(*facs, N=p), f, N=p)

    logger.debug("#3 completed: factors:")
    for h in facs:
        logger.debug(h)

    # 4
    f_tree = make_tree(h_list)
    print("s:", f_tree.s)
    print("t:", f_tree.t)
    b, f_tree = Multifactor_Hensel_Lifting(p, f, mul_inv(b, N=p), l, f_tree)

    return f


if __name__ == "__main__":

    # create logger
    logger = logging.getLogger("polynomial_fact_q")
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

    p = polynomial((1, 3))
    p.set_coef(2, 4)
    p.set_coef(0, 2)
    print("p: ", p)
    print("return: ", poly_fact_z(p))
