from polynomial import *
from gen_primes import *
from polynomial_fact_finite import *
from factor_tree import *
from hensel import Multifactor_Hensel_Lifting
from short_vector import *
import logging


def deriv_polynomial(p):
    assert isinstance(p, polynomial)
    p_out = p.copy_meta()

    p_out.set_coef(p.degree, 0)
    for i in range(p.degree):
        p_out.set_coef(i, p.get_coef(i + 1) * (i + 1))
    p_out.update_degree()
    return p_out


def poly_2_vec(p, degree):
    assert isinstance(p, polynomial) and float(degree).is_integer()
    return vector(p.coef[0:degree])


def vec_2_poly(v, N=N):
    assert isinstance(v, vector)
    p_ = polynomial(N=N)
    for i in range(v.shape[0]):
        p_.set_coef(i, v.val[i])
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

        i = 0
        while i < 100:
            h_list = factorise_polynomial_int_finite(f_lc, N=p)
            cond = True
            for h in h_list:
                if h.max_norm() > (p / 2):
                    cond = False
            if cond:
                break
            i += 1
        if cond:
            break

    b_poly = polynomial((0, b), N=p)
    facs = [b_poly] + h_list
    assert equal_polynomial(mul_polynomial(*facs, N=p), f, N=p)

    logger.debug("#3 completed: factors:")
    for h in facs:
        logger.debug(h)

    # 4
    print("l: ", l)
    print("f: ", f)
    print("mul_inv(b, N=p): ", mul_inv(b, N=p))
    print("p: ", p)

    f_tree = make_tree(h_list, N=p)
    _, f_tree = Multifactor_Hensel_Lifting(p, f, mul_inv(b, N=p), l, f_tree)
    g = [t.value for t in f_tree.get_leaves()]
    for gi, hi in zip(g, h_list):
        print("p: ", p)
        print(gi)
        print(hi)
        assert gi.max_norm() < p ** l / 2
        assert equal_polynomial(gi, hi, N=p)

    logger.debug("#4 completed")

    # 5
    T = set(range(len(h_list)))
    G = set()
    f_ = copy.deepcopy(f)

    logger.debug("#5 completed")
    logger.debug("#6 completed")

    # 6
    while len(T) > 0:

        # 7
        u = h_list[max(T, key=lambda p: h_list[p].degree)]
        d = u.degree
        n_ = f_.degree

        for i in T:
            print(h_list[i])
        logger.debug("#7 completed")

        for j in range(d + 1, n_ + 1):
            # 8
            basis_ = basis(N=j)
            c = 0
            x_poly = polynomial((0, 1), N=p ** (l + 1))
            for k in range(0, (j - d)):
                new_poly = mul_polynomial(u, x_poly, N=p ** (l + 1))
                new_vec = poly_2_vec(new_poly, j)
                basis_.set_vec(k, new_vec)
                x_poly = mul_polynomial(
                    x_poly, polynomial((1, 1), N=p ** (l + 1)), N=p ** (l + 1)
                )
            p_poly = polynomial((0, p ** l), N=p ** (l + 1))
            for k in range((j - d), j):
                new_vec = poly_2_vec(p_poly, j)
                basis_.set_vec(k, new_vec)
                p_poly = mul_polynomial(
                    p_poly, polynomial((1, 1), N=p ** (l + 1)), N=p ** (l + 1)
                )

            g_vec = basis_reduction(basis_).get_vec(0)
            g_ = vec_2_poly(g_vec)

            logger.debug(
                "#8 completed: {} of {}: g*={}".format((j - d), (n_ - d), g_)
            )

            # 9
            h_ = polynomial((0, b), N=p ** l)
            for i in T:
                S = set()
                if equal_polynomial(
                    mod_polynomial(g_, h_list[i], N=p)[1], polynomial(N=p), N=p
                ):
                    S.add(i)
                else:

                    h_ = mul_polynomial(h_, g[i], N=p ** l)

            assert h_.max_norm() < ((p ** l) / 2)
            if g_.pp().one_norm() * h_.pp().one_norm() <= B:
                for i in S:
                    print("S: ", h_list[i])
                T = T - S
                G.add(g_.pp())
                f_ = h_.pp()
                b = f_.lc()
                logger.debug(
                    "#9 completed: {} of {}: g*={}".format(
                        (j - d), (n_ - d), g_
                    )
                )
                break
    G.add(f_)

    return G


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
    r = poly_fact_z(p)
    print("return: ")
    for f in r:
        print(f)
