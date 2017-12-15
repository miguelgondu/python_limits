import math
import random
import re
import sympy
from functools import reduce

x = sympy.Symbol('x')
y = sympy.Symbol('y')

# Auxiliar functions

def substitute_in_poly(p, list_of_subst):
    '''
    Substitutes substitutions (in the form of tuples) from list_of_subst
    in the sympy polynomial p. 
    '''
    return p.as_expr().subs(list_of_subst, simultaneous=True)

def h_getter(f, g):
    return (-x*(g*sympy.diff(f, y) - f*sympy.diff(g, y))
            + y*(g*sympy.diff(f, x) - f*sympy.diff(g, x))).simplify()

def monic_maker(h, n=None):
    '''
    This function peforms a rotation on a polynomial of two variables
    h(x,y). That is, it maps h(x,y) to h(x+ny, -nx + y)) for some n
    such that the resulting polynomial is monic in y.

    TO-DO:
        - implement a smarter choosing of n, perhaps considering the
          roots of the polynomial.
    '''
    new_h = h
    while sympy.LC(new_h, y) not in sympy.CC:
        if n == None:
            n = random.randint(1, 10)
        # print('n is {}'.format(n))
        new_h = h.as_expr().subs([(x, x + n*y), (y, -n*x + y)],
                                 simultaneous=True)
    new_h = new_h * (sympy.Rational(1, sympy.LC(new_h, y)))
    return new_h

def newton_automorphism(h, q, p):
    '''
    This function performs the Newton Automorphism (x maps to x ** q
    and y maps to y(x ** p)) to a polynomial h, returning a sympy expr
    (not a polynomial). Here, p and q are expected to be rational
    numbers (i.e. sympy.Rational objects).
    '''
    return substitute_in_poly(h, [(x, x ** q), (y, y*x**p)])

def components_of_y(F):
    '''
    This function returns a list with the homogenous
    components of h as a polynomial in y in increasing order.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    list_of_Fprimes = [F]
    list_of_fs = [F.subs(y, 0)]
    for j in range(1, sympy.degree(F, y)+1):
        nextFprime = sympy.diff(list_of_Fprimes[j-1], y)
        list_of_Fprimes.append(nextFprime)
        list_of_fs.append(sympy.Rational(1, math.factorial(j))*nextFprime.subs(y, 0))
    return list_of_fs

def components_of_x(F):
    '''
    This function returns a list with the homogenous
    components of h as a polynomial in y in increasing order.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    list_of_Fprimes = [F]
    list_of_fs = [F.subs(x, 0)]
    for j in range(1, sympy.degree(F, x)+1):
        # print('the past prime is: ' + str(list_of_Fprimes[j-1]))
        nextFprime = sympy.diff(list_of_Fprimes[j-1], x)
        list_of_Fprimes.append(nextFprime)
        list_of_fs.append(sympy.Rational(1, math.factorial(j))*nextFprime.subs(x, 0))
    return list_of_fs

def phi_automorphism(h):
    '''
    This function performs the phi automorphism to a polynomial that's monic in y,
    so that the term that's with degree(h, y) - 1 banishes. Also returns the
    term that was substracted for the construction of the inverse.
    '''
    y = sympy.Symbol('y')
    list_of_bs = components_of_y(h)
    b1 = list_of_bs[-2]
    #print(sympy.degree(h, y))
    #print(h)
    #print(list_of_bs)
    return substitute_in_poly(h, [(y, y- sympy.Rational(1,sympy.degree(h, y))*b1.as_expr())]), b1

def least_degree_in_x(b):
    '''
    Returns the least degree in x of a polynomial b(x) (assumed to only be a
    polynomial in x, not x and y).
    '''
    x = sympy.Symbol('x')
    if b == 0:
        return sympy.simplify('oo')
    if b.subs(x, 0) == b:
        return 0
    k = 0
    # print(b)
    while sympy.simplify(b/x**k).subs(x, 0) == 0:
        k += 1
    return k

def ur_getter(h):
    '''
    Given a polynomial
        h(x,y) = y ** d + b_2(x)y**(d-2) + ... + b_d(x)
    this function returns the tuple (r, u_r) where u_r is the
    degree of b_r such that u_r/r is minimal.
    '''
    list_of_bs = components_of_y(h)
    list_of_bs = list_of_bs[::-1] # To follow the convention of the proof.
    # Here, position 0 is really from b2.
    list_of_least_degrees = [least_degree_in_x(b) for b in list_of_bs[2:]] 
    # print(list_of_least_degrees)
    list_of_quotients = [ui/(i+2) for (i, ui) in enumerate(list_of_least_degrees)]
    r_indexer = list_of_quotients.index(min(list_of_quotients))
    return r_indexer + 2, list_of_least_degrees[r_indexer]

def poly_mod(h, n):
    '''
    This function takes a polynomial h in one variable and returns what rests
    after eliminating every monomial whose degree is bigger than n (i.e.
    after erasing every x**k with k >= n).
    '''
    if len(h.gens) > 1:
        raise ValueError('Polynomial should only have one generator')
    if n > sympy.degree(h):
        '''
        If n is bigger than the degree, there's nothing to erase.
        '''
        return h

    x = h.gen
    list_of_coeffs = h.all_coeffs()[::-1] # i.e. list[k] is the coeff of x**k.
    list_of_modulus_coeffs = list_of_coeffs[:n]
    return sympy.Poly(list_of_modulus_coeffs[::-1], x)

def get_denominators_of_x(sympy_expression):
    '''
    this function takes a sympy EXPRESSION f of the form
    $f(x) = \sum c_i x^{p_i}$, where $p_i$s are rational, and returns
    a list with all the denominators of said rationals (1 is not conunted).

    To-Do:
        - Implement the possibility of another generator instead of x.
    '''
    expr_string = sympy.srepr(sympy_expression)
    list_of_matches = re.findall(r"Pow\(Symbol\('x'\), Rational\(\d+, \d+\)",
                                 expr_string)
    list_of_exponents = [int(string[29:-1]) for string in list_of_matches]
    return list_of_exponents

def gcd_for_lists(list_of_denominators):
    '''
    This function returns the gcd of a list of integers.
    '''
    return reduce(math.gcd, list_of_denominators)

def upper_bound_getter(h):
    '''
    This function takes a polynomial h(x,y) = y^d + a_1(x)y^{d-1} + ... + a_d(x)
    and returns N = ceil(e^{d/e}u_d)
    '''
    list_of_as = components_of_y(h)
    u_d = least_degree_in_x(list_of_as[0])
    return math.ceil(math.e ** (sympy.degree(h, y)/math.e) * u_d)

# Hensel's lemma lemmas.
def lemma_1(f,p,q):
    '''
    f, p, q are polynomials in one variable (say K[y]).

    Necesitamos implementar el caso f == 0 aparte, ¿no?
    '''
    r = sympy.degree(p)
    s = sympy.degree(q)
    if sympy.degree(f.as_poly(y), y) >= r+s:
        raise ValueError('polynomial f must have degree less that deg(p)+deg(q)')
    (phi, psi, gcd) = sympy.gcdex(p.as_poly(y),q.as_poly(y))
    phi = phi.as_expr()
    psi = psi.as_expr()
    if gcd != 1:
        raise ValueError('polynomials p and q are not relatively prime')
    l, h = sympy.div(f*psi, p)
    g = f*phi + l*q
    return g, h

def fgetter(F):
    '''
    This functions gets the homogenous components with respect to x. They are
    interpreted as polynomials in y.
    
    F: a polynomial in two variables x and y.
    returns: a list of homogenous components.
    
    TO-DO:
        - find a way to get the generators of a polynomial, so that I
         can implement this in a more general fashion.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    list_of_Fprimes = [F]
    list_of_fs = [F.eval(x, 0)]
    for j in range(1, sympy.degree(F, x)+1):
        nextFprime = sympy.diff(list_of_Fprimes[j-1], x)
        list_of_Fprimes.append(nextFprime)
        list_of_fs.append((1/math.factorial(j))*nextFprime.eval(x, 0))
    return list_of_fs

def pqgetter(f):
    '''
    This function takes a polynomial f = pq with gcd(p,q) = 1 and returns the tuple (p, q).
    '''
    # if f.is_irreducible:
    #     raise ValueError('f has no factorization')
    # _list_of_factors = list(sympy.Mul.make_args(f.factor()))
    # if len(_list_of_factors) < 2:
    #     raise ValueError('f might have no factorization f=pq with gcd(p,q) = 1')
    # if len(_list_of_factors) == 2:
    #     if sympy.gcd(_list_of_factors[0], _list_of_factors[1]) == 1:
    #         return _list_of_factors[0], _list_of_factors[1]
    #     else:
    #         raise ValueError('f has no factorization f=pq with gcd(p,q) = 1')
    # else:
    #     first_factor = _list_of_factors[0]
    #     second_factor = 1
    #     for index in range(1, len(_list_of_factors)):
    #         second_factor = second_factor * _list_of_factors[index]
    #     if sympy.gcd(first_factor, second_factor) == 1:
    #         return first_factor, second_factor
    #     else:
    #         raise ValueError('f has no factorization f=pq with gcd(p,q) = 1')
    list_of_factors = sympy.factor_list(f)[1]
    p, mult_p = list_of_factors[0]
    if sympy.degree(p, y) > 1:
        raise ValueError('There are no real roots')
    q = f/(p**mult_p)
    return p**mult_p, q

def Hensels_lemma(F, n):
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    if sympy.LC(F, y) != 1:
        raise ValueError('F must be monic in y')

    list_of_ps = [0 for k in range(0,n+1)]
    list_of_qs = [0 for k in range(0,n+1)]

    degree_of_F_x = sympy.degree(F, x)

    # list_of_fs = [0 for k in range(0,n+1)]
    # for k in range(0, degree_of_F_x+1):
    #     try:
    #         list_of_fs[k] = fgetter(F)[k]
    #     except:
    #         pass

    list_of_fs = components_of_x(F)
    print('list of fs: ')
    print(list_of_fs)


    f0 = F.subs(x, 0)
    p, q = pqgetter(f0)
    q1, p1 = lemma_1(list_of_fs[1], p, q.simplify())
    
    list_of_ps[0] = p
    list_of_ps[1] = p1
    list_of_qs[0] = q
    list_of_qs[1] = q1

    list_of_ms = [list_of_fs[1], list_of_fs[2] - p1*q1]
    for i in range(2,n+1):
        print('for i={}'.format(i))
        print('pi={}'.format(list_of_ps[i]))
        print('qi={}'.format(list_of_qs[i]))
        ri = 0
        for j in range(1,i):
            ri = ri + list_of_qs[i-j]*list_of_ps[j]
        if i < degree_of_F_x:
            mi = list_of_fs[i] - ri
        elif i >= degree_of_F_x:
            mi = (-1)*ri
        list_of_ms.append(mi)
        print('mi={}'.format(list_of_ms[i-1]))
        qi, pi = lemma_1(mi, p, q.simplify())
        list_of_ps[i] = pi
        list_of_qs[i] = qi
    P = 0
    Q = 0

    for k in range(0,n+1):
        Q = Q + x**k*list_of_qs[k]
        P = P + x**k*list_of_ps[k]
    
    print('list of fs: ' + str(list_of_fs))
    print('list of ps: ' + str(list_of_ps))
    print('list of qs: ' + str(list_of_qs))
    print('list of ms: ' + str(list_of_ms))

    return P, Q

def agrupador(g):
    list_of_factors = sympy.factor_list(g)[1]
    dict_of_imaginary_roots = {}
    dict_of_real_roots = {}
    for _tuple in list_of_factors:
        p, mult_p = _tuple
        list_of_roots_of_p = sympy.solve(p)
        for root in list_of_roots_of_p:
            if root in sympy.CC and root not in sympy.RR:
                dict_of_imaginary_roots[root] = mult_p
            else:
                dict_of_real_roots[root] = mult_p

    real_pol = 1
    imaginary_pol = 1

    for real_root in dict_of_real_roots:
        real_pol *= (y - real_root) ** dict_of_real_roots[real_root]

    imaginary_root_list = []

    for im_root in dict_of_imaginary_roots:
        if sympy.conjugate(im_root) in imaginary_root_list:
            continue
        else:
            imaginary_root_list.append(im_root)
            imaginary_pol *= ((y - im_root)*(y-sympy.conjugate(im_root).simplify())).expand() ** dict_of_imaginary_roots[im_root]
    print(imaginary_root_list)
    return real_pol, imaginary_pol
            


def limit_poly(f, g):
    # First, we get the quotient's variety divider h
    h = h_getter(f,g)
    print('h: ' + str(h.expand()))

    # Rotating it so that it's monic
    new_h = monic_maker(h, 1)
    print('h after rotating: ' + str(new_h.expand()))

    # Removing the d-1th term in y
    new_h, b1 = phi_automorphism(new_h)
    print('h after phi: ' + str(new_h.simplify().expand())) # Change everything to sympy expressions.

    # Applying Newton's

    # Añadir un if que verifique si ya se factoriza.
    r, ur = ur_getter(new_h)
    print('r is: ' + str(r))
    print('ur is: ' + str(ur))
    h_newton = newton_automorphism(new_h, r, ur)
    print('h after Newton\'s: ' + str(h_newton.simplify().expand()))

    h_newton = h_newton.simplify()
    # h_newton = h_newton.as_poly()
    d = sympy.degree(h_newton, y)
    print('This is h: ' + str(h_newton.factor()/x**(d*ur)))
    h_newton = h_newton.factor()/x**(d*ur)
    new_h = h_newton.subs(x, 0)
    list_of_factors = sympy.factor_list(new_h)[1]
    print(list_of_factors)

    p, q = pqgetter(new_h)
    # print('p is: ' + str(p))
    # print('q is: ' + str(q.simplify()))

    print(components_of_x(h_newton))

    # P, Q = Hensels_lemma(h_newton, 3)
    # print('P is: {}'.format(P))
    # print('Q is: {}'.format(Q))
    # # Reimplementar el lema de Hensel para que acepte expresiones en vez de polinomios.
    # return P, Q, h_newton
