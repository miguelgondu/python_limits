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
    F = F.as_poly()
    list_of_Fprimes = [F]
    list_of_fs = [F.eval(x, 0)]
    for j in range(1, sympy.degree(F, x)+1):
        next_F_prime = sympy.diff(list_of_Fprimes[j-1], x)
        list_of_Fprimes.append(next_F_prime)
        list_of_fs.append((1/math.factorial(j))*next_F_prime.eval(x, 0).as_poly(y))
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

    Question:
        -must it be strictly positive?
    '''
    x = sympy.Symbol('x')
    k = 0
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

def limit_poly(f, g):
    # First, we get the quotient's variety divider h
    h = h_getter(f, g)
    h = monic_maker(h)
