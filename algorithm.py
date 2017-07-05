import math
import random
import sympy

x = sympy.Symbol('x')
y = sympy.Symbol('y')

# Auxiliar functions

def substitute_in_poly(p, list_of_subst):
    '''
    Substitutes substitutions (in the form of tuples) from list_of_subst
    in the sympy polynomial p. 
    '''
    return p.as_expr().subs(list_of_subst, simultaneous=True).as_poly()

def h_getter(f, g):
    return (x.as_poly()*(g*sympy.diff(f, x) - f*sympy.diff(g, x))
            - y.as_poly()*(g*sympy.diff(f, y) - f*sympy.diff(g, y)))

def monic_maker(h):
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
        n = random.randint(1, 10)
        new_h = h.as_expr().subs([(x, x + n*y), (y, -n*x + y)], simultaneous=True)
    new_h = new_h * (1/sympy.LC(new_h, y))
    return new_h.as_poly()

def newton_automorphism(h, q, p):
    '''
    This function performs the Newton Automorphism (x maps to x ** q
    and y maps to y(x ** p)) to a polynomial h, returning a sympy expr
    (not a polynomial). Here, p and q are expected to be rational
    numbers.
    '''
    return substitute_in_poly(h, [(x, x ** q), (y, y*x**p)])

def components_of_y(F):
    '''
    This function returns a list with the homogenous
    components of h as a polynomial in y in increasing order.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    F = F.as_poly()
    list_of_Fprimes = [F]
    list_of_fs = [F.eval(y, 0)]
    for j in range(1, sympy.degree(F, y)+1):
        nextFprime = sympy.diff(list_of_Fprimes[j-1], y)
        list_of_Fprimes.append(nextFprime)
        list_of_fs.append((1/math.factorial(j))*nextFprime.eval(y, 0).as_poly(x))
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
    return substitute_in_poly(h, [(y, y- (1/sympy.degree(h, y))*b1.as_expr())]), b1

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
    list_of_quotients = [k/(i+2) for (i, k) in enumerate(list_of_least_degrees)]
    r_indexer = list_of_quotients.index(min(list_of_quotients))
    return r_indexer + 2, list_of_least_degrees[r_indexer]
