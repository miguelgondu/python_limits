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
    k = 0
    while sympy.LC(new_h, y) not in sympy.CC:
        # if n == None:
        #     n = random.randint(1, 10)
        # # print('n is {}'.format(n))
        new_h = h.as_expr().subs([(x, x + k*y), (y, -k*x + y)],
                                 simultaneous=True)
        k += 1
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

def newton_inverse(h, q, p):
    '''
    This function performs the Newton Automorphism (x maps to x ** q
    and y maps to y(x ** p)) to a polynomial h, returning a sympy expr
    (not a polynomial). Here, p and q are expected to be rational
    numbers (i.e. sympy.Rational objects).
    '''
    return newton_automorphism(h, sympy.Rational(1,q), sympy.Rational(-p, q))

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
    b1 = sympy.Rational(1,sympy.degree(h, y))*list_of_bs[-2]
    #print(sympy.degree(h, y))
    #print(h)
    #print(list_of_bs)
    return substitute_in_poly(h, [(y, y- b1)]), b1

def phi_inverse(h, b1):
    '''
    This function performs the inverse phi automorphism to a polynomial that's monic in y,
    so that the term that's with degree(h, y) - 1 banishes.
    '''
    y = sympy.Symbol('y')
    #print(sympy.degree(h, y))
    #print(h)
    #print(list_of_bs)
    return substitute_in_poly(h, [(y, y + b1)])

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
    and returns N = ceil(e^{(d/e)}u_d)
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
    if sympy.degree(f.simplify(), y) >= r+s:
        print('p: ' + str(p))
        print('q: ' + str(q))
        raise ValueError('polynomial {} must have degree less that deg(p)+deg(q)'.format(f.expand().simplify()))
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

    -Modificar este para que siempre saque un par, así sea irreducible el polinomio.
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
    dict_of_roots = sympy.roots(f)
    list_of_real_roots = [key for key in dict_of_roots if key in sympy.RR]
    if list_of_real_roots == []:
        raise ValueError('There are no real roots')
    r = list_of_real_roots[0]
    p = (y - r) ** dict_of_roots[r]
    q, res = sympy.div(f, p)
    print('res: ' + str(res.simplify()))
    return p, q
    

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
            ri = ri + (list_of_qs[i-j]*list_of_ps[j]).expand().simplify()
        if i < degree_of_F_x:
            mi = list_of_fs[i] - ri
        elif i >= degree_of_F_x:
            mi = (-1)*ri
        list_of_ms.append(mi.expand().simplify().evalf())
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

from queue import LifoQueue
from operator import mul

def lcm(a,b):
    return a*b//math.gcd(a,b)

def rotate(f, n):
    return f.subs([(x, x + n*y), (y, -x + n*y)], simultaneous=True)

def limit_poly2(f, g, N):
    list_of_tray = []
    h = h_getter(f, g)
    print('h: ' + str(h.expand().simplify()))
    h = monic_maker(h)
    print('monic h: ' + str(h.expand().simplify()))
    # try:
    #     N = upper_bound_getter(h)
    #     print('N: ' + str(N))
    # except OverflowError:
    #     list_of_factors = sympy.factor_list(h)[1]
    #     for _tuple in list_of_factors:
    #         a, b = _tuple
    #         if a == y:
    #             h = (h / (y**b)).simplify()
    #             N = upper_bound_getter(h)
    #             print('N: ' + str(N))
    #             list_of_tray.append(y)
    #             break
    print('h: ' + str(h))
    h, b1 = phi_automorphism(h)
    print('h after phi: ' + str(h.expand().simplify()) + ' with inverse ' + str(b1))
    r, ur = ur_getter(h)
    print('r is: ' + str(r))
    print('ur is: ' + str(ur))
    if ur == sympy.simplify('oo'):
        return []
    d = sympy.degree(h, y)
    h = newton_automorphism(h, r, ur)/(x**(d*ur))
    print('h after newton: ' + str(h.expand().simplify()))

    if sympy.degree(h, x) == 0:
        try:
            P, Q = pqgetter(h)
        except ValueError:
            return []
        s = sympy.degree(P, y)
        t = sympy.degree(Q, y)
        P = (phi_inverse(newton_inverse(P, r, ur), b1).subs(x, x**r)) * x**(s*ur)
        Q = (phi_inverse(newton_inverse(Q, r, ur), b1).subs(x, x**r)) * x**(t*ur)
    else:
        P, Q = Hensels_lemma(h.simplify(), N)
        s = sympy.degree(P, y)
        t = sympy.degree(Q, y)
        P = (phi_inverse(newton_inverse(P, r, ur), b1).subs(x, x**r)) * x**(s*ur)
        Q = (phi_inverse(newton_inverse(Q, r, ur), b1).subs(x, x**r)) * x**(t*ur)
    stack = LifoQueue()
    stack.put(P.expand().simplify().evalf())
    stack.put(Q.expand().simplify().evalf())
    list_of_rs = [r]

    # We find the trayectiories
    while not stack.empty():
        H = stack.get().simplify()
        print('From the stack I got: ' + str(H.simplify()))
        H = H.simplify().expand()
        # if H.subs([(x, 0), (y, 0)], simultaneous=True) != 0:
        #     print('H doesn\'t go through 0!, because H(0,0) = {}'.format(
        #         H.subs([(x, 0), (y, 0)], simultaneous=True)))
        #     continue
        if sympy.degree(H, y) == 1:
            print('H works!')
            list_of_tray.append(H)
            continue
        
        print('H needs more work')

        H = monic_maker(H)
        print('monic H: ' + str(H.expand().simplify()))
        # N = upper_bound_getter(H)
        H, b1 = phi_automorphism(H)
        print('H after phi: ' + str(H.expand().simplify()))
        r, ur = ur_getter(H)
        list_of_rs.append(r)
        # print('r is: ' + str(r))
        # print('ur is: ' + str(ur))
        if ur == sympy.simplify('oo'):
            print('H is of the form y^d!')
            list_of_tray.append(sympy.simplify(0))
            continue
        d = sympy.degree(H, y)
        H = newton_automorphism(H, r, ur)/(x**(d*ur))
        print('H after newton: ' + str(H.expand().simplify()))
        if sympy.degree(H, x) == 0:
            try:
                P, Q = pqgetter(H)
            except ValueError:
                print('H has no real roots!')
                continue
            s = sympy.degree(P, y)
            t = sympy.degree(Q, y)
            P = (phi_inverse(newton_inverse(P, r, ur), b1).subs(x, x**r)) * x**(s*ur)
            Q = (phi_inverse(newton_inverse(Q, r, ur), b1).subs(x, x**r)) * x**(t*ur)
        else:
            P, Q = Hensels_lemma(h.simplify(), N)
            s = sympy.degree(P, y)
            t = sympy.degree(Q, y)
            P = (phi_inverse(newton_inverse(P, r, ur), b1).subs(x, x**r)) * x**(s*ur)
            Q = (phi_inverse(newton_inverse(Q, r, ur), b1).subs(x, x**r)) * x**(t*ur)
            try:
                P, Q = Hensels_lemma(H.simplify(), N)
            except ValueError:
                print('H has no real roots!')
                continue
        print('Putting {} in the stack'.format(P))
        stack.put(P.expand().simplify().evalf())
        print('Putting {} in the stack'.format(Q))
        stack.put(Q.expand().simplify().evalf())

    # We evaluate the one-variable limit.
    new_list_of_tray = [-(tray - y) for tray in list_of_tray]
    print(new_list_of_tray)
    R = reduce(lcm, list_of_rs)
    print('R: ' + str(R))
    F = rotate(f, 1)
    G = rotate(g, 1)
    print('F: ' + str(F))
    print('G: ' + str(G))
    list_of_quotients = [(F.subs([(x, x), (y, new_list_of_tray[i])])
                          /G.subs([(x, x), (y, new_list_of_tray[i])]))
                          for i in range(len(list_of_tray))]
    print(list_of_quotients)
    list_of_limits = [sympy.limit(q, x, 0) for q in list_of_quotients]
    print(list_of_limits)
    return list_of_limits


