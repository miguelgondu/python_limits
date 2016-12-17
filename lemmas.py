import sympy
import math

def lemma1(f,p,q):
    r = sympy.degree(p)
    s = sympy.degree(q)
    if sympy.degree(f) >= r+s:
        raise ValueError('polynomial f must have degree less ' +
            ' that deg(p)+deg(q)')
    (phi, psi, gcd) = sympy.gcdex(p,q)
    if gcd != 1:
        raise ValueError('polynomials p and q are not relatively prime')
    l, h = sympy.div(f*psi, p)
    g = f*phi + l*q
    return g, h

def fgetter(F):
    '''
    This functions gets the homogenous components with respect to x.
    They are interpreted as polynomials in y.
    
    F: a polynomial in two variables x and y.
    returns: a list of homogenous components.
    
    TO-DO:
        - find a way to get the generators of a polynomial, so that
         I can implement this in a more general fashion.
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
    This function takes a polynomial f = pq with gcd(p,q) = 1 and
    returns the tuple (p, q).
    '''
    if f.is_irreducible:
        raise ValueError('f has no factorization')
    _list_of_factors = list(sympy.Mul.make_args(f.factor()))
    if len(_list_of_factors) < 2:
        raise ValueError('f might have no factorization ' + 
         'f=pq with gcd(p,q) = 1')
    if len(_list_of_factors) == 2:
        if sympy.gcd(_list_of_factors[0], _list_of_factors[1]) == 1:
            return _list_of_factors[0], _list_of_factors[1]
        else:
            raise ValueError('f has no factorization f=pq '+
                ' with gcd(p,q) = 1')
    else:
        first_factor = _list_of_factors[0]
        second_factor = 1
        for index in range(1, len(_list_of_factors)):
            second_factor = second_factor * _list_of_factors[index]
        if sympy.gcd(first_factor, second_factor) == 1:
            return first_factor, second_factor
        else:
            raise ValueError('f has no factorization f=pq ' + 
                ' with gcd(p,q) = 1')

def HenselsLemma(F, n):
    '''
    This function takes a polynomial in two variables F(x,y) monic in
    y such that F(0,y) = p(y)q(y), and finds two polynomials in
    two variables P and Q such that P(0,y) = p(y), Q(0,y) = q(y) and
    F = PQ up to deg(F). That is, this function lifts the
    factorization F(0,y) = p(y)q(y).
    
    F: a polynomial in two variables x and y.
    returns: a tuple (P,Q) such that F = PQ up to deg(F).
    
    TO-DO:
        - this function only accepts polynomials in explicitly x and y, we need one that accepts polys in 
          arbitrary generators.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    if sympy.LC(F, y) != 1:
        raise ValueError('F must be monic in y')

    list_of_ps = [0 for k in range(0,n+1)]
    list_of_qs = [0 for k in range(0,n+1)]

    degree_of_F_x = sympy.degree(F, x)

    list_of_fs = [0 for k in range(0,n+1)]
    for k in range(0, degree_of_F_x + 1):
        try:
            list_of_fs[k] = fgetter(F)[k]
        except:
            pass


    f0 = F.subs(x, 0)
    p, q = pqgetter(f0)
    q1, p1 = lemma1(list_of_fs[1], p, q)
    
    list_of_ps[0] = p
    list_of_ps[1] = p1
    list_of_qs[0] = q
    list_of_qs[1] = q1

    list_of_ms = [0]
    for i in range(2,n+1):
        ri = 0
        for j in range(1,i):
            ri = ri + list_of_qs[i-j]*list_of_ps[j]
        mi = list_of_fs[i] - ri
        list_of_ms.append(mi)
        qi, pi = lemma1(mi, p, q)
        list_of_ps[i] = sympy.Poly(pi,x,y)
        list_of_qs[i] = sympy.Poly(qi,x,y)
    P = 0
    Q = 0

    for k in range(0,n+1):
        Q = Q + sympy.Poly(x**k, x, y)*_list_of_qs[k]
        P = P + sympy.Poly(x**k, x, y)*_list_of_ps[k]

    return P, Q
