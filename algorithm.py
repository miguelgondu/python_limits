import sympy
import random
import math

x = sympy.Symbol('x')
y = sympy.Symbol('y')

# Auxiliar functions

def hgetter(f,g):
    return x.as_poly()*(g*sympy.diff(f, x) - f*sympy.diff(g, x)) - y.as_poly()*(g*sympy.diff(f, y) - f*sympy.diff(g, y))

def monicMaker(h):
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
        n = random.randint(1,10)
        new_h = h.as_expr().subs([(x, x + n*y), (y, -n*x + y)], simultaneous=True)
    new_h = new_h * (1/sympy.LC(new_h, y))
    return new_h.as_poly()

def newtonAutomorphism(h, q, p):
    '''
    This function performs the Newton Automorphism (x maps to x ** q
    and y maps to y(x ** p)) to a polynomial h, returning a sympy expr
    (not a polynomial). Here, p and q are expected to be rational
    numbers.
    '''
    return h.as_expr().subs([(x, x ** q), (y, y*x**p)], simultaneous=True)

def componentsOfy(F):
    '''
    This function returns a list with the homogenous
    components of h as a polynomial in y in increasing order.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    list_of_Fprimes = [F]
    list_of_fs = [F.eval(y, 0)]
    for j in range(1, sympy.degree(F, y)+1):
        nextFprime = sympy.diff(list_of_Fprimes[j-1], y)
        list_of_Fprimes.append(nextFprime)
        list_of_fs.append((1/math.factorial(j))*nextFprime.eval(y, 0).as_poly(x))
    return list_of_fs

def componentsOfx(F):
    '''
    This function returns a list with the homogenous
    components of h as a polynomial in y in increasing order.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    list_of_Fprimes = [F]
    list_of_fs = [F.eval(x, 0)]
    for j in range(1, sympy.degree(F, x)+1):
        nextFprime = sympy.diff(list_of_Fprimes[j-1], x)
        list_of_Fprimes.append(nextFprime)
        list_of_fs.append((1/math.factorial(j))*nextFprime.eval(x, 0).as_poly(y))
    return list_of_fs

def phiAutomorphism(h):
    '''
    This function performs the phi automorphism to a polynomial that's monic in y,
    so that the term that's with degree(h, y) - 1 banishes. Also returns the 
    term that was substracted for the construction of the inverse.
    '''
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    list_of_bs = componentsOfy(h)
    b1 = list_of_bs[-2]
    #print(sympy.degree(h, y))
    #print(h)
    #print(list_of_bs)
    return h.as_expr().subs([(y, y-b1/sympy.degree(h, y).as_expr())], simultaenous=True).as_poly(x,y), b1

def leastDegreeInx(b):
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    k = 0
    while(sympy.simplify(b/x**k).subs(x, 0) == 0):
        k += 1
    return k

def urgetter(h):
    '''
    Given a polynomial 
        h(x,y) = y ** d + b_2(x)y**(d-2) + ... + b_d(x)
    this function returns the tuple r, u_r where u_r is the
    degree of b_r such that u_r/r is minimal.
    '''
    list_of_bs = componentsOfy(h)
    list_of_bs = list_of_bs[::-1] # To follow the convention of the proof.
    list_of_least_degrees = [leastDegreeInx(b) for b in list_of_bs[2:]] #Here, position 0 is really from b2.
    list_of_quotients = [k/(i+2) for (i, k) in enumerate(list_of_least_degrees)]
    r_indexer = list_of_quotients.index(min(list_of_quotients))
    return r_indexer + 2, list_of_least_degrees[r_indexer]
