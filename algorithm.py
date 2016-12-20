import sympy
import random
from lemmas import HenselsLemma

x = sympy.Symbol('x')
y = sympy.Symbol('y')

def monicMaker(h n_limit = 10):
	'''
	This function peforms a rotation on a polynomial of two variables
	h(x,y). That is, it maps h(x,y) to h(x+ny, -nx + y)) for some n
	such that the resulting polynomial is monic in y.
	'''
	new_h = h
	while sympy.LC(new_h, y) not in sympy.CC:
		n = random.randint(1, n_limit)
		new_h = h.as_expr().subs([(x, x + n*y), (y, -n*x + y)])
	new_h = new_h * (1 / sympy.LC(new_h, y))
	return new_h.as_poly()

def newtonAutomorphism(h, q, p):
	'''
	This function performs the Newton Automorphism (x maps to x ** q
	and y maps to y(x ** p)) to a polynomial h, returning a sympy expr
	(not a polynomial). Here, p and q are expected to be rational
	numbers.
	'''
	return h.as_expr().subs([(x, x ** q), (y, y*x**p)])

def componentsOfy(F):
	'''
	This function returns a list with the homogenous
	components of h as a polynomial in y in increasing order.
	'''
	list_of_Fprimes = [F]
	list_of_fs = [F.eval(x, 0)]
	for j in range(1, sympy.degree(F, x)+1):
		nextFprime = sympy.diff(list_of_Fprimes[j-1], x)
		list_of_Fprimes.append(nextFprime)
		list_of_fs.append((1/math.factorial(j))*nextFprime.eval(x, 0))
	return list_of_fs

def phiAutomorphism(h):
	list_of_bs = componentsOfy(h)
	b1 = list_of_bs[-2]
	return h.as_expr().subs([(y, y-b1/sympy.degree(h, y))]).as_poly()
