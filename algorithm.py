import sympy
import random
from lemmas import HenselsLemma

x = sympy.Symbol('x')
y = sympy.Symbol('y')

def monicMaker(h):
	'''
	This function peforms a rotation on a polynomial of two variables
	h(x,y). That is, it maps h(x,y) to h(x+ny, -nx + y)) for some n
	such that the resulting polynomial is monic in y.
	'''
	new_h = h
	while sympy.LC(new_h, y) != 1:
		n = random.randint(1,1000)
		new_h.subs([(x, x+n*y), (y, -n*x + y)])

	return new_h
