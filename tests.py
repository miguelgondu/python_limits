import sympy
from algorithm import *

x = sympy.Symbol('x')
y = sympy.Symbol('y')

def test_hgetter1():
	f = sympy.Poly(x, x, y)
	g = sympy.Poly(y, x, y)

	h = hgetter(f, g)
	assert h == 2*x*y

def test_hgetter2():
	f = sympy.Poly(x**2 + y**2)
	g = sympy.Poly(2*y**3 + 4*x**2*y)

	h = hgetter(f, g)
	assert h == 4*x**4*y - 2*x**2*y**3 + 2*y**5

def test_monicMaker():
	h = 4*x**4*y - 2*x**2*y**3 + 2*y**5
	h1 = monicMaker(h)
	assert h1.as_poly(y).LC() == 1
	