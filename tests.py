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

def test_phi_automorphism():
	h = (y**2 - sympy.Rational(80,9)*x*y - x**2).as_poly(y)
	h1, a = phiAutomorphism(h)
	bs = componentsOfy(h1)
	assert bs[-2] == 0

def test_newton_automorphism():
	h = 2*x*y
	new_h = newtonAutomorphism(h, sympy.Rational(1,2), sympy.Rational(1,3)) 
	assert new_h == 2*x**(sympy.Rational(5,6))*y

def test_least_degree_in_x():
	fx = x**2 + x
	assert leastDegreeInx(fx) == 1

def test_urgetter():
	h = y**5 + (x**2+1)*y**3 + (2*x**3)*y**2 + (x+1)*y + x**3 + 2
	r, ur = urgetter(h)
	assert ur == 0 and r==2