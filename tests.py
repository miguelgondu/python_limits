import sympy
from algorithm import *

x = sympy.Symbol('x')
y = sympy.Symbol('y')

def test_h_getter1():
	f = sympy.Poly(x, x, y)
	g = sympy.Poly(y, x, y)

	h = h_getter(f, g)
	assert h == 2*x*y

def test_h_getter2():
	f = sympy.Poly(x**2 + y**2)
	g = sympy.Poly(2*y**3 + 4*x**2*y)

	h = h_getter(f, g)
	assert h == 4*x**4*y - 2*x**2*y**3 + 2*y**5

def test_monic_maker():
	h = 4*x**4*y - 2*x**2*y**3 + 2*y**5
	h1 = monic_maker(h)
	assert h1.as_poly(y).LC() == 1

def test_phi_automorphism():
	h = (y**2 - sympy.Rational(80,9)*x*y - x**2).as_poly(y)
	h1, a = phi_automorphism(h)
	bs = components_of_y(h1)
	assert bs[-2] == 0

def test_newton_automorphism():
	h = 2*x*y
	new_h = newton_automorphism(h, sympy.Rational(1,2), sympy.Rational(1,3)) 
	assert new_h == 2*x**(sympy.Rational(5,6))*y

def test_least_degree_in_x():
	fx = x**2 + x
	assert least_degree_in_x(fx) == 1

def test_urgetter():
	h = y**5 + (x**2+1)*y**3 + (2*x**3)*y**2 + (x+1)*y + x**3 + 2
	r, ur = ur_getter(h)
	assert ur == 0 and r==2

def test_poly_mod():
	h = sympy.Poly(y**5 + y**4 + y + 3)
	assert poly_mod(h, 3) == y + 3

def test_poly_mod2():
	h = sympy.Poly(y**5 + y**4 + y + 3)
	assert poly_mod(h, 6) == h

def test_get_denominators1():
	h = (x**sympy.Rational(1,5) + y*x**sympy.Rational(3,4)
		+ y**2*x**sympy.Rational(3,10))
	assert set(get_denominators_of_x(h)) == set([5,4,10])

def test_get_denominators2():
	h = (x**sympy.Rational(1,1024) + y*x**sympy.Rational(1,2022)
		+ y**2*x**sympy.Rational(4,123))
	assert set(get_denominators_of_x(h)) == set([1024, 2022, 123])
