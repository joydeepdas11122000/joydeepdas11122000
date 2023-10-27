import sympy as smp
import scipy.integrate as spi

x = smp.symbols('x')
a, b = list(map(int, input('enter Semi-Major & Semi-Minor Axis Length : ').split(',')))
a, b = max(a, b), min(a, b)

g = (b / a) * ((a ** 2 - x ** 2) ** 0.5)
gd1 = smp.diff(g, x)
# def f(value):                         # used for lambdifying the sympy expressions
#     return g.subs(x, value)
# def f1(value):
#     return g1.subs(x, value)
f = smp.lambdify(x, g)  # Returns the value at the point, also called lamdifying
fd1 = smp.lambdify(x, gd1)  # Returns the Slope at the point
p = lambda x: (1 + (fd1(x)) ** 2) ** 0.5
############################################
c = (a ** 2 - b ** 2) ** 0.5
e = (1 - (b ** 2 / a ** 2)) ** 0.5
Asp = 4 * spi.quad(f, 0, a)[0]
Psp = 4 * spi.quad(p, 0, a)[0]
###########################################
print(f'Eccentricity = {e:.3f}')
print(f'Focal Length = {c:.3f} unit')
print(f'Area by formula     = {3.141592653589793 * a * b:.3f} unit\u00b2')
print(f'Area by integration = {Asp:.3f} unit\u00b2')
print(f'Perimeter by Integration= {Psp:.3f} unit')
