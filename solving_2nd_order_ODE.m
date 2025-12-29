clear all
close all
clc

syms y(t) f(t)

dy = diff(y, t, 1)
d2y = diff(y, t, 2)

y0 = 2
dy0 = 5

a = 2
b = 3
c = 1
f(t) = exp(3*t)

sol = simplify(dsolve(a*d2y + b*dy + c*y == f(t), y(0)==y0, dy(0)==dy0))

answer = vpa(subs(sol, t, log(3)))