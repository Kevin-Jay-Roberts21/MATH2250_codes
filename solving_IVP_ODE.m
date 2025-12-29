% The following 3 lines are clearing all scripts and runnning program to 
% ensure that when you hit the "run" button, the only code running is from 
% this script.
clear all
close all
clc

% here we are defining a function y(t). You must define a function such as
% this if you wish to do differential operations on it.
syms y(t)

% here I am defining y' as dy. In english, I am saying that I want to take
% the derivative of y with respect to t one time (only one derivative, can 
% be replaced with a 2 or 3 if I want to define y'' or y''' for example).
dy = diff(y, t, 1)

% here I am asking matlab to do all the computation. I am giving it an
% equation (y' - 3y = 2) and I'm telling it to solve for the exact solution
% if the initial condidtion is y(0) = 1. (Note that we can solve for the
% general solution just as well by leaving out the y(0) == 1).
sol = dsolve(dy - 3*y == 2, y(0) == 1)


% At this point in the program, I have my exact solution which I have
% defined as sol. Now I want to plug in a value for t. You can do this by
% using the subs() function. Here I say "for sol, plug in log(2) for the t
% variable". And finally I wrap it with vpa() so that the answer, if
% simplified to a digit value (as opposed to somehting like exp(5) * 4 for 
% example).
answ = vpa(subs(sol, t, log(2)))


A = [-22 84; 
     -6 23]
 
[S, V, W] = eig(A)


 


