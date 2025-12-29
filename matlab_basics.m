clear all
close all
clc


%%%%%%%%%%%%%%%
% READ  MEEEE %
%%%%%%%%%%%%%%%
% To run a certain section, click on a section of code and it'll be
% highlihgted yellow. Then above, instead of clicking the big green arrow
% with "Run" underneath it. Click on the "Run Section" arrow instead.




%%%%%%%%%%%%%%%%%%%
% Chapter 1 and 2 %
%%%%%%%%%%%%%%%%%%%
% Mostly we deal with differential equations here. Let's jump right into an
% example.

%%
syms y(t) f(t) p(t)% defining aribitrary functions and variables

% dy is essentially defining y'. The first parameter of diff() is the
% function you want to take the derivative of, the second parameter is what
% you'd like to take the derivative with respect to, and the third
% parameter is denoting how many times you'd like to take the derivative of
% the function
dy = diff(y, t, 1)

% defining an initial condition
y0 = 2

% defining the forcing function and p(t) function
f(t) = 2*t % (make sure you add the * symbol for multiplication)
p(t) = 4*t

% solving the differential equation (ignoring initial condition, the C1 in 
% the output is simply the constant left over) be sure to add the double 
% equals sign
solution1 = dsolve(dy + p(t)*y == f(t))

% solving the differential equation (with initial condition, we don't always
% have to do y(0), we could do y(-2) for example)
solution2 = dsolve(dy + p(t)*y == f(t), y(0) == y0)
%%
% You can apply the exact same process for second order differential
% equations here's another example:

syms y(t) f(t) p(t)% defining aribitrary functions and variables
a = 5
b = 10
c = 4
f(t) = cos(3*t)
d2y = diff(y, t, 2)
dy = diff(y, t, 1)
y0 = 10
dy0 = 20

solution3 = dsolve(a*d2y + b*dy + c*y == f(t), y(0)==y0, dy(0)==dy0)
%%
%%%%%%%%%%%%%%%
% Some basics %
%%%%%%%%%%%%%%%

syms x a b % defining arbitrary variables

% Sometimes we're given a large polynormial or a complicated equation where
% we want to solve for a certain variable. Here's a few examples of how to 
% do so:

% The second parameter is which variable I'm asking matlab to solve for and
% the first parameter is the equation
solution1 = solve(x^2 - 9 == 0, x) % remember to use two equals signs
solutoin2 = solve(x^4 - 2*x^3 - 5*x^2 + 8*x + 4 == 0, x) % Remember to add * when multiplying
solution3 = solve(2*b + a*sqrt(b) == a^2, a)

%%

% Integration
% Integrating a function isn't too complex, first of course we have to
% define some arbitrary variables
syms t s

% new define a function here
f = t*cos(5*(t-s))*cos(s)

% Now to integrate, the first parameter is the function you want to
% integrate, the second parameter is what you'd like to integrate with
% respect to, the third parameter is the lower bound and the fourth
% parameter is the upper bound. In the following I'll integrate f with
% respect to s, from 0 to t

answer = int(f, s, 0, t)

% here is a simpler one (where vpa approximates to decimal points)
g = exp(5 - 6*t)

answer1 = vpa(int(g, t, 0, log(5)))


%%

% Substition! Given a function y(t), if the function is very complicated,
% it can be laborsome to plug in a value t and evaluate by hand. Here's how
% to do is in matlab!

syms y(t)

y(t) = exp(3*t)*t + cos(4*t) + exp(sqrt(t) - 5*t)

% plugging in t = 4*ln(2). The vpa() function simplifies the answer to a
% decimal. For the subs() function, the first parameter is the function,
% the second parameter is the variable we want to plug in a value for, and 
% third parameter is the value for t that you plug in.
solution = vpa(subs(y, t, 4*log(2)))

% ANOTHER NOTE: if you'd like to round to a certain decimal place using
% vpa(), then you can you the round() function. Notice the following rounds
% the variable to 2 decimal places:
z = 57.8932
rounded_z = round(z, 2)

%%
%%%%%%%%%%%%%%%%%%
% LINEAR ALGEBRA %
%%%%%%%%%%%%%%%%%%

% Here's a few examples of how to define a matrix (the brackets enclose the 
% entire matrix or vector, and the semi colons start a new row)

% The following is a 3x3 matrix
A = [1 2 3; 
     4 5 6; 
     7 8 9]

% This is a 1x4 vector
v1 = [8 9 10 11]

% This is a 5x1 vector
v2 = [3; 4; 2; 7; 9]

% Say I want to get just the first row, or perhaps the second row of A:
row1_of_A = A(1,:)
row2_of_A = A(2,:)

% Say I want to get the third column of A: 
col3_of_A = A(:,3)

% Finally, to get any element of A, you simply need to type in it's
% position A(a,b) (ath row and bth col element)
element_1_1_of_A = A(1,1)
element_2_3_of_A = A(2,3)

%%
% Now let's get into some matrix operations. Remember the dimensions and
% positions of matrix operations matter. Notice:

A = [2 3 4; 
     5 4 1]
b = [1; 1; 1]

Ab = A*b % this possible, BUT b*A is not

% Now what about the transpose of a matrix. Two ways: transpose(A) or A'
AT1 = transpose(A)
AT2 = A'

% What about the inverse of a matrix? Just use the inv() operation.
W = [1 5; 
     2 9]

W_inverse = inv(W)
%%
% What if the inverse isn't so nice and we're left with decimals when we
% want fractions? Include "format rat" at the top of your file to force the
% document to have only factions and no decimals. Notice the difference:
A = [1 2; 
     5 4]
A_inv = inv(A)

%%
% (once this is ran, the fractions seem to stay. If you want 
% to get rid of the fractions, try commenting this out, closing the
% software and opening it again. OR just add vpa() to get decimals)
% format rat
A = [1 2; 
     5 4]
A_inv = inv(A)

%%
% Now what about eigen values and eigenvectors? This is tricky because
% matlab will default to have all of the eigevectors have a magnitude of 1.
% This makes the eigenvectors looks ugly, but nevertheless matlab can do it

A = [2 3 4; 
     1 5 3; 
     4 2 1]
 
% To find just the eigenvalues of A, use the eig() function
eigen_values_of_A = eig(A)

% To find the eigenvalues and eigenvectors of A, we must do a funky
% operation. Below, the matrix V is the eigenvector matrix and E is the
% eigenvalue matrix. (Again, use vpa is you want decimals instead)
[V, E] = eig(A)

% To perform complete diagonalization (i.e. to get the inverse of the
% eigenvector matrix), add another variable:
[V, E, Vinv] = eig(A)

%%
% I think it's important to note that given the following matrix A, when we
% want to perform exponentials to matrices, matlab automatically performs 
% the matrix diagonalization to reduce computation time:
A = [1 2; 
     3 2]
 
A_to_the_power_of_10 = A^10
% special note: sometimes matlab will output * if the value is too large or
% too small (try computing A^100 and you will see this)

%%
% SVD (Single Value Decomposition)
% We can do this with a whole bunch of small operations, but matlab have a
% very nice feature: the svd() operation. U is the left singular vector
% matrix, E is the singlular value matrix, and V is the right singular matrix.
A = [-2 9 7 -7 7; 
     -6 -10 -7 -3 10; 
     7 9 -9 -1 -3; 
     4 5 8 9 8]
[U E V] = svd(A)

% SPECIAL NOTE: We have to transpose V to get the complete svd (so redefine
% V):
V = V'







