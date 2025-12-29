clear all
close all
clc

syms y(t) p r x(t) f(t) a b

% APPLICATIONS

%%%%%%%%%%%
% NEWTONS %
%%%%%%%%%%%
% 1) 

% A = 76.80
% dy = diff(y, t, 1)
% y0 = 44.80
% t0 = 0
% t1 = 16 
% T1 = 53.93
% T2 = 59.40
% sol = dsolve(dy == r*(A - y), y(t0)==y0)
% 
% rate_r = solve(subs(sol, t, t1)==T1, r)
% new_sol = dsolve(dy == rate_r*(A - y), y(t0)==y0)
% time = vpa(solve(new_sol == T2, t))

% 2) 
% A = 80.00
% T1 = 95.00
% T0 = 98.6 
% t0 = 0
% t1 = 14
% T2 = 89.86
% dy = diff(y, t, 1)
% 
% sol = dsolve(dy == r*(A - y), y(t0)==T1)
% 
% rate_r = vpa(solve(subs(sol, t, t1)==T2, r))
% 
% time = vpa(solve(subs(sol, r, rate_r)==T0, t))

%%%%%%%%%%%%%%
% RETIREMENT %
%%%%%%%%%%%%%%
% 1)
% dy = diff(y, t, 1)
% y0 = -5000
% r = 0.059
% monthly_pay = 82.46
% 
% sol = simplify(dsolve(dy == r*y + monthly_pay*12, y(0)==y0))
% 
% solve_for_t = vpa(solve(sol==0, t)) * 12

% 2)
% dy = diff(y, t, 1)
% y0 = -5750
% r = 0.1176
% monthly_pay = 77.65
% sol = simplify(dsolve(dy == r*y + monthly_pay*12, y(0)==y0))
% 
% solve_for_t = vpa(solve(sol==0, t))
% 
% money = solve_for_t*12*monthly_pay
% 
% answ = vpa(money+y0)

% 3)
% dy = diff(y, t, 1)
% years = 41
% todays_dollars = 74000
% inflation_over_years = 0.045
% interest = .1309
% 
% sol = simplify(dsolve(dy == inflation_over_years*y, y(0)==todays_dollars))
% money_in_years = vpa(subs(sol, t, years))
% 
% sol1 = simplify(dsolve(dy == interest*y + p*12, y(0)==0))
% sol2 = subs(sol1, t, years)
% 
% answ = vpa(solve(sol2==money_in_years/interest, p))

%%%%%%%%%%%
% SPRINGS %
%%%%%%%%%%%

% 1)
% dy = diff(y, t, 1)
% d2y = diff(y, t, 2)
% y0 = 0
% dy0 = 28
% a = 101
% c = 396
% b = sqrt(4*a*c)
% p(t) = 0
% sol = simplify(dsolve(a*d2y + b*dy + c*y == p(t), y(0)==y0, dy(0)==dy0))
% 
% dsol = diff(sol, t, 1)
% 
% time_t = vpa(solve(dsol == 0, t))
% peak = vpa(subs(sol, t, time_t))

% 2)
% dy = diff(y, t, 1)
% d2y = diff(y, t, 2)
% y0 = -2
% dy0 = 24
% a = 1
% b = 4
% c = 13
% f(t) = 0
% sol = simplify(dsolve(a*d2y + b*dy + c*y == f(t), y(0)==y0, dy(0)==dy0))
% 
% dsol = diff(sol, t, 1)
% 
% time_t = vpa(solve(dsol == 0, t))
% peak = vpa(subs(sol, t, time_t))

%%%%%%%%%%%
% BUNNIES %
%%%%%%%%%%%

% eqn1 = (450 - 6*a - 3*b)
% eqn2 = (800 -4*a - 12*b)
% 
% this_b = solve(eqn1 == eqn2, b)
% x = solve(subs(eqn1, b, this_b)==0, a)




%%%%%%%%%%%%%
% ANIMATION %
%%%%%%%%%%%%%

% 1)
% point1 = [-1; 3]
% point2 = [3; 4]
% point3 = [4; 8]
% avg_x = (point1(1,:) + point2(1,:) + point3(1,:))/3
% avg_y = (point1(2,:) + point2(2,:) + point3(2,:))/3
% theta = -pi/2
% matrix1 = [1 0 -avg_x;
%            0 1 -avg_y;
%            0 0 1]
% matrix2 = [cos(theta) -sin(theta) 0; 
%            sin(theta) cos(theta) 0; 
%            0 0 1]
% matrix3 = [1 0 avg_x;
%            0 1 avg_y;
%            0 0 1]
% final_matrix = matrix3*matrix2*matrix1

% 2) (RECALL that clockwise and counterclockwise are different matrices!)
% point1 = [6; 2]
% point2 = [3; 0]
% point3 = [3; 4]
% avg_x = (point1(1,:) + point2(1,:) + point3(1,:))/3
% avg_y = (point1(2,:) + point2(2,:) + point3(2,:))/3
% theta = pi/3
% matrix1 = [1 0 -avg_x;
%            0 1 -avg_y;
%            0 0 1]
% matrix2 = [cos(theta) sin(theta) 0; 
%            -sin(theta) cos(theta) 0; 
%            0 0 1]
% matrix3 = [1 0 avg_x;
%            0 1 avg_y;
%            0 0 1]
% final_matrix = round(vpa(matrix3*matrix2*matrix1), 2)


%%%%%%%%%%%
% TURTLES %
%%%%%%%%%%%
% 1)
% M = [1 0 0; 
%      0.16 0.92 0; 
%      0 0.08 0.96]
%  
% u = [73; 0; 0]
% 
% yrs = 2023 - 1127
% 
% ans = sum(M^yrs*u)

% 2) 
% M = [1 0 0 0; 
%      0.3 0.7 0 0; 
%      0 0.3 0.7 0; 
%      0 0 0.3 0.75]
% 
% ans = inv(0.3) + inv(0.3) + inv(0.3) + inv(0.25)

% 3) 
% M = [1.01 0.06 0; 
%      0 0.94 0.08; 
%      0.03 0 0.91]
%  
% u = [1; 99; 321]
% SRI = M^999 * u
% I = SRI(3,:)
% answ = vpa(SRI(3,:)/sum(SRI))

%%%%%%%%%%%%%%
% REGRESSION %
%%%%%%%%%%%%%%

% 1)
A = [(-2)^2 exp(-2);
     (0)^2 exp(0);
     (1)^2 exp(1);
     (3)^2 exp(3)]
y = [16; 0; 1; 21]

x = inv(A'*A)*A'*y

%%%%%%%
% SVD %
%%%%%%%

% 1) 
% v1 = [sqrt(2)/2;
%       sqrt(2)/2]
% w1 = [sqrt(6)/6 sqrt(2)*sqrt(3)/3 sqrt(6)/6]
% l1 = 3*sqrt(6)
% 
% R1 = l1*v1*w1

% 2)
% A = [7 3 -7 -2 6;
%      4 4 -6 -5 -3;
%      -5 -4 -7 4 9;
%      -1 -3 6 7 -4];
% 
% [U E V] = svd(A)
% V = V'
% 
% % find out sum of singular values that make up at least 75%
% e1 = E(1,1);
% e2 = E(2,2);
% e3 = E(3,3);
% e4 = E(4,4);
% sum = e1+e2+e3+e4;
% 
% percentage = (e1 + e2 + e3)/sum
% 
% 
% u1 = U(:,1)
% v1 = V(1,:)
% u2 = U(:,2)
% v2 = V(2,:)
% u3 = U(:,3)
% v3 = V(3,:)
% u4 = U(:,4)
% v4 = V(4,:)
% 
% R1 = e1*u1*v1
% R2 = e2*u2*v2
% R3 = e3*u3*v3
% R4 = e4*u4*v4
% 
% R = R2 + R1 + R4

% 3)
% Define matrices
% A = [-10 4; 2 9];
% B = [-7 10; 0 -7];
% C = [1 -8; -1 7];
% D = [8 7; -4 0];
% O = [9 2; -8 10];
% 
% % Step 1: Flatten each matrix column-wise and store as rows
% W = [A(:)'; B(:)'; C(:)'; D(:)'];  % 4x4 matrix, each row = one sample
% 
% % Step 2: Transpose W so that each column is a sample
% W = W';  % Now 4x4: each column is a sample
% 
% % Step 3: Subtract the mean from each row (i.e., mean of each feature)
% mu = mean(W, 2);             % 4x1 mean vector (feature-wise)
% W_zero_mean = W - mu;        % Subtract mean column-wise
% 
% % Step 4: Compute covariance of the centered data
% covW = cov(W_zero_mean');    % Samples are rows → transpose first
% 
% % Step 5: Perform SVD and extract dominant direction
% [U, S, V] = svd(covW);
% u1 = U(:,1);  % First left singular vector
% 
% % Step 6: Project each column (i.e., each matrix) into 1D
% projA = u1' * W_zero_mean(:,1);
% projB = u1' * W_zero_mean(:,2);
% projC = u1' * W_zero_mean(:,3);
% projD = u1' * W_zero_mean(:,4);
% 
% % Step 7: Project matrix O
% O_vec = O(:);                  % 4x1 vector
% O_centered = O_vec - mu;       % Subtract same mean
% projO = u1' * O_centered;
% 
% % Step 8: Compute distances in 1D
% distances = [abs(projA - projO);
%              abs(projB - projO);
%              abs(projC - projO);
%              abs(projD - projO)];
% 
% [min_dist, idx] = min(distances);
% 
% fprintf('Minimum distance is %.4f from matrix %c\n', min_dist, 'A' + idx - 1);



% 4)
A = [-10 4; 2 9];
B = [-7 10; 0 -7];
C = [1 -8; -1 7];
D = [8 7; -4 0];
O = [9 2; -8 10];

W = [A(:,1)' A(:,2)';
     B(:,1)' B(:,2)';
     C(:,1)' C(:,2)';
     D(:,1)' D(:,2)']

mu = mean(W) % Calculating the mean of each column of data (the dimensions) 
 
Ones = ones(1,4) % Subtract the mean from each column of data
A1 = W - Ones'*mu

covW = cov(W)

[U E V] = svd(covW)
V= V'

A1 = W(1,:) % row 1 of W
B1 = W(2,:)
C1 = W(3,:)
D1 = W(4,:)

u1 = U(:,1) % Extracting the first vector of U
O = [O(:,1)' O(:,2)']

distA = u1'*A1'
distB = u1'*B1'
distC = u1'*C1'
distD = u1'*D1'
distO = u1'*O'

all_distances = [abs(distA - distO);
                 abs(distB - distO);
                 abs(distC - distO);
                 abs(distD - distO)]

ans = min(all_distances)













