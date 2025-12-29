clear all 
close all 
clc

image = imread('flower_image.jpg');
imshow(image);

gray_image = im2gray(image);
imshow(gray_image);

[U,S,V] = svdsketch(double(gray_image),1e-2);
A = uint8(U*S*V');
imshow(A);
u = U;
s = S;
v = V';

col1_u = u(:,1);
row1_v = v(1,:);
s11 = s(1, 1);
first_rank1_matrix = s11*col1_u*row1_v;

operator_norm = s(1, 1)

rank_approx_choice = 80;

for i = 2:rank_approx_choice
    col1_u = u(:,i);
    row1_v = v(i,:);
    s11 = s(i, i);
    
    this_matrix = s11*col1_u*row1_v;
    
    first_rank1_matrix = first_rank1_matrix + this_matrix;
end

imshow(uint8(first_rank1_matrix))


% singular_values = zeros(50, 1);
% for i = 1:50
%    singular_values(i, 1) = s(i, i);
% end
% 
% t = 1:50;
% plot(t, singular_values)
% legend('Singular Values')
% xlabel('Column of Singular Value'), ylabel('Value of Singular Values')
% title('Singular Value Plot')





