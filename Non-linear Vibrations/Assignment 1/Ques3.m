clc;
clear;

A = zeros(5,5); % initializing the A matrix

for j=1:5
    A(:,j) = sin(1:5).*cos(1:5); % substituting the values
end

syms x;
polyA = charpoly(A,x); % charactyeristic polynomial of the matrix

[v,d] = eig(A); 
% columns of v are the eigenvectors and the digonal elements of d are the eigwevalues

v(:,1)

for i=1:5
fprintf('\\begin{bmatrix} \n');
fprintf('%0.4f\\\\ ', v(:,i));
fprintf('\n \\end{bmatrix},\n');
end

i =1
