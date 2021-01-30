clc;
clear;

A = zeros(5,5); % initializing the A matrix

for i = 1:5
    for j=1:5
        A(i,j) = sin(i)*cos(j); % substituting the values
    end
end

syms x;
polyA = charpoly(A,x); % charactyeristic polynomial of the matrix

[v,d] = eig(A); 
% columns of v are the eigenvectors and the digonal elements of d are the eigwevalues

v(:,1)

for j=1:5
    fprintf('\\begin{bmatrix} \n');
    for i=1:5
        fprintf('%0.4f + %0.4fi\\\\ ', real(v(i,j)), imag(v(i,j)));
    end
    fprintf('\n \\end{bmatrix},\n');
end

i =1
