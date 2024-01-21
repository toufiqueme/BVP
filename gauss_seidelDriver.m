
A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1;0 3 -1 8];
b=[6;25;-11;15];

max_iter = 1000;
tol = 1e-4;

x_solution = gauss_seidel(A, b, max_iter, tol);
disp('Solution:');
disp(x_solution);
