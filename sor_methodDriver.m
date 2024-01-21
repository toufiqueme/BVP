A = [4 3 0; 3 4 -1; 0 -1 4];
b = [24; 30; -24];
omega = 1.25;  
max_iter = 1000;
tol = 1e-4;

initial_solution = [1; 1; 1]; 

result_table = sor_method(A, b, omega, initial_solution, max_iter, tol);
disp('Table of Iteration Values:');
disp(result_table);
