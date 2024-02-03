function x = gauss_seidel(A, b, x_initial, max_iter, tol)
    % A: Coefficient matrix
    % b: Right-hand side vector
    % x_initial: Initial solution vector
    % max_iter: Maximum number of iterations
    % tol: Tolerance for convergence
    
    n = length(b);
    x = x_initial;
    
    for k = 1:max_iter
        x_old = x;
        
        for i = 1:n
            x(i) = (b(i) - A(i, 1:i-1)*x(1:i-1) - A(i, i+1:end)*x_old(i+1:end)) / A(i, i);
        end
        
        % Check for convergence
        if norm(x - x_old, inf) < tol
            break;
        end
    end
    
    if k == max_iter
        fprintf('Gauss-Seidel did not converge within the specified number of iterations.\n');
    else
        fprintf('Gauss-Seidel converged in %d iterations.\n', k);
    end
    disp(x)
end
