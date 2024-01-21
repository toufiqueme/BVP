function x = gauss_seidel(A, b, max_iter, tol)
    % A: Coefficient matrix
    % b: Right-hand side vector
    % max_iter: Maximum number of iterations
    % tol: Tolerance for convergence
    
    n = length(b);
    x = zeros(n, 1);
    
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
        disp('Gauss-Seidel did not converge within the specified number of iterations.');
    else
        disp(['Gauss-Seidel converged in ', num2str(k), ' iterations.']);
    end
end
