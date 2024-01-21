function result = sor_method(A, b, omega, initial_solution, max_iter, tol)
  
     n = length(b);
    x = initial_solution;
    iter_values = zeros(max_iter, n + 1);
    
    for k = 1:max_iter
        x_old = x;
        
        for i = 1:n
            sum1 = A(i, 1:i-1) * x(1:i-1);
            sum2 = A(i, i+1:end) * x_old(i+1:end);
            x(i) = (1 - omega) * x_old(i) + (omega / A(i, i)) * (b(i) - sum1 - sum2);
        end
        
        
        iter_values(k, :) = [x', norm(x - x_old, inf)];
        
        
        if norm(x - x_old, inf) < tol
            break;
        end
    end
    
    if k == max_iter
        disp('SOR method did not converge within the specified number of iterations.');
    else
        disp(['SOR method converged in ', num2str(k), ' iterations.']);
    end
    
   
    header = cell(1, n + 1);
    for i = 1:n
        header{i} = ['x', num2str(i)];
    end
    header{n + 1} = 'Error';
    
    result = array2table(iter_values(1:k, :), 'VariableNames', header);
end
