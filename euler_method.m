function y = euler_method(n, t0, t1, y0)
    h = (t1 - t0) / n;
    t(1) = t0;
    y(1) = y0;

    for i = 1:n
        t(i + 1) = t(i) + h;
        y(i + 1) = y(i) + h * (y(i) - t(i)^2 + 1);  
    end

    exact_solution = (t + 1).^2 - 0.5 * exp(t);  
    
    error = abs(y - exact_solution);  
    
    T = table(t', exact_solution', y', error', 'VariableNames', {'t', 'ExactSolution', 'NumericalSolution', 'Error'});
    
    disp('Table of Values:');
    disp(T);
    
    plot(t, y, 'b', t, exact_solution, 'r--');
    legend('Euler method', 'Exact Solution');
    title('Euler method vs Exact Solution');
    xlabel('t');
    ylabel('y');
end
% function y = test(n, t0, t1, y0, y_prime0)
%     h = (t1 - t0) / n;
%     t(1) = t0;
%     y(1) = y0;
%     y_prime(1) = y_prime0;
% 
%     for i = 1:n
%         t(i + 1) = t(i) + h;
%         y(i + 1) = y(i) + h * y_prime(i);
%         y_prime(i + 1) = y_prime(i) + h * (-y(i) - y_prime(i) - t(i));
%     end
% 
%     T = table(t', y', y_prime', 'VariableNames', {'t', 'NumericalSolution', 'NumericalSolutionPrime'});
%     
%     disp('Table of Values:');
%     disp(T);
%     
%     plot(t, y, 'b', t, y_prime, 'g--');
%     legend('y(t)', "y'(t)");
%     title('Numerical Solution using Euler method');
%     xlabel('t');
%     ylabel('y');
% end
