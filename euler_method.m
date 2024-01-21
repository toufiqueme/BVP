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
