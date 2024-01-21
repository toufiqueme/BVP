function y = adam_bashforth_adam_moulton(n, t0, t1, y0)
    h = (t1 - t0) / n;
    t(1) = t0;
    y(1) = y0;

    for i = 1:3  
        t(i + 1) = t(i) + h;
        k1 = h * (y(i) - t(i)^2 + 1);
        k2 = h * ((y(i) + 0.5 * k1) - (t(i) + 0.5 * h)^2 + 1);
        k3 = h * ((y(i) + 0.5 * k2) - (t(i) + 0.5 * h)^2 + 1);
        k4 = h * ((y(i) + k3) - (t(i) + h)^2 + 1);
        y(i + 1) = y(i) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    for i = 4:n
        t(i + 1) = t(i) + h;
        
        % Adams-Bashforth predictor
        predictor = y(i) + h * (55 * (y(i) - t(i)^2 + 1) - 59 * (y(i - 1) - t(i - 1)^2 + 1) + 37 * (y(i - 2) - t(i - 2)^2 + 1) - 9 * (y(i - 3) - t(i - 3)^2 + 1)) / 24;
        
        % Adams-Moulton corrector
        corrector = y(i) + h * (9 * (predictor - (t(i) + h)^2 + 1) + 19 * (y(i) - t(i)^2 + 1) - 5 * (y(i - 1) - t(i - 1)^2 + 1) + (y(i - 2) - t(i - 2)^2 + 1)) / 24;
        
        y(i + 1) = corrector;
    end

    exact_solution = (t + 1).^2 - 0.5 * exp(t);  
    
    error = abs(y - exact_solution);  % Error calculation
    
    
    T = table(t', exact_solution', y', error', 'VariableNames', {'t', 'ExactSolution', 'NumericalSolution', 'Error'});
    
    disp('Table of Values:');
    disp(T);
    
    plot(t, y, 'b', t, exact_solution, 'r--');
    legend('Adams-Bashforth-Adam-Moulton', 'Exact Solution');
    title('Adams-Bashforth-Adam-Moulton vs Exact Solution');
    xlabel('t');
    ylabel('y');
end
