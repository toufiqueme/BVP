function y = predictor_corrector_method(n, t0, t1, y0)
    h = (t1 - t0) / n;
    t(1) = t0;
    y(1) = y0;

    for i = 1:n
        t(i + 1) = t(i) + h;
        
        % Predictor step (using Adams-Bashforth 2nd order)
        f1 = h * (y(i) - t(i)^2 + 1);
        f2 = h * (y(i) + f1 - (t(i) + h)^2 + 1);
        predictor = y(i) + (3/2) * f2 - (1/2) * f1;
        
        % Corrector step (using Adams-Moulton 2nd order)
        corrector = y(i) + (h/2) * ((y(i) - t(i)^2 + 1) + (predictor - (t(i) + h)^2 + 1));
        
        y(i + 1) = corrector;
    end

    exact_solution = (t + 1).^2 - 0.5 * exp(t);  
    
    error = abs(y - exact_solution);  
    
    T = table(t', exact_solution', y', error', 'VariableNames', {'t', 'ExactSolution', 'NumericalSolution', 'Error'});
    
    disp('Table of Values:');
    disp(T);
    
    plot(t, y, 'b', t, exact_solution, 'r--');
    legend('Predictor-Corrector method', 'Exact Solution');
    title('Predictor-Corrector method vs Exact Solution');
    xlabel('t');
    ylabel('y');
end
