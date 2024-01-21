function y = runge_kutta_4th_order(n, t0, t1, y0)
    h = (t1 - t0) / n;
    t(1) = t0;
    y(1) = y0;

    for i = 1:n
        t(i + 1) = t(i) + h;
        
        k1 = h * (y(i) - t(i)^2 + 1);
        k2 = h * ((y(i) + 0.5 * k1) - (t(i) + 0.5 * h)^2 + 1);
        k3 = h * ((y(i) + 0.5 * k2) - (t(i) + 0.5 * h)^2 + 1);
        k4 = h * ((y(i) + k3) - (t(i) + h)^2 + 1);
        
        y(i + 1) = y(i) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    exact_solution = (t + 1).^2 - 0.5 * exp(t);  
    
    error = abs(y - exact_solution);  
    
    
    T = table(t', exact_solution', y', error', 'VariableNames', {'t', 'ExactSolution', 'NumericalSolution', 'Error'});
    
    disp('Table of Values:');
    disp(T);
    
    plot(t, y, 'b', t, exact_solution, 'r--');
    legend('Runge-Kutta 4th order', 'Exact Solution');
    title('Runge-Kutta 4th order vs Exact Solution');
    xlabel('t');
    ylabel('y');
end
