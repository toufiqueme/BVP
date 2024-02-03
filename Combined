function merge_methods_with_common_h_and_error()
    % Prompt user for input
    n = input('Enter the number of steps (n): ');
    t0 = input('Enter the initial time (t0): ');
    t1 = input('Enter the final time (t1): ');
    y0 = input('Enter the initial value (y0): ');

    % Common step size for all methods
    h = (t1-t0)/n;

    % Euler method
    [t_euler, y_euler] = euler_method(n, t0, t1, y0, h);
    exact_solution_euler = (t_euler + 1).^2 - 0.5 * exp(t_euler);
    error_euler = abs(y_euler - exact_solution_euler);

    % Modified Euler method
    [t_modified_euler, y_modified_euler] = modified_euler_method(n, t0, t1, y0, h);
    exact_solution_modified_euler = (t_modified_euler + 1).^2 - 0.5 * exp(t_modified_euler);
    error_modified_euler = abs(y_modified_euler - exact_solution_modified_euler);

    % Runge-Kutta 4th order
    [t_rk4, y_rk4] = runge_kutta_4th_order(n, t0, t1, y0, h);
    exact_solution_rk4 = (t_rk4 + 1).^2 - 0.5 * exp(t_rk4);
    error_rk4 = abs(y_rk4 - exact_solution_rk4);

    % Concatenate results and errors into a single table
    T = table(t_euler', y_euler', exact_solution_euler', error_euler', ...
              y_modified_euler', exact_solution_modified_euler', error_modified_euler', ...
              y_rk4', exact_solution_rk4', error_rk4', ...
              'VariableNames', {'t_Euler', 'NumericalSolution_Euler', 'ExactSolution_Euler', 'Error_Euler', ...
                                'NumericalSolution_ModifiedEuler', 'ExactSolution_ModifiedEuler', 'Error_ModifiedEuler', ...
                                'NumericalSolution_RK4', 'ExactSolution_RK4', 'Error_RK4'});

    % Displaying the combined table
    disp('Table of Values and Errors:');
    disp(T);

    % Plotting results
    figure;

    subplot(2,1,1);
    plot(t_euler, y_euler, 'b', t_modified_euler, y_modified_euler, 'g--', t_rk4, y_rk4, 'r-.', t_euler, exact_solution_euler, 'k:');
    legend('Euler method', 'Modified Euler method', 'Runge-Kutta 4th order', 'Exact Solution (Euler)');
    title('Numerical Solutions');
    xlabel('t');
    ylabel('y');

    subplot(2,1,2);
    plot(t_euler, error_euler, 'b', t_modified_euler, error_modified_euler, 'g--', t_rk4, error_rk4, 'r-.');
    legend('Error (Euler method)', 'Error (Modified Euler method)', 'Error (Runge-Kutta 4th order)');
    title('Errors');
    xlabel('t');
    ylabel('Absolute Error');
end

% Function definitions for each method

function [t, y] = euler_method(n, t0, t1, y0, h)
    t = zeros(1, n + 1);
    y = zeros(1, n + 1);
    t(1) = t0;
    y(1) = y0;

    for i = 1:n
        t(i + 1) = t(i) + h;
        y(i + 1) = y(i) + h * (y(i) - t(i)^2 + 1);  
    end
end

function [t, y] = modified_euler_method(n, t0, t1, y0, h)
    t = zeros(1, n + 1);
    y = zeros(1, n + 1);
    t(1) = t0;
    y(1) = y0;

    for i = 1:n
        t(i + 1) = t(i) + h;
        y_pred = y(i) + h * (y(i) - t(i)^2 + 1);
        y(i + 1) = y(i) + 0.5 * h * ((y(i) - t(i)^2 + 1) + (y_pred - (t(i) + h)^2 + 1));
    end
end

function [t, y] = runge_kutta_4th_order(n, t0, t1, y0, h)
    t = zeros(1, n + 1);
    y = zeros(1, n + 1);
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
end
