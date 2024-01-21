
N = 100; 
L = 10;
h = L/N;
x = linspace(0, L, N+1); 
y = zeros(N+1, 1);

a = ones(N+1, 1);
b = -2*ones(N+1, 1);
c = ones(N+1, 1);
d = -h^2*(x.^2+1);


y(1) = 0;
y(N+1) = 0;


for i = 2:N
    y(i) = (d(i) - a(i)*y(i-1) - c(i)*y(i+1)) / b(i)
end

plot(x, y);
xlabel('x');
ylabel('y');
title('Solution of y''''+y''+y+1=0 using finite difference method');
